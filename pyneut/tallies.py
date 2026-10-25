"""Event-based tallies shared by both transport engines.

A `Scorer` holds every tally of a run. The engines report four kinds of event
to it and never consult it otherwise; it draws no random numbers, so passing
one cannot change a transport result.

    collision(hist, cell, x, y, z, E, w, sigma_t)   at each collision site
    track(hist, cell, length, E, w)                  for each flight segment
    absorb(hist, cell, x, y, z, E, w)                weight removed by absorption
    leak(hist, x, y, z, E, w)                        weight leaving the model

Arguments are arrays (vector engine) or scalars (scalar engine). `hist` maps
each event to its source history, which `set_batches` turns into a batch
index; the scalar engine passes `hist=None` and sets `scorer.batch` before
each history instead. E is the energy before the event and w the weight
before it, so a collision is scored ahead of roulette and implicit capture.

Each tally is a score with optional filters:

    score      flux (track-length or collision estimator), absorption, leakage
    filters    cells (medium indices), mesh (regular grid), energy_bins

Values accumulate per batch; `finalize` divides by the source neutrons of
each batch and reports the batch mean and its standard error. Flux is
volume-integrated (track length per source neutron, cm), as OpenMC reports it.
"""
import numpy as np

from .statistics import mean_sem

SCORES = ("flux", "absorption", "leakage")


class TallyDef:
    """One tally: a score plus optional filters.

    cell_map:   array over medium index -> output bin (-1 = not tallied), or None
    mesh:       {"min": [..], "max": [..], "dims": [nx, ny, nz]} or None
    energy_bins: increasing bin edges (eV) or None
    estimator:  "track-length" or "collision" (flux only)
    """

    def __init__(self, name, score, cell_map=None, cell_names=None, mesh=None,
                 energy_bins=None, estimator=None):
        if score not in SCORES:
            raise ValueError(f"unknown score {score!r}")
        if score == "flux":
            estimator = estimator or ("collision" if mesh else "track-length")
            if estimator not in ("track-length", "collision"):
                raise ValueError(f"unknown flux estimator {estimator!r}")
            if estimator == "track-length" and mesh is not None:
                raise ValueError("a mesh flux needs the collision estimator")
        elif estimator is not None:
            raise ValueError("only flux tallies take an estimator")
        if score == "leakage" and (cell_map is not None or mesh is not None):
            raise ValueError("leakage happens outside every cell; it takes energy_bins only")
        self.name = name
        self.score = score
        self.estimator = estimator
        self.cell_map = None if cell_map is None else np.asarray(cell_map, dtype=np.int64)
        self.cell_names = list(cell_names) if cell_names is not None else None
        self.energy_bins = None if energy_bins is None else np.asarray(energy_bins, dtype=float)
        self.mesh = None
        shape = []
        if self.cell_map is not None:
            shape.append(int(self.cell_map.max()) + 1)
        if mesh is not None:
            self.mesh = {"min": np.asarray(mesh["min"], dtype=float),
                         "max": np.asarray(mesh["max"], dtype=float),
                         "dims": np.asarray(mesh["dims"], dtype=np.int64)}
            shape.extend(int(d) for d in self.mesh["dims"])
        if self.energy_bins is not None:
            shape.append(self.energy_bins.size - 1)
        self.shape = tuple(shape)

    @property
    def axes(self):
        axes = []
        if self.cell_map is not None:
            axes.append("cell")
        if self.mesh is not None:
            axes += ["x", "y", "z"]
        if self.energy_bins is not None:
            axes.append("energy")
        return axes


class Scorer:
    def __init__(self, defs, n_batches):
        self.defs = list(defs)
        self.n_batches = int(n_batches)
        self.values = {d.name: np.zeros((self.n_batches,) + d.shape) for d in self.defs}
        self.totals = {d.name: np.zeros(self.n_batches) for d in self.defs}
        # sum of w*E for leakage tallies: the leaked-weight mean energy
        self.energy_sums = {d.name: np.zeros(self.n_batches) for d in self.defs
                            if d.score == "leakage"}
        self.batch = 0               # scalar engine: batch of the current history
        self.hist_batch = None       # vector engine: history index -> batch
        self.active = True           # eigenvalue: False during inactive generations
        self._by_event = {
            "collision": [d for d in self.defs if d.score == "flux" and d.estimator == "collision"],
            "track": [d for d in self.defs if d.score == "flux" and d.estimator == "track-length"],
            "absorb": [d for d in self.defs if d.score == "absorption"],
            "leak": [d for d in self.defs if d.score == "leakage"],
        }

    # ----------------------------------------------------------- engine side
    def wants(self, event):
        return self.active and bool(self._by_event[event])

    def set_batches(self, hist_batch):
        """Vector engine: `hist_batch[h]` is the batch of history h."""
        self.hist_batch = np.asarray(hist_batch, dtype=np.int64)

    def collision(self, hist, cell, x, y, z, E, w, sigma_t):
        if not self.wants("collision"):
            return
        w, sigma_t = np.atleast_1d(w).astype(float), np.atleast_1d(sigma_t).astype(float)
        value = np.where(sigma_t > 0, w / np.where(sigma_t > 0, sigma_t, 1.0), 0.0)
        for d in self._by_event["collision"]:
            self._add(d, hist, cell, (x, y, z), E, value)

    def track(self, hist, cell, length, E, w):
        if not self.wants("track"):
            return
        value = np.atleast_1d(w).astype(float) * np.atleast_1d(length).astype(float)
        for d in self._by_event["track"]:
            self._add(d, hist, cell, None, E, value)

    def absorb(self, hist, cell, x, y, z, E, w):
        if not self.wants("absorb"):
            return
        value = np.atleast_1d(w).astype(float)
        for d in self._by_event["absorb"]:
            self._add(d, hist, cell, (x, y, z), E, value)

    def leak(self, hist, x, y, z, E, w):
        if not self.wants("leak"):
            return
        value = np.atleast_1d(w).astype(float)
        energy = np.atleast_1d(E).astype(float)
        for d in self._by_event["leak"]:
            batch = self._batch(hist, value.size)
            np.add.at(self.energy_sums[d.name], batch, value * energy)
            self._add(d, hist, None, (x, y, z), E, value, batch=batch)

    # ------------------------------------------------------------- internals
    def _batch(self, hist, n):
        if hist is None:
            return np.full(n, self.batch, dtype=np.int64)
        return self.hist_batch[np.atleast_1d(hist)]

    def _add(self, d, hist, cell, xyz, E, value, batch=None):
        n = value.size
        if batch is None:
            batch = self._batch(hist, n)
        np.add.at(self.totals[d.name], batch, value)
        if not d.shape:
            np.add.at(self.values[d.name], batch, value)
            return
        keep = np.ones(n, dtype=bool)
        index = [batch]
        if d.cell_map is not None:
            cb = np.broadcast_to(d.cell_map[np.atleast_1d(cell)], (n,))
            keep &= cb >= 0
            index.append(cb)
        if d.mesh is not None:
            lo, hi, dims = d.mesh["min"], d.mesh["max"], d.mesh["dims"]
            for axis, coord in enumerate(xyz):
                c = np.broadcast_to(np.atleast_1d(coord).astype(float), (n,))
                i = np.floor((c - lo[axis]) / (hi[axis] - lo[axis]) * dims[axis]).astype(np.int64)
                keep &= (i >= 0) & (i < dims[axis])
                index.append(i)
        if d.energy_bins is not None:
            eb = np.searchsorted(d.energy_bins, np.broadcast_to(np.atleast_1d(E), (n,)),
                                 side="right") - 1
            keep &= (eb >= 0) & (eb < d.energy_bins.size - 1)
            index.append(eb)
        if keep.all():
            np.add.at(self.values[d.name], tuple(index), value)
        elif keep.any():
            np.add.at(self.values[d.name], tuple(i[keep] for i in index), value[keep])

    # ----------------------------------------------------------- bookkeeping
    def state(self):
        """Picklable accumulators, for returning from a worker process."""
        return {"values": self.values, "totals": self.totals,
                "energy_sums": self.energy_sums}

    def merge(self, state):
        for key in ("values", "totals", "energy_sums"):
            mine = getattr(self, key)
            for name, arr in state[key].items():
                mine[name] += arr

    def finalize(self, per_batch_sources):
        """Per-source-neutron batch means and standard errors for every tally.

        `per_batch_sources[b]` is the number of source neutrons in batch b.
        """
        n_b = np.asarray(per_batch_sources, dtype=float)
        out = {}
        for d in self.defs:
            norm = n_b.reshape((-1,) + (1,) * len(d.shape))
            per = self.values[d.name] / norm
            total_m, total_s = mean_sem(self.totals[d.name] / n_b)
            res = {"type": d.score, "total": total_m, "total_sem": total_s}
            if d.score == "flux":
                res["estimator"] = d.estimator
            if d.shape:
                mean = per.mean(axis=0)
                sem = (per.std(axis=0, ddof=1) / np.sqrt(per.shape[0])
                       if per.shape[0] > 1 else np.full_like(mean, np.nan))
                res.update(axes=d.axes, mean=mean, sem=sem)
                if d.cell_names is not None:
                    res["cells"] = d.cell_names
                if d.energy_bins is not None:
                    res["energy_bins"] = d.energy_bins
                if d.mesh is not None:
                    res["mesh"] = {k: v.tolist() for k, v in d.mesh.items()}
            if d.score == "leakage":
                t = self.totals[d.name]
                ratio = np.where(t > 0, self.energy_sums[d.name] / np.where(t > 0, t, 1.0), np.nan)
                res["mean_energy"], res["mean_energy_sem"] = mean_sem(ratio)
            out[d.name] = res
        return out
