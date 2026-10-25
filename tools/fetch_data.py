"""Download the ENDF/B-VIII.0 nuclear data PyNeut's tests and validation use.

The files come from OpenMC's official ENDF/B-VIII.0 HDF5 library
(https://openmc.org/data-libraries/). The archive is streamed and only the
requested nuclides are written, to <dest>/neutron/<nuclide>.h5, so nothing
else touches the disk. Each file is then checked against the SHA-256 in
tools/endfb-viii.0.sha256: the snapshot the committed validation results were
produced with. A file that differs is reported, because a cross-code
comparison only holds for one data snapshot.

Usage
-----
  python tools/fetch_data.py                      # every nuclide the repo uses
  python tools/fetch_data.py --set fixed-source   # only the fixed-source suite
  python tools/fetch_data.py U235 Pu239           # named nuclides
  python tools/fetch_data.py --verify             # check what is already there

The archive is 3.4 GB and has to be read up to the last requested file, so
a fetch takes as long as downloading most of it, even though only ~0.8 GB is
kept.
"""
import argparse
import hashlib
import os
import sys
import tarfile
import time
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CHECKSUMS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "endfb-viii.0.sha256")
URL = "https://anl.box.com/shared/static/uhbxlrx7hvxqw27psymfbhi7bx7s6u6a.xz"

SETS = {
    # Validation/OpenMC_Comparison/run_all.py, the mixture and the tally checks
    "fixed-source": ["Al27", "B10", "Be9", "C12", "Cd112", "Fe56", "Ni58",
                     "O16", "Pb208"],
    # validate_keff.py: the U235 spheres and --benchmarks
    "keff": ["Ga69", "Ga71", "Pu239", "Pu240", "Pu241", "U233", "U234",
             "U235", "U238"],
    # validate_pincell.py: the BEAVRS pin cell
    "pincell": ["B10", "B11", "Cr50", "Cr52", "Cr53", "Cr54", "Fe54", "Fe56",
                "Fe57", "Fe58", "H1", "H2", "He3", "He4", "O16", "O17", "O18",
                "Sn112", "Sn114", "Sn115", "Sn116", "Sn117", "Sn118", "Sn119",
                "Sn120", "Sn122", "Sn124", "U234", "U235", "U236", "U238",
                "Zr90", "Zr91", "Zr92", "Zr94", "Zr96"],
    # the examples and tests add these
    "extra": ["Mn55", "Pb207"],
}


def load_checksums():
    sums = {}
    with open(CHECKSUMS) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                digest, name = line.split()
                sums[name[:-len(".h5")]] = digest
    return sums


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify(names, neutron_dir, sums):
    """Print one line per nuclide; return the names that are missing or differ."""
    bad = []
    for name in names:
        path = os.path.join(neutron_dir, f"{name}.h5")
        if not os.path.exists(path):
            print(f"  {name:<7} missing")
            bad.append(name)
        elif name not in sums:
            print(f"  {name:<7} present (no reference checksum)")
        elif sha256(path) != sums[name]:
            print(f"  {name:<7} DIFFERS from the validated snapshot")
            bad.append(name)
        else:
            print(f"  {name:<7} ok")
    return bad


class _Counter:
    """File-like wrapper that reports how much of the download has been read."""

    def __init__(self, raw):
        self.raw, self.n, self.t0, self.last = raw, 0, time.time(), 0.0

    def read(self, size=-1):
        data = self.raw.read(size)
        self.n += len(data)
        if time.time() - self.last > 5:
            self.last = time.time()
            rate = self.n / max(self.last - self.t0, 1e-9) / 1e6
            print(f"    ... {self.n / 1e9:.2f} GB read ({rate:.0f} MB/s)", flush=True)
        return data


def fetch(names, neutron_dir, url=URL):
    """Stream the archive and write each wanted <name>.h5 into neutron_dir."""
    wanted = set(names)
    os.makedirs(neutron_dir, exist_ok=True)
    print(f"Streaming {url}")
    with urllib.request.urlopen(url) as response:
        stream = _Counter(response)
        with tarfile.open(fileobj=stream, mode="r|xz") as archive:
            for member in archive:
                parts = member.name.split("/")
                if (not member.isfile() or len(parts) < 2
                        or parts[-2] != "neutron" or not parts[-1].endswith(".h5")):
                    continue
                name = parts[-1][:-len(".h5")]
                if name not in wanted:
                    continue
                # only the base name is used, so no path from the archive
                # can write outside neutron_dir
                final = os.path.join(neutron_dir, f"{name}.h5")
                partial = final + ".part"
                with archive.extractfile(member) as src, open(partial, "wb") as dst:
                    for block in iter(lambda: src.read(1 << 20), b""):
                        dst.write(block)
                os.replace(partial, final)
                wanted.discard(name)
                print(f"  wrote {name}.h5 ({member.size / 1e6:.1f} MB); "
                      f"{len(wanted)} to go", flush=True)
                if not wanted:
                    break
    return sorted(wanted)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("nuclides", nargs="*",
                        help="nuclides to fetch, e.g. U235 Pb208 (default: --set all)")
    parser.add_argument("--set", choices=[*SETS, "all"], default=None,
                        help="a predefined group; 'all' is every nuclide the repo uses")
    parser.add_argument("--dest", default=os.path.join(ROOT, "endfb"),
                        help="data directory (default: endfb/ in the repository)")
    parser.add_argument("--force", action="store_true",
                        help="download files that are already present and match")
    parser.add_argument("--verify", action="store_true",
                        help="only check the files already present; download nothing")
    args = parser.parse_args(argv)

    if args.nuclides and args.set:
        parser.error("give nuclide names or --set, not both")
    if args.nuclides:
        names = sorted(set(args.nuclides))
    elif args.set and args.set != "all":
        names = sorted(set(SETS[args.set]))
    else:
        names = sorted({n for group in SETS.values() for n in group})

    neutron_dir = os.path.join(args.dest, "neutron")
    sums = load_checksums()

    if args.verify:
        print(f"Checking {len(names)} nuclides in {neutron_dir}")
        return 1 if verify(names, neutron_dir, sums) else 0

    todo = names if args.force else [
        n for n in names
        if not (os.path.exists(os.path.join(neutron_dir, f"{n}.h5"))
                and (n not in sums
                     or sha256(os.path.join(neutron_dir, f"{n}.h5")) == sums[n]))]
    if todo:
        print(f"Fetching {len(todo)} of {len(names)} nuclides into {neutron_dir}")
        not_found = fetch(todo, neutron_dir)
        if not_found:
            print(f"Not in the library: {', '.join(not_found)}")
    else:
        print(f"All {len(names)} nuclides are already present.")

    print("Checking against the validated snapshot")
    bad = verify(names, neutron_dir, sums)
    if bad:
        print(f"{len(bad)} file(s) missing or different: {', '.join(bad)}")
        return 1
    print(f"All {len(names)} files match. Results will compare with the "
          "committed validation tables.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
