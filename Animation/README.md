# Trajectory animation

A [Manim](https://www.manim.community/) scene that replays recorded neutron
paths in 2-D.

1. Record trajectories: in `main.py` set `track = True`, then run it from the
   repository root. It writes `all_trajectories.json` (one polyline per source
   neutron, descendants included). Keep the particle count small, since each
   history stores every vertex.
2. Copy or move `all_trajectories.json` into this folder.
3. Render: `manim -pqh manim_scene.py ParticleSimulation` (or `run.bat` on
   Windows).

The JSON is a generated file and is git-ignored. The shield outline (a
rectangle), its bounds and the scale factor are hard-coded in
`manim_scene.py`. Edit them to match the geometry you simulated. Making this self-describing is item 1 of the roadmap in
the main README.
