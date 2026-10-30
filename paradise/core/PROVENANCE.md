# Provenance of `paradise/core`

These modules are Poly-Time's `poly/` package, copied from
[CaliTarheel/Poly-Time](https://github.com/CaliTarheel/Poly-Time) at commit
`a5a3d2b372801eb66dca024aefe166d08c822951`. The import paths were changed from `poly` to
`paradise.core`, and nothing else, apart from the changes listed below.

| Module | What it is |
| --- | --- |
| `sphere.py` | Points, rotations and oriented S2 polygons; polygon repair; shared-edge simplification |
| `regions.py` | Plates as regions stepped under any rotations: crust ledger, boundary typing, sutures |
| `world.py`, `topology.py`, `build.py` | Crust pieces, the boundary network the regions world grew out of, test worlds |
| `mesh.py` | Continental thickness on a mesh that moves with its plate |
| `render.py` | Equirectangular map of crust by kind and ocean age |
| `dynamics.py` | The torque balance. Paradise does not use it to move plates; it is kept as a possible advisor |

Poly-Time's tests for these modules are in `tests/test_core_*.py`.

## Changes since the copy

- `sphere.loop_area_km2` fans each loop from the mean direction of its vertices instead of its
  first vertex. The fan from a vertex is degenerate when the loop passes through the point opposite
  it, which is what a plate painted from pole to pole does; such plates failed to rotate.
- `sphere.tidy_ring` ends on a ring that turns back at every vertex (a zig-zag along one line,
  enclosing nothing). Taking one of each two neighbouring spikes left nothing to take when every
  vertex was one, so the loop repeated unchanged for ever; a four-vertex sliver of accreted crust
  hung a history run inside `clean`. Every other spike goes instead.
- `ContinentMesh.spread` takes one implicit (backward Euler) step instead of sub-cycling an explicit
  one. On a world with collisions, small mesh triangles made the explicit scheme need so many
  cycles that it took most of each step's time. The implicit step is stable for any step and keeps
  the volume exactly, but it diffuses a little differently from the explicit one over the same time.
- `sphere.part` handles a plate covering most of the sphere (one ocean plate round a few
  continental ones). S2 writes such a polygon's shell as one hole's boundary, running clockwise.
  The area check on S2's own result now covers a clockwise shell, which encloses the sphere less
  its area; before, it was skipped, and a result missing a whole hole's worth was accepted. And the
  loop-by-loop rebuild counts the polygon's holes when it decides whether a shell loop is really an
  inner loop, since such a shell is larger than the polygon by its other holes.
- `RegionWorld._initiate` lets a world with `continents_override` set put a continental plate
  (`CONTINENTAL_SHARE` continent or more) over an oceanic one when a trench starts between them. The
  core's own rule looks 150 km either side, which on a continent's passive margin finds ocean on
  both and often sent the continental plate down.
- `RegionWorld._sutures` welds only two continental plates when the world sets
  `weld_continents_only`. An ocean plate that had picked up an accreted terrane otherwise welded to
  every continental plate it met, and one run welded the whole world into a single plate.
- `RegionWorld` reads two more world attributes when they are set, for rigid continents as the
  hand-built Artifexia project has them (its continental area never shrinks, and colliding
  continents meet without overlapping): `suture_km`, the closing after which a collision welds
  (and after which converging continents become a collision), and `escape_fraction`, the share of
  shortened continent that escapes sideways as bands along the coast (0: none).
- `RegionWorld.step` merges ocean into age bands of the world's `age_band_myr` when it has one
  (default `AGE_BAND_MYR`, 10 Myr). Bands 10 Myr apart made a thousand-Myr history carry five
  times the vertices of the hand-built Artifexia project, whose isochrons are its 50 Myr stages.
- `render.draw` writes to a file object as well as a path.
- `RegionWorld._initiate` reads the closing needed to start a trench from the world's
  `initiation_km` attribute when it has one (default `INITIATION_KM`), so a project can leave every
  trench to be declared.
- `RegionWorld._sutures` leaves alone any pair of plates in the world's `kept_apart` set, when it
  has one, so a person can let two plates collide without welding them.
