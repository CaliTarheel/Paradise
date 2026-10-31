# GPlates to goSPL

A hand-built GPlates reconstruction, made the way Artifexian made Artifexia, drives a goSPL
landscape-evolution run. GPlates holds every decision (where plates rift, what collides, which
terranes accrete); this package only reads the project and writes goSPL's inputs, following
goSPL's own global plate-motion workflow (`Global-examples/plate_mvt` in goSPL-examples).

## What it makes

For a project directory and a span of ages (`python -m gplates_gospl.build --help`):

| File | Holds |
| --- | --- |
| `input/mesh.npz` | icosphere nodes `v` (m), triangles `c`, the first stage's elevation `z` |
| `input/forcing_NNN.npz` | per interval: `hdisp`, plate displacement (m/yr); `t`, vertical rate (m/yr) |
| `input/rain.npz` | `r`, rainfall (m/yr), 1 m/yr everywhere for now |
| `input/maps.npz` | every stage's map `z`, the surface goSPL's fast mode reaches (`expected`), plate IDs |
| `input.yml` | goSPL with erosion and deposition |
| `forcing-check.yml` | goSPL's `fast` mode: the forcing alone, no surface processes |

A landscape run goes through `run.py`, interval by interval, not through `input.yml` in one go
(see *Keeping the landscape on the reconstruction* below):

```bash
/opt/gospl-env/bin/python -m gplates_gospl.run BUILD_DIR RUN_DIR --np 4 --K 2e-7
```

goSPL time starts at 0 at the oldest stage and counts years forward.

- **Elevation** (`elevation.py`) comes from what each reconstructed feature is, by rules of thumb:
  - ocean depth from crust age (GDH1), the age being when the crust's strip formed;
  - continents, cratons, shelf and slope;
  - orogenies by age class (active, former, old);
  - island arcs, LIPs, trenches, rifts.
- **Plates** (`motion.py`) come from the crust under each node: the crust features carry the plate
  IDs the rotation file moves. Artifexia's plate outlines are labelled 2001, 2002... and never
  move, so the topologies are not used.
- **Vertical forcing** is what turns goSPL's own surface, carried by the plates, into the next
  map. goSPL re-interpolates the surface each time it carries it (`build.advected` copies how), and
  forcing worked out from the last map instead let that difference build up stage on stage.

`ARTIFEXIA` in `project.py` says which of Artifexia's layer files plays which part; another project
passes its own role map.

## Running goSPL in this container

goSPL (the `CaliTarheel/gospl` fork, 2026.7.14) runs from a Python 3.12 virtual environment over
Ubuntu's PETSc, apart from the system Python:

```bash
apt-get install gfortran libopenmpi-dev openmpi-bin python3-petsc4py python3-mpi4py python3-numpy \
  python3-scipy python3-h5py python3-pandas python3-ruamel.yaml python3-vtk9 meson ninja-build \
  python3.12-dev python3.12-venv
python3.12 -m venv --system-site-packages /opt/gospl-env
/opt/gospl-env/bin/python -m pip install "meson-python>=0.15" build numpy-indexed \
  "numpy<2" "astropy<6.1" "pyshtools<4.13"
/opt/gospl-env/bin/python -m pip install --no-build-isolation --no-deps path/to/gospl
```

numpy stays below 2: Ubuntu's petsc4py is built against 1.26. Ubuntu's PETSc has no ParMETIS, so
a parallel run names PETSc's own partitioner:

```bash
PETSC_OPTIONS="-petscpartitioner_type simple" mpirun -x PETSC_OPTIONS -np 4 /opt/gospl-env/bin/gospl -i input.yml
```

## Artifexia, 2000 to 1150 Ma

Built at mesh level 6 (40,962 nodes, about 120 km apart), with stages 10 Myr apart:

- **Build.** 86 maps and 85 intervals in about 8 minutes.
- **Plate speeds.** Median 2.5 cm/yr; 99% under 13 cm/yr.
- **Vertical rates.** 99% between -0.15 and +0.6 mm/yr.
- **goSPL's fast mode.** It runs the whole 850 Myr in 20 s on one process.
  - Its surface matches the emulated one (`expected`) to 0.24 mm at every stage.
  - It stays within 150-520 m RMS of the maps throughout, with no drift.
  - 1.3-2% of nodes fall on the other side of sea level. This is coastlines and ranges carried
    across a 120 km mesh; a finer mesh narrows it.

## Keeping the landscape on the reconstruction

goSPL run in one go over the build's forcing erodes the land and nothing gives it back. In the
first 10 Myr of Artifexia:

| | goSPL in one go | Map |
| --- | --- | --- |
| Land | 21% of the sphere | 26% |
| Highest point | 150 m | 1,490 m |

`run.py` therefore runs one interval at a time from where the last one ended. The vertical
forcing is the build's plus two corrections, each smoothed over 400 km with land and sea kept apart:

- what erosion and deposition have changed so far (goSPL's surface against the emulated
  surface without them);
- what they will change over the coming interval, going by the last interval's rates.

The regional height returns to the reconstruction's; the valleys, basins and shelves goSPL makes
stay. goSPL's `zfit` cannot do this: it resets only the sea floor.

goSPL's example erodibility, `K = 1e-6`, is for a 30 km mesh. On the 120 km mesh it wore down
2 km of land in 30 Myr, faster each interval. `K = 2e-7` settles.

First 50 Myr of Artifexia (2000 to 1950 Ma), mesh level 6, 4 processes, about 8 minutes an
interval:

| Stage | Land (map) | Highest point (map) | RMS from map | Median eroded on land |
| --- | --- | --- | --- | --- |
| 1990 Ma | 25.8% (26.3%) | 808 m (1,492 m) | 186 m | 104 m |
| 1980 Ma | 25.8% (26.3%) | 645 m (1,566 m) | 325 m | 167 m |
| 1970 Ma | 24.9% (26.3%) | 693 m (1,442 m) | 271 m | 180 m |
| 1960 Ma | 25.1% (26.3%) | 705 m (1,715 m) | 254 m | 200 m |
| 1950 Ma | 25.4% (26.4%) | 776 m (1,766 m) | 279 m | 209 m |

- **What goSPL makes:** branching drainage across the supercontinent, sediment shelves along its
  coasts, and the rift seaway the reconstruction opens.
- **Peaks:** they stay well under the maps'. A range narrower than the 400 km smoothing is not
  given back what it loses, and at 120 km the mesh holds few nodes across a range.
- **Full span:** 850 Myr at this level is about 85 intervals, roughly 12 hours.

