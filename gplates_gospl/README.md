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

Smoothed corrections give a range back its average loss, and crests lose least. Without a limit,
the highest point climbed about 450 m an interval past the map's, from 5,864 m to 6,708 m against
5,600 m. So no node is lifted past the map's height where its plate takes it; erosion may still
take it below.

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

## Rainfall from ROCKE-3D

`rocke3d.py` takes a world to ROCKE-3D and its climate back, in the formats of Aerie's
`map_definitions/rocke3d_4x5`:

- **Out.** `grid_npz` writes a stage's elevation as Aerie's 72 x 46 grid NPZ: land fraction,
  mean, land and ocean elevations, and cell areas. Each cell is sampled on its own fine grid,
  since 5 degrees of longitude at 86 N is under 40 km.
- **Back.** `read_precipitation` reads the extractor's climate CSV (`pNN_mm_day` for twelve
  months), weights each month by its length in days, and gives annual rainfall in m/yr.
  `rain_on_mesh` interpolates that onto goSPL's mesh.
- **Into a run.** `run.py --rain CLIMATE.csv` uses it in every interval.
- **Stand-in.** `random_climate_csv` writes a made-up climate in the same CSV format, to prove the
  route before a real ROCKE-3D run exists.

ROCKE-3D's cells are about 450 km across, so the rain carries no rain shadows at goSPL's scale.
goSPL's own orographic precipitation could add them.

## Artifexia's last 50 Myr at 60 km

Built at level 7 (163,842 nodes, about 60 km apart), 1200 to 1150 Ma, with a random climate
read through the ROCKE-3D CSV route. Land rainfall is 0.05-1.8 m/yr, 0.45 on average. The run used
`K = 2e-7`, steps of 250,000 years and 4 processes: about 10-14 s a step, 8 minutes an interval.

| Stage | Land (map) | Highest point (map) | RMS from map | Median eroded on land |
| --- | --- | --- | --- | --- |
| 1190 Ma | 30.3% (30.3%) | 5,328 m (5,799 m) | 514 m | 96 m |
| 1180 Ma | 30.4% (30.5%) | 5,876 m (5,993 m) | 445 m | 98 m |
| 1170 Ma | 30.5% (30.5%) | 5,194 m (5,805 m) | 397 m | 89 m |
| 1160 Ma | 30.6% (30.5%) | 5,733 m (6,166 m) | 419 m | 88 m |
| 1150 Ma | 30.7% (30.5%) | 5,159 m (5,808 m) | 472 m | 89 m |

- **Land and peaks.** Land stays within a fifth of a point of the reconstruction. Peaks stay under
  the map's height without climbing past it.
- **Erosion.** It holds near 90 m per 10 Myr.
- **Rivers.** They run off the ranges, with catchments following the relief and the rain bands.
- **Sediment.** Shelves build along the coasts and inland basins fill.
- **Limit.** Where a coastal range sheds into the sea, deposits reach several km over the run.
  The loop restores the land each interval, so the sediment supply never runs out.

## On a bigger machine

goSPL is built on PETSc and MPI, and the mesh is the limit, not goSPL. Measured here (4 cores,
15 GB):

| Level | Spacing | Nodes | Per goSPL step, 4 processes | Memory |
| --- | --- | --- | --- | --- |
| 6 | 120 km | 40,962 | about 5 s | small |
| 7 | 60 km | 163,842 | about 10-14 s | under 4 GB |
| 8 | 30 km | 655,362 | about 90 s | about 15 GB: past this container's limit after two steps |

Memory runs to about 22 KB a node.

| Level | Spacing | Nodes | Memory needed |
| --- | --- | --- | --- |
| 9 | 15 km | 2.6 million | about 60 GB: a 24-core, 128 GB machine takes it |
| 10 | 7.5 km | 10.5 million | about 230 GB: past that machine |

Three things make a fine mesh run at all:

- **`oFill`, in the `diffusion` section, set to -500 m.** goSPL fills depressions above this depth.
  Its default, -6,000 m, takes in nearly the whole ocean floor, whose age bands and trenches make
  countless pits. At 30 km its spill graph asked for 33 GB and the run died.
- **Land that drains.** Each map's land is priority-flooded from the sea, continents rise
  gently inland, and ranges rise towards their middles. 60% of land edges had been dead flat.
- **Mesh numbering.** The nodes are numbered along a Morton curve, so that PETSc's simple
  partitioner gives each process one patch of the sphere. Numbered as subdivision made them, each
  process held most of the mesh as halo, and at 60 km a step took twice as long.

1. Install pygplates in the Python that builds, and goSPL with its MPI in the one that runs (see
   *Running goSPL in this container*).
2. Build the window:

   ```bash
   python -m gplates_gospl.build PROJECT OUT --from 1200 --to 1150 --every 10 --level 9
   ```

3. Run it:

   ```bash
   python -m gplates_gospl.run OUT RUN --np 22 --K 2e-7 --dt 2.5e5 --rain CLIMATE.csv
   ```

   `run.py` falls back to PETSc's simple partitioner, which Ubuntu's PETSc needs (it has no
   ParMETIS). On a PETSc with ParMETIS, set `PETSC_OPTIONS=-petscpartitioner_type parmetis` (or
   any other options) before running, and `run.py` keeps them.

`K` needs tuning for each resolution. Watch the median erosion that `run.py` logs: it should level
off, not climb interval on interval.

