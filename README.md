# Paradise

A guided tectonic history for worldbuilding. Paint a starting world, then step it forward 10 Myr
at a time, choosing how each plate moves. Paradise keeps the plates, crust and boundaries
consistent as they move, and exports the history to GPlates.

Paradise is not a physics solver. You choose the plate motions; Paradise does the bookkeeping
that a hand-built GPlates reconstruction otherwise needs:
- new ocean at ridges, and crust consumed at trenches;
- boundaries typed by their motion;
- trenches that form as plates converge, and continents that collide and weld.

See [the design](docs/DESIGN.md) for how it works and what comes next.

## Run it

Use Python 3.12 or later. From the repository:

```sh
python -m venv .venv
```

```powershell
# Windows PowerShell
.\.venv\Scripts\Activate.ps1
```

```sh
# Linux or macOS
source .venv/bin/activate
```

```sh
python -m pip install -r requirements.txt
```

```sh
python -m paradise.server --project my-world.paradise.json --open
```

The app opens at [127.0.0.1:8790](http://127.0.0.1:8790). The project file is created when you
build a world, and saved after every step. A `.cache` file beside it saves replaying the
history when you open the project again; you can delete it at any time.

## Paint the start

1. **Plates.** Under **Start from**, make random plates (a count and a seed), or paint your own.
   Choose a plate swatch, or **New plate**, and paint with the brush. **Fill** floods a region,
   which is quick for redrawing a plate whose outline you have painted.
2. **Crust.** Switch the layer to **Crust** and paint continents and cratons onto the plates.
   Everything not painted is ocean.
3. **World.** Set the start age (the history counts down from it, as GPlates does) and the starting
   ocean age, then **Build world**.

## Step the history

Each plate has a stage pole: latitude, longitude and a rate in degrees per Myr, positive
anticlockwise seen from above the pole. This is the GPlates convention.
- A plate with a rate of 0 stays still.
- Click a plate on the map to select its row.
- Shift-click the map to put the selected plate's pole there.
- **Max cm/yr** is the fastest speed anywhere on the plate. It turns red above 20, about the speed
  of Earth's fastest plates.

**Preview** shows where the plates would be after the next 10 Myr, as dashed outlines, without
keeping the step. **Commit step** keeps it, with your note. Under **History**, **Go back here**
drops every later step.

The notes under the buttons report what happened in the step, such as trenches that formed and
continents that collided.

## Export to GPlates

**Export to GPlates** writes a folder beside the project file containing:
- `rotation.rot`;
- `Plates`, `Continents`, `Cratons`, `Ocean crust`, `Divergent`, `Convergent` and `Transform`
  `.gpml` layers.

Open them all in GPlates. Each step's features show from that step's age until the next step's.
Edits made in GPlates do not come back into Paradise.

## Limitations

- **Trenches.** A converging boundary becomes a trench only after it has closed 100 km. You cannot
  yet declare a trench, its polarity, or a rift yourself.
- **Time resolution.** Sub-steps are a fixed 2.5 Myr, so a fast plate moves a few hundred
  kilometres per sub-step. A 10 Myr step takes about ten seconds for a dozen plates.
- **Export.**
  - Every plate's rotation is relative to the anchor plate 000; there is no plate hierarchy yet.
  - Plates are exported as polygons for each step, not as GPlates topologies.
- **Painting resolution.** The painting is 0.5 degrees per cell, and the cells' staircase is
  smoothed away. Detailed coastlines are left for later map rendering.

## Tests

```sh
python -m pytest -q tests
```

`tests/test_core_*.py` are Poly-Time's tests for the copied core; `tests/test_paradise.py` covers
painting, the step history and the GPlates export.
