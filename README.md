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
   ocean age, then **Build world**. Clear **Trenches form by themselves** if you want to declare
   every trench yourself. Otherwise a converging boundary also becomes a trench once it has closed
   100 km, with continent overriding ocean and younger ocean overriding older. **Largest terrane
   that accretes** sets the size above which continent reaching a trench collides instead (0: it
   always collides).

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

## Declare a trench

1. Under the map, set **Clicking the map** to **picks a boundary**.
2. Click a boundary. That picks the whole stretch of it between the same two plates and of one
   type. To declare only part of it, click it a second time: the pick shrinks to the stretch
   between your two clicks.
3. Choose which plate goes under. The declaration appears under **Declared for the next step**,
   drawn on the map with teeth on the overriding plate's side. You can remove it before
   committing.
4. **Preview** or **Commit step** as usual.

A declared trench lasts as the plates move, consuming the plate that goes under, until you end it
or a continent reaches it and it becomes a collision. To flip its polarity, pick it (on the map or
in the **Trenches** table) and declare the other plate under. **End subduction** turns it back into
an ordinary boundary, typed from the motion.

The **Trenches** table measures every trench after each step:
- its length;
- how fast it closes, negative where the plates pull apart;
- the age of the ocean going down;
- the area it consumes per step at that rate.

Warnings point out continent arriving at a trench, plates pulling away from their trench, and a
long trench whose downgoing plate barely moves. Scotese's first rule of thumb: plates carried by
a slab move fast.

## Weld two plates, or keep them apart

When continents collide, the two plates weld into one once the collision has closed 100 km. To
decide that yourself, pick the boundary between them:
- **Weld these plates now** joins them at the start of the step, whatever their boundary.
- **Keep these plates apart** lets them collide without welding, from this step on. The continent
  between them goes on shortening and thickening, and a warning says so each step.
- **Let them weld** undoes that.

## Rift a plate, or open a back-arc

**Rifts.**
1. Set **Clicking the map** to **draws a rift**.
2. Click points along the line, starting and ending beyond the edges of the plate it splits.
3. Check the plate, which Paradise fills in from the middle of the line, and choose which side of
   the line becomes the new plate. Left and right are seen walking the line in the order you
   clicked it.
4. Optionally give the new plate a pole for this step. Without one it moves with the plate it came
   from until you give it one.
5. **Declare the rift**, then preview or commit.

Moved apart, the new boundary opens as a ridge. Moved together, it shortens and thickens any
continent along it. The split takes the plate's crust, cratons, arcs, terranes and old sutures on
that side with it.

**Back-arcs.** Pick a trench, set a pole for the new plate if you want, and choose **Open a
back-arc behind it**. The overriding plate's strip from the trench to 250 km inboard, its forearc
and arc, becomes a new plate that overrides the trench. Rolled back towards the trench, it opens a
back-arc basin behind it.

## Island arcs and terranes

These follow from the trenches; there is nothing to declare.

- **Island arcs.** A trench that has lasted a whole step and closes at least 1 cm/yr raises a
  volcanic arc about 150 km inboard, on the overriding plate.
  - Over ocean the arc is new, thin crust: an island arc.
  - Over continent it is a continental arc, as in the Andes. It is recorded, but adds no crust.
  - An arc stays where it formed when its trench moves away, as a remnant arc.
  - Island arc crust starts 25 km thick and thickens 2 km every 10 Myr its trench still feeds
    it, up to 35 km.
- **Terranes.** A piece of continental crust of up to a million km² (an island arc, or a
  microcontinent you painted) that reaches a trench its plate goes down is not subducted. It
  accretes to the overriding plate, along with the ocean around it out to the trench. The trench
  steps outboard of it and carries on. Larger continents collide. The size limit is a project
  setting.

The map outlines island arcs, continental arcs and terranes. **Arcs and terranes** lists each one
with its story: where it formed, and when and to which plate it accreted.

## Export to GPlates

**Export to GPlates** writes a folder beside the project file containing:
- `rotation.rot`;
- `Plates`, `Continents`, `Cratons`, `Island arcs`, `Continental arcs`, `Accreted terranes`,
  `Rifts`, `Failed rifts`, `Ocean crust`, `Divergent`, `Convergent` and `Transform` `.gpml`
  layers. A rift whose two sides have welded back into one plate is a failed rift.

A plate made by a rift or a back-arc rotates relative to the plate it came from, as in a
hand-built rotation tree. Every other plate rotates relative to the anchor plate 000.

Each arc and terrane carries its story in its description, in ages.

Open them all in GPlates. Each step's features show from that step's age until the next step's.
Trenches are in `Convergent` as subduction zones, carried by the overriding plate, with
the polarity GPlates reads as the overriding plate's side. Edits made in GPlates do not come back into Paradise.

## Limitations

- **Trenches.**
  - A declaration applies to whole boundary segments within 60 km of the stretch you pick.
    Segments are at most about 110 km long at the start, and up to 150 km once the core has
    simplified moved outlines, so a declared trench can run up to about that much past either end
    of the pick.
  - A stretch shorter than three segments takes the type of its neighbours; Paradise warns when
    that happens.
  - Where the plates pull apart across a trench, the new ocean is split between the two plates,
    unless you open a back-arc there.
- **Arcs and terranes.**
  - Arcs grow only as a band 30 km either side of the volcanic front, and render as plain
    lowland, without volcanic relief.
  - Whether a piece of continent accretes or collides is decided by its size alone; you cannot yet
    choose for one piece.
- **Rifts.**
  - A rift splits a plate at once; there is no slow continental rifting stage before the ocean
    opens.
  - A failed rift is only one whose sides have welded back together. A rift that stops opening
    while its plates stay separate still counts as a rift.
- **Time resolution.** Sub-steps are a fixed 2.5 Myr, so a fast plate moves a few hundred
  kilometres per sub-step. A 10 Myr step takes about ten seconds for a dozen plates.
- **Export.**
  - Only plates made by rifts and back-arcs have a parent in the rotation file. The painted
    plates all hang from 000.
  - Plates are exported as polygons for each step, not as GPlates topologies.
- **Painting resolution.** The painting is 0.5 degrees per cell, and the cells' staircase is
  smoothed away. Detailed coastlines are left for later map rendering.

## Tests

```sh
python -m pytest -q tests
```

`tests/test_core_*.py` are Poly-Time's tests for the copied core; `tests/test_paradise.py` covers
painting, the step history and the GPlates export.
