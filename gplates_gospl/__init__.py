"""From a GPlates reconstruction to a goSPL landscape-evolution run.

A hand-built GPlates project (rotation file and feature layers, as in Artifexian's Artifexia)
supplies the plate motions and the features; this package turns them, stage by stage, into what
goSPL's global plate-motion workflow needs:

- a spherical mesh (`mesh`);
- an elevation map at every stage boundary, drawn from the reconstructed layers (`elevation`);
- each interval's plate displacement, from the rotation model and the plates (`motion`);
- each interval's vertical forcing: the rate that turns one stage's map, carried by the plates,
  into the next one, as goSPL's own examples work it out (`forcing`);
- goSPL's input files (`build`).
"""
