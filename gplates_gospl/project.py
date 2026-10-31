"""A GPlates project on disk: its rotation model and its layers, sorted by what they mean for relief.

Layer files are named by the person who made the project, so which file plays which part is
given by a role map. `ARTIFEXIA` is the map for Artifexian's Artifexia; another project passes
its own.
"""
from __future__ import annotations

import glob
import os

import pygplates

ARTIFEXIA = {
    # Crust: ocean strips carry their age as their begin time (the stage they formed in).
    'ocean': ['ocean crust.gpml', 'extant ocean crust.gpml'],
    'continent': ['continents.gpml', 'microcontinents.gpml', 'accreted terranes.gpml'],
    'craton': ['cratons.gpml'],
    'orogeny_active': ['active orogenies.gpml'],
    'orogeny_former': ['former orogenies.gpml', 'former orogenies culled.gpml'],
    'orogeny_old': ['old orogenies.gpml'],
    'arc': ['island arcs.gpml'],
    'lip_active': ['active LIP.gpml'],
    'lip_former': ['former LIP.gpml'],
    'trench': ['subduction zones.gpml'],
    'rift': ['rifts.gpml'],
}


class Project:
    def __init__(self, directory: str, roles: dict = ARTIFEXIA, rotation: str = 'rotation.rot'):
        self.directory = directory
        self.rotation_model = pygplates.RotationModel(os.path.join(directory, rotation))
        self.layers = {}
        for role, names in roles.items():
            self.layers[role] = [pygplates.FeatureCollection(os.path.join(directory, n)) for n in names
                                 if os.path.exists(os.path.join(directory, n))]
        # The plates are closed topologies whose sections are features in any layer (ridges,
        # subduction zones, plain lines), so every layer is loaded to resolve them.
        self.layers['topology'] = [pygplates.FeatureCollection(f)
                                   for f in sorted(glob.glob(os.path.join(directory, '*.gpml')))]

    def reconstructed(self, role: str, age: float) -> list:
        """(feature, reconstructed geometry) for every feature of `role` that exists at `age`."""
        out = []
        for collection in self.layers.get(role, []):
            found = []
            pygplates.reconstruct(collection, self.rotation_model, found, age)
            out += [(r.get_feature(), r.get_reconstructed_geometry()) for r in found]
        return out

    def topologies(self, age: float) -> list:
        resolved = []
        pygplates.resolve_topologies(self.layers['topology'], self.rotation_model, resolved, age)
        return resolved
