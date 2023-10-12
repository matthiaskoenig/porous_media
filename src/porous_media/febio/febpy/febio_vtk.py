"""
Enums for the reading for XPLT.
This file has been adapted from Jacob Sturdy https://gitlab.com/utsekaj/febpy.
"""
import numpy as np
import pyvista as pv

from .xplt import domainClass as Region
from .xplt_enums import Elem_Type, nodesPerElementClass


vtk = pv._vtk

NN_PER_ELEMENT = {
    "TRI": 3,
    "TRI3": 3,
    "TRI6": 6,
    "TRI10": 10,
    "QUAD4": 4,
    "QUAD8": 8,
    "QUAD9": 9,
    "HEX8": 8,
    "HEX20": 20,
    "HEX27": 27,
    "PENTA6": 6,
    "PENTA15": 15,
    "PENTA18": 18,
    "TET4": 4,
    "TET10": 10,
}

ELEMENT_DIMENSION = {
    "TRI": 2,
    "TRI3": 2,
    "TRI6": 2,
    "TRI10": 2,
    "QUAD4": 2,
    "QUAD8": 2,
    "QUAD9": 2,
    "HEX8": 3,
    "HEX20": 3,
    "HEX27": 3,
    "PENTA6": 3,
    "PENTA15": 3,
    "PENTA18": 3,
    "TET4": 3,
    "TET10": 3,
}

VTK_TYPES = {  # How to merge with Enums.py?
    "TRI": vtk.VTK_TRIANGLE,
    "TRI3": vtk.VTK_TRIANGLE,
    "TRI6": vtk.VTK_QUADRATIC_TRIANGLE,
    "TRI10": vtk.VTK_BIQUADRATIC_TRIANGLE,
    "QUAD4": vtk.VTK_QUAD,
    "QUAD8": vtk.VTK_QUADRATIC_QUAD,
    "QUAD9": vtk.VTK_BIQUADRATIC_QUAD,
    "HEX8": vtk.VTK_HEXAHEDRON,
    "HEX": vtk.VTK_HEXAHEDRON,  # TODO why is this needed?
    "HEX20": vtk.VTK_QUADRATIC_HEXAHEDRON,
    "HEX27": vtk.VTK_TRIQUADRATIC_HEXAHEDRON,
    "PENTA": vtk.VTK_WEDGE,
    "PENTA6": vtk.VTK_WEDGE,
    "PENTA15": vtk.VTK_QUADRATIC_WEDGE,
    "PENTA18": 32,
    "TET4": vtk.VTK_TETRA,
    "TET10": vtk.VTK_QUADRATIC_TETRA,
}


def reverse_penta6(element):
    """swap connectivity from FEBIO sense to VTK sense this is an involution"""
    return element[3], element[4], element[5], element[0], element[1], element[2]


# VTK_TYPES = {
#     Elem_Type.ELEM_TRI: vtk.VTK_TRIANGLE,
#     # 'TRI3': vtk.VTK_TRIANGLE,
#     Elem_Type.ELEM_TRI6: vtk.VTK_QUADRATIC_TRIANGLE,
#     Elem_Type.ELEM_TRI10: vtk.VTK_BIQUADRATIC_TRIANGLE,
#     Elem_Type.ELEM_QUAD: vtk.VTK_QUAD,
#     Elem_Type.ELEM_QUAD8: vtk.VTK_QUADRATIC_QUAD,
#     Elem_Type.ELEM_QUAD9: vtk.VTK_BIQUADRATIC_QUAD,
#     Elem_Type.ELEM_HEX: vtk.VTK_HEXAHEDRON,
#     Elem_Type.ELEM_HEX20: vtk.VTK_QUADRATIC_HEXAHEDRON,
#     Elem_Type.ELEM_HEX27: vtk.VTK_TRIQUADRATIC_HEXAHEDRON,
#     Elem_Type.ELEM_PENTA: vtk.VTK_WEDGE,
#     Elem_Type.ELEM_PENTA15: vtk.VTK_QUADRATIC_WEDGE,
#     'PENTA18': 32,
#     Elem_Type.ELEM_TET4: vtk.VTK_TETRA,
#     Elem_Type.ELEM_TET10: vtk.VTK_QUADRATIC_TETRA,
# }


class RegionToVTKConverter:
    # # TODO: Need to look into numpy and nptyping
    # point: np.ndarray
    # cell_data: np.ndarray[int]
    # types: np.array[bool]
    # regionID: np.array[int]
    # partID: np.array[int]
    # point_ids: np.array[int]
    # points_mask: np.array(bool)
    def __init__(self, points: np.ndarray, region: Region):
        self.points = None
        self.cell_data = None
        self.point_ids = None
        self.types = None
        self.regionID = None
        self.partID = None

        if region is None:
            return

        self.points = (
            points  # TODO should this be a copy to prevent accidental modification?
        )
        etype = region.elemType
        # enodes = NN_PER_ELEMENT[etype]
        enodes = nodesPerElementClass[etype]
        vtk_type = VTK_TYPES[etype.split("ELEM_")[-1]]
        elements = region.elements
        if vtk_type == VTK_TYPES["PENTA6"]:
            cells = np.fromiter(region.elements.values(), dtype=np.dtype((int, enodes)))
            cells = np.apply_along_axis(reverse_penta6, 1, cells)
        else:
            cells = np.fromiter(region.elements.values(), dtype=np.dtype((int, enodes)))
        self.point_ids = np.unique(cells)
        self.cell_data = np.empty((len(elements), enodes + 1), dtype=int)
        self.cell_data[:, 0] = enodes
        self.cell_data[:, 1:] = cells
        self.cell_data = self.cell_data.reshape(self.cell_data.size)
        self.types = np.empty(len(elements), dtype=int)
        self.types[:] = vtk_type

        self.regionID = np.empty(len(elements), dtype=int)
        self.regionID[:] = region.domainID

        self.partID = np.empty(len(elements), dtype=int)
        self.partID[:] = region.partID

    def generate_unstructured_grid(self):
        grid = pv.UnstructuredGrid(self.cell_data, self.types, self.points)
        grid.cell_data.set_scalars(self.regionID, "Region")
        grid.cell_data.set_scalars(self.partID, "Part")
        return grid

    def __add__(self, other):
        assert np.allclose(self.points, other.points)
        new = RegionToVTKConverter(None, None)
        new.points = self.points.copy()
        new.cell_data = np.concatenate((self.cell_data, other.cell_data))
        new.types = np.concatenate((self.types, other.types))
        new.regionID = np.concatenate((self.regionID, other.regionID))
        new.partID = np.concatenate((self.partID, other.partID))
        new.point_ids = np.unique(np.concatenate((self.point_ids, other.point_ids)))
        return new

    def set_data_from_xplt(self, xplt_reader, step):
        # if xplt_reader.readMode != 'readAllStates':
        #     # 'Assuming  no data has been read in yet' # TODO: fix xplt_reader for caching/rereading steps
        #     xplt_reader.readSteps([step, ])
        regionID = self.regionID[0]
        grid = self.generate_unstructured_grid()
        data_arrays = dict()
        # TODO implement direct from domain??

        region_key = xplt_reader.mesh.regionName(regionID)
        for variable, dataObj in xplt_reader.results.items():
            if "global" in dataObj.data_dict:
                data = dataObj.data_dict["global"]
                n_items = grid.n_points
                region_start = 0
                region_end = n_items
            elif region_key in dataObj.data_dict:
                data = dataObj.data_dict[region_key]
                n_items = grid.n_cells
                region_start = 0
                region_end = n_items
            else:
                continue

            if len(data) == 0:
                continue

            data = data[step]
            if not variable in data_arrays:
                n_dims = data[0].shape[-1]
                var_array = np.nan * np.zeros((n_items, n_dims))
                data_arrays[variable] = var_array

            data_arrays[variable][region_start:region_end] = data[:]

        for variable, values in data_arrays.items():
            if values.shape[0] == grid.n_cells:
                grid_data = grid.cell_data
            else:
                grid_data = grid.point_data
            grid_data.set_array(values, variable)
        return grid

    def removeUnusedPoints(self):
        # https://vtk.org/doc/nightly/html/classvtkRemoveUnusedPoints.html
        pass


class DomainToVTKConverter:
    def __init__(self, domain_name, xplt_reader):
        domain_index = xplt_reader.mesh.parts[domain_name]
        xplt_regions = xplt_reader.mesh.parts[domain_name].domains
        domain = None
        regions = []
        for region in xplt_regions:
            part = RegionToVTKConverter(xplt_reader.mesh.nodes.copy(), region)
            regions.append(part)
            if domain is None:
                domain = part
            else:
                domain += part
        self.points = domain.points
        self.domain = domain
        self.regions = regions

    def generate_unstructured_grid(self):
        grid = None
        for part in self.regions:
            region_grid = part.generate_unstructured_grid()
            if grid is None:
                grid = region_grid
            else:
                grid += region_grid
        return grid

    def set_data_from_xplt(self, xplt_reader, step):
        grid = None
        for part in self.regions:
            region_grid = part.set_data_from_xplt(xplt_reader, step)
            if grid is None:
                grid = region_grid
            else:
                grid += region_grid
        return grid
