from pyearth.system.define_global_variables import *
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_file_area


sFilename_mesh_boudary = '/qfs/people/liao313/data/hexwatershed/conus/vector/conus_boundary.geojson'
dArea = calculate_polygon_file_area(sFilename_mesh_boudary)


dAarea_cell = 3000 * 3000

nCells = dArea / dAarea_cell

print(f"Area of the mesh boundary: {dArea} square meters")
print(f"Number of cells in the mesh: {nCells}")
print(f"Area of each cell: {dAarea_cell} square meters")