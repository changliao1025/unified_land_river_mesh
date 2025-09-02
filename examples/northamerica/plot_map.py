#a special function to plot a customized map
import sys, os
import cartopy.crs as ccrs
from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent

sPath = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/'
#add the path to the python path
sys.path.append(sPath)

from codes.map.map_multiple_vector_files import map_multiple_vector_files
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/northamerica/pyhexwatershed20250703011'
#feature to be plotter
sFilename_mesh = os.path.join(sWorkspace_output, 'pyflowline', 'mpas.geojson')


#list of basin
aFilename_river_network_real=list()
aFilename_river_conceptual=list()

aFilename_flow_direction = list()
aFilename_watershed_boundary = list()

sFilename_flow_direction = os.path.join(sWorkspace_output, 'hexwatershed', 'mpas_flow_direction.geojson')
nBasin = 10

for iBasin in range(nBasin):
    sBasin = str(iBasin + 1).zfill(8)
    aFilename_watershed_boundary.append( os.path.join(sWorkspace_output, 'pyflowline', sBasin, 'watershed_boundary.geojson') )
    aFilename_river_network_real.append( os.path.join(sWorkspace_output, 'pyflowline' , sBasin , 'flowline_filter.geojson') )
    aFilename_river_conceptual.append( os.path.join(sWorkspace_output, 'pyflowline' , sBasin , 'flowline_conceptual.geojson') )
    aFilename_flow_direction.append( os.path.join(sWorkspace_output, 'hexwatershed' , sBasin , 'flowline_direction.geojson') )



aFiletype_in =[]
aFilename_in=list()
aFlag_thickness_in = []
aFlag_color_in = []
aVariable_in = list()
aColor_in = []

#aFiletype_in.append(3)
#aFilename_in.append(sFilename_mesh)
#aFlag_thickness_in.append(0)
#aFlag_color_in.append(1)  #use color for watershed boundary
#aVariable_in.append('')
#aColor_in.append('black')

#add watershed boundary first
for sFilename in aFilename_watershed_boundary:
    aFiletype_in.append(3)
    aFilename_in.append(sFilename)
    aFlag_thickness_in.append(0)
    aFlag_color_in.append(1)  #use color for watershed boundary
    aVariable_in.append('')
    aColor_in.append('red')

#add flow direction
aFiletype_in.append(2)
aFilename_in.append(sFilename_flow_direction)
aFlag_thickness_in.append(1)
aFlag_color_in.append(1)
aVariable_in.append('drainage_area')
aColor_in.append('blue')
#for sFilename in aFilename_river_network_real:
#    aFiletype_in.append(2)
#    aFilename_in.append(sFilename)
#    aFlag_thickness_in.append(0)
#    aFlag_color_in.append(0)
#    aVariable_in.append('')
#    aColor_in.append('black')

#for sFilename in aFilename_river_conceptual:
#    aFiletype_in.append(2)
#    aFilename_in.append(sFilename)
#    aFlag_thickness_in.append(0)
#    aFlag_color_in.append(1)
#    aVariable_in.append('')
#    aColor_in.append('red')

sFilename_output = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/figures/flow_direction_saint_2.png'

sFilename_mesh_boudary = '/qfs/people/liao313/data/hexwatershed/conus/vector/conus_boundary.geojson'
aExtent_conus = gdal_get_vector_extent(sFilename_mesh_boudary)

#use saint louis extent
aExtent_conus = [-93.0, -87.0, 35.0, 40.0]  # [dLon_min, dLon_max, dLat_min, dLat_max]

dLon_min = aExtent_conus[0]
dLon_max = aExtent_conus[1]
dLat_min = aExtent_conus[2]
dLat_max = aExtent_conus[3]


pProjection_map = ccrs.Orthographic(central_longitude =  0.50*(dLon_max+dLon_min),
                                            central_latitude = 0.50*(dLat_max+dLat_min),
                                            globe=None)

map_multiple_vector_files(aFiletype_in,
                             aFilename_in,
                             iFlag_title_in = None,
                             iFlag_zebra_in=1,
                             iFlag_filter_in = None,
                             iFlag_arrow_in = None,
                             aFlag_thickness_in = aFlag_thickness_in,
                             aFlag_color_in = aFlag_color_in,
                             aFlag_colorbar_in = None,
                             aFlag_discrete_in = None,
                             aFlag_fill_in = None,
                             aVariable_in = aVariable_in,
                             sFilename_output_in=sFilename_output,
                             iFlag_scientific_notation_colorbar_in = None,
                             iFlag_openstreetmap_in = None,
                             iFlag_terrain_image_in = None,
                             iFlag_esri_hydro_image_in = 1,
                             iFont_size_in=None,
                             iBasemap_zoom_level_in = None,
                             sColormap_in = None,
                             sTitle_in = None,
                             iDPI_in = None,
                             iSize_x_in = None,
                             iSize_y_in = None,
                             aThickness_in = None,
                             aColor_in = aColor_in,
                             aMissing_value_in = None,
                             aData_max_in = None,
                             aData_min_in = None,
                             sExtend_in =None,
                             sFont_in = None,
                             sUnit_in=None,
                             aLegend_in = None,
                             aExtent_in = aExtent_conus,
                             pProjection_map_in=pProjection_map,
                             pProjection_data_in=None

)

