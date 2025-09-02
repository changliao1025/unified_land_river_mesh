
import cartopy.crs as ccrs
from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.visual.map.vector.map_vector_polyline_file import map_vector_polyline_file
from codes.map.map_multiple_vector_files import map_multiple_vector_files

sWorkspace_simulation = '/compyfs/liao313/00raw/mosart'
sFilename_parameter_in = sWorkspace_simulation + '/mosart_parameter_16th.nc'
sFilename_geojson_out = sWorkspace_simulation + '/DRT_2_FDR_globe.geojson'
sFilename_png = sWorkspace_simulation + '/DRT_2_FDR_globe.png'

sFilename_mesh_boudary = '/qfs/people/liao313/data/hexwatershed/conus/vector/conus_boundary.geojson'
aExtent_conus = gdal_get_vector_extent(sFilename_mesh_boudary)



dLon_min = aExtent_conus[0]
dLon_max = aExtent_conus[1]
dLat_min = aExtent_conus[2]
dLat_max = aExtent_conus[3]


pProjection_map = ccrs.Orthographic(central_longitude =  0.50*(dLon_max+dLon_min),
                                            central_latitude = 0.50*(dLat_max+dLat_min),
                                            globe=None)
#map_vector_polyline_file(1,
#                             sFilename_geojson_out,
#                             sFilename_output_in = sFilename_png,
#                             iFlag_thickness_in =1,
#                             iFlag_esri_hydro_image_in=1,
#                             sField_thickness_in='Drainage',
#                             pProjection_map_in = pProjection_map,
#                             aExtent_in=aExtent_conus
#                          )
aFiletype_in = []
aFilename_in = []
aFlag_thickness_in = []
aFlag_color_in = []

#add watershed boundary first
aFilename_watershed_boundary = list()
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/northamerica/pyhexwatershed20250703011'
sFilename_flow_direction = os.path.join(sWorkspace_output, 'hexwatershed', 'mpas_flow_direction.geojson')
nBasin = 10

for iBasin in range(nBasin):
    sBasin = str(iBasin + 1).zfill(8)
    aFilename_watershed_boundary.append( os.path.join(sWorkspace_output, 'pyflowline', sBasin, 'watershed_boundary.geojson') )
aVariable_in = []
aColor_in = []
for sFilename in aFilename_watershed_boundary:
    aFiletype_in.append(3)
    aFilename_in.append(sFilename)
    aFlag_thickness_in.append(0)
    aFlag_color_in.append(1)  #use color for watershed boundary
    aVariable_in.append('')
    aColor_in.append('red')


aFiletype_in.append(2)



aFilename_in.append(sFilename_geojson_out)


aVariable_in.append('Drainage')
sFilename_output = sWorkspace_simulation + '/DRT_2_FDR_globe_w_boundary.png'


aColor_in.append('blue')

aFlag_thickness_in.append(1)

aFlag_color_in.append(1)


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
print('finished')