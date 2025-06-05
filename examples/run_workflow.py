

#refer to the e3sm confluence page for more details of the workflow

#each step may require different input datasets
#however, all the output will be saved in the same output directory for easy access

#you can change this to your preferred output directory
import os
sWorkspace_output = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/conus/output'

#Step 0 set up the resolution

dResolution_land = 10.0  #unit in km
dResolution_coastal = 5.0  #unit in km
sResolution_coastal = "{:d}".format(dResolution_coastal) + 'km'
dThreshold_area_island = 1.0E5 #unit km2
dDistance_tolerance_in=dResolution_land
dDrainage_area_threshold_in= dResolution_land * dResolution_land *5  #at least five grid cells of drainage area

#Step 1
#prepare the river network and coastal line dataset
from pyflowline.algorithms.simplification.simplify_hydrosheds import simplify_hydrosheds_river_network
sFilename_flowline_hydrosheds_in=''
sFilename_flowline_hydrosheds_out=''

simplify_hydrosheds_river_network(sFilename_flowline_hydrosheds_in,
                       sFilename_flowline_hydrosheds_out,
                       dDistance_tolerance_in,
                        dDrainage_area_threshold_in)

sFilename_coastal_in = ''
sFilename_geojson_out  = os.path.join(sWorkspace_output, 'coastal_line_wo_island.geojson')
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
remove_small_polygon(sFilename_coastal_in, sFilename_geojson_out, dThreshold_area_island )
from pyearth.toolbox.geometry.create_gcs_buffer_zone import create_gcs_buffer_zone_polygon
sFilename_geojson_buffer = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_buffer_' + sResolution_coastal + '.geojson')
create_gcs_buffer_zone_polygon(sFilename_geojson_out, sFilename_geojson_buffer,
                                      dBuffer_distance_in = dResolution_coastal * 1000)

#Step 2 - 4
#run the hexwatershed model, this step include three steps merged together.
#for debug purpose, you can also run then one by one, using the iFlag_debug flag to control

iFlag_debug = 1  #0 for normal run, 1 for debug
if iFlag_debug == 1:

    pass
else:
    pass

print('Congratulations! The workflow has been completed successfully!')
