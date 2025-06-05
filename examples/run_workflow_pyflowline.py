

#refer to the e3sm confluence page for more details of the workflow

#each step may require different input datasets
#however, all the output will be saved in the same output directory for easy access

#you can change this to your preferred output directory
import os
from datetime import datetime
sWorkspace_input = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/conus/input'
sWorkspace_output = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/conus/output'

sDate_today = datetime.now().strftime('%Y%m%d')
#Step 0 set up the resolution
iCase_index = 1
sMesh_type = 'mpas'  #
dResolution_land = 10.0  #unit in km
dResolution_coastal = 5.0  #unit in km
sResolution_coastal = "{:d}".format(dResolution_coastal) + 'km'
dThreshold_area_island = 1.0E5 #unit km2
dDistance_tolerance_in=dResolution_land
dDrainage_area_threshold_in= dResolution_land * dResolution_land *5  #at least five grid cells of drainage area

#Step 1
#prepare the river network and coastal line dataset
from pyearth.toolbox.management.vector.fields import get_field_and_value, add_field_to_vector_file
from pyflowline.algorithms.simplification.simplify_hydrosheds import simplify_hydrosheds_river_network
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.geometry.create_gcs_buffer_zone import create_gcs_buffer_zone_polygon
sFilename_flowline_hydrosheds_in=''
sFilename_flowline_hydrosheds_out=''
sFilename_geojson_geometery_feature =''

#step 1: record attribute from the MPAS tools
aField, aValue = get_field_and_value(sFilename_geojson_geometery_feature)

simplify_hydrosheds_river_network(sFilename_flowline_hydrosheds_in,
                       sFilename_flowline_hydrosheds_out,
                       dDistance_tolerance_in,
                        dDrainage_area_threshold_in)

sFilename_coastal_in = ''
sFilename_geojson_out  = os.path.join(sWorkspace_output, 'coastal_line_wo_island.geojson')

remove_small_polygon(sFilename_coastal_in, sFilename_geojson_out, dThreshold_area_island )

sFilename_geojson_buffer = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_buffer_' + sResolution_coastal + '.geojson')
create_gcs_buffer_zone_polygon(sFilename_geojson_out, sFilename_geojson_buffer,
                                      dBuffer_distance_in = dResolution_coastal * 1000)

sFilename_vector_geojson_merged = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_merged_' + sResolution_coastal + '.geojson')
merge_features(sFilename_geojson_buffer, sFilename_vector_geojson_merged)
add_field_to_vector_file(sFilename_vector_geojson_merged, aField, aValue)

#Step 2 - 4
#run the hexwatershed model, this step include three steps merged together.
#for debug purpose, you can also run then one by one, using the iFlag_debug flag to control

iFlag_debug = 1  #0 for normal run, 1 for debug
from pyflowline.configuration.read_configuration_file import pyflowline_read_configuration_file
from pyflowline.configuration.config_manager import create_template_configuration_file, create_template_jigsaw_configuration_file
from pyflowline.configuration.change_json_key_value import change_json_key_value
if iFlag_debug == 1:
    #create the mesh using pyflowline
    #create a template configuration file
    sFilename_configuration_json = os.path.join(sWorkspace_output, 'pyflowline_configuration.json')
    create_template_configuration_file(sFilename_configuration_json,
        sWorkspace_input = sWorkspace_input,
        sWorkspace_output = sWorkspace_output,
        iFlag_standalone_in=1,
        sMesh_type_in=sMesh_type,
        sModel_in='pyflowline')

    #if the stream burning is turned on, a basin configuration file is also automatically created under the same output directory
    sFilename_basins_json = os.path.join(sWorkspace_output, 'pyflowline_configuration_basins.json')

    #now we will update the configuration file with our own settings
    change_json_key_value(sFilename_configuration_json, "iFlag_flowline", 0) #disable the flowline simplification
    change_json_key_value(sFilename_configuration_json, "iFlag_run_jigsaw", 1) #turn on the jigsaw mesh generation
    change_json_key_value(sFilename_configuration_json, "iFlag_intersect", 0) #turn off the intersection

    change_json_key_value(sFilename_configuration_json, "sFilename_mesh_boundary", 0) #set the target mesh boundary
    change_json_key_value(sFilename_configuration_json, "sFilename_coastal_boundary", 0) #set the coastline boundary

    #now we also need to set the jigsaw configuration file
    sFilename_jigsaw_configuration_json = os.path.join(sWorkspace_output, 'pyflowline_configuration_jigsaw.json')
    create_template_jigsaw_configuration_file(sFilename_jigsaw_configuration_json)
    #update the jigsaw configuration file
    change_json_key_value(sFilename_jigsaw_configuration_json, "iFlag_spac", "true") #enable resolution control
    change_json_key_value(sFilename_jigsaw_configuration_json, "iFlag_geom", "true") # enable geometry control

    #now we can run the pyflowline to create the mesh
    oPyflowline = pyflowline_read_configuration_file(sFilename_configuration_json,
                    iCase_index_in=iCase_index,
                    sDate_in= sDate_today,
                    sMesh_type_in = sMesh_type)

    oPyflowline._pyflowline_create_hpc_job(sSlurm_in = 'slurm', hours_in = 30 )
    #now you should manually submit the job


else:
    pass

print('Congratulations! The workflow has been completed successfully!')
