

#refer to the e3sm confluence page for more details of the workflow

#each step may require different input datasets
#however, all the output will be saved in the same output directory for easy access

#you can change this to your preferred output directory
import os
from datetime import datetime
sWorkspace_input = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/global/input'

sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/global/river_network'
if os.path.exists(sWorkspace_output) is False:
    os.makedirs(sWorkspace_output)

sDate_today = datetime.now().strftime('%Y%m%d')
#Step 0 set up the resolution
iCase_index = 1
sMesh_type = 'mpas'  #
dResolution_land = 25.0  #unit in km
dResolution_coastal = 10.0  #unit in km
dResolution_x_in = 30.0/3600 * dResolution_coastal
dResolution_y_in = dResolution_x_in
nrow = int(180 / dResolution_y_in)
ncolumn = int(360 / dResolution_x_in)
sResolution_coastal = "{:.1f}".format(dResolution_coastal) + 'km'
dThreshold_area_island = 1.0E5 #unit km2
dDistance_tolerance_in = dResolution_land * 0.5 * 1.0E3 #how far away two river need to be for mesh generation
dDrainage_area_threshold_in= dResolution_land * dResolution_land *5 * 1.0E6 #at least five grid cells of drainage area
sDistance_tolerance = "{:.2E}".format(dDistance_tolerance_in)
sDrainage_area_threshold = "{:.2E}".format(dDrainage_area_threshold_in)

#setup flags for debugging
iFlag_simplify_hydrosheds_river_network = 1
iFlag_process_watershed_boundary = 1
iFlag_process_coastal_line = 0

#Step 1
#prepare the river network and coastal line dataset
from pyearth.toolbox.management.vector.fields import get_field_and_value, add_field_to_vector_file
from pyflowline.algorithms.simplification.simplify_hydrosheds import simplify_hydrosheds_river_network
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.geometry.create_gcs_buffer_zone import create_buffer_zone_polygon_file
from pyearth.toolbox.conversion.convert_vector_to_geojson import convert_vector_to_geojson
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster

sFilename_flowline_hydrosheds_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp'
sFilename_flowline_hydroshed_tmp = 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '.geojson'
sFilename_flowline_hydrosheds_out = os.path.join(sWorkspace_output, sFilename_flowline_hydroshed_tmp)
sFilename_geojson_geometery_feature = '/qfs/people/liao313/data/hexwatershed/global/vector/region.geojson'

#step 1: record attribute from the MPAS tools
aField, aValue = get_field_and_value(sFilename_geojson_geometery_feature)

#the river flowline simplficiation process already generated the basin configuration file
sFilename_pyflowline_configuration = os.path.join(sWorkspace_output, 'pyflowline_configuration.json')
sFilename_pyflowline_configuration_basins = os.path.join(sWorkspace_output, 'pyflowline_configuration_basins.json')
if iFlag_simplify_hydrosheds_river_network == 1:
    simplify_hydrosheds_river_network(sFilename_flowline_hydrosheds_in,
                       sFilename_flowline_hydrosheds_out,
                       dDistance_tolerance_in,
                        dDrainage_area_threshold_in,
                        iFlag_pyflowline_configuration_in=1,
                        nOutlet_largest=50)

exit()

sFilename_coastal_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrobasin/geojson/hybas_lake_all_lev03.geojson'


sFilename_geojson_basin_global = sFilename_coastal_in
sFilename_geojson_wo_island  = os.path.join(sWorkspace_output, 'coastal_line_wo_island.geojson')
sFilename_tif_wo_island = os.path.join(sWorkspace_output, 'coastal_line_wo_island.tif')
sFilename_geojson_buffer = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_buffer_' + sResolution_coastal + '.geojson')
sFilename_vector_geojson_merged = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_merged_' + sResolution_coastal + '.geojson')

if iFlag_process_coastal_line == 1:
    #convert_vector_to_geojson(sFilename_coastal_in, sFilename_geojson_basin_global)
    remove_small_polygon(sFilename_geojson_basin_global, sFilename_geojson_wo_island, dThreshold_area_island )
    #convert the vector to raster, this is a global raster, so we need to set the resolution
    convert_vector_to_global_raster(sFilename_geojson_wo_island,
                                    sFilename_tif_wo_island,
                                    dResolution_x_in, dResolution_y_in,
                                    iFlag_boundary_only_in = 1,
                                    dFill_value_in = 2)
    create_buffer_zone_polygon_file(sFilename_geojson_wo_island, sFilename_geojson_buffer,
                                          dBuffer_distance_in = dResolution_coastal ) #* 1000)
    merge_features(sFilename_geojson_buffer, sFilename_vector_geojson_merged)
    add_field_to_vector_file(sFilename_vector_geojson_merged, aField, aValue)

#Step 2 - 4
#run the hexwatershed model, this step include three steps merged together.
#for debug purpose, you can also run then one by one, using the iFlag_debug flag to control

iFlag_debug = 1  #0 for normal run, 1 for debug
from shutil import copy2
from pyflowline.configuration.read_configuration_file import pyflowline_read_configuration_file
from pyflowline.configuration.config_manager import create_template_configuration_file, create_template_jigsaw_configuration_file
from pyflowline.configuration.change_json_key_value import change_json_key_value
if iFlag_debug == 1:
    #create the mesh using pyflowline
    #create a template configuration file
    #sFilename_pyflowline_configuration = os.path.join(sWorkspace_output, 'pyflowline_configuration.json')
    #oPyflowline_template = create_template_configuration_file(sFilename_configuration_json)
    #sFilename_basins_json = sFilename_pyflowline_configuration_basins #oPyflowline_template.sFilename_basins

    #copy the template configuration file to the output directory
    change_json_key_value(sFilename_pyflowline_configuration, "sWorkspace_output", sWorkspace_output)
    oPyflowline = pyflowline_read_configuration_file(sFilename_pyflowline_configuration, \
    iCase_index_in=iCase_index, sDate_in=sDate_today, iFlag_create_directory_in=None)
    sWorkspace_output_case = oPyflowline.sWorkspace_output
    sFilename_configuration_copy = os.path.join(sWorkspace_output_case, 'pyflowline_configuration_copy.json')
    copy2(sFilename_pyflowline_configuration, sFilename_configuration_copy)
    #if the stream burning is turned on, a basin configuration file is also automatically created under the same output directory
    #copy the basin configuration file to the output directory as well
    sFilename_configuration_basins_copy = os.path.join( sWorkspace_output_case, 'pyflowline_configuration_basins_copy.json' )
    copy2(sFilename_pyflowline_configuration_basins, sFilename_configuration_basins_copy)


    #now we also need to set the jigsaw configuration file
    sFilename_jigsaw_configuration_json = os.path.join(sWorkspace_output, 'pyflowline_configuration_jigsaw.json')
    create_template_jigsaw_configuration_file(sFilename_jigsaw_configuration_json)
    sFilename_jigsaw_configuration_copy = os.path.join( sWorkspace_output_case, 'jigsaw_configuration_copy.json' )
    copy2(sFilename_jigsaw_configuration_json, sFilename_jigsaw_configuration_copy)

    #now we will update the configuration file with our own settings
    change_json_key_value(sFilename_configuration_copy, "iFlag_simplification", 0) #disable the flowline simplification
    change_json_key_value(sFilename_configuration_copy, "iFlag_create_mesh", 1)
    change_json_key_value(sFilename_configuration_copy, "iFlag_run_jigsaw", 1) #turn on the jigsaw mesh generation
    change_json_key_value(sFilename_configuration_copy, "iFlag_intersect", 0) #turn off the intersection
    change_json_key_value(sFilename_configuration_copy, "sFilename_mesh_boundary", sFilename_vector_geojson_merged) #set the target mesh boundary
    change_json_key_value(sFilename_configuration_copy, "sFilename_coastal_boundary", sFilename_vector_geojson_merged) #set the coastline boundary
    change_json_key_value(sFilename_configuration_copy, "sFilename_jigsaw_configuration", sFilename_jigsaw_configuration_copy)
    change_json_key_value(sFilename_configuration_copy, "sFilename_basins", sFilename_configuration_basins_copy) #set the river network
    #update basin configuration file below
    #the basin configuration is already edited during the simplification process

    #update the jigsaw configuration file below
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac", "true") #enable resolution control
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom", "true") # enable geometry control
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom_river_network", "true") #set the resolution
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_coastline", "true") #set the resolution for coastal line
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_coastal", dResolution_coastal) #set the resolution for coastal line
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_land", dResolution_land) #set the resolution for land
    change_json_key_value(sFilename_jigsaw_configuration_copy, "ncolumn_space", ncolumn) #set the resolution for x direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "nrow_space", nrow) #set the resolution for y direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_river_network_vector", sFilename_flowline_hydrosheds_out) #set the resolution for x direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_coastline_raster", sFilename_tif_wo_island) #set the resolution for y direction

    #now we can set up the actual pyflowline to create the mesh
    oPyflowline = pyflowline_read_configuration_file(sFilename_configuration_copy,
                    iCase_index_in=iCase_index,
                    sDate_in= sDate_today,
                    sMesh_type_in = sMesh_type)

    oPyflowline._pyflowline_create_hpc_job(sSlurm_in = 'slurm', hours_in = 30 )
    #now you should manually submit the job


else:
    pass

print('Congratulations! The workflow has been completed successfully!')
