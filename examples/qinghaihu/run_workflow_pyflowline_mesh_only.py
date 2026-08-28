#refer to the e3sm confluence page for more details of the workflow

#each step may require different input datasets
#however, all the output will be saved in the same output directory for easy access

#you can change this to your preferred output directory
import os
import glob
from datetime import datetime
from shutil import copy2

from pyearth.toolbox.management.vector.fields import get_field_and_value, add_field_to_vector_file
from pyearth.gis.gdal.write.vector.gdal_write_wkt_to_vector_file import gdal_write_wkt_to_vector_file
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
from pyearth.toolbox.analysis.image.raster_process import create_raster_buffer_zone, fix_raster_antimeridian_issue

from hexwatershed_utility.preprocess.feature.river.simplify_hydrorivers_networks import simplify_hydrorivers_networks
from hexwatershed_utility.preprocess.feature.river.tag_river_outlet import tag_river_outlet
from hexwatershed_utility.preprocess.feature.watershed_boundary.find_minimal_hydrobasins_watershed_boundary import find_minimal_hydrobasins_watershed_boundary
from hexwatershed_utility.preprocess.feature.river.get_outlet_location import get_outlet_location
from hexwatershed_utility.preprocess.feature.coastline.create_land_ocean_mask_from_hydrobasin import create_land_ocean_mask_from_hydrobasin
from hexwatershed_utility.preprocess.feature.coastline.create_land_ocean_mask_from_naturalearth import create_land_ocean_mask_from_naturalearth
from hexwatershed_utility.preprocess.feature.coastline.fix_naturalearth_hydrosheds_incompatibility import fix_naturalearth_hydrosheds_incompatibility

from pyflowline.configuration.config_manager import create_pyflowline_template_configuration_file
from pyflowline.configuration.read_configuration_file import pyflowline_read_configuration_file
from pyflowline.configuration.config_manager import create_jigsaw_template_configuration_file
from pyflowline.configuration.change_json_key_value import change_json_key_value

##======================================================================
# The only thing you need to change for different runs is the settings below,
# which is used to control the resolution and other settings for the mesh generation. You can also set up different flags to turn on/off certain process for debugging purpose. The output will be saved in the same output directory, which is defined below as well.
# start of common user settings
##======================================================================
#date time for simulation
sDate_today = '20260802'  #use a fixed date for easy repeatability
sDate_today = datetime.now().strftime('%Y%m%d')
#index for different runs
iCase_index = 1

#flag for component
iFlag_run_jigsaw = 1
#setup flags for debugging
iFlag_simplify_hydrosheds_river_network = 0
iFlag_flexible_river_mouth = 0
iFlag_process_watershed_boundary = 0
iFlag_process_coastline = 0

iFlag_reprocess_config = 0

#number of largest outlet to be processed
nOutlet_largest = 100

#resolution settings
dResolution_ocean = 100
dResolution_land = 12  #unit in km
dResolution_river_network = 3
dResolution_river_outlet = 3
dResolution_coastline = 10  #unit in km
dResolution_lake_boundary = 3
dResolution_watershed_boundary =3

if platform == 'Linux':
    sWorkspace_input = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/global/input'
    sWorkspace_output = '/data2/share/liaochang/04model/jigsaw/qinghaihu'
    #define global output directory
    sWorkspace_river_network_output = '/data2/share/liaochang/04model/jigsaw/qinghaihu/river_network'   
    sWorkspace_coastline_output = '/data2/share/liaochang/04model/jigsaw/global/coastline'    
    sWorkspace_data = '/public/home/liaochang/data/hexwatershed/qinghaihu/'
    sFilename_flowline_hydrosheds_in = '/data2/share/liaochang/data/raw/hydrology/hydrosheds/hydroriver/asian/HydroRIVERS_v10_as_shp/HydroRIVERS_v10_as.shp'
    sFilename_flowline_hydroshed_tmp = 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '.geojson'
    sFilename_geojson_geometery_feature = '/public/home/liaochang/data/hexwatershed/global/vector/region.geojson'
    sWorkspace_watershed_boundary_in = '/data2/share/liaochang/data/raw/hydrology/hydrosheds/hydrobasin'
else:
    if platform == 'Windows':
        sWorkspace_input = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/global/input'
        sWorkspace_output = 'D:\\scratch\\04model\\jigsaw\\qinghaihu'
        sWorkspace_river_network_output = 'D:\\scratch\\04model\\jigsaw\\qinghaihu\\river_network'   
        sWorkspace_coastline_output = 'D:\\scratch\\04model\\jigsaw\\global\\coastline'   
        sWorkspace_data = 'D:\\data\\modeldata\\hexwatershed\\qinghaihu'
        sFilename_flowline_hydroshed_tmp = 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '.geojson'
        sFilename_geojson_geometery_feature = 'D:\\data\\modeldata\\hexwatershed\\global\\vector\\region.geojson'
        sWorkspace_watershed_boundary_in = '/data2/share/liaochang/data/raw/hydrology/hydrosheds/hydrobasin'


#sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/global/'
#sWorkspace_watershed_boundary_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrobasin'
#sFilename_flowline_hydrosheds_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp'
#sFilename_geojson_geometery_feature = '/qfs/people/liao313/data/hexwatershed/global/vector/region.geojson'

##======================================================================
# end of common user settings
##======================================================================

##======================================================================
# start uncommon user settings
##======================================================================
iFlag_dam = 0
sFilename_dam = '/compyfs/liao313/00raw/dam/GRanD_Version_1_3/GRanD_dams_v1_3_merged.geojson' #should consider both on and snapped dams in this dataset

##======================================================================
# end of uncommon user settings
##======================================================================

#for coastline
dThreshold_area_island = dResolution_ocean * dResolution_ocean * 10 * 1.0E6  #unit m2, this one may need to be adjusted based on the resolution
dResolution_coastline_buffer = dResolution_coastline * 1.0E3  #buffer zone for coastline line
#small island removal threshold
dDrainage_area_threshold= dResolution_land * dResolution_land * 100 * 1.0E6  #at least ten grid cells of drainage area, this may be adjusted as well


if os.path.exists(sWorkspace_river_network_output) is False:
    os.makedirs(sWorkspace_river_network_output)    
if os.path.exists(sWorkspace_coastline_output) is False:
    os.makedirs(sWorkspace_coastline_output)

#add the threshol into the output folder
sWorkspace_river_network_output = os.path.join(sWorkspace_river_network_output,  sDistance_tolerance + '_' + sDrainage_area_threshold)
if os.path.exists(sWorkspace_river_network_output) is False:
    os.makedirs(sWorkspace_river_network_output)

sWorkspace_coastline_output = os.path.join(sWorkspace_coastline_output,  sCoastline_buffer + '_' + sThreshold_area_island )
if os.path.exists(sWorkspace_coastline_output) is False:
    os.makedirs(sWorkspace_coastline_output)

sMesh_type = 'mpas'  #

#for jigsaw resolution control
dResolution_x_in = 30.0/3600 * dResolution_coastline
dResolution_y_in = dResolution_x_in
nrow = int(180 / dResolution_y_in)
ncolumn = int(360 / dResolution_x_in)

#string format for file names
dDistance_tolerance = dResolution_river_network * 1.0E3 #how far away two river need to be for mesh generation
sDistance_tolerance = "{:.2E}".format(dDistance_tolerance)
sDrainage_area_threshold = "{:.2E}".format(dDrainage_area_threshold) # m2
#coastline second
sCoastline_buffer = "{:.1E}".format(dResolution_coastline_buffer  ) # to m
sThreshold_area_island = "{:.1E}".format(dThreshold_area_island ) # to m2



print(sWorkspace_river_network_output)
print(sWorkspace_watershed_boundary_output)
print(sWorkspace_coastline_output)

#Step 1
#prepare the river network and coastline line dataset
sFilename_flowline_hydroshed_tmp = 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '.geojson'
sFilename_flowline_hydroshed_outlet = 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '_outlet.geojson'
sFilename_flowline_hydrosheds_out = os.path.join(sWorkspace_river_network_output, sFilename_flowline_hydroshed_tmp)

#step 1: record attribute from the MPAS tools
aField, aValue = get_field_and_value(sFilename_geojson_geometery_feature)

#the river flowline simplficiation process already generated the basin configuration file
sFilename_pyflowline_configuration = os.path.join(sWorkspace_river_network_output, 'pyflowline_configuration.json')
sFilename_pyflowline_configuration_basins = os.path.join(sWorkspace_river_network_output, 'pyflowline_configuration_basins.json')
sFilename_river_network_raster = os.path.join(sWorkspace_river_network_output, 'river_network_raster.tif')

if iFlag_simplify_hydrosheds_river_network == 1:
    simplify_hydrorivers_networks(sFilename_flowline_hydrosheds_in,
                       sFilename_flowline_hydrosheds_out,
                       dDistance_tolerance,
                        dDrainage_area_threshold,
                        iFlag_pyflowline_configuration_in=1,
                        nOutlet_largest=nOutlet_largest)
else:
    sFilename_flowline_hydrosheds_out = '/public/home/liaochang/data/hexwatershed/qinghaihu/vector/flowline_hydroshed_simplified_clipped_clean.geojson'
    if iFlag_reprocess_config == 1:
        #if flowline is pre-processed, we just need to create the configuration file
        create_pyflowline_template_configuration_file(sFilename_pyflowline_configuration,
                sWorkspace_river_network_output = sWorkspace_river_network_output,
                iFlag_standalone_in=1,
                nOutlet = nOutlet_largest,
                sMesh_type_in='mpas',
                sModel_in='pyflowline')

        for i in range(1, nOutlet_largest+1):
            sBasin_id = '{:04d}'.format(i)
            sFilename_flowline_simplified_basin = os.path.join(sWorkspace_river_network_output, 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold +'_'+ sBasin_id + '.geojson')
            #read the geojson file using gdal to get the outlet location, which is the end point of the first feature
            dLongitude_outlet, dLatitude_outlet = get_outlet_location(sFilename_flowline_simplified_basin)
            change_json_key_value(sFilename_pyflowline_configuration_basins, 'dLatitude_outlet_degree', dLatitude_outlet, iFlag_basin_in=1, iBasin_index_in=i-1)
            change_json_key_value(sFilename_pyflowline_configuration_basins, 'dLongitude_outlet_degree', dLongitude_outlet, iFlag_basin_in=1, iBasin_index_in=i-1)
            change_json_key_value(sFilename_pyflowline_configuration_basins, 'sFilename_flowline_filter', sFilename_flowline_simplified_basin, iFlag_basin_in=1, iBasin_index_in=i-1)




if iFlag_process_watershed_boundary == 1:
    for i in range(1, nOutlet_largest+1):
        sBasin_id  = f'{i:04d}'
        aFile = glob.glob(os.path.join(sWorkspace_river_network_output, f'*{sBasin_id}.geojson'))
        if len(aFile) == 0:
            print(f"No file found for basin {sBasin_id}")
            continue
        sFilename_river_network_in = aFile[0]
        wkt = find_minimal_hydrobasins_watershed_boundary(sFilename_river_network_in, sWorkspace_watershed_boundary_in)
        if wkt is not None:
            sFilename_out = os.path.join(sWorkspace_watershed_boundary_output, f"watershed_boundary_{sDistance_tolerance}_{sDrainage_area_threshold}_{sBasin_id}.geojson")
            gdal_write_wkt_to_vector_file(wkt, sFilename_out)
            change_json_key_value(sFilename_pyflowline_configuration_basins, 'sFilename_watershed_boundary', sFilename_out, iFlag_basin_in=1, iBasin_index_in=i-1)

else:
    #only need to update the basin boundary file path in the configuration file
    for i in range(1, nOutlet_largest+1):
        sBasin_id = '{:04d}'.format(i)
        sFilename_watershed_boundary_basin = os.path.join(sWorkspace_watershed_boundary_output, f"watershed_boundary_{sDistance_tolerance}_{sDrainage_area_threshold}_{sBasin_id}.geojson")
        change_json_key_value(sFilename_pyflowline_configuration_basins, 'sFilename_watershed_boundary', sFilename_watershed_boundary_basin, iFlag_basin_in=1, iBasin_index_in=i-1)
    pass


sFilename_vector_coastline_merged = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_merged.geojson')
sFilename_tif_wo_island = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.tif')
sFilename_tif_wo_island_buffered = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_buffered.tif')
sFilename_tif_wo_island_buffered_fixed = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_buffered_fixed.tif')
if iFlag_process_coastline == 1:
    sFilename_tif_wo_island, sFilename_vector_coastline = create_land_ocean_mask_from_naturalearth(sWorkspace_coastline_output,
                                                                             dResolution_x_in, dResolution_y_in,
                                                                             dThreshold_area_island,
                                                                             dResolution_coastline_buffer,
                                                                             iRaster_buffer_pixel=2)

    ##fix the incompatibilty between hydrosheds and naturalearth
    aFilename_flowline = list()
    for i in range(1, nOutlet_largest+1):
        sBasin_id = '{:04d}'.format(i)
        sFilename_flowline_simplified_basin = os.path.join(sWorkspace_river_network_output, 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold +'_'+ sBasin_id + '.geojson')
        aFilename_flowline.append(sFilename_flowline_simplified_basin)

    sFilename_vector_coastline_updated = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_fixed.geojson')
    fix_naturalearth_hydrosheds_incompatibility(aFilename_flowline, sFilename_vector_coastline, sFilename_vector_coastline_updated )
    #should be merged into one single function
    merge_features(sFilename_vector_coastline_updated, sFilename_vector_coastline_merged, iFlag_force= True)
    add_field_to_vector_file(sFilename_vector_coastline_merged, aField, aValue)
else:
    sFilename_wo_island = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.geojson')
    convert_vector_to_global_raster(sFilename_wo_island,
                                    sFilename_tif_wo_island,
                                    dResolution_x_in,
                                    dResolution_y_in,
                                    iFlag_boundary_only_in = 0,
                                    dFill_value_in = 2)
    create_raster_buffer_zone(sFilename_tif_wo_island, sFilename_tif_wo_island_buffered, 1, 2)
    fix_raster_antimeridian_issue(sFilename_tif_wo_island_buffered, sFilename_tif_wo_island_buffered_fixed,1, 2, iRaster_buffer_pixel=2)

sFilename_tif_wo_island = sFilename_tif_wo_island_buffered_fixed

#Step 2 - 4
#run the hexwatershed model, this step include three steps merged together.
try:
    #copy the template configuration file to the output directory
    change_json_key_value(sFilename_pyflowline_configuration, "sWorkspace_output", sWorkspace_output)
    oPyflowline = pyflowline_read_configuration_file(sFilename_pyflowline_configuration, \
    iCase_index_in=iCase_index, sDate_in=sDate_today)
    sWorkspace_output_case = oPyflowline.sWorkspace_output
    sFilename_configuration_copy = os.path.join(sWorkspace_output_case, 'pyflowline_configuration_copy.json')
    copy2(sFilename_pyflowline_configuration, sFilename_configuration_copy)
    #copy the basin configuration file to the output directory as well
    sFilename_configuration_basins_copy = os.path.join( sWorkspace_output_case, 'pyflowline_configuration_basins_copy.json' )
    copy2(sFilename_pyflowline_configuration_basins, sFilename_configuration_basins_copy)
    sFilename_jigsaw_configuration_json = os.path.join(sWorkspace_river_network_output, 'pyflowline_configuration_jigsaw.json')
    create_jigsaw_template_configuration_file(sFilename_jigsaw_configuration_json)
    sFilename_jigsaw_configuration_copy = os.path.join( sWorkspace_output_case, 'jigsaw_configuration_copy.json' )
    change_json_key_value(sFilename_configuration_copy, "iFlag_run_jigsaw", iFlag_run_jigsaw) #turn off the jigsaw mesh generation

    #now we will update the configuration file with our own settings
    change_json_key_value(sFilename_configuration_copy, "iFlag_simplification", 0) #disable the flowline simplification
    change_json_key_value(sFilename_configuration_copy, "iFlag_create_mesh", 1)
    change_json_key_value(sFilename_configuration_copy, "iFlag_global", 1)
    change_json_key_value(sFilename_configuration_copy, "iFlag_intersect", 0) #turn off the intersection
    change_json_key_value(sFilename_configuration_copy, "sFilename_coastline_boundary", sFilename_vector_coastline_merged) #set the coastline boundary
    change_json_key_value(sFilename_configuration_copy, "sFilename_jigsaw_configuration", sFilename_jigsaw_configuration_copy)
    change_json_key_value(sFilename_configuration_copy, "iFlag_force_watershed_boundary", 1) #turn on the intersection
    change_json_key_value(sFilename_configuration_copy, "sFilename_basins", sFilename_configuration_basins_copy) #set the river network

    if iFlag_run_jigsaw == 1:
        copy2(sFilename_jigsaw_configuration_json, sFilename_jigsaw_configuration_copy)
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom", "true") # enable geometry control
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom_river_network", "true") #set the resolution
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom_dam", "true")
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac", "true") #enable resolution control
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_ocean", "true")
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_river_network", "true")
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_RRS18to6_ocean", "true")
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_land", "true")
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_coastline", "true") #set the resolution for coastline line
        change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_ocean", dResolution_ocean) #set the resolution for ocean
        change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_coastline", dResolution_coastline) #set the resolution for coastline line
        change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_land", dResolution_land) #set the resolution for land
        change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_river_network", dResolution_river_network) #set the resolution for river network
        change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_river_outlet", dResolution_river_outlet) #set the small island removal threshold
        change_json_key_value(sFilename_jigsaw_configuration_copy, "ncolumn_space", ncolumn) #set the resolution for x direction
        change_json_key_value(sFilename_jigsaw_configuration_copy, "nrow_space", nrow) #set the resolution for y direction
        change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_dam_vector", sFilename_dam) #set the dam file
        change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_river_network_vector", sFilename_flowline_hydrosheds_out) #set the resolution for x direction
        change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_coastline_raster", sFilename_tif_wo_island) #set the resolution for y direction
        change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_river_network_raster", sFilename_river_network_raster) #set the small island removal threshold

    #now we can set up the actual pyflowline to create the mesh
    oPyflowline = pyflowline_read_configuration_file(sFilename_configuration_copy,
                    iCase_index_in=iCase_index,
                    sDate_in= sDate_today,
                    sMesh_type_in = sMesh_type)

    oPyflowline._pyflowline_create_hpc_job(sSlurm_in = 'slurm', hours_in = 8 )
    #now you should manually submit the job

except Exception as e:
    print(f"An error occurred: {e}")

print('Congratulations! The workflow has been completed successfully!', sWorkspace_output_case )
