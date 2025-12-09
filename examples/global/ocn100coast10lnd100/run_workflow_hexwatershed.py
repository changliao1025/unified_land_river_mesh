

#refer to the e3sm confluence page for more details of the workflow

#each step may require different input datasets
#however, all the output will be saved in the same output directory for easy access

#you can change this to your preferred output directory
import os
import glob
from datetime import datetime
from shutil import copy2
#prepare the river network and coastal line dataset
from pyearth.toolbox.management.vector.fields import get_field_and_value, add_field_to_vector_file
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.geometry.create_gcs_buffer_zone import create_buffer_zone_polygon_file
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
from pyearth.gis.gdal.write.vector.gdal_write_wkt_to_vector_file import gdal_write_wkt_to_vector_file

from hexwatershed_utility.preprocess.features.rivers.simplify_hydrorivers_networks import simplify_hydrorivers_networks
from hexwatershed_utility.preprocess.features.rivers.get_outlet_location import get_outlet_location
from hexwatershed_utility.preprocess.features.watershed_boundary.find_minimal_hydrobasins_watershed_boundary import find_minimal_hydrobasins_watershed_boundary
from hexwatershed_utility.preprocess.features.coastline.create_land_ocean_mask_from_hydrobasin import create_land_ocean_mask_from_hydrobasin
from hexwatershed_utility.preprocess.features.coastline.create_land_ocean_mask_from_naturalearth import create_land_ocean_mask_from_naturalearth
from hexwatershed_utility.preprocess.features.coastline.fix_naturalearth_hydrosheds_incompatibility import fix_naturalearth_hydrosheds_incompatibility

from pyflowline.configuration.config_manager import create_jigsaw_template_configuration_file

from pyhexwatershed.configuration.config_manager import create_pyhexwatershed_template_configuration_file
from pyhexwatershed.configuration.read_configuration_file import pyhexwatershed_read_configuration_file
from pyhexwatershed.configuration.change_json_key_value import change_json_key_value

sDate_today = datetime.now().strftime('%Y%m%d')

sDate_today = '20250928'  #use a fixed date for easy repeatability

#Step 0 set up the resolution
iCase_index = 1


#for coastline
#resolution settings
#for rivers and watershed
dResolution_land = 100  #unit in km
dResolution_coastline = 10.0  #unit in km

#for coastline
dThreshold_area_island = 1.0E3  #unit km2, this one may need to be adjusted based on the resolution
dResolution_coastline_buffer = dResolution_coastline   #buffer zone for coastline line
#small island removal threshold

dDrainage_area_threshold= dResolution_land * dResolution_land *10  #at least ten grid cells of drainage area, this may be adjusted as well

#setup flags for debugging
iFlag_simplify_hydrosheds_river_network = 0
iFlag_process_watershed_boundary = 0
iFlag_process_coastline = 0

nOutlet_largest = 5


sWorkspace_input = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/global/input'
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/global/'
sWorkspace_river_network_output = '/compyfs/liao313/04model/pyhexwatershed/global/river_network'
if os.path.exists(sWorkspace_river_network_output) is False:
    os.makedirs(sWorkspace_river_network_output)

sWorkspace_watershed_boundary_output = '/compyfs/liao313/04model/pyhexwatershed/global/watershed_boundary'
if os.path.exists(sWorkspace_watershed_boundary_output) is False:
    os.makedirs(sWorkspace_watershed_boundary_output)

sWorkspace_coastline_output = '/compyfs/liao313/04model/pyhexwatershed/global/coastline'
if os.path.exists(sWorkspace_coastline_output) is False:
    os.makedirs(sWorkspace_coastline_output)



sMesh_type = 'mpas'  #
#for jigsaw resolution control
dResolution_x_in = 30.0/3600 * dResolution_coastline
dResolution_y_in = dResolution_x_in
nrow = int(180 / dResolution_y_in)
ncolumn = int(360 / dResolution_x_in)

#string format for file names
#river first
dDistance_tolerance = dResolution_land * 1.0E3 #how far away two river need to be for mesh generation
sDistance_tolerance = "{:.2E}".format(dDistance_tolerance)
sDrainage_area_threshold = "{:.2E}".format(dDrainage_area_threshold * 1.0E6) #convert to m2

#coastline second
sCoastline_buffer = "{:.1E}".format(dResolution_coastline_buffer * 1.0E3) #convert to m
sThreshold_area_island = "{:.1E}".format(dThreshold_area_island * 1.0E6) #convert to m2

sWorkspace_river_network_output = os.path.join(sWorkspace_river_network_output,  sDistance_tolerance + '_' + sDrainage_area_threshold)
if os.path.exists(sWorkspace_river_network_output) is False:
    os.makedirs(sWorkspace_river_network_output)
sWorkspace_watershed_boundary_output = os.path.join(sWorkspace_watershed_boundary_output,  sDistance_tolerance + '_' + sDrainage_area_threshold)
if os.path.exists(sWorkspace_watershed_boundary_output) is False:
    os.makedirs(sWorkspace_watershed_boundary_output)


#Step 1

sFilename_flowline_hydrosheds_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp'
sFilename_flowline_hydroshed_tmp = 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '.geojson'
sFilename_flowline_hydrosheds_out = os.path.join(sWorkspace_river_network_output, sFilename_flowline_hydroshed_tmp)
sFilename_geojson_geometery_feature = '/qfs/people/liao313/data/hexwatershed/global/vector/region.geojson'
sFolder_watershed_boundary_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrobasin'
#step 1: record attribute from the MPAS tools
aField, aValue = get_field_and_value(sFilename_geojson_geometery_feature)

#the river flowline simplficiation process already generated the basin configuration file
sFilename_pyhexwatershed_configuration = os.path.join(sWorkspace_river_network_output, 'pyhexwatershed_configuration.json')
sFilename_pyflowline_configuration_basins = os.path.join(sWorkspace_river_network_output, 'pyflowline_configuration_basins.json')

if iFlag_simplify_hydrosheds_river_network == 1:
    simplify_hydrorivers_networks(sFilename_flowline_hydrosheds_in,
                       sFilename_flowline_hydrosheds_out,
                       dDistance_tolerance,
                        dDrainage_area_threshold,
                        iFlag_pyflowline_configuration_in=1,
                        nOutlet_largest=nOutlet_largest)
else:
    #if flowline is pre-processed, we just need to create the configuration file
    sFilename_configuration_json = os.path.join(sWorkspace_river_network_output, 'pyhexwatershed_configuration.json')
    create_pyhexwatershed_template_configuration_file(sFilename_configuration_json,
            sWorkspace_river_network_output = sWorkspace_river_network_output,
            nOutlet = nOutlet_largest,
            sMesh_type_in='mpas',
            sModel_in='pyhexwatershed')
    sFilename_configuration_basin_json = os.path.join(sWorkspace_river_network_output, 'pyflowline_configuration_basins.json')
    for i in range(1, nOutlet_largest+1):
        sBasin_id = '{:04d}'.format(i)
        sFilename_flowline_simplified_basin = os.path.join(sWorkspace_river_network_output, 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold +'_'+ sBasin_id + '.geojson')
        dLongitude_outlet, dLatitude_outlet = get_outlet_location(sFilename_flowline_simplified_basin)
        change_json_key_value(sFilename_configuration_basin_json, 'dLatitude_outlet_degree', dLatitude_outlet, iFlag_basin_in=1, iBasin_index_in=i-1)
        change_json_key_value(sFilename_configuration_basin_json, 'dLongitude_outlet_degree', dLongitude_outlet, iFlag_basin_in=1, iBasin_index_in=i-1)
        change_json_key_value(sFilename_configuration_basin_json, 'sFilename_flowline_filter', sFilename_flowline_simplified_basin, iFlag_basin_in=1, iBasin_index_in=i-1)

if iFlag_process_watershed_boundary ==1:
     for i in range(1, nOutlet_largest+1):
        sBasin_id  = f'{i:04d}'
        aFile = glob.glob(os.path.join(sWorkspace_river_network_output, f'*{sBasin_id}.geojson'))
        if len(aFile) == 0:
            print(f"No file found for basin {sBasin_id}")
            continue
        sFilename_river_network_in = aFile[0]
        wkt = find_minimal_hydrobasins_watershed_boundary(sFilename_river_network_in, sFolder_watershed_boundary_in)
        if wkt is not None:
            sFilename_out = os.path.join(sWorkspace_watershed_boundary_output, f"watershed_boundary_{sDistance_tolerance}_{sDrainage_area_threshold}_{sBasin_id}.geojson")
            gdal_write_wkt_to_vector_file(wkt, sFilename_out)
            change_json_key_value(sFilename_configuration_basin_json, 'sFilename_watershed_boundary', sFilename_out, iFlag_basin_in=1, iBasin_index_in=i-1)
else:
    #only need to update the basin boundary file path in the configuration file
    sFilename_configuration_json = os.path.join(sWorkspace_river_network_output, 'pyhexwatershed_configuration.json')
    sFilename_configuration_basin_json = os.path.join(sWorkspace_river_network_output, 'pyflowline_configuration_basins.json')
    for i in range(1, nOutlet_largest+1):
        sBasin_id = '{:04d}'.format(i)
        sFilename_watershed_boundary_basin = os.path.join(sWorkspace_watershed_boundary_output, f"watershed_boundary_{sDistance_tolerance}_{sDrainage_area_threshold}_{sBasin_id}.geojson")
        change_json_key_value(sFilename_configuration_basin_json, 'sFilename_watershed_boundary', sFilename_watershed_boundary_basin, iFlag_basin_in=1, iBasin_index_in=i-1)


sFilename_mpas_mesh_netcdf = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250926008/jigsaw/out/invert_mesh.nc'
sFilename_vector_coastline_updated = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_fixed.geojson')
sFilename_vector_coastline_merged = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_merged.geojson')
if iFlag_process_coastline == 1:
    sFilename_tif_wo_island, sFilename_vector_coastline = create_land_ocean_mask_from_naturalearth(sWorkspace_coastline_output,
                                                                             dResolution_x_in, dResolution_y_in,
                                                                             dThreshold_area_island,
                                                                             dResolution_coastline_buffer)

    aFilename_flowline = list()
    for i in range(1, nOutlet_largest+1):
        sBasin_id = '{:04d}'.format(i)
        sFilename_flowline_simplified_basin = os.path.join(sWorkspace_river_network_output, 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold +'_'+ sBasin_id + '.geojson')
        aFilename_flowline.append(sFilename_flowline_simplified_basin)
    fix_naturalearth_hydrosheds_incompatibility(aFilename_flowline, sFilename_vector_coastline, sFilename_vector_coastline_updated )
    #should be merged into one single function
    merge_features(sFilename_vector_coastline_updated, sFilename_vector_coastline_merged, iFlag_force= True)
    add_field_to_vector_file(sFilename_vector_coastline_merged, aField, aValue)
else:
    sFilename_tif_wo_island = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.tif')
    pass

#Step 2 - 4
iFlag_debug = 1  #0 for normal run, 1 for debug

if iFlag_debug == 1:
    #copy the template configuration file to the output directory
    change_json_key_value(sFilename_pyhexwatershed_configuration, "sWorkspace_output", sWorkspace_output)
    oPyhexwatershed = pyhexwatershed_read_configuration_file(sFilename_pyhexwatershed_configuration, \
    iCase_index_in=iCase_index, sDate_in=sDate_today)
    sWorkspace_output_case = oPyhexwatershed.sWorkspace_output
    sFilename_configuration_copy = os.path.join(sWorkspace_output_case, 'pyhexwatershed_configuration_copy.json')
    copy2(sFilename_pyhexwatershed_configuration, sFilename_configuration_copy)
    #if the stream burning is turned on, a basin configuration file is also automatically created under the same output directory
    #copy the basin configuration file to the output directory as well
    sFilename_configuration_basins_copy = os.path.join( sWorkspace_output_case, 'pyflowline_configuration_basins_copy.json' )
    copy2(sFilename_pyflowline_configuration_basins, sFilename_configuration_basins_copy)

    #now we also need to set the jigsaw configuration file
    sFilename_jigsaw_configuration_json = os.path.join(sWorkspace_river_network_output, 'pyflowline_configuration_jigsaw.json')
    create_jigsaw_template_configuration_file(sFilename_jigsaw_configuration_json)
    sFilename_jigsaw_configuration_copy = os.path.join( sWorkspace_output_case, 'jigsaw_configuration_copy.json' )
    copy2(sFilename_jigsaw_configuration_json, sFilename_jigsaw_configuration_copy)

    #now we will update the configuration file with our own settings
    change_json_key_value(sFilename_configuration_copy, "iFlag_simplification", 1) #disable the flowline simplification
    change_json_key_value(sFilename_configuration_copy, "iFlag_create_mesh", 1)
    change_json_key_value(sFilename_configuration_copy, "iFlag_global", 1)
    change_json_key_value(sFilename_configuration_copy, "iFlag_run_jigsaw", 0) #turn off the jigsaw mesh generation
    change_json_key_value(sFilename_configuration_copy, "iFlag_intersect", 1) #turn off the intersection
    #change_json_key_value(sFilename_configuration_copy, "sFilename_mesh_boundary", sFilename_vector_geojson_merged) #set the target mesh boundary
    change_json_key_value(sFilename_configuration_copy, "sFilename_coastline_boundary", sFilename_vector_coastline_merged) #set the coastline boundary
    change_json_key_value(sFilename_configuration_copy, "sFilename_jigsaw_configuration", sFilename_jigsaw_configuration_copy)
    change_json_key_value(sFilename_configuration_copy, "sFilename_mpas_mesh_netcdf", sFilename_mpas_mesh_netcdf) #set the target mesh boundary
    change_json_key_value(sFilename_configuration_copy, "iFlag_force_watershed_boundary", 1) #turn on the intersection
    change_json_key_value(sFilename_configuration_copy, "sFilename_basins", sFilename_configuration_basins_copy) #set the river network
    #update basin configuration file below
    #the basin configuration is already edited during the simplification process

    #update the jigsaw configuration file below
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac", "true") #enable resolution control
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom", "true") # enable geometry control
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom_river_network", "true") #set the resolution
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_coastline", "true") #set the resolution for coastal line
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_coastline", dResolution_coastline) #set the resolution for coastal line
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_land", dResolution_land) #set the resolution for land
    change_json_key_value(sFilename_jigsaw_configuration_copy, "ncolumn_space", ncolumn) #set the resolution for x direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "nrow_space", nrow) #set the resolution for y direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_river_network_vector", sFilename_flowline_hydrosheds_out) #set the resolution for x direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_coastline_raster", sFilename_tif_wo_island) #set the resolution for y direction

    #now we can set up the actual pyhexwatershed to create the mesh
    oPyhexwatershed = pyhexwatershed_read_configuration_file(sFilename_configuration_copy,
                    iCase_index_in=iCase_index,
                    sDate_in= sDate_today,
                    sMesh_type_in = sMesh_type)

    oPyhexwatershed._pyhexwatershed_create_hpc_job(sSlurm_in = 'slurm', hours_in = 30 )
    #now you should manually submit the job


else:
    pass

print('Congratulations! The workflow has been completed successfully!')
