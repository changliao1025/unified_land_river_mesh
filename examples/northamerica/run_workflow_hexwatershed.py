

#refer to the e3sm confluence page for more details of the workflow

#each step may require different input datasets
#however, all the output will be saved in the same output directory for easy access

#you can change this to your preferred output directory
import os
from shutil import copy2
from datetime import datetime
sWorkspace_input = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/conus/input'
sWorkspace_output = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/conus/output'
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/northamerica'

sDate_today = datetime.now().strftime('%Y%m%d')
sDate_today = '20250801'  #for debugging purpose, you can set this to a fixed date
#Step 0 set up the resolution
iCase_index = 2
iFlag_visualization = 0
iFlag_hpc_job = 1  #set to 1 if you want to run this on HPC
iFlag_run_hexwatershed_utility = 1
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
iFlag_simplify_hydrosheds_river_network = 0
iFlag_process_coastal_line = 0
iFlag_create_buffer_zone = 0
#Step 1
#prepare the river network and coastal line dataset
from pyearth.toolbox.management.vector.fields import get_field_and_value, add_field_to_vector_file
from pyflowline.algorithms.simplification.simplify_hydrosheds import simplify_hydrosheds_river_network
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.geometry.create_gcs_buffer_zone import create_buffer_zone_polygon_file
from pyearth.toolbox.conversion.convert_vector_to_geojson import convert_vector_to_geojson
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
from pyearth.toolbox.data.geojson.copy_geometry_without_attributes import copy_geometery_without_attributes
from pyhexwatershed.configuration.change_json_key_value import change_json_key_value
from hexwatershed_utility.mosart.convert_hexwatershed_output_to_mosart import convert_hexwatershed_json_to_mosart_netcdf
from pye3sm.mosart.map.unstructured.mosart_map_unstructured_parameters import mosart_map_unstructured_parameters
from pye3sm.mosart.map.unstructured.mosart_map_unstructured_flow_direction import mosart_map_unstructured_flow_direction
sFilename_flowline_hydrosheds_in=  '/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na.shp'
sFilename_flowline_hydroshed_tmp = 'HydroRIVERS_v10_na_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '.geojson'
sWorkspace_river_network = os.path.join(sWorkspace_output, 'river_network')
sWorkspace_watershed_boundary = os.path.join(sWorkspace_output, 'watershed_boundary')
sFilename_flowline_hydrosheds_out = os.path.join(sWorkspace_river_network, sFilename_flowline_hydroshed_tmp)
#get basename without folder and extension
sFilename_flowline_hydrosheds = os.path.splitext(os.path.basename(sFilename_flowline_hydrosheds_out))[0]


sFilename_geojson_geometery_feature = '/qfs/people/liao313/data/hexwatershed/global/vector/region.geojson'

#step 1: record attribute from the MPAS tools
aField, aValue = get_field_and_value(sFilename_geojson_geometery_feature)

#the river flowline simplficiation process already generated the basin configuration file
sFilename_pyhexwatershed_configuration = os.path.join(sWorkspace_output, 'pyhexwatershed_configuration.json')
sFilename_pyflowline_configuration= os.path.join(sWorkspace_output, 'pyflowline_configuration.json')
sFilename_pyhexwatershed_configuration_basins = os.path.join(sWorkspace_output, 'pyflowline_configuration_basins.json')

nOutlet = 12
if iFlag_simplify_hydrosheds_river_network == 1:
    simplify_hydrosheds_river_network(sFilename_flowline_hydrosheds_in,
                       sFilename_flowline_hydrosheds_out,
                       dDistance_tolerance_in,
                        dDrainage_area_threshold_in,
                        iFlag_pyflowline_configuration_in=1,
                        nOutlet_largest=nOutlet)

copy2(sFilename_pyflowline_configuration, sFilename_pyhexwatershed_configuration)
change_json_key_value(sFilename_pyhexwatershed_configuration, "sModel", 'pyhexwatershed')


sFilename_coastal_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrobasin/hybas_lake_na_lev01-12_v1c/hybas_lake_na_lev01_v1c.shp'
sFilename_geojson_na = os.path.join(sWorkspace_output, 'na_hydrobasin.geojson')
sFilename_geojson_wo_island  = os.path.join(sWorkspace_output, 'coastal_line_wo_island.geojson')
sFilename_tif_wo_island = os.path.join(sWorkspace_output, 'coastal_line_wo_island.tif')
sFilename_geojson_buffer = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_buffer_' + sResolution_coastal + '.geojson')
sFilename_vector_geojson_merged = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_merged_' + sResolution_coastal + '.geojson')
sFilename_dem = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrosheds/dem/na_dem_3s.tif'
sFilename_mpas_mesh_netcdf = '/compyfs/liao313/04model/pyhexwatershed/northamerica/pyflowline20250702011/jigsaw/out/invert_mesh.nc'

if iFlag_process_coastal_line == 1:
    convert_vector_to_geojson(sFilename_coastal_in, sFilename_geojson_na)
    remove_small_polygon(sFilename_geojson_na, sFilename_geojson_wo_island, dThreshold_area_island )
    #convert the vector to raster, this is a global raster, so we need to set the resolution
    convert_vector_to_global_raster(sFilename_geojson_wo_island,
                                    sFilename_tif_wo_island,
                                    dResolution_x_in, dResolution_y_in,
                                    iFlag_boundary_only_in = 0,
                                    dFill_value_in = 2)
    if iFlag_create_buffer_zone == 1:
        create_buffer_zone_polygon_file(sFilename_geojson_wo_island, sFilename_geojson_buffer,
                                          dBuffer_distance_in = dResolution_coastal * 1000)
    else:
        copy_geometery_without_attributes(sFilename_geojson_wo_island, sFilename_geojson_buffer)

    merge_features(sFilename_geojson_buffer, sFilename_vector_geojson_merged)
    add_field_to_vector_file(sFilename_vector_geojson_merged, aField, aValue)

#Step 2 - 4
#run the hexwatershed model, this step include three steps merged together.
#for debug purpose, you can also run then one by one, using the iFlag_debug flag to control

iFlag_debug = 1  #0 for normal run, 1 for debug

from pyhexwatershed.configuration.read_configuration_file import pyhexwatershed_read_configuration_file
from pyflowline.configuration.config_manager import create_template_jigsaw_configuration_file


from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
aExtent_na = gdal_get_vector_extent(sFilename_vector_geojson_merged)

aLegend = list()

dLongitude_delaware = -75.51052 # -75.51052,39.70195
dLatitude_delaware = 39.70195

dLongitude_susquehanna = -76.09755 #-76.09755,39.56820
dLatitude_susquehanna = 39.56820

if iFlag_debug == 1:
    #create the mesh using pyflowline

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
    copy2(sFilename_pyhexwatershed_configuration_basins, sFilename_configuration_basins_copy)

    for i in range(0, nOutlet):
        sBasin  = f'{i+1:04d}'
        #get basename without folder of sFilename_flowline_hydrosheds_out
        sFilename_flowline_filter_tmp = sFilename_flowline_hydrosheds_out.replace('.geojson', f'_{i+1:04d}_outlet_simplified.geojson')
        sFilename_flowline_filter = os.path.join(sWorkspace_river_network, sFilename_flowline_filter_tmp)
        sFilename_watershed_boundary = os.path.join(sWorkspace_watershed_boundary, f"watershed_boundary_{sBasin}.geojson")
        if i ==10:
            sFilename_flowline_filter = os.path.join(sWorkspace_river_network, 'delaware.geojson')
            sFilename_watershed_boundary = os.path.join(sWorkspace_watershed_boundary, f"watershed_boundary_delaware.geojson")
            #also change the outlet location
            change_json_key_value(sFilename_configuration_basins_copy, "dLongitude_outlet_degree", dLongitude_delaware, iFlag_basin_in=1, iBasin_index_in=i)
            change_json_key_value(sFilename_configuration_basins_copy, "dLatitude_outlet_degree", dLatitude_delaware, iFlag_basin_in=1, iBasin_index_in=i)
        else:
            if i ==11:
                sFilename_flowline_filter = os.path.join(sWorkspace_river_network, 'susquehanna.geojson')
                sFilename_watershed_boundary = os.path.join(sWorkspace_watershed_boundary, f"watershed_boundary_susquehanna.geojson")
                #also change the outlet location
                change_json_key_value(sFilename_configuration_basins_copy, "dLongitude_outlet_degree", dLongitude_susquehanna, iFlag_basin_in=1, iBasin_index_in=i)
                change_json_key_value(sFilename_configuration_basins_copy, "dLatitude_outlet_degree", dLatitude_susquehanna, iFlag_basin_in=1, iBasin_index_in=i)

        change_json_key_value(sFilename_configuration_basins_copy, "sFilename_flowline_filter", sFilename_flowline_filter, iFlag_basin_in=1, iBasin_index_in=i)
        change_json_key_value(sFilename_configuration_basins_copy, "sFilename_watershed_boundary", sFilename_watershed_boundary, iFlag_basin_in=1, iBasin_index_in=i)


    #now we also need to set the jigsaw configuration file
    sFilename_jigsaw_configuration_json = os.path.join(sWorkspace_output, 'pyflowline_configuration_jigsaw.json')
    create_template_jigsaw_configuration_file(sFilename_jigsaw_configuration_json)
    sFilename_jigsaw_configuration_copy = os.path.join( sWorkspace_output_case, 'jigsaw_configuration_copy.json' )
    copy2(sFilename_jigsaw_configuration_json, sFilename_jigsaw_configuration_copy)

    #now we will update the configuration file with our own settings
    change_json_key_value(sFilename_configuration_copy, "iFlag_simplification", 1) #disable the flowline simplification
    change_json_key_value(sFilename_configuration_copy, "iFlag_create_mesh", 1)
    change_json_key_value(sFilename_configuration_copy, "iFlag_mesh_boundary", 1) #enable the river network generation
    change_json_key_value(sFilename_configuration_copy, "iFlag_run_jigsaw", 1) #turn on the jigsaw mesh generation
    change_json_key_value(sFilename_configuration_copy, "iFlag_multiple_outlet", 1) #turn on the jigsaw mesh generation
    change_json_key_value(sFilename_configuration_copy, "iFlag_intersect", 1) #turn off the intersection
    change_json_key_value(sFilename_configuration_copy, "iFlag_stream_burning_topology", 1) #turn off the intersection
    change_json_key_value(sFilename_configuration_copy, "iFlag_export_individual_watershed", 1) #turn off the intersection
    change_json_key_value(sFilename_configuration_copy, "iFlag_force_watershed_boundary", 1) #turn on the intersection
    change_json_key_value(sFilename_configuration_copy, "sFilename_mesh_boundary", sFilename_vector_geojson_merged) #set the target mesh boundary
    change_json_key_value(sFilename_configuration_copy, "sFilename_mpas_mesh_netcdf", sFilename_mpas_mesh_netcdf) #set the target mesh boundary
    change_json_key_value(sFilename_configuration_copy, "sFilename_coastal_boundary", sFilename_vector_geojson_merged) #set the coastline boundary
    change_json_key_value(sFilename_configuration_copy, "sFilename_jigsaw_configuration", sFilename_jigsaw_configuration_copy)
    change_json_key_value(sFilename_configuration_copy, "sFilename_basins", sFilename_configuration_basins_copy) #set the river network
    #dem
    change_json_key_value(sFilename_configuration_copy, "iFlag_use_mesh_dem", 0) #enable the dem file
    change_json_key_value(sFilename_configuration_copy, "sFilename_dem", sFilename_dem) #set the dem file

    #update basin configuration file below
    #the basin configuration is already edited during the simplification process

    #update the jigsaw configuration file below
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac", True) #enable resolution control
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom", True) # enable geometry control
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_land", True) #set the resolution for x direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom_river_network", True) #set the resolution
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_coastline", True) #set the resolution for coastal line
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_coastline", dResolution_coastal) #set the resolution for coastal line
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_land", dResolution_land) #set the resolution for land
    change_json_key_value(sFilename_jigsaw_configuration_copy, "ncolumn_space", ncolumn) #set the resolution for x direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "nrow_space", nrow) #set the resolution for y direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_river_network_vector", sFilename_flowline_hydrosheds_out) #set the resolution for x direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_coastline_raster", sFilename_tif_wo_island) #set the resolution for y direction

    #now we can set up the actual pyflowline to create the mesh
    oPyhexwatershed = pyhexwatershed_read_configuration_file(sFilename_configuration_copy,
                    iCase_index_in=iCase_index,
                    sDate_in= sDate_today,
                    sMesh_type_in = sMesh_type)

    if iFlag_hpc_job == 1:
        oPyhexwatershed._pyhexwatershed_create_hpc_job(sSlurm_in = 'slurm', hours_in = 10 )
        #now you should manually submit the job
    else:
        #oPyhexwatershed.pyhexwatershed_export()
        pass

    if iFlag_visualization ==1:
        sVariable_in = 'flow_direction'

        sFilename = os.path.join( oPyhexwatershed.sWorkspace_output_hexwatershed, 'flow_direction.png' )
        oPyhexwatershed.plot(sFilename_output_in=sFilename, iFlag_title_in =1,
                             iFlag_esri_hydro_image_in=1,
                             sVariable_in = sVariable_in,
                             aExtent_in = aExtent_na,
                             iDPI_in = 600,
                             sFilename_boundary_in=sFilename_vector_geojson_merged)
        pass

    sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
    sFilename_mosart_parameter_out = os.path.join(oPyhexwatershed.sWorkspace_output_hexwatershed, 'mosart_northamerica_parameter.nc')
    sFilename_mosart_unstructured_domain= os.path.join(oPyhexwatershed.sWorkspace_output_hexwatershed, 'mosart_northamerica_domain.nc')
    sFilename_mosart_unstructured_script = os.path.join(oPyhexwatershed.sWorkspace_output_hexwatershed, 'mosart_northamerica_scriptgrid.nc')
    if iFlag_run_hexwatershed_utility == 1:
        sFilename_json_in = oPyhexwatershed.sFilename_hexwatershed_json
        convert_hexwatershed_json_to_mosart_netcdf(sFilename_json_in,
                sFilename_mosart_parameter_in,
                sFilename_mosart_parameter_out,
                sFilename_mosart_unstructured_domain)

    sFilename_geojson_out = os.path.join(oPyhexwatershed.sWorkspace_output_hexwatershed, 'mosart_northamerica_parameter.geojson')
    aVariable_parameter= ['rwid', 'rdep']
    aVariable_short= ['rwid', 'rdep']
    aTitle = ['Main channel width', 'Main channel depth']
    aColormap = ['YlGnBu', 'Blues']
    aFlag_colorbar = [1, 1]
    aFlag_scientific_notation_colorbar = [1, 1]
    aUnit = ['Unit: m', 'Unit: m']
    aData_max = [None, None]
    aData_min = [0, 0]
    iSize_x = 8
    iFlag_map_parameter = 1
    if iFlag_map_parameter == 1:
        mosart_map_unstructured_parameters(sFilename_mosart_unstructured_domain, sFilename_mosart_parameter_out,
                                       sFilename_geojson_out,
                                       aVariable_parameter,
                                       aVariable_short,
                                       aTitle,
                                        aFlag_colorbar_in = aFlag_colorbar,
                                       aFlag_scientific_notation_colorbar_in = aFlag_scientific_notation_colorbar,
                                       aUnit_in=aUnit,
                                        aColormap_in=aColormap,
                                       aData_max_in=aData_max,
                                       aData_min_in=aData_min,
                                          iSize_x_in=8,
                                            iSize_y_in=8)

    iFlag_map_flow_direction = 1
    sFilename_geojson_out = os.path.join(oPyhexwatershed.sWorkspace_output_hexwatershed, 'mosart_northamerica_flow_direction.geojson')
    if iFlag_map_flow_direction == 1:
        mosart_map_unstructured_flow_direction(sFilename_mosart_parameter_out,
                                           sFilename_geojson_out,
                                           sLengend_in=None,
                                           iSize_x_in = 8,
                                           iSize_y_in = 8,
                                            aLegend_in=None,
                                           aExtent_in=None)


else:
    pass

print('Congratulations! The workflow has been completed successfully!')
