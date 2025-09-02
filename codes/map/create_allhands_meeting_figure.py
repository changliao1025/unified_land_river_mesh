#i need to create a figure for the allhands meeting poster

#in this poster figure, I want to present a few key points:
#1. The global mesh  that include both land and ocean
#2. the land mesh, but that is already in the global mesh
#3. The watershed boudnary files
#4. The simplied river network that is produced from the hydrosheds
#5. The river network that is produced from the hexwatershed
#6. The coastaline,

import os
import datetime
import textwrap
import numpy as np
from urllib.error import URLError
from osgeo import  osr, gdal, ogr
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.cm as cm
from matplotlib.collections import PathCollection
from matplotlib.path import Path
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon as mpolygon
from matplotlib.patches import Polygon
import shapely.geometry as sgeom
import cartopy as cpl
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.system.define_global_variables import *
from pyearth.toolbox.math.stat.remap import remap
from pyearth.visual.map.zebra_frame import zebra_frame
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.visual.formatter import OOMFormatter
from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
from pyearth.visual.map.map_servers import StadiaStamen, EsriTerrain, EsriRelief, EsriHydro
from pyearth.visual.map.map_servers import Stadia_terrain_images, Esri_terrain_images, Esri_relief_images, Esri_hydro_images


#step 1, convert the global mesh to a geojson file

sFilename_geojson_base = '/compyfs/liao313/04model/pyhexwatershed/northamerica/pyflowline20250702011/jigsaw/out/base_mesh.geojson'

from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
dResolution_land = 25.0  #unit in km
dResolution_coastal = 10.0  #unit in km
dResolution_x_in = 30.0/3600 * dResolution_coastal
dResolution_y_in = dResolution_x_in
nrow = int(180 / dResolution_y_in)
ncolumn = int(360 / dResolution_x_in)
sResolution_coastal = "{:.1f}".format(dResolution_coastal) + 'km'


sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/northamerica'
sFilename_vector_geojson_merged = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_merged_' + sResolution_coastal + '.geojson')

aExtent_na = gdal_get_vector_extent(sFilename_vector_geojson_merged)
sFilename_boundary = '/qfs/people/liao313/data/hexwatershed/conus/vector/conus_boundary.geojson'
aExtent_conus = gdal_get_vector_extent(sFilename_boundary)
#give land mesh a different color
sFilename_cull_mesh = '/compyfs/liao313/04model/pyhexwatershed/northamerica/pyhexwatershed20250801001/pyflowline/mpas.geojson'
sFilename_coastal = '/compyfs/liao313/04model/pyhexwatershed/northamerica/pyhexwatershed20250801001/pyflowline/coastal_boundary.geojson'
pDriver = ogr.GetDriverByName('GeoJSON')
pSRS_wgs84 = ccrs.PlateCarree()  # for latlon data only
pSRS_geodetic = ccrs.Geodetic()
pProjection_data = pSRS_wgs84
dLon_max = aExtent_conus[1]
dLon_min = aExtent_conus[0]
dLat_max = aExtent_conus[3]
dLat_min = aExtent_conus[2]
aExtent = [dLon_min, dLon_max, dLat_min, dLat_max]
pProjection_map = ccrs.Orthographic(central_longitude =  0.50*(dLon_max+dLon_min),
                                            central_latitude = 0.50*(dLat_max+dLat_min),
                                            globe=None)
iDPI = 300
fig = plt.figure( dpi = iDPI )
fig_width_inches = 12
fig_height_inches = 8
fig.set_figwidth( fig_width_inches )
fig.set_figheight( fig_height_inches )
plot_width_inch = fig.get_size_inches()[0] * fig.dpi
char_width_inch = 0.1 * fig.dpi
cwidth = int(plot_width_inch / char_width_inch)
ax = fig.add_axes([0.08, 0.1, 0.62, 0.7], projection= pProjection_map  ) #projection=ccrs.PlateCarree()
ax.set_global()
ax.coastlines(color='black', linewidth=1,resolution='10m')
ax.set_extent(aExtent, crs = pSRS_wgs84)
minx, maxx, miny,  maxy = aExtent

image_size = [2000, 1000]
image_size = [int(fig_width_inches * iDPI), int(fig_height_inches * iDPI)]  # [7200, 4800]
scale_denominator = calculate_scale_denominator(aExtent, image_size)
pSrc = osr.SpatialReference()
pSrc.ImportFromEPSG(3857) # mercator
pProjection = pSrc.ExportToWkt()
iBasemap_zoom_level = calculate_zoom_level(scale_denominator, pProjection, dpi=int(iDPI))
iFlag_basemap_hydro = 1

if iFlag_basemap_hydro == 1:
    esri_terrain = EsriTerrain()
    esri_hydro = EsriHydro()
    ll_target_domain = sgeom.box(minx, miny, maxx, maxy)
    multi_poly = esri_hydro.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
    target_domain = multi_poly.geoms[0]
    _, aExtent_terrain, _ = esri_terrain.image_for_domain(target_domain, iBasemap_zoom_level)
    _, aExtent_hydro, _ = esri_hydro.image_for_domain(target_domain, iBasemap_zoom_level)
    img_esri_terrain  = Esri_terrain_images(aExtent, iBasemap_zoom_level)
    img_eari_hydro  = Esri_hydro_images(aExtent, iBasemap_zoom_level)
    ax.imshow(img_esri_terrain,  extent=aExtent_terrain, transform=esri_terrain.crs)
    ax.imshow(img_eari_hydro,  extent=aExtent_hydro, transform=esri_hydro.crs)
    #add the license information
    sLicense_info = "© Esri Hydro Reference Overlay"
    sLicense_info_wrapped = "\n".join(textwrap.wrap(sLicense_info, width=60))
    ax.text(0.5, 0.05, sLicense_info_wrapped, transform=ax.transAxes, ha='center', va='center', fontsize=6,
            color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

iFlag_basemap_terrain = 0
if iFlag_basemap_terrain == 1:
    stamen_terrain = StadiaStamen()
    ll_target_domain = sgeom.box(minx, miny, maxx, maxy)
    multi_poly = stamen_terrain.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
    target_domain = multi_poly.geoms[0]
    _, aExtent_terrain, _ = stamen_terrain.image_for_domain(target_domain, iBasemap_zoom_level)
    img_stadia_terrain = Stadia_terrain_images(aExtent, iBasemap_zoom_level)
    ax.imshow(img_stadia_terrain,  extent=aExtent_terrain, transform=stamen_terrain.crs)
    #add the license information
    sLicense_info = "© Stamen Design, under a Creative Commons Attribution (CC BY 3.0) license."
    sLicense_info_wrapped = "\n".join(textwrap.wrap(sLicense_info, width=cwidth))
    ax.text(0.5, 0.05, sLicense_info_wrapped, transform=ax.transAxes, ha='center', va='center', fontsize=6,
            color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

dAlpha = 0.8

#open the global mesh

#sColor_base_mesh = 'lightsteelblue'
#sColor_land_mesh = 'lightgreen'
#sColor_watershed_boundary = 'lightblue'
#sColor_hydrosheds_flowline = 'black'
#sColor_hexwatershed_flowline = 'blue'
# Improved color scheme with better contrast
sColor_base_mesh = "#2B41E2"  # Very light gray for background mesh
sColor_land_mesh = '#2E7D32'  # Dark green for land mesh
sColor_watershed_boundary = '#1976D2'  # Blue for watersheds
sColor_hydrosheds_flowline = '#00838F'  # Dark teal
sColor_hexwatershed_flowline = '#0D47A1'  # Dark blue for your method


iFlag_plot_base_mesh = 1
if iFlag_plot_base_mesh == 1:
    #open the global mesh
    pDataset = ogr.Open(sFilename_geojson_base, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    #pLayer.SetSpatialFilterRect(minx, miny, maxx, maxy)
    aPolygon = []
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if sGeometry_type == 'POLYGON':
            aCoords_gcs = get_geometry_coordinates(pGeometry_in)
            # Ensure the polygon is valid and has enough points
            if len(aCoords_gcs) >= 3:
                x_coords = aCoords_gcs[:, 0]
                y_coords = aCoords_gcs[:, 1]
                #use limit check to filter out the polygons that are not in the extent
                if (np.max(x_coords) - np.min(x_coords))> 10 or (np.max(y_coords) - np.min(y_coords)) > 10:
                    continue
                else:
                    aPolygon.append(aCoords_gcs[:, 0:2])
                # Plot the polygon outline
                #ax.plot(x_coords, y_coords, color='grey', alpha=dAlpha, linewidth=0.5, transform=pProjection_data)

    # Filter out invalid polygons
    aPolygon = [poly for poly in aPolygon if len(poly) >= 3]
    if aPolygon:  # Only create collection if we have valid polygons
        aPatch = [Polygon(poly, closed=True) for poly in aPolygon]
        pPC = PatchCollection(aPatch, alpha=0.9, edgecolor=sColor_base_mesh,
                                  facecolor='none', linewidths=0.25,
                            linestyles='dotted', transform=pProjection_data)
        ax.add_collection(pPC)



#next the land mesh
iFlag_plot_land_mesh = 1
if iFlag_plot_land_mesh == 1:
    pDataset = ogr.Open(sFilename_cull_mesh, gdal.GA_ReadOnly)
    pLayer_land = pDataset.GetLayer(0)
    aPolygon_land = []
    for pFeature in pLayer_land:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if sGeometry_type == 'POLYGON':
            aCoords_gcs = get_geometry_coordinates(pGeometry_in)
            if len(aCoords_gcs) >= 3:
                aPolygon_land.append(aCoords_gcs[:, 0:2])

    # Filter out invalid polygons
    aPolygon_land = [poly for poly in aPolygon_land if len(poly) >= 3]

    if aPolygon_land:  # Only create collection if we have valid polygons
        aPatch_land = [Polygon(poly, closed=True) for poly in aPolygon_land]
        pPC_land = PatchCollection(aPatch_land, alpha=0.6, edgecolor=sColor_land_mesh,
                                  facecolor='none', linewidths=0.25,
                                  transform=pProjection_data)
        ax.add_collection(pPC_land)
#save the mesh as a pdf

#now watershed boundary

nOutlet = 12
sWorkspace_watershed_boundary  = '/compyfs/liao313/04model/pyhexwatershed/northamerica/watershed_boundary'
sWorkspace_river_network = '/compyfs/liao313/04model/pyhexwatershed/northamerica/river_network'
sFilename_flowline_hydrosheds_out=  os.path.join(sWorkspace_river_network,  'HydroRIVERS_v10_na_simplified_1.25E+04_3.12E+09.geojson')
iFlag_plot_hydrosheds_flowline = 1
if iFlag_plot_hydrosheds_flowline == 1:
    aPolyline = []
    pDataset = ogr.Open(sFilename_flowline_hydrosheds_out, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if sGeometry_type =='LINESTRING':
            aCoords_gcs = get_geometry_coordinates(pGeometry_in)
            aCoords_gcs = aCoords_gcs[:,0:2]
            nvertex = len(aCoords_gcs)
            codes = np.full(nvertex, mpl.path.Path.LINETO, dtype=int )
            codes[0] = mpl.path.Path.MOVETO
            path = mpl.path.Path(aCoords_gcs, codes)
            x, y = zip(*path.vertices)
            aPolyline.append(list(zip(x, y)))
    pLC = LineCollection(aPolyline,  alpha=dAlpha, edgecolor=sColor_hydrosheds_flowline,
                         facecolor='none', linewidths=0.25, transform=pProjection_data)

    ax.add_collection(pLC)

iFlag_plot_hydrosheds_flowline_each = 0
for i in range(nOutlet):
    sBasin  = f'{i+1:04d}'
    sFilename_flowline_filter_tmp = sFilename_flowline_hydrosheds_out.replace('.geojson', f'_{i+1:04d}_outlet_simplified.geojson')
    sFilename_flowline_filter = os.path.join(sWorkspace_river_network, sFilename_flowline_filter_tmp)
    sFilename_watershed_boundary = os.path.join(sWorkspace_watershed_boundary, f"watershed_boundary_{sBasin}.geojson")
    if i ==10:
        sFilename_flowline_filter = os.path.join(sWorkspace_river_network, 'delaware.geojson')
        sFilename_watershed_boundary = os.path.join(sWorkspace_watershed_boundary, f"watershed_boundary_delaware.geojson")
    else:
        if i ==11:
            sFilename_flowline_filter = os.path.join(sWorkspace_river_network, 'susquehanna.geojson')
            sFilename_watershed_boundary = os.path.join(sWorkspace_watershed_boundary, f"watershed_boundary_susquehanna.geojson")
    #draw the watershed boundary
    pDataset = ogr.Open(sFilename_watershed_boundary, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if sGeometry_type == 'POLYGON':
            aCoords_gcs = get_geometry_coordinates(pGeometry_in)
            mPolygon = mpolygon(aCoords_gcs, closed=True, edgecolor=sColor_watershed_boundary, facecolor='none',
                                      alpha=dAlpha,  transform=pProjection_data)
            # Add the polygon to the map
            ax.add_patch(mPolygon)
        else:
            if sGeometry_type == 'MULTIPOLYGON':
                for j in range(pGeometry_in.GetGeometryCount()):
                    pPolygon = pGeometry_in.GetGeometryRef(j)
                    aCoords_gcs = get_geometry_coordinates(pPolygon)
                    mPolygon = mpolygon(aCoords_gcs, closed=True, edgecolor=sColor_watershed_boundary, facecolor='none',
                                              alpha=dAlpha,  transform=pProjection_data)
                    # Add the polygon to the map
                    ax.add_patch(mPolygon)
            else:
                print('Not a polygon geometry type:', sGeometry_type)

    #draw the river network
    if iFlag_plot_hydrosheds_flowline_each == 1:
        pDataset = ogr.Open(sFilename_flowline_filter, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        aPolyline = []

        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='LINESTRING':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                aCoords_gcs = aCoords_gcs[:,0:2]
                nvertex = len(aCoords_gcs)
                codes = np.full(nvertex, mpl.path.Path.LINETO, dtype=int )
                codes[0] = mpl.path.Path.MOVETO
                path = mpl.path.Path(aCoords_gcs, codes)
                x, y = zip(*path.vertices)
                aPolyline.append(list(zip(x, y)))
        pLC = LineCollection(aPolyline,  alpha=dAlpha, edgecolor='black',
                             facecolor='none', linewidths=0.25, transform=pProjection_data)

        ax.add_collection(pLC)

    pass

#draw the hydroshed flowline


#draw the flow direction
sFilename_flow_direction= '/compyfs/liao313/04model/pyhexwatershed/northamerica/pyhexwatershed20250801001/hexwatershed/mpas_flow_direction.geojson'
pDataset = ogr.Open(sFilename_flow_direction, gdal.GA_ReadOnly)
pLayer = pDataset.GetLayer(0)
iThickness_min = 0.25
iThickness_max = 2.5

aPolyline = []
aThickness = []
aColors = []
#get the min and max drainage area
aDrainage_area = []
for pFeature in pLayer:
    pGeometry_in = pFeature.GetGeometryRef()
    sGeometry_type = pGeometry_in.GetGeometryName()
    dValue_thickness = pFeature.GetField('drainage_area')
    if sGeometry_type == 'LINESTRING':
        aDrainage_area.append(dValue_thickness)
aDrainage_area = np.array(aDrainage_area)
dValue_thickness_max = np.max(aDrainage_area)
dValue_thickness_min = np.min(aDrainage_area)
for pFeature in pLayer:
    pGeometry_in = pFeature.GetGeometryRef()
    sGeometry_type = pGeometry_in.GetGeometryName()
    dValue_thickness = pFeature.GetField('drainage_area')
    if sGeometry_type =='LINESTRING':
        aCoords_gcs = get_geometry_coordinates(pGeometry_in)
        aCoords_gcs = aCoords_gcs[:,0:2]
        nvertex = len(aCoords_gcs)
        codes = np.full(nvertex, mpl.path.Path.LINETO, dtype=int )
        codes[0] = mpl.path.Path.MOVETO
        path = mpl.path.Path(aCoords_gcs, codes)
        x, y = zip(*path.vertices)
        iThickness = remap( dValue_thickness, dValue_thickness_min, dValue_thickness_max, iThickness_min, iThickness_max )

        # Add color variation based on drainage area
        if dValue_thickness > dValue_thickness_max * 0.1:  # Major rivers
            color = '#0D47A1'  # Dark blue
        elif dValue_thickness > dValue_thickness_max * 0.01:  # Medium rivers
            color = '#1976D2'  # Medium blue
        else:  # Small streams
            color = '#42A5F5'  # Light blue
        aColors.append(color)
        aThickness.append(iThickness)
        aPolyline.append(list(zip(x, y)))
pLC = LineCollection(aPolyline,  alpha=dAlpha, edgecolor=aColors,
                     facecolor='none', linewidths=aThickness, transform=pProjection_data)
ax.add_collection(pLC)


#add the title
sTitle = 'Unified Land-River-Ocean Mesh for North America:\nIntegrated Flow Direction Network'
sTitle_wrapped = "\n".join(textwrap.wrap(sTitle, width=cwidth))
ax.set_title(sTitle_wrapped, fontsize=10, fontweight='bold', pad=10, loc='left', color='black')

# Replace text labels with proper legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Create legend elements
legend_elements = [
    Patch(facecolor='none', alpha=0.3, edgecolor=sColor_base_mesh,
          label='Ocean Mesh'),
    Patch(facecolor='none', alpha=0.3, edgecolor=sColor_land_mesh,
          label='Land Mesh'),
    Line2D([0], [0], color=sColor_watershed_boundary, linewidth=2,
           label='Watershed Boundary'),
    Line2D([0], [0], color=sColor_hydrosheds_flowline, linewidth=1,
           label='HydroSHEDS River Network'),
    Line2D([0], [0], color=sColor_hexwatershed_flowline, linewidth=2,
           label='Hexwatershed Flow Direction')
]

# Position legend better
ax.legend(handles=legend_elements, loc='upper left',
         bbox_to_anchor=(0.05, 0.95), fontsize=6,
         frameon=True, fancybox=True, shadow=True)


gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--',
                      xlocs=np.arange(minx, maxx+(maxx-minx)/9, (maxx-minx)/8),
                      ylocs=np.arange(miny, maxy+(maxy-miny)/9, (maxy-miny)/8))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mpl.ticker.MaxNLocator(4)
gl.ylocator = mpl.ticker.MaxNLocator(4)
gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation': 0, 'ha': 'right'}
gl.ylabel_style = {'size': 10, 'color': 'k',
                       'rotation': 90, 'weight': 'normal'}
ax.set_axis_off()
ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)

ax.set_extent(aExtent, crs = pSRS_wgs84)
pDataset = pLayer = pFeature  = None


sFilename_output_in = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/figures/allhands_meeting_poster_figure.png'

sDirname = os.path.dirname(sFilename_output_in)
sFilename = os.path.basename(sFilename_output_in)
sFilename_out = os.path.join(sDirname, sFilename)
sExtension = os.path.splitext(sFilename)[1]
if sExtension == '.png':
    plt.savefig(sFilename_out, bbox_inches='tight')
else:
    if sExtension == '.pdf':
        plt.savefig(sFilename_out, bbox_inches='tight')
    else:
        plt.savefig(sFilename_out, bbox_inches='tight', format ='ps')
#clean cache
plt.close('all')
plt.clf()
print('The plot is saved as ', sFilename_out)