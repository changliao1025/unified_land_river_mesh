"""
Visualization module for uraster class - Refactored with pyearthviz3d.

This module contains visualization methods refactored to use the pyearthviz3d library
for improved code organization and functionality.

Features:
- 3D mesh visualization using GeoVista via pyearthviz3d
- Interactive and static rendering modes
- Comprehensive error handling and validation
- Support for multiple output formats
"""

import os
import logging
import traceback
from typing import Optional, Dict, Any, Tuple
import numpy as np
from osgeo import gdal, ogr, osr
gdal.UseExceptions()

from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import extract_unique_vertices_and_connectivity
from pyearthviz3d.geovista.utility import (
    VisualizationConfig,
    ScalarBarConfig,
    PlotterManager,
    configure_camera_enhanced,
    add_geographic_context_enhanced,
    add_mesh_to_plotter,
    VALID_ANIMATION_FORMATS
)

# Set up logging
logger = logging.getLogger(__name__)
CRS = "EPSG:4326"

# Constants for visualization
DEFAULT_ZOOM_FACTOR = 0.7
VALID_IMAGE_FORMATS = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']


def load_mpas_mesh_optimized(filename: str, verbose: bool = False) -> Optional[Dict[str, np.ndarray]]:
    """
    Optimized MPAS mesh loading with chunked reading and memory efficiency.

    Args:
        filename: Path to MPAS mesh NetCDF file
        verbose: Enable verbose logging

    Returns:
        Dictionary containing mesh data or None if failed
    """
    try:
        import netCDF4 as nc

        with nc.Dataset(filename, 'r') as dataset:
            # Read coordinates with optimized chunking
            lon_vertex = dataset.variables['lonVertex'][:]
            lat_vertex = dataset.variables['latVertex'][:]
            # Convert to degrees efficiently using vectorized operations
            lon_vertex = np.rad2deg(lon_vertex)
            lat_vertex = np.rad2deg(lat_vertex)
            # Fix longitude range [0, 360] to [-180, 180] efficiently
            lon_vertex = np.where(lon_vertex > 180, lon_vertex - 360, lon_vertex)
            connectivity = dataset.variables['verticesOnCell'][:]
            cell_ids = dataset.variables['indexToCellID'][:]

        return {
            'vertices_longitude': lon_vertex,
            'vertices_latitude': lat_vertex,
            'connectivity': connectivity,
            'cell_ids': cell_ids
        }

    except Exception as e:
        logger.error(f"Failed to load MPAS mesh from {filename}: {e}")
        return None


def optimize_connectivity_processing(connectivity: np.ndarray, num_vertices: int, verbose: bool = False) -> np.ndarray:
    """
    Optimized connectivity processing with vectorized operations.

    Args:
        connectivity: Raw connectivity array (1-based indexing)
        num_vertices: Number of vertices in mesh
        verbose: Enable verbose logging

    Returns:
        Processed connectivity array with masking
    """
    # Convert 1-based to 0-based indexing efficiently
    fill_value = 0
    connectivity_0based = np.where(connectivity > fill_value, connectivity - 1, -1)
    # Vectorized validation
    valid_mask = (connectivity_0based >= 0) & (connectivity_0based < num_vertices)
    # Count valid vertices per cell efficiently
    valid_count = np.sum(valid_mask, axis=1)
    if verbose:
        logger.info(f"Connectivity stats - Valid vertices per cell: "
                   f"min={np.min(valid_count)}, max={np.max(valid_count)}, "
                   f"mean={np.mean(valid_count):.1f}")
    # Create masked array for invalid indices
    connectivity_masked = np.ma.masked_where(~valid_mask, connectivity_0based)
    return connectivity_masked


def process_polylines_optimized(filename: str, line_width: float = 1.0, verbose: bool = False,
                               width_attribute: str = None, width_min: float = 0.5,
                               width_max: float = 3.0) -> Optional[Dict[str, Any]]:
    """
    Optimized polyline processing with batch operations and attribute-based line width scaling.

    Args:
        filename: Path to polyline vector file
        line_width: Default line width (used when width_attribute is None)
        verbose: Enable verbose logging
        width_attribute: Attribute name to use for line width scaling (e.g., 'discharge', 'flow_rate')
        width_min: Minimum line width for scaling (default: 0.5)
        width_max: Maximum line width for scaling (default: 3.0)

    Returns:
        Dictionary containing processed polyline data with scaled line widths or None if failed
    """
    try:
        # Open dataset
        dataset = ogr.Open(filename, gdal.GA_ReadOnly)
        if dataset is None:
            logger.error(f'Could not open vector file: {filename}')
            return None
        layer = dataset.GetLayer(0)
        if layer is None:
            logger.error('Could not get layer from dataset')
            return None
        feature_count = layer.GetFeatureCount()
        if verbose:
            logger.info(f'Processing {feature_count} polyline features')
        if feature_count == 0:
            logger.error('No features found in the input file')
            return None

        # Get layer definition for attribute extraction
        layer_defn = layer.GetLayerDefn()
        field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]

        # Check if width attribute exists
        use_attribute_width = False
        if width_attribute is not None:
            if width_attribute in field_names:
                use_attribute_width = True
                if verbose:
                    logger.info(f'Using attribute "{width_attribute}" for line width scaling')
            else:
                logger.warning(f'Width attribute "{width_attribute}" not found in file. Available: {field_names}')
                logger.warning(f'Using default line width: {line_width}')

        # Batch process geometries with attribute extraction
        all_points = []
        all_lines = []
        attribute_values = []
        point_offset = 0
        valid_polylines = 0
        layer.ResetReading()
        for feature in layer:
            geometry = feature.GetGeometryRef()
            if geometry is None:
                continue
            geom_name = geometry.GetGeometryName()
            if geom_name not in ['LINESTRING', 'MULTILINESTRING']:
                continue
            # Clone and transform geometry if needed
            geom_clone = geometry.Clone()
            # Extract points efficiently
            points = geom_clone.GetPoints()
            if points is None or len(points) < 2:
                continue
            # Convert to 3D coordinates on unit sphere (vectorized)
            coords_array = np.array(points)
            lon_rad = np.radians(coords_array[:, 0])
            lat_rad = np.radians(coords_array[:, 1])
            x = np.cos(lat_rad) * np.cos(lon_rad)
            y = np.cos(lat_rad) * np.sin(lon_rad)
            z = np.sin(lat_rad)
            # Add points to batch
            sphere_points = np.column_stack([x, y, z])
            all_points.extend(sphere_points.tolist())
            # Create line connectivity
            num_points_in_line = len(points)
            line = [num_points_in_line]
            line.extend(range(point_offset, point_offset + num_points_in_line))
            all_lines.extend(line)
            point_offset += num_points_in_line
            valid_polylines += 1
            # Extract attribute value for line width scaling
            if use_attribute_width:
                try:
                    attr_value = feature.GetField(width_attribute)
                    if attr_value is not None:
                        attribute_values.append(float(attr_value))
                    else:
                        attribute_values.append(0.0)  # Use 0 for missing values
                except Exception as e:
                    if verbose:
                        logger.warning(f'Could not read attribute {width_attribute}: {e}')
                    attribute_values.append(0.0)
            else:
                attribute_values.append(line_width)

        if valid_polylines == 0:
            logger.error('No valid polylines found for processing')
            return None
        if verbose:
            logger.info(f'Successfully processed {valid_polylines} polylines with {len(all_points)} total points')

        # Calculate scaled line widths based on attribute values
        scaled_widths = []
        if use_attribute_width and len(attribute_values) > 0:
            # Calculate min/max range of attribute values (excluding zeros for missing data)
            valid_values = [v for v in attribute_values if v > 0]
            if len(valid_values) > 0:
                data_min, data_max = min(valid_values), max(valid_values)
                if data_max > data_min:
                    # Scale attribute values to [width_min, width_max] range
                    for val in attribute_values:
                        if val > 0:
                            scaled_width = width_min + (val - data_min) / (data_max - data_min) * (width_max - width_min)
                            scaled_widths.append(scaled_width)
                        else:
                            scaled_widths.append(width_min)  # Use minimum width for missing/zero values
                    if verbose:
                        logger.info(f'Scaled line widths using attribute "{width_attribute}":')
                        logger.info(f'  - Attribute range: {data_min:.2f} to {data_max:.2f}')
                        logger.info(f'  - Width range: {width_min:.2f} to {width_max:.2f}')
                        logger.info(f'  - Average width: {np.mean(scaled_widths):.2f}')
                else:
                    # All values are the same, use default width
                    scaled_widths = [line_width] * valid_polylines
                    if verbose:
                        logger.info(f'All attribute values are identical ({data_min}), using default width')
            else:
                # No valid values found, use default width
                scaled_widths = [line_width] * valid_polylines
                if verbose:
                    logger.warning('No valid attribute values found, using default width')
        else:
            # Use uniform line width
            scaled_widths = [line_width] * valid_polylines

        return {
            'points': all_points,
            'lines': all_lines,
            'polyline_count': valid_polylines,
            'line_width': line_width,  # Default/fallback width
            'line_widths': scaled_widths,  # Per-polyline scaled widths
            'attribute_values': attribute_values,  # Raw attribute values
            'use_variable_width': use_attribute_width and len(valid_values) > 0,
            'width_attribute': width_attribute,
            'width_range': (width_min, width_max)
        }

    except Exception as e:
        logger.error(f'Error processing polylines from {filename}: {e}')
        return None


# Caching mechanism for expensive operations
_mesh_cache = {}


def get_cached_mesh(filename: str, cache_key: str = None) -> Optional[Dict[str, np.ndarray]]:
    """Get cached mesh data if available."""
    key = cache_key or filename
    if key in _mesh_cache:
        logger.info(f"Using cached mesh data for {filename}")
        return _mesh_cache[key]
    return None


def cache_mesh_data(filename: str, data: Dict[str, np.ndarray], cache_key: str = None):
    """Cache mesh data for future use."""
    key = cache_key or filename
    _mesh_cache[key] = data
    logger.info(f"Cached mesh data for {filename}")


def clear_mesh_cache():
    """Clear mesh cache to free memory."""
    global _mesh_cache
    _mesh_cache.clear()
    logger.info("Mesh cache cleared")


def process_elevation_data_optimized(filename: str, variable: str, verbose: bool = False) -> Optional[Dict[str, Any]]:
    """
    Optimized elevation data processing with caching and efficient operations.

    Args:
        filename: Path to elevation raster file
        variable: Variable name to extract
        verbose: Enable verbose logging

    Returns:
        Dictionary containing elevation data or None if failed
    """
    try:
        # Check cache first
        cache_key = f"{filename}_{variable}"
        cached_data = get_cached_mesh(filename, cache_key)
        if cached_data is not None:
            return cached_data

        # Rebuild mesh topology efficiently
        field_unique_id = 'cellid'
        raster_mesh_info = rebuild_mesh_topology(
            filename,
            sField_unique_id=field_unique_id,
            iFlag_verbose_in=verbose
        )
        if raster_mesh_info is None:
            logger.warning(f'Failed to rebuild mesh topology for raster: {filename}')
            return None

        # Extract mesh data
        vertex_longitude = raster_mesh_info['vertices_longitude']
        vertex_latitude = raster_mesh_info['vertices_latitude']
        connectivity = raster_mesh_info['connectivity']

        # Process connectivity efficiently
        connectivity_masked = np.ma.masked_where(connectivity == -1, connectivity)

        # Validate connectivity indices
        valid_connectivity = connectivity[connectivity >= 0]
        if len(valid_connectivity) > 0:
            max_vertex_idx = len(vertex_longitude) - 1
            if np.max(valid_connectivity) > max_vertex_idx:
                logger.error(f'Connectivity contains invalid vertex index: '
                           f'max={np.max(valid_connectivity)}, vertices={len(vertex_longitude)}')
                return None

        # Extract elevation data
        elevation_result = _extract_target_mesh_data(filename, variable, verbose)
        if elevation_result is None:
            logger.error(f'Failed to extract elevation data from {filename}')
            return None

        elevation_data, num_features = elevation_result

        result = {
            'vertices_longitude': vertex_longitude,
            'vertices_latitude': vertex_latitude,
            'connectivity_masked': connectivity_masked,
            'elevation_data': elevation_data,
            'variable_name': variable
        }

        return result

    except Exception as e:
        logger.error(f'Error processing elevation data from {filename}: {e}')
        return None


def rebuild_mesh_topology(sFilename_mesh_in: str, iFlag_verbose_in: bool = False, sField_unique_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """
    Rebuild mesh topology from source mesh file by extracting vertices,
    connectivity, and centroids for unstructured mesh processing.

    Args:
        sFilename_mesh_in (str): Path to the source mesh file (GeoJSON, Shapefile, etc.)
        iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
        sField_unique_id (str, optional): Field name for unique cell IDs.

    Returns:
        dict: Comprehensive mesh topology information
        Returns None on failure.
    """
    try:
        # Open the input data source
        pDataset = ogr.Open(sFilename_mesh_in, 0)  # Read-only
        if pDataset is None:
            logger.error(f'Failed to open file: {sFilename_mesh_in}')
            return None
        if iFlag_verbose_in:
            logger.info(f'Successfully opened mesh file: {sFilename_mesh_in}')

        # Get the first layer
        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            logger.error('Failed to get layer from the dataset.')
            pDataset = None
            return None

        # Get layer information
        pLayerDefn = pLayer.GetLayerDefn()
        if pLayerDefn is None:
            logger.error('Failed to get layer definition.')
            pDataset = None
            return None

        nFeatures = pLayer.GetFeatureCount()
        iFieldCount = pLayerDefn.GetFieldCount()
        if nFeatures == 0:
            logger.warning('Layer contains no features.')
            pDataset = None
            return None

        aCellID = []
        if iFlag_verbose_in:
            logger.info(f'Processing {nFeatures} features with {iFieldCount} fields')

        # Get the field name
        if sField_unique_id is None:
            sVariable = pLayerDefn.GetFieldDefn(0).GetName() if iFieldCount > 0 else None
        else:
            sVariable = sField_unique_id

        # Initialize lists for storing geometry data
        lons_list = []
        lats_list = []
        aArea_list = []

        # Process features
        pLayer.ResetReading()
        iFeature_index = 0
        invalid_geometry_count = 0

        for pFeature in pLayer:
            if pFeature is None:
                iFeature_index += 1
                continue

            pGeometry = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry.GetGeometryName()

            if sGeometry_type == 'POLYGON':
                try:
                    aCoord = get_geometry_coordinates(pGeometry)
                    if aCoord is not None and len(aCoord) >= 3:
                        lons = aCoord[:, 0]
                        lats = aCoord[:, 1]
                        lons_list.append(lons)
                        lats_list.append(lats)

                        try:
                            dArea = calculate_polygon_area(lons, lats)
                            aArea_list.append(dArea)
                        except Exception as area_error:
                            logger.warning(f'Could not calculate area for feature {iFeature_index}: {area_error}')
                            aArea_list.append(0.0)

                        if sVariable:
                            try:
                                field_value = pFeature.GetFieldAsInteger(sVariable)
                                aCellID.append(int(field_value))
                            except (ValueError, TypeError, AttributeError) as e:
                                logger.warning(f'Could not read integer field value for feature {iFeature_index}: {e}')
                                aCellID.append(len(aCellID))
                        else:
                            aCellID.append(len(aCellID))
                    else:
                        invalid_geometry_count += 1
                except Exception as e:
                    logger.warning(f'Error processing feature {iFeature_index}: {str(e)}')
                    invalid_geometry_count += 1

            elif sGeometry_type == 'MULTIPOLYGON':
                try:
                    if iFlag_verbose_in:
                        logger.info(f'Processing multipolygon feature {iFeature_index}')

                    current_cellid = None
                    if sVariable:
                        try:
                            field_value = pFeature.GetFieldAsInteger(sVariable)
                            current_cellid = int(field_value)
                        except (ValueError, TypeError, AttributeError) as e:
                            logger.warning(f'Could not read integer field value for multipolygon feature {iFeature_index}: {e}')
                            current_cellid = len(aCellID)
                    else:
                        current_cellid = len(aCellID)

                    for iPart in range(pGeometry.GetGeometryCount()):
                        pPolygon_part = pGeometry.GetGeometryRef(iPart)
                        if pPolygon_part is None:
                            continue

                        aCoord_part = get_geometry_coordinates(pPolygon_part)
                        if aCoord_part is not None and len(aCoord_part) >= 3:
                            lons_part = aCoord_part[:, 0]
                            lats_part = aCoord_part[:, 1]
                            lons_list.append(lons_part)
                            lats_list.append(lats_part)

                            try:
                                dArea_part = calculate_polygon_area(lons_part, lats_part)
                                aArea_list.append(dArea_part)
                            except Exception as area_error:
                                logger.warning(f'Could not calculate area for multipolygon part {iPart}: {area_error}')
                                aArea_list.append(0.0)

                            aCellID.append(current_cellid)

                except Exception as e:
                    logger.warning(f'Error processing multipolygon feature {iFeature_index}: {str(e)}')
                    invalid_geometry_count += 1

            iFeature_index += 1

        # Clean up dataset
        pDataset = None

        if not lons_list:
            logger.error('No valid polygon features found in mesh file')
            return None

        # Calculate maximum vertices and pad coordinates
        max_vertices = max(len(coord) for coord in lons_list)
        nVertex_max = max_vertices

        # Pre-allocate arrays for better memory efficiency
        num_polygons = len(lons_list)
        lons_padded = np.full((num_polygons, max_vertices), np.nan, dtype=np.float64)
        lats_padded = np.full((num_polygons, max_vertices), np.nan, dtype=np.float64)

        # Fill padded arrays efficiently
        for i, (lon_coords, lat_coords) in enumerate(zip(lons_list, lats_list)):
            lon_coords = np.asarray(lon_coords, dtype=np.float64)
            lat_coords = np.asarray(lat_coords, dtype=np.float64)
            coord_len = min(len(lon_coords), len(lat_coords))
            if coord_len > 0:
                lons_padded[i, :coord_len] = lon_coords[:coord_len]
                lats_padded[i, :coord_len] = lat_coords[:coord_len]

        lons = lons_padded
        lats = lats_padded

        # Calculate centroids
        cell_lons_1d = np.zeros(len(lons_list), dtype=np.float64)
        cell_lats_1d = np.zeros(len(lats_list), dtype=np.float64)

        for i in range(len(lons_list)):
            valid_mask = ~np.isnan(lons[i])
            if np.any(valid_mask):
                valid_lons = lons[i][valid_mask]
                valid_lats = lats[i][valid_mask]
                cell_lons_1d[i] = np.mean(valid_lons)
                cell_lats_1d[i] = np.mean(valid_lats)

        # Extract unique vertices and connectivity
        xv, yv, connectivity, vertex_to_index = extract_unique_vertices_and_connectivity(lons_list, lats_list)

        if xv is None or yv is None or connectivity is None:
            logger.error('Failed to extract unique vertices and connectivity')
            return None

        # Ensure aCellID matches the number of valid mesh cells
        if len(aCellID) != len(cell_lons_1d):
            if len(aCellID) > len(cell_lons_1d):
                aCellID = aCellID[:len(cell_lons_1d)]
            else:
                missing_count = len(cell_lons_1d) - len(aCellID)
                aCellID.extend(range(len(aCellID), len(aCellID) + missing_count))

        aCellID = np.array(aCellID)

        # Calculate area statistics
        if aArea_list:
            area_array = np.array(aArea_list)
            valid_areas = area_array[area_array > 0]
            if len(valid_areas) > 0:
                dArea_min = float(np.min(valid_areas))
                dArea_max = float(np.max(valid_areas))
                dArea_mean = float(np.mean(valid_areas))
            else:
                dArea_min = 0.0
                dArea_max = 0.0
                dArea_mean = 0.0

        # Return comprehensive mesh topology information
        mesh_info = {
            'vertices_longitude': xv,
            'vertices_latitude': yv,
            'connectivity': connectivity,
            'cell_centroids_longitude': cell_lons_1d,
            'cell_centroids_latitude': cell_lats_1d,
            'cell_ids': aCellID,
            'area_min': dArea_min,
            'area_max': dArea_max,
            'area_mean': dArea_mean,
            'max_vertices_per_cell': nVertex_max,
            'num_cells': nFeatures,
            'num_polygns': len(cell_lons_1d),
            'num_vertices': len(xv),
            'success': True
        }

        return mesh_info

    except Exception as e:
        logger.error(f'Unexpected error in rebuild_mesh_topology: {str(e)}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return None


def _extract_target_mesh_data(sFilename: str, sVariable: str, iFlag_verbose_in: bool = False) -> Optional[Tuple[np.ndarray, int]]:
    """
    Extract variable data from target mesh file.

    Args:
        sFilename: Path to target mesh file
        sVariable: Variable name to extract
        iFlag_verbose_in: Enable verbose logging

    Returns:
        Tuple of (data_array, feature_count) or None if failed
    """
    try:
        pDataset = ogr.Open(sFilename, 0)
        if pDataset is None:
            logger.error(f'Failed to open target mesh file: {sFilename}')
            return None

        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            logger.error('Failed to get layer from target mesh dataset')
            pDataset = None
            return None

        pLayerDefn = pLayer.GetLayerDefn()
        if pLayerDefn is None:
            logger.error('Failed to get layer definition from target mesh')
            pDataset = None
            return None

        iFieldCount = pLayerDefn.GetFieldCount()
        nFeatures = pLayer.GetFeatureCount()

        if iFlag_verbose_in:
            logger.info(f'Target mesh contains {nFeatures} features with {iFieldCount} fields')

        # Check if variable field exists
        aField_names = [pLayerDefn.GetFieldDefn(i).GetName() for i in range(iFieldCount)]
        if sVariable not in aField_names:
            logger.error(f'Variable "{sVariable}" not found in target mesh')
            logger.error(f'Available fields: {", ".join(aField_names)}')
            pDataset = None
            return None

        if iFlag_verbose_in:
            logger.info(f'Extracting variable: {sVariable}')

        # Extract variable data from features
        aData_list = []
        pLayer.ResetReading()
        pFeature = pLayer.GetNextFeature()
        iFeature_count = 0

        while pFeature is not None:
            pGeometry = pFeature.GetGeometryRef()
            if pGeometry is not None:
                sGeometry_type = pGeometry.GetGeometryName()
                if sGeometry_type == 'POLYGON':
                    try:
                        dField_value = pFeature.GetField(sVariable)
                        aData_list.append(dField_value if dField_value is not None else np.nan)
                    except Exception as e:
                        logger.warning(f'Error reading field {sVariable} from feature {iFeature_count}: {e}')
                        aData_list.append(np.nan)

                elif sGeometry_type == 'MULTIPOLYGON':
                    try:
                        dField_value = pFeature.GetField(sVariable)
                        dData_value = dField_value if dField_value is not None else np.nan

                        nGeometryParts = pGeometry.GetGeometryCount()
                        for iPart in range(nGeometryParts):
                            pPolygon_part = pGeometry.GetGeometryRef(iPart)
                            if pPolygon_part is not None and pPolygon_part.IsValid():
                                aData_list.append(dData_value)
                            else:
                                aData_list.append(np.nan)

                    except Exception as e:
                        logger.warning(f'Error reading field {sVariable} from multipolygon feature {iFeature_count}: {e}')
                        aData_list.append(np.nan)

            iFeature_count += 1
            pFeature = pLayer.GetNextFeature()

        pDataset = None

        if not aData_list:
            logger.error('No data extracted from target mesh')
            return None

        aData = np.array(aData_list, dtype=np.float64)
        return aData, nFeatures

    except Exception as e:
        logger.error(f'Error extracting target mesh data: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return None


def _handle_visualization_output(pPlotter, sFilename: Optional[str], iFlag_verbose_in: bool = False, window_size: Optional[list] = None) -> bool:
    """
    Handle visualization output (save file or show interactive).

    Args:
        pPlotter: GeoVista plotter instance
        sFilename: Output filename or None for interactive
        iFlag_verbose_in: Enable verbose logging

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        if window_size is None:
            window_size = [1200, 800]

        if sFilename is not None:
            # Save screenshot
            sExt = os.path.splitext(sFilename)[1].lower()
            if sExt in VALID_IMAGE_FORMATS:
                pPlotter.screenshot(sFilename, window_size=window_size)
            else:
                pPlotter.save_graphic(sFilename)  # pdf format

            if iFlag_verbose_in:
                logger.info(f'✓ Visualization saved to: {sFilename}')

                if os.path.exists(sFilename):
                    iFile_size = os.path.getsize(sFilename)
                    logger.info(f'  File size: {iFile_size / 1024:.1f} KB')
                else:
                    logger.warning(f'Screenshot command executed but file not found: {sFilename}')

            pPlotter.close()
            return True
        else:
            # Interactive display
            if iFlag_verbose_in:
                logger.info('Opening interactive visualization window...')
            pPlotter.show()
            return True

    except Exception as e:
        logger.error(f'Failed to handle visualization output: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        try:
            pPlotter.close()
        except Exception:
            pass
        return False


def visualize_unified_mesh_dataset(sFilename_mpas_mesh_in: str,
                      sFilename_uraster_in: str,
                        sFilename_flow_directions_in: str,
                        iFlag_full_mesh: int = 1,
                         iFlag_elevation: int = 0,
                         iFlag_flow_directions: int = 0,
                          sFilename_out: Optional[str] = None,
                          dLongitude_focus_in: Optional[float] = 0.0,
                          dLatitude_focus_in: Optional[float] = 0.0,
                          dZoom_factor: float = 0.7,
                          iFlag_show_coastlines: bool = True,
                          iFlag_show_graticule: bool = False,
                          sCoastline_color: str = 'black',
                          dCoastline_width: float = 1.0,
                          sEdge_color: str = 'black',
                          sElevation_colormap: str = 'terrain',
                          dElevation_opacity: float = 0.85,
                          sPolyline_width_attribute: Optional[str] = None,
                          dPolyline_width_min: float = 0.5,
                          dPolyline_width_max: float = 3.0,
                         iFlag_verbose_in: Optional[bool] = False,
                         window_size: Tuple[int, int] = (800, 600)
                         ) -> bool:

    """
    Visualize the source mesh topology using GeoVista 3D globe rendering with pyearthviz3d.

    Creates an interactive or saved 3D visualization of the unstructured mesh
    with proper geographic context including coastlines and coordinate grid.

    Args:
        sFilename_mpas_mesh_in: Path to MPAS mesh NetCDF file
        sFilename_uraster_in: Path to elevation raster file (GeoJSON/Shapefile)
        sFilename_flow_directions_in: Path to flow directions vector file
        sFilename_out: Output screenshot file path. If None, displays interactive viewer.
        dLongitude_focus_in: Camera focal point longitude in degrees (-180 to 180).
        dLatitude_focus_in: Camera focal point latitude in degrees (-90 to 90).
        dZoom_factor: Camera zoom level. Higher values zoom in. Default is 0.7.
        iFlag_show_coastlines: Show coastline overlay. Default is True.
        iFlag_show_graticule: Show coordinate grid with labels. Default is False.
        sCoastline_color: Color for coastlines. Default is 'black'.
        dCoastline_width: Line width for coastlines. Default is 1.0.
        sEdge_color: Color of mesh edges. Default is 'black'.
        sElevation_colormap: Colormap for elevation data. Default is 'terrain'.
        dElevation_opacity: Opacity of elevation mesh (0.0 to 1.0). Default is 0.85.
        sPolyline_width_attribute: Attribute name for variable polyline widths.
        dPolyline_width_min: Minimum line width for attribute-based scaling. Default is 0.5.
        dPolyline_width_max: Maximum line width for attribute-based scaling. Default is 3.0.
        iFlag_verbose_in: If True, print detailed progress messages. Default is False.
        iFlag_full_mesh: If 1, display the full MPAS mesh wireframe. Default is 1.
        iFlag_elevation: If 1, display elevation data as colored mesh. Default is 0.
        iFlag_flow_directions: If 1, display flow direction polylines. Default is 0.
        window_size: Window size for rendering (width, height). Default is (800, 600).

    Returns:
        True if visualization successful, False otherwise
    """
    import geovista as gv

    # Define a base window size and font size
    base_window_size = 1200
    base_font_size = 8
    base_line_width = 0.1

    # Calculate scaling factor
    scaling_factor = max(window_size) / base_window_size

    # Adjust font size
    scaled_font_size = int(base_font_size * scaling_factor)
    scaled_line_width = base_line_width * scaling_factor

    # Create configuration object
    pConfig = VisualizationConfig(
        longitude_focus=dLongitude_focus_in,
        latitude_focus=dLatitude_focus_in,
        zoom_factor=dZoom_factor,
        show_coastlines=iFlag_show_coastlines,
        show_graticule=iFlag_show_graticule,
        coastline_color=sCoastline_color,
        coastline_width=dCoastline_width,
        window_size=window_size,
        verbose=iFlag_verbose_in
    )

    try:
        # Setup plotter using pyearthviz3d
        pPlotter = PlotterManager.setup_geovista_plotter(
            off_screen=(sFilename_out is not None),
            verbose=pConfig.verbose,
            window_size=pConfig.window_size,
            use_xvfb=pConfig.use_xvfb,
            force_xvfb=pConfig.force_xvfb
        )

        if iFlag_full_mesh:
            # Load MPAS mesh
            mesh_data = load_mpas_mesh_optimized(sFilename_mpas_mesh_in, pConfig.verbose)
            if mesh_data is None:
                logger.error(f'Failed to load MPAS mesh from {sFilename_mpas_mesh_in}')
                return False

            aVertex_longititude = mesh_data['vertices_longitude']
            aVertex_latitude = mesh_data['vertices_latitude']
            aConnectivity = mesh_data['connectivity']

            if pConfig.verbose:
                logger.info('Creating mesh visualization...')
                logger.info(f'  - Vertices: {len(aVertex_longititude)}')
                logger.info(f'  - Connectivity shape: {aConnectivity.shape}')

            # Process connectivity
            connectivity_masked = optimize_connectivity_processing(
                aConnectivity, len(aVertex_longititude), pConfig.verbose
            )

            # Transform to GeoVista mesh
            pMesh_full = gv.Transform.from_unstructured(
                aVertex_longititude,
                aVertex_latitude,
                connectivity=connectivity_masked,
                crs=CRS
            )

            # Create valid cell indices
            aValid_cell_indices = np.arange(connectivity_masked.shape[0])
            if hasattr(connectivity_masked, 'mask'):
                valid_vertices_per_cell = np.sum(~connectivity_masked.mask, axis=1)
            else:
                valid_vertices_per_cell = np.sum(connectivity_masked >= 0, axis=1)
            aValid_cell_indices = aValid_cell_indices[valid_vertices_per_cell >= 3]

            pPlotter.add_base_layer(mesh=pMesh_full, texture=gv.natural_earth_hypsometric())

            mesh_success = add_mesh_to_plotter(
                plotter=pPlotter,
                mesh=pMesh_full,
                style='wireframe',
                show_edges=True,
                edge_color=sEdge_color,
                valid_indices=aValid_cell_indices,
                colormap=None,
            )

        if iFlag_elevation == 1:
            # Process elevation data
            elevation_data = process_elevation_data_optimized(
                sFilename_uraster_in,
                'mean',
                pConfig.verbose
            )

            if elevation_data is None:
                logger.warning(f'Failed to process elevation data from: {sFilename_uraster_in}')
            else:
                pMesh_elevation = gv.Transform.from_unstructured(
                    elevation_data['vertices_longitude'],
                    elevation_data['vertices_latitude'],
                    connectivity=elevation_data['connectivity_masked'],
                    crs=CRS
                )

                sScalars = 'Elevation (m)'
                pMesh_elevation.cell_data[sScalars] = elevation_data['elevation_data']

                # Create valid cell indices for elevation mesh
                aValid_cell_indices_elevation = np.arange(elevation_data['connectivity_masked'].shape[0])
                if hasattr(elevation_data['connectivity_masked'], 'mask'):
                    valid_vertices_per_cell = np.sum(~elevation_data['connectivity_masked'].mask, axis=1)
                else:
                    valid_vertices_per_cell = np.sum(elevation_data['connectivity_masked'] >= 0, axis=1)
                aValid_cell_indices_elevation = aValid_cell_indices_elevation[valid_vertices_per_cell >= 3]

                mesh_success = add_mesh_to_plotter(
                    plotter=pPlotter,
                    mesh=pMesh_elevation,
                    scalar_name=sScalars,
                    valid_indices=aValid_cell_indices_elevation,
                    scalar_config=ScalarBarConfig(
                        title=sScalars,
                        title_font_size=scaled_font_size,
                        label_font_size=scaled_font_size,
                        position_x=0.85,
                        position_y=0.1,
                        orientation='vertical',
                        n_labels=5,
                    ),
                    colormap=sElevation_colormap,
                    unit="m",
                    opacity=dElevation_opacity,
                    show_edges=True,
                    edge_color=sEdge_color
                )

                if pConfig.verbose:
                    logger.info('✓ Successfully added elevation mesh')

        if iFlag_flow_directions == 1:
            import pyvista as pv

            # Process polylines
            polyline_data = process_polylines_optimized(
                sFilename_flow_directions_in,
                scaled_line_width * 1.5,
                pConfig.verbose,
                width_attribute=sPolyline_width_attribute,
                width_min=dPolyline_width_min,
                width_max=dPolyline_width_max
            )

            if polyline_data is None:
                logger.error('Failed to process polylines')
                return False

            polydata = pv.PolyData(polyline_data['points'], lines=polyline_data['lines'])

            sColor_polyline_in = 'blue'
            pPlotter.add_mesh(
                polydata,
                color=sColor_polyline_in,
                line_width=polyline_data['line_width'],
                name="combined_polylines",
                label='Flow directions'
            )

        # Configure camera using pyearthviz3d
        camera_success = configure_camera_enhanced(
            plotter=pPlotter,
            config=pConfig,
            use_enhanced_controller=True
        )

        # Add geographic context using pyearthviz3d
        context_results = add_geographic_context_enhanced(
            plotter=pPlotter,
            config=pConfig
        )

        # Output or display
        return _handle_visualization_output(pPlotter, sFilename_out, pConfig.verbose, window_size=window_size)

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False


if __name__ == "__main__":
    sFilename_mpas_mesh_in = r'/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20251122001/jigsaw/out/base_mesh.nc'
    sFilename_uraster_in = r'/qfs/people/liao313/workspace/python/uraster/data/example_12/output/mpas_uraster.geojson'
    sFilename_flow_directions_in = r'/compyfs/liao313/04model/pyhexwatershed/global/pyhexwatershed20251208002/hexwatershed/mpas_flow_direction.geojson'
    sFilename_out = r'/qfs/people/liao313/workspace/python/unified_land_river_mesh/figures/global_mesh_visualization.jpg'

    visualize_unified_mesh_dataset(sFilename_mpas_mesh_in,
                      sFilename_uraster_in,
                      sFilename_flow_directions_in,
                      sFilename_out=sFilename_out,
                         iFlag_full_mesh=1,
                         iFlag_elevation=1,
                         iFlag_flow_directions=1,
                      dLongitude_focus_in=-97,
                      dLatitude_focus_in=43,
                      dZoom_factor=0.7,
                      sPolyline_width_attribute=None,
                      window_size=(8000, 7000))

    pass
