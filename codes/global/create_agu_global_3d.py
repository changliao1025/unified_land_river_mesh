"""
Visualization module for uraster class.

This module contains all visualization-related methods that were moved from the main uraster class
to reduce the size of the main uraster.py file and improve code organization.

Features:
- 3D mesh visualization using GeoVista
- Interactive and static rendering modes
- Animation support with rotation and camera movement
- Comprehensive error handling and validation
- Support for multiple output formats
"""

import os
import logging
import traceback
import math
import time
from functools import wraps
import matplotlib.pyplot as plt
from typing import Optional, List, Tuple, Union, Dict, Any
import numpy as np
from osgeo import gdal, ogr, osr
gdal.UseExceptions()

from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import extract_unique_vertices_and_connectivity
from pyearth.visual.animate.animate_polyline_file_on_sphere import animate_polyline_file_on_sphere
from geovista.geodesic import line as gv_line
# Set up logging
logger = logging.getLogger(__name__)
CRS = "EPSG:4326"

# Constants for visualization
DEFAULT_EARTH_RADIUS = 1.0
DEFAULT_CAMERA_DISTANCE_MULTIPLIER = 3.0
DEFAULT_ZOOM_FACTOR = 0.7
VALID_ANIMATION_FORMATS = ['mp4', 'gif', 'avi']
VALID_IMAGE_FORMATS = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']
COORDINATE_BOUNDS = {'longitude': (-180, 180), 'latitude': (-90, 90)}

# Performance monitoring utilities
def profile_function(func):
    """Decorator to profile function execution time."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        logger.info(f"⏱️ {func.__name__} executed in {end_time - start_time:.3f} seconds")
        return result
    return wrapper

class PerformanceMonitor:
    """Context manager for monitoring performance of code blocks."""

    def __init__(self, operation_name: str, verbose: bool = False):
        self.operation_name = operation_name
        self.verbose = verbose
        self.start_time = None

    def __enter__(self):
        self.start_time = time.time()
        if self.verbose:
            logger.info(f"🚀 Starting {self.operation_name}...")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        end_time = time.time()
        duration = end_time - self.start_time
        if self.verbose:
            logger.info(f"✅ {self.operation_name} completed in {duration:.3f} seconds")

@profile_function
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

        with PerformanceMonitor("MPAS mesh loading", verbose):
            with nc.Dataset(filename, 'r') as dataset:
                # Read coordinates with optimized chunking
                with PerformanceMonitor("Reading vertex coordinates", verbose):
                    lon_vertex = dataset.variables['lonVertex'][:]
                    lat_vertex = dataset.variables['latVertex'][:]

                    # Convert to degrees efficiently using vectorized operations
                    lon_vertex = np.rad2deg(lon_vertex)
                    lat_vertex = np.rad2deg(lat_vertex)

                    # Fix longitude range [0, 360] to [-180, 180] efficiently
                    lon_vertex = np.where(lon_vertex > 180, lon_vertex - 360, lon_vertex)

                # Read connectivity with validation
                with PerformanceMonitor("Reading connectivity", verbose):
                    connectivity = dataset.variables['verticesOnCell'][:]
                    cell_ids = dataset.variables['indexToCellID'][:]

                if verbose:
                    logger.info(f"Loaded mesh: {len(lon_vertex)} vertices, {connectivity.shape[0]} cells")

        return {
            'vertices_longitude': lon_vertex,
            'vertices_latitude': lat_vertex,
            'connectivity': connectivity,
            'cell_ids': cell_ids
        }

    except Exception as e:
        logger.error(f"Failed to load MPAS mesh from {filename}: {e}")
        return None

@profile_function
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
    with PerformanceMonitor("Connectivity optimization", verbose):
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

@profile_function
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
        with PerformanceMonitor("Polyline processing", verbose):
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

            # Setup coordinate transformation if needed
            spatial_ref = layer.GetSpatialRef()
            wgs84_srs = osr.SpatialReference()
            wgs84_srs.ImportFromEPSG(4326)

            transform = None
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
            with PerformanceMonitor("Geometry extraction", verbose):
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
_polyline_cache = {}

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

@profile_function
def process_elevation_data_optimized(filename: str, variable: str,  verbose: bool = False) -> Optional[Dict[str, Any]]:
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
        with PerformanceMonitor("Elevation data processing", verbose):
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

            # Cache the result
            #cache_mesh_data(filename, result, cache_key)

            return result

    except Exception as e:
        logger.error(f'Error processing elevation data from {filename}: {e}')
        return None

class VisualizationConfig:
    """Configuration class for visualization parameters."""

    def __init__(self,
                 longitude_focus: float = 0.0,
                 latitude_focus: float = 0.0,
                 zoom_factor: float = DEFAULT_ZOOM_FACTOR,
                 show_coastlines: bool = True,
                 show_graticule: bool = True,
                 colormap: str = 'terrain',
                 coastline_color: str = 'black',
                 coastline_width: float = 1.0,
                 verbose: bool = False):
        self.longitude_focus = self._validate_longitude(longitude_focus)
        self.latitude_focus = self._validate_latitude(latitude_focus)
        self.zoom_factor = self._validate_zoom_factor(zoom_factor)
        self.show_coastlines = show_coastlines
        self.show_graticule = show_graticule
        self.colormap = colormap
        self.coastline_color = coastline_color
        self.coastline_width = coastline_width
        self.verbose = verbose

    def _validate_longitude(self, lon: float) -> float:
        """Validate and clamp longitude to valid range."""
        if not (-180 <= lon <= 180):
            logger.warning(f'Longitude {lon} out of range [-180, 180], clamping')
            return np.clip(lon, -180, 180)
        return lon

    def _validate_latitude(self, lat: float) -> float:
        """Validate and clamp latitude to valid range."""
        if not (-90 <= lat <= 90):
            logger.warning(f'Latitude {lat} out of range [-90, 90], clamping')
            return np.clip(lat, -90, 90)
        return lat

    def _validate_zoom_factor(self, zoom: float) -> float:
        """Validate zoom factor."""
        if zoom <= 0:
            logger.warning(f'Invalid zoom factor {zoom}, using default {DEFAULT_ZOOM_FACTOR}')
            return DEFAULT_ZOOM_FACTOR
        return zoom

class CameraController:
    """Handles camera positioning and movement calculations."""

    @staticmethod
    def calculate_camera_position(longitude: float, latitude: float,
                                zoom_factor: float = DEFAULT_ZOOM_FACTOR) -> Tuple[List[float], List[float]]:
        """
        Calculate camera position and focal point from geographic coordinates.

        Args:
            longitude: Longitude in degrees
            latitude: Latitude in degrees
            zoom_factor: Camera zoom level

        Returns:
            Tuple of (focal_point, camera_position) as [x, y, z] lists
        """
        # Convert to radians
        lon_rad = math.radians(longitude)
        lat_rad = math.radians(latitude)

        # Calculate positions
        earth_radius = DEFAULT_EARTH_RADIUS
        camera_distance = earth_radius * DEFAULT_CAMERA_DISTANCE_MULTIPLIER

        # Focal point on Earth surface
        x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
        y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
        z_focal = earth_radius * math.sin(lat_rad)

        # Camera position away from Earth
        x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
        y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
        z_camera = camera_distance * math.sin(lat_rad)

        focal_point = [x_focal, y_focal, z_focal]
        camera_position = [x_camera, y_camera, z_camera]

        return focal_point, camera_position

    @staticmethod
    def validate_camera_setup(focal_point: List[float], camera_position: List[float]) -> bool:
        """Validate camera setup to ensure proper positioning."""
        distance = math.sqrt(
            sum((c - f)**2 for c, f in zip(camera_position, focal_point))
        )
        return distance >= 0.1  # Minimum distance threshold

class AnimationConfig:
    """Configuration class for animation parameters."""

    def __init__(self,
                 frames: int = 36,
                 speed: float = 1.0,
                 format: str = 'mp4',
                 amplitude_deg: float = 20.0,
                 cycles: float = 1.0,
                 phase: float = 0.0):
        # Convert to int if string and validate
        try:
            frames_int = int(frames) if isinstance(frames, str) else frames
            self.frames = max(1, frames_int)  # Ensure at least 1 frame
        except (ValueError, TypeError):
            logger.warning(f'Invalid frames value {frames}, using default 36')
            self.frames = 36

        # Convert to float if string and validate
        try:
            speed_float = float(speed) if isinstance(speed, str) else speed
            self.speed = max(0.1, speed_float)  # Ensure positive speed
        except (ValueError, TypeError):
            logger.warning(f'Invalid speed value {speed}, using default 1.0')
            self.speed = 1.0

        self.format = format.lower()
        self.amplitude_deg = amplitude_deg
        self.cycles = cycles
        self.phase = phase

        # Validate format
        if self.format not in VALID_ANIMATION_FORMATS:
            logger.warning(f'Invalid animation format {format}, using mp4')
            self.format = 'mp4'

def rebuild_mesh_topology(sFilename_mesh_in: str, iFlag_verbose_in: bool = False, sField_unique_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """
    Rebuild mesh topology from source mesh file by extracting vertices,
    connectivity, and centroids for unstructured mesh processing.

    Args:
        sFilename_mesh_in (str): Path to the source mesh file (GeoJSON, Shapefile, etc.)
        iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.
        sField_unique_id (str, optional): Field name for unique cell IDs. If None, uses first field or feature index.
            Note: Field is always treated as integer type since setup_mesh_cellid() enforces this.

    Returns:
        dict: Comprehensive mesh topology information with keys:
            - 'vertices_longitude': np.ndarray of unique vertex longitudes
            - 'vertices_latitude': np.ndarray of unique vertex latitudes
            - 'connectivity': np.ndarray connectivity matrix
            - 'cell_centroids_longitude': np.ndarray of cell centroid longitudes
            - 'cell_centroids_latitude': np.ndarray of cell centroid latitudes
            - 'cell_ids': np.ndarray of cell IDs
            - 'area_min': float minimum cell area
            - 'area_max': float maximum cell area
            - 'area_mean': float mean cell area
            - 'max_vertices_per_cell': int maximum vertices per cell
            - 'num_cells': int total number of cells
            - 'num_vertices': int total number of unique vertices
            - 'success': bool whether processing was successful
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
        aCellID = []  # Will be populated dynamically as features are processed
        if iFlag_verbose_in:
            logger.info(
                f'Processing {nFeatures} features with {iFieldCount} fields')
        # Get the first field name (assuming it contains the data variable)
        if sField_unique_id is None:
            sVariable = pLayerDefn.GetFieldDefn(
                0).GetName() if iFieldCount > 0 else None
        else:
            sVariable = sField_unique_id
        # Initialize lists for storing geometry data
        lons_list = []
        lats_list = []
        aArea_list = []
        # Process features with enhanced error handling
        pLayer.ResetReading()
        iFeature_index = 0
        invalid_geometry_count = 0
        iCount_multipolygon_cells = 0
        for pFeature in pLayer:
            if pFeature is None:
                iFeature_index += 1
                continue
            pGeometry = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry.GetGeometryName()
            if sGeometry_type == 'POLYGON':
                try:
                    # Get coordinates of the polygon (validation already done by check_geometry_validity)
                    aCoord = get_geometry_coordinates(pGeometry)
                    if aCoord is not None and len(aCoord) >= 3:
                        lons = aCoord[:, 0]
                        lats = aCoord[:, 1]
                        lons_list.append(lons)
                        lats_list.append(lats)
                        # Calculate polygon area
                        try:
                            dArea = calculate_polygon_area(lons, lats)
                            aArea_list.append(dArea)
                        except Exception as area_error:
                            logger.warning(
                                f'Could not calculate area for feature {iFeature_index}: {area_error}')
                            aArea_list.append(0.0)
                        # Get field data (always integer since setup_mesh_cellid enforces it)
                        if sVariable:
                            try:
                                field_value = pFeature.GetFieldAsInteger(
                                    sVariable)
                                aCellID.append(int(field_value))
                            except (ValueError, TypeError, AttributeError) as e:
                                logger.warning(
                                    f'Could not read integer field value for feature {iFeature_index}: {e}')
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
                    # Process multipolygon by extracting each constituent polygon
                    if iFlag_verbose_in:
                        logger.info(f'Processing multipolygon feature {iFeature_index}')
                    multipolygon_processed = False
                    iCount_multipolygon_cells += 1

                    # Get field data for the multipolygon feature (same for all parts)
                    current_cellid = None
                    if sVariable:
                        try:
                            field_value = pFeature.GetFieldAsInteger(sVariable)
                            current_cellid = int(field_value)
                        except (ValueError, TypeError, AttributeError) as e:
                            logger.warning(
                                f'Could not read integer field value for multipolygon feature {iFeature_index}: {e}')
                            current_cellid = len(aCellID)
                    else:
                        current_cellid = len(aCellID)

                    for iPart in range(pGeometry.GetGeometryCount()):
                        pPolygon_part = pGeometry.GetGeometryRef(iPart)
                        if pPolygon_part is None:
                            logger.warning(f'Multipolygon part {iPart} is None in feature {iFeature_index}')
                            continue
                        # Get coordinates of the polygon part (validation already done by check_geometry_validity)
                        aCoord_part = get_geometry_coordinates(pPolygon_part)
                        if aCoord_part is not None and len(aCoord_part) >= 3:
                            lons_part = aCoord_part[:, 0]
                            lats_part = aCoord_part[:, 1]
                            lons_list.append(lons_part)
                            lats_list.append(lats_part)
                            # Calculate polygon area for this part
                            try:
                                dArea_part = calculate_polygon_area(
                                    lons_part, lats_part)
                                aArea_list.append(dArea_part)
                            except Exception as area_error:
                                logger.warning(
                                    f'Could not calculate area for multipolygon part {iPart} in feature {iFeature_index}: {area_error}')
                                aArea_list.append(0.0)
                            # Add cell ID for this part
                            aCellID.append(current_cellid)
                            multipolygon_processed = True
                        else:
                            logger.warning(
                                f'Failed to extract coordinates from multipolygon part {iPart} in feature {iFeature_index}')
                    if not multipolygon_processed:
                        logger.warning(
                            f'No valid parts found in multipolygon feature {iFeature_index}')
                        invalid_geometry_count += 1
                except Exception as e:
                    logger.warning(
                        f'Error processing multipolygon feature {iFeature_index}: {str(e)}')
                    invalid_geometry_count += 1
            elif sGeometry_type in ['POINT', 'LINESTRING']:
                logger.warning(
                    f'Geometry type {sGeometry_type} not supported in feature {iFeature_index}, skipping')
                invalid_geometry_count += 1
            else:
                logger.warning(
                    f'Unknown geometry type {sGeometry_type} in feature {iFeature_index}, skipping')
                invalid_geometry_count += 1

            iFeature_index += 1
        # Report processing statistics
        valid_mesh_cells = len(lons_list)
        if iFlag_verbose_in:
            logger.info(f'Feature processing summary:')
            logger.info(f'  - Total input features: {iFeature_index}')
            logger.info(f'  - Valid mesh cells created: {valid_mesh_cells}')
            logger.info(
                f'  - Invalid/skipped features: {invalid_geometry_count}')
            logger.info(
                f'  - Success rate: {((iFeature_index-invalid_geometry_count)/iFeature_index*100):.1f}%' if iFeature_index > 0 else '  - Success rate: 0%')
            # Report multipolygon handling statistics
            multipolygon_cells = valid_mesh_cells - \
                (iFeature_index - invalid_geometry_count)
            if multipolygon_cells > 0:
                logger.info(
                    f'  - Additional cells from multipolygons: {multipolygon_cells}')
                logger.info(
                    f'  - Total mesh cells (including multipolygon parts): {valid_mesh_cells}')
        # Clean up dataset
        pDataset = None
        if not lons_list:
            logger.error('No valid polygon features found in mesh file')
            return None
        if iFlag_verbose_in:
            logger.info(
                f'Successfully processed {len(lons_list)} polygon features')
        # Calculate maximum vertices and pad coordinates efficiently
        try:
            if not lons_list:
                logger.error('No coordinate data found')
                return None
            max_vertices = max(len(coord) for coord in lons_list)
            if max_vertices == 0:
                logger.error('No vertices found in any polygon')
                return None
            nVertex_max = max_vertices
            if iFlag_verbose_in:
                logger.info(f'Maximum vertices per polygon: {max_vertices}')
            # Pre-allocate arrays for better memory efficiency
            num_polygons = len(lons_list)
            lons_padded = np.full(
                (num_polygons, max_vertices), np.nan, dtype=np.float64)
            lats_padded = np.full(
                (num_polygons, max_vertices), np.nan, dtype=np.float64)
            # Fill padded arrays efficiently
            for i, (lon_coords, lat_coords) in enumerate(zip(lons_list, lats_list)):
                # Ensure coordinates are numpy arrays with proper dtype
                lon_coords = np.asarray(lon_coords, dtype=np.float64)
                lat_coords = np.asarray(lat_coords, dtype=np.float64)
                # Validate coordinate data
                if len(lon_coords) != len(lat_coords):
                    logger.warning(
                        f'Coordinate length mismatch in polygon {i}: lon={len(lon_coords)}, lat={len(lat_coords)}')
                    min_len = min(len(lon_coords), len(lat_coords))
                    lon_coords = lon_coords[:min_len]
                    lat_coords = lat_coords[:min_len]
                # Check for valid coordinate values
                if not (np.all(np.isfinite(lon_coords)) and np.all(np.isfinite(lat_coords))):
                    logger.warning(f'Invalid coordinates found in polygon {i}')
                    # Remove invalid coordinates
                    valid_mask = np.isfinite(
                        lon_coords) & np.isfinite(lat_coords)
                    lon_coords = lon_coords[valid_mask]
                    lat_coords = lat_coords[valid_mask]
                coord_len = len(lon_coords)
                if coord_len > 0:
                    lons_padded[i, :coord_len] = lon_coords
                    lats_padded[i, :coord_len] = lat_coords
                else:
                    logger.warning(
                        f'No valid coordinates remaining for polygon {i}')
            # Convert to the expected format for backward compatibility
            lons = lons_padded
            lats = lats_padded
        except Exception as e:
            logger.error(f'Error during coordinate padding: {str(e)}')
            logger.error(f'Traceback: {traceback.format_exc()}')
            return None
        # Calculate centroids efficiently using vectorized operations
        try:
            cell_lons_1d = []
            cell_lats_1d = []
            # Pre-allocate arrays for better performance
            cell_lons_1d = np.zeros(len(lons_list), dtype=np.float64)
            cell_lats_1d = np.zeros(len(lats_list), dtype=np.float64)
            for i in range(len(lons_list)):
                # Calculate centroid of each cell (ignoring NaN values)
                valid_mask = ~np.isnan(lons[i])
                if np.any(valid_mask):
                    valid_lons = lons[i][valid_mask]
                    valid_lats = lats[i][valid_mask]
                    # Use vectorized operations for better performance
                    centroid_lon = np.mean(valid_lons)
                    centroid_lat = np.mean(valid_lats)
                    # Validate centroid coordinates
                    if np.isfinite(centroid_lon) and np.isfinite(centroid_lat):
                        cell_lons_1d[i] = centroid_lon
                        cell_lats_1d[i] = centroid_lat
                    else:
                        logger.warning(
                            f'Invalid centroid calculated for cell {i}: lon={centroid_lon}, lat={centroid_lat}')
                        # Use geometric center of bounding box as fallback
                        if len(valid_lons) > 0 and len(valid_lats) > 0:
                            cell_lons_1d[i] = (
                                np.min(valid_lons) + np.max(valid_lons)) / 2.0
                            cell_lats_1d[i] = (
                                np.min(valid_lats) + np.max(valid_lats)) / 2.0
                        else:
                            cell_lons_1d[i] = 0.0
                            cell_lats_1d[i] = 0.0
                else:
                    logger.warning(f'No valid coordinates found for cell {i}')
                    cell_lons_1d[i] = 0.0
                    cell_lats_1d[i] = 0.0
            if iFlag_verbose_in:
                logger.info(
                    f'Calculated centroids for {len(cell_lons_1d)} cells')
            # Validate centroid ranges
            lon_range = (np.min(cell_lons_1d), np.max(cell_lons_1d))
            lat_range = (np.min(cell_lats_1d), np.max(cell_lats_1d))
            if not (-180 <= lon_range[0] <= 180 and -180 <= lon_range[1] <= 180):
                logger.warning(
                    f'Longitude centroids outside valid range: {lon_range}')
            if not (-90 <= lat_range[0] <= 90 and -90 <= lat_range[1] <= 90):
                logger.warning(
                    f'Latitude centroids outside valid range: {lat_range}')
        except Exception as e:
            logger.error(f'Error during centroid calculation: {str(e)}')
            return None
        # Extract unique vertices and connectivity
        try:
            if iFlag_verbose_in:
                logger.info('Extracting unique vertices and connectivity...')
            xv, yv, connectivity, vertex_to_index = extract_unique_vertices_and_connectivity(
                lons_list, lats_list
            )
            if xv is None or yv is None or connectivity is None:
                logger.error(
                    'Failed to extract unique vertices and connectivity')
                return None
            if iFlag_verbose_in:
                logger.info(f'Extracted {len(xv)} unique vertices')
                logger.info(
                    f'Created connectivity matrix with shape: {connectivity.shape}')
        except Exception as e:
            logger.error(
                f'Error during vertex/connectivity extraction: {str(e)}')
            return None
        # Store results in class attributes
        aVertex_longititude = xv
        aVertex_latitude = yv
        aCenter_longititude = cell_lons_1d
        aCenter_latitude = cell_lats_1d
        aConnectivity = connectivity
        # Ensure aCellID matches the number of valid mesh cells
        if len(aCellID) != len(cell_lons_1d):
            logger.warning(
                f"aCellID length ({len(aCellID)}) doesn't match mesh cells ({len(cell_lons_1d)})")
            if len(aCellID) > len(cell_lons_1d):
                # Truncate aCellID to match mesh cells
                logger.warning("Truncating aCellID to match mesh cell count")
                aCellID = aCellID[:len(cell_lons_1d)]
            else:
                # Extend aCellID with sequential indices
                logger.warning(
                    "Extending aCellID with sequential indices to match mesh cell count")
                missing_count = len(cell_lons_1d) - len(aCellID)
                aCellID.extend(
                    range(len(aCellID), len(aCellID) + missing_count))
        aCellID = np.array(aCellID)
        if iFlag_verbose_in:
            logger.info(f'Final aCellID array length: {len(aCellID)}')
        # Calculate and store area statistics
        if aArea_list:
            area_array = np.array(aArea_list)
            # Exclude zero areas from statistics
            valid_areas = area_array[area_array > 0]
            if len(valid_areas) > 0:
                dArea_min = float(np.min(valid_areas))
                dArea_max = float(np.max(valid_areas))
                dArea_mean = float(np.mean(valid_areas))
                dArea_max = float(np.max(valid_areas))
                dArea_min = float(np.min(valid_areas))
                if iFlag_verbose_in:
                    logger.info(f'Mesh area statistics:')
                    logger.info(f'  - Min area: {dArea_min:.6f}')
                    logger.info(f'  - Max area: {dArea_max:.6f}')
                    logger.info(f'  - Mean area: {dArea_mean:.6f}')
            else:
                logger.warning('No valid polygon areas calculated')
                dArea_min = 0.0
                dArea_max = 0.0
                dArea_mean = 0.0
        # Enhanced validation of final results
        validation_passed = True
        if len(aVertex_longititude) == 0:
            logger.error('No unique vertices extracted')
            validation_passed = False
        if len(aCenter_longititude) != len(lons_list):
            logger.error(
                f'Centroid count mismatch: expected {len(lons_list)}, got {len(aCenter_longititude)}')
            validation_passed = False
        if aConnectivity is None or aConnectivity.size == 0:
            logger.error('Empty connectivity matrix')
            validation_passed = False
        # Validate connectivity indices
        if aConnectivity is not None:
            max_vertex_index = len(aVertex_longititude) - 1
            valid_connectivity = aConnectivity[aConnectivity >= 0]
            if len(valid_connectivity) > 0 and np.max(valid_connectivity) > max_vertex_index:
                logger.error(
                    'Connectivity matrix contains invalid vertex indices')
                validation_passed = False
        # Check for reasonable mesh bounds
        if len(aVertex_longititude) > 0:
            vertex_lon_range = (np.min(aVertex_longititude),
                                np.max(aVertex_longititude))
            vertex_lat_range = (np.min(aVertex_latitude),
                                np.max(aVertex_latitude))
            # Basic range reporting (detailed validation done by check_geometry_validity)
            if iFlag_verbose_in:
                logger.debug(f'Vertex longitude range: {vertex_lon_range}')
                logger.debug(f'Vertex latitude range: {vertex_lat_range}')
        if not validation_passed:
            logger.error('Mesh topology rebuild failed validation')
            return None
        if iFlag_verbose_in:
            logger.info('Mesh topology successfully rebuilt')
            logger.info(f'Final mesh statistics:')
            logger.info(f'  - Unique vertices: {len(aVertex_longititude)}')
            logger.info(f'  - Mesh cells: {len(aCenter_longititude)}')
            logger.info(f'  - Max vertices per cell: {nVertex_max}')
            logger.info(f'  - Connectivity shape: {aConnectivity.shape}')
            logger.info(
                f'  - Vertex longitude range: [{np.min(aVertex_longititude):.3f}, {np.max(aVertex_longititude):.3f}]')
            logger.info(
                f'  - Vertex latitude range: [{np.min(aVertex_latitude):.3f}, {np.max(aVertex_latitude):.3f}]')

        # Return comprehensive mesh topology information
        mesh_info = {
            'vertices_longitude': aVertex_longititude,
            'vertices_latitude': aVertex_latitude,
            'connectivity': aConnectivity,
            'cell_centroids_longitude': aCenter_longititude,
            'cell_centroids_latitude': aCenter_latitude,
            'cell_ids': aCellID,
            'area_min': dArea_min,
            'area_max': dArea_max,
            'area_mean': dArea_mean,
            'max_vertices_per_cell': nVertex_max,
            'num_cells': nFeatures,
            'num_polygns': len(aCenter_longititude),
            'num_vertices': len(aVertex_longititude),
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
        # Open target mesh file
        pDataset = ogr.Open(sFilename, 0)  # Read-only
        if pDataset is None:
            logger.error(f'Failed to open target mesh file: {sFilename}')
            return None

        # Get first layer
        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            logger.error('Failed to get layer from target mesh dataset')
            pDataset = None
            return None

        # Get layer definition
        pLayerDefn = pLayer.GetLayerDefn()
        if pLayerDefn is None:
            logger.error('Failed to get layer definition from target mesh')
            pDataset = None
            return None

        # Get field information
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

        # Extract variable data from features, handling multipolygons correctly
        aData_list = []
        pLayer.ResetReading()
        pFeature = pLayer.GetNextFeature()
        iFeature_count = 0

        iCount_multipolygons = 0

        while pFeature is not None:
            pGeometry = pFeature.GetGeometryRef()
            if pGeometry is not None:
                sGeometry_type = pGeometry.GetGeometryName()
                if sGeometry_type == 'POLYGON':
                    # Single polygon - add one data value
                    try:
                        dField_value = pFeature.GetField(sVariable)
                        aData_list.append(dField_value if dField_value is not None else np.nan)
                    except Exception as e:
                        logger.warning(f'Error reading field {sVariable} from feature {iFeature_count}: {e}')
                        aData_list.append(np.nan)

                elif sGeometry_type == 'MULTIPOLYGON':
                    iCount_multipolygons += 1
                    # Multipolygon - add the same data value for each polygon part
                    try:
                        dField_value = pFeature.GetField(sVariable)
                        dData_value = dField_value if dField_value is not None else np.nan

                        # Add the same data value for each polygon part in the multipolygon
                        nGeometryParts = pGeometry.GetGeometryCount()
                        iValid_parts = 0

                        for iPart in range(nGeometryParts):
                            pPolygon_part = pGeometry.GetGeometryRef(iPart)
                            if pPolygon_part is not None and pPolygon_part.IsValid():
                                aData_list.append(dData_value)
                                iValid_parts += 1
                            else:
                                pPolygon_part.FlattenTo2D()
                                sWkt = pPolygon_part.ExportToWkt()
                                logger.warning(f'Invalid polygon part {iPart} in multipolygon feature {iFeature_count}: {sWkt} ')
                                aData_list.append(np.nan)

                        if iValid_parts == 0:
                            logger.warning(f'No valid parts found in multipolygon feature {iFeature_count}')
                            #aData_list.append(np.nan)

                    except Exception as e:
                        logger.warning(f'Error reading field {sVariable} from multipolygon feature {iFeature_count}: {e}')
                        aData_list.append(np.nan)
                else:
                    logger.warning(f'Feature {iFeature_count} has unsupported geometry type: {sGeometry_type}')
                    aData_list.append(np.nan)
            else:
                logger.warning(f'Feature {iFeature_count} has no geometry')
                aData_list.append(np.nan)

            iFeature_count += 1
            pFeature = pLayer.GetNextFeature()

        # Close dataset
        pDataset = None

        if not aData_list:
            logger.error('No data extracted from target mesh')
            return None

        # Convert to numpy array
        aData = np.array(aData_list, dtype=np.float64)

        return aData, nFeatures

    except Exception as e:
        logger.error(f'Error extracting target mesh data: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return None

def _validate_output_path(sFilename: Optional[str]) -> bool:
    """
    Validate output file path and create directories if needed.

    Args:
        sFilename: Output file path to validate

    Returns:
        bool: True if path is valid and accessible, False otherwise
    """
    if sFilename is None:
        return True  # Interactive mode, no file validation needed

    if not isinstance(sFilename, str) or not sFilename.strip():
        logger.error('Output filename must be a non-empty string')
        return False

    # Check output directory exists
    sOutput_dir = os.path.dirname(sFilename)
    if sOutput_dir and not os.path.exists(sOutput_dir):
        try:
            os.makedirs(sOutput_dir, exist_ok=True)
            logger.info(f'Created output directory: {sOutput_dir}')
        except Exception as e:
            logger.error(f'Cannot create output directory {sOutput_dir}: {e}')
            return False

    # Check supported file extensions
    sFile_ext = os.path.splitext(sFilename)[1].lower()
    aAll_valid_extensions = VALID_IMAGE_FORMATS + [f'.{fmt}' for fmt in VALID_ANIMATION_FORMATS]
    if sFile_ext not in aAll_valid_extensions:
        logger.warning(f'File extension {sFile_ext} may not be supported. '
                      f'Recommended: {", ".join(VALID_IMAGE_FORMATS)}')

    return True

def _setup_geovista_plotter(iFlag_off_screen: bool = False, iFlag_verbose_in: bool = False):
    """
    Set up GeoVista plotter with error handling.

    Args:
        iFlag_off_screen: Whether to create off-screen plotter
        iFlag_verbose_in: Enable verbose logging

    Returns:
        GeoVista plotter instance or None if failed
    """
    try:
        import geovista as gv
        if iFlag_verbose_in:
            logger.info('GeoVista library imported successfully')
    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return None

    try:
        if iFlag_off_screen:
            pPlotter = gv.GeoPlotter(off_screen=True)
            if iFlag_verbose_in:
                logger.debug('Created off-screen plotter')
        else:
            pPlotter = gv.GeoPlotter()
            if iFlag_verbose_in:
                logger.debug('Created interactive plotter')
        return pPlotter
    except Exception as e:
        logger.error(f'Failed to create GeoVista plotter: {e}')
        logger.error('This may be due to missing graphics context or display')
        if iFlag_off_screen:
            logger.error('For headless systems, ensure proper OpenGL/Mesa setup')
        else:
            logger.error('For interactive mode, ensure display environment (X11/Wayland) is available')
        return None

def _add_geographic_context(pPlotter, pConfig: VisualizationConfig):
    """
    Add geographic context (coastlines, graticule, axes) to plotter.

    Args:
        pPlotter: GeoVista plotter instance
        pConfig: Visualization configuration
    """
    # Add coastlines
    if pConfig.show_coastlines:
        try:
            # You can set coastline color using the 'color' parameter
            # Common options: 'black', 'white', 'red', 'blue', 'gray', etc.
            # You can also use RGB tuples like (1.0, 0.0, 0.0) for red
            pPlotter.add_coastlines(color=pConfig.coastline_color, line_width=pConfig.coastline_width)
            if pConfig.verbose:
                logger.debug(f'Added coastlines overlay (color: {pConfig.coastline_color}, width: {pConfig.coastline_width})')
        except Exception as e:
            logger.warning(f'Could not add coastlines: {e}')

    # Add coordinate axes
    try:
        pPlotter.add_axes()
        if pConfig.verbose:
            logger.debug('Added coordinate axes')
    except Exception as e:
        logger.warning(f'Could not add axes: {e}')

    # Add graticule (coordinate grid)
    if pConfig.show_graticule:
        try:
            pPlotter.add_graticule(show_labels=True)
            if pConfig.verbose:
                logger.debug('Added coordinate graticule with labels')
        except Exception as e:
            logger.warning(f'Could not add graticule: {e}')

def _configure_camera(pPlotter, pConfig: VisualizationConfig) -> bool:
    """
    Configure camera position and orientation.

    Args:
        pPlotter: GeoVista plotter instance
        pConfig: Visualization configuration

    Returns:
        bool: True if camera configured successfully, False otherwise
    """
    try:
        aFocal_point, aCamera_position = CameraController.calculate_camera_position(
            pConfig.longitude_focus, pConfig.latitude_focus, pConfig.zoom_factor
        )

        # Validate camera setup
        if not CameraController.validate_camera_setup(aFocal_point, aCamera_position):
            logger.warning('Camera and focal point are too close, using default view')
            raise ValueError('Invalid camera positioning')

        pPlotter.camera.focal_point = aFocal_point
        pPlotter.camera.position = aCamera_position
        pPlotter.camera.zoom(pConfig.zoom_factor)
        pPlotter.camera.up = [0, 0, 1]  # Z-up orientation

        if pConfig.verbose:
            logger.debug(f'Camera configured successfully:')
            logger.debug(f'  Focal point: {aFocal_point}')
            logger.debug(f'  Camera position: {aCamera_position}')

        return True

    except Exception as e:
        logger.warning(f'Error setting camera position: {e}. Using default view.')
        try:
            pPlotter.reset_camera()
        except Exception:
            pass
        return False

def add_camera_aware_arrow(pPlotter, normalized_x: float, normalized_y: float,
                          surface_lon: float, surface_lat: float,
                          arrow_length: float = 0.5, arrowhead_length: float = 0.08,
                          arrowhead_radius: float = 0.02, shaft_width: float = 1.0,
                          distance_from_camera: float = 2.0, color: str = 'red'):
    """
    Add a camera-aware arrow pointing from screen-space coordinates to Earth surface.

    Args:
        pPlotter: GeoVista plotter instance
        normalized_x: Normalized x coordinate [0,1] in screen space
        normalized_y: Normalized y coordinate [0,1] in screen space
        surface_lon: Target longitude on Earth surface (degrees)
        surface_lat: Target latitude on Earth surface (degrees)
        arrow_length: Total arrow length (default: 0.5)
        arrowhead_length: Length of arrowhead cone (default: 0.08)
        arrowhead_radius: Radius of arrowhead cone (default: 0.02)
        shaft_width: Width of arrow shaft line (default: 1.0)
        distance_from_camera: Distance from camera for space point (default: 2.0)
        color: Arrow color (default: 'red')
    """
    import numpy as np
    import pyvista as pv

    # Convert surface point to 3D coordinates
    surface_lon_rad = np.radians(surface_lon)
    surface_lat_rad = np.radians(surface_lat)
    surface_point = np.array([
        1.0 * np.cos(surface_lat_rad) * np.cos(surface_lon_rad),  # Earth radius = 1.0
        1.0 * np.cos(surface_lat_rad) * np.sin(surface_lon_rad),
        1.0 * np.sin(surface_lat_rad)
    ])

    # Get current camera information for screen-space to world-space conversion
    camera = pPlotter.camera
    camera_pos = np.array(camera.position)
    camera_focal = np.array(camera.focal_point)
    camera_up = np.array(camera.up)

    # Calculate camera coordinate system
    camera_forward = camera_focal - camera_pos
    camera_forward = camera_forward / np.linalg.norm(camera_forward)

    camera_right = np.cross(camera_forward, camera_up)
    camera_right = camera_right / np.linalg.norm(camera_right)

    camera_up_corrected = np.cross(camera_right, camera_forward)
    camera_up_corrected = camera_up_corrected / np.linalg.norm(camera_up_corrected)

    # Map normalized coordinates to screen space [-1, 1]
    screen_x = normalized_x * 2 - 1
    screen_y = normalized_y * 2 - 1

    # Calculate space point in camera coordinate system
    space_point = (camera_pos +
                   camera_forward * distance_from_camera +
                   camera_right * screen_x * distance_from_camera * 0.5 +
                   camera_up_corrected * screen_y * distance_from_camera * 0.5)

    # Calculate direction vector (from space to surface)
    direction = surface_point - space_point
    direction = direction / np.linalg.norm(direction)  # Normalize

    # Use the camera-relative space point as arrow start
    arrow_start = space_point

    # Calculate shaft end point (where arrowhead base starts)
    shaft_end = surface_point - direction * arrowhead_length

    # Position cone so its tip touches the surface point
    arrowhead_center = surface_point - direction * (arrowhead_length / 2)

    # Create arrow shaft (line from space point to shaft end)
    shaft_points = np.array([arrow_start, shaft_end])
    pPlotter.add_lines(
        shaft_points,
        color=color,
        width=shaft_width
    )

    # Create arrowhead (cone with tip at surface point)
    cone = pv.Cone(
        center=arrowhead_center,
        direction=direction,      # Point toward surface
        height=arrowhead_length,  # Independent arrowhead length
        radius=arrowhead_radius,  # Independent arrowhead width
        resolution=8
    )
    pPlotter.add_mesh(cone, color=color)

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
                          iFlag_create_animation: Optional[bool] = False,
                         iAnimation_frames: Optional[int] = 36,
                         dAnimation_speed: Optional[float] = 1.0,
                         sAnimation_format: Optional[str] = 'mp4',
                          sEdge_color: str = 'black',
                          sElevation_colormap: str = 'terrain',
                          dElevation_opacity: float = 0.85,
                          sPolyline_width_attribute: Optional[str] = None,
                          dPolyline_width_min: float = 0.5,
                          dPolyline_width_max: float = 3.0,
                         iFlag_verbose_in: Optional[bool] = False,
                         iFlag_enable_profiling: Optional[bool] = False,
                         window_size: Tuple[int, int] = (800, 600)
                         ) -> bool:

    """
    Visualize the source mesh topology using GeoVista 3D globe rendering.

    Creates an interactive or saved 3D visualization of the unstructured mesh
    with proper geographic context including coastlines and coordinate grid.

    Args:
        sFilename_mpas_mesh_in: Path to MPAS mesh NetCDF file
        sFilename_uraster_in: Path to elevation raster file (GeoJSON/Shapefile)
        sFilename_flow_directions_in: Path to flow directions vector file
        sFilename_out: Output screenshot file path. If None, displays interactive viewer.
            Supports formats: .png, .jpg, .svg
        dLongitude_focus_in: Camera focal point longitude in degrees (-180 to 180).
            Default is 0.0 (prime meridian).
        dLatitude_focus_in: Camera focal point latitude in degrees (-90 to 90).
            Default is 0.0 (equator).
        dZoom_factor: Camera zoom level. Higher values zoom in. Default is 0.7.
        iFlag_show_coastlines: Show coastline overlay. Default is True.
        iFlag_show_graticule: Show coordinate grid with labels. Default is False.
        sCoastline_color: Color for coastlines. Default is 'black'.
            Examples: 'white', 'red', 'blue', 'gray', or RGB tuples like (1.0, 0.0, 0.0).
        dCoastline_width: Line width for coastlines. Default is 1.0.
        iFlag_create_animation: Create animation instead of static image. Default is False.
        iAnimation_frames: Number of animation frames. Default is 36.
        dAnimation_speed: Animation speed multiplier. Default is 1.0.
        sAnimation_format: Animation format ('mp4', 'gif', 'avi'). Default is 'mp4'.
        sEdge_color: Color of mesh edges. Default is 'black'.
        sElevation_colormap: Colormap for elevation data. Default is 'terrain'.
            Popular options: 'terrain', 'viridis', 'plasma', 'coolwarm', 'jet', 'rainbow',
            'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'bone', 'copper',
            'gray', 'pink', 'flag', 'prism', 'ocean', 'gist_earth', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv', 'tab10', 'tab20'
        dElevation_opacity: Opacity of elevation mesh (0.0 to 1.0). Default is 1.0.
        sPolyline_width_attribute: Attribute name for variable polyline widths (e.g., 'discharge', 'flow_rate').
            If None, uses uniform width. Default is None.
        dPolyline_width_min: Minimum line width for attribute-based scaling. Default is 0.5.
        dPolyline_width_max: Maximum line width for attribute-based scaling. Default is 3.0.
        iFlag_verbose_in: If True, print detailed progress messages. Default is False.
        iFlag_enable_profiling: If True, enable performance profiling. Default is False.
        iFlag_full_mesh: If 1, display the full MPAS mesh wireframe. Default is 1.
        iFlag_elevation: If 1, display elevation data as colored mesh. Default is 1.
        iFlag_flow_directions: If 1, display flow direction polylines. Default is 0.

    Returns:
        True if visualization successful, False otherwise

    Note:
        - Requires 'geovista' package: pip install geovista
        - Interactive mode requires display environment
        - Mesh topology must be built before visualization (call rebuild_mesh_topology first)

    Examples:
        # Basic usage with terrain colormap
        visualize_unified_mesh_dataset(mesh_file, raster_file, flow_file)

        # Custom colormap for elevation
        visualize_unified_mesh_dataset(mesh_file, raster_file, flow_file,
                                     sElevation_colormap='viridis')

        # Semi-transparent elevation with custom colormap
        visualize_unified_mesh_dataset(mesh_file, raster_file, flow_file,
                                     sElevation_colormap='coolwarm',
                                     dElevation_opacity=0.7)
    """
    # read the mpas mesh file
    import netCDF4 as nc
    import geovista as gv

    # Define a base window size and font size
    base_window_size = 1200  # Reference window size
    base_font_size = 8     # Reference font size
    base_line_width = 0.1   # Reference line width

    # Calculate scaling factor
    scaling_factor = max(window_size) / base_window_size

    # Adjust font size
    scaled_font_size = int(base_font_size * scaling_factor)
    scaled_line_width = base_line_width * scaling_factor

    # Flag parameters are now passed as function arguments

    if not _validate_output_path(sFilename_out):
        return False

    # Create configuration object
    config = VisualizationConfig(
        longitude_focus=dLongitude_focus_in,
        latitude_focus=dLatitude_focus_in,
        zoom_factor=dZoom_factor,
        show_coastlines=iFlag_show_coastlines,
        show_graticule=iFlag_show_graticule,
        coastline_color=sCoastline_color,
        coastline_width=dCoastline_width,
        verbose=iFlag_verbose_in
    )

    animation_config = AnimationConfig(
        frames=iAnimation_frames,
        speed=dAnimation_speed,
        format=sAnimation_format
    ) if iFlag_create_animation else None
    import numpy as np
    try:

        # Setup plotter
        pPlotter = _setup_geovista_plotter(iFlag_off_screen=(sFilename_out is not None), iFlag_verbose_in=config.verbose)
        if pPlotter is None:
            return False

        if config.verbose:
            logger.info(f'Created GeoVista mesh with {pMesh_full.n_cells} cells and {pMesh_full.n_points} points')

        if iFlag_full_mesh:
            # Use optimized mesh loading
            mesh_data = load_mpas_mesh_optimized(sFilename_mpas_mesh_in, config.verbose or iFlag_enable_profiling)
            if mesh_data is None:
                logger.error(f'Failed to load MPAS mesh from {sFilename_mpas_mesh_in}')
                return False

            aVertex_longititude = mesh_data['vertices_longitude']
            aVertex_latitude = mesh_data['vertices_latitude']
            aConnectivity = mesh_data['connectivity']
            aCellID = mesh_data['cell_ids']

            if config.verbose:
                logger.info('Creating mesh visualization...')
                logger.info(f'  - Vertices: {len(aVertex_longititude)}')
                logger.info(f'  - Connectivity shape: {aConnectivity.shape}')
                logger.info(f'  - Focus: ({config.longitude_focus:.2f}°, {config.latitude_focus:.2f}°)')
                logger.info(f'  - Zoom factor: {config.zoom_factor}')

            # Validate connectivity array structure
            if aConnectivity.ndim != 2:
                logger.error(f'Connectivity array must be 2D, got {aConnectivity.ndim}D')
                return False

            # Use optimized connectivity processing
            connectivity_masked = optimize_connectivity_processing(
                aConnectivity, len(aVertex_longititude), config.verbose or iFlag_enable_profiling
            )

            # Transform to GeoVista unstructured mesh with performance monitoring
            with PerformanceMonitor("GeoVista mesh creation", config.verbose or iFlag_enable_profiling):
                pMesh_full = gv.Transform.from_unstructured(
                    aVertex_longititude,
                    aVertex_latitude,
                    connectivity=connectivity_masked,
                    crs=CRS
                )

            sScalars = 'Cell ID'
            #pMesh_full.cell_data[sScalars] = aCellID
            sargs = {
                "title": sScalars,
                "shadow": True,
                "title_font_size": scaled_font_size,
                "label_font_size": scaled_font_size,
                "fmt": "%.0f",  # Integer formatting for cell IDs
                "position_x": 0.9,  # Move to the right (closer to 1.0)
                "position_y": 0.10,  # Adjust vertical position (closer to 0.0 for bottom)
                "height": 0.8,
                "vertical": True,    # Make the colorbar vertical

            }
            pPlotter.add_base_layer(mesh=pMesh_full, texture= gv.natural_earth_hypsometric())
            pPlotter.add_mesh(
                pMesh_full,
                style='wireframe',
                line_width=scaled_line_width,
                color=sEdge_color,
                show_edges=True,
                #scalars=sScalars,
                #scalar_bar_args=sargs,
                show_scalar_bar=False,
                label='Unified MPAS global mesh'
            )


        if iFlag_elevation == 1:
            # Use optimized elevation processing
            elevation_data = process_elevation_data_optimized(
                sFilename_uraster_in,
                'mean',
                config.verbose or iFlag_enable_profiling
            )

            if elevation_data is None:
                logger.warning(f'Failed to process elevation data from: {sFilename_uraster_in}')
            else:
                # Create GeoVista mesh with optimized data
                with PerformanceMonitor("Elevation mesh creation", config.verbose or iFlag_enable_profiling):
                    pMesh_elevation = gv.Transform.from_unstructured(
                        elevation_data['vertices_longitude'],
                        elevation_data['vertices_latitude'],
                        connectivity=elevation_data['connectivity_masked'],
                        crs=CRS
                    )

                    sScalars = 'Elevation (m)'
                    pMesh_elevation.cell_data[sScalars] = elevation_data['elevation_data']

                    sargs = {
                        "title": sScalars,
                        "shadow": True,
                        "title_font_size": scaled_font_size * 2,
                        "label_font_size": int(scaled_font_size * 1.5),
                        "height": 0.8,
                        "fmt": "%.0f",
                        "position_x": 0.9,
                        "position_y": 0.1,
                        "vertical": True

                    }

                    pPlotter.add_mesh(
                        pMesh_elevation,
                        scalars=sScalars,
                        scalar_bar_args=sargs,
                        cmap=sElevation_colormap,
                        opacity=dElevation_opacity,
                        show_edges=True,
                        edge_color=sEdge_color,
                        line_width=scaled_line_width * 0.8,
                        label='Surface elevation'
                    )

                if config.verbose:
                    logger.info('✓ Successfully added elevation mesh using optimized processing')


        if iFlag_flow_directions == 1:
            import pyvista as pv

            # Use optimized polyline processing with attribute-based width scaling
            polyline_data = process_polylines_optimized(
                sFilename_flow_directions_in,
                scaled_line_width * 1.5,
                config.verbose or iFlag_enable_profiling,
                width_attribute=sPolyline_width_attribute,
                width_min=dPolyline_width_min,
                width_max=dPolyline_width_max
            )

            if polyline_data is None:
                logger.error('Failed to process polylines')
                return False

            # Handle variable line widths vs uniform line widths
            if polyline_data.get('use_variable_width', False):
                # Variable line widths: process each polyline individually
                with PerformanceMonitor("Variable width polyline creation", config.verbose or iFlag_enable_profiling):
                    sColor_polyline_in = 'blue'

                    # Process polylines individually for variable widths
                    points = polyline_data['points']
                    lines = polyline_data['lines']
                    line_widths = polyline_data['line_widths']

                    # Parse lines array to extract individual polylines
                    line_idx = 0
                    polyline_idx = 0

                    while line_idx < len(lines):
                        num_points = lines[line_idx]
                        line_indices = lines[line_idx + 1:line_idx + 1 + num_points]

                        # Extract points for this polyline
                        polyline_points = [points[i] for i in line_indices]

                        # Create individual polyline
                        individual_lines = [num_points] + list(range(num_points))
                        individual_polydata = pv.PolyData(polyline_points, lines=individual_lines)

                        # Add with individual line width
                        current_width = line_widths[polyline_idx] if polyline_idx < len(line_widths) else scaled_line_width
                        pPlotter.add_mesh(
                            individual_polydata,
                            color=sColor_polyline_in,
                            line_width=current_width,
                            name=f"polyline_{polyline_idx}",
                            label='Flow directions'
                        )

                        line_idx += num_points + 1
                        polyline_idx += 1

                    if config.verbose:
                        logger.info(f'✓ Successfully added {polyline_idx} polylines with variable widths')
                        logger.info(f'  Width range: {min(line_widths):.2f} to {max(line_widths):.2f}')
            else:
                # Uniform line width: use efficient batch processing
                with PerformanceMonitor("Uniform width polyline creation", config.verbose or iFlag_enable_profiling):
                    polydata = pv.PolyData(polyline_data['points'], lines=polyline_data['lines'])

                    sColor_polyline_in = 'blue'
                    pPlotter.add_mesh(
                        polydata,
                        color=sColor_polyline_in,
                        line_width=polyline_data['line_width'],
                        name="combined_polylines",
                        label='Flow directions'
                    )

                    if config.verbose:
                        logger.info(f'✓ Successfully added {polyline_data["polyline_count"]} polylines with uniform width {polyline_data["line_width"]}')

        # Configure camera
        _configure_camera(pPlotter, config)

        # Add geographic context
        _add_geographic_context(pPlotter, config)

        pPlotter.add_legend(loc='center left')

        #pPlotter.add_text("MPAS mesh", position=(0.05, 0.05), font_size=scaled_font_size * 1.5)

        iFlag_arrow = 1
        if iFlag_arrow:
            # Add arrows based on the number of legend entries added
            legend_entries = []
            legend_colors = []
            legend_targets = []

            # Collect legend information from enabled features
            if iFlag_full_mesh == 1:
                legend_entries.append('Unified MPAS global mesh')
                legend_colors.append('black')
                legend_targets.append((-143.0, 52.0))  # ocean

            if iFlag_elevation == 1:
                legend_entries.append('Surface elevation')
                legend_colors.append('black')
                legend_targets.append((-111.0, 40.0))  # North Americar

            if iFlag_flow_directions == 1:
                legend_entries.append('Flow directions')
                legend_colors.append('black')
                legend_targets.append((-100.0, 30.0))  #

            # Calculate legend position based on 'upper left' location
            # For 'upper left': typically starts around (0.02, 0.98) and extends down
            legend_base_x = 0.02  # Left margin
            legend_base_y = 0.58  # Top margin
            legend_item_height = 0.05  # Approximate height per legend item
            legend_width = 0.15   # Approximate legend width

            # Add arrows for each legend entry
            for i, (entry, color, target) in enumerate(zip(legend_entries, legend_colors, legend_targets)):
                # Calculate arrow start position relative to legend position
                # Start arrows from the right edge of the legend, at each item's vertical position
                normalized_x = legend_base_x + legend_width + 0.02  # Right edge + small offset
                normalized_y = legend_base_y - (i * legend_item_height) - (legend_item_height / 2)  # Center of legend item
                target_lon, target_lat = target

                add_camera_aware_arrow(
                    pPlotter,
                    normalized_x=normalized_x, normalized_y=normalized_y,
                    surface_lon=target_lon, surface_lat=target_lat,
                    arrow_length=0.4 + (i * 0.1),  # Vary arrow length
                    arrowhead_length=0.06, arrowhead_radius=0.015,
                    shaft_width=scaled_line_width * 5,
                    color=color  # Keep existing black color
                )


        # Output or display
        return _handle_visualization_output(pPlotter, sFilename_out, config.verbose, window_size=window_size)

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

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
            #check extension
            sExt = os.path.splitext(sFilename)[1].lower()
            #use different method for vector and raster formats
            if sExt in VALID_IMAGE_FORMATS:
                pPlotter.screenshot(sFilename, window_size=window_size)
            else:
                pPlotter.save_graphic(sFilename) #pdf format
                #pPlotter.show()
                #plt.imshow(pPlotter.image)
                #plt.axis("off")
                #plt.savefig(sFilename, dpi=500, bbox_inches="tight", pad_inches=0)
            if iFlag_verbose_in:
                logger.info(f'✓ Visualization saved to: {sFilename}')

                # Verify file was created
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


if __name__ == "__main__":
    sFilename_mpas_mesh_in = r'C:\workspace\python\unified_land_river_mesh\data\global\base_mesh.nc'  # specify input file
    sFilename_mpas_mesh_invert_in = r'C:\workspace\python\unified_land_river_mesh\data\global\invert_mesh.nc'  # specify inverted input file
    sFilename_uraster_in = r'C:\workspace\python\uraster\data\example_12\output\mpas_uraster.geojson'  # specify output file
    sFilename_flow_directions_in = r'C:\workspace\python\unified_land_river_mesh\data\global\mpas_flow_direction.geojson'  # specify flow directions file

    sFilename_out =  r'C:\workspace\python\unified_land_river_mesh\figures\global_mesh_visualization.png'  # specify output visualization file

    visualize_unified_mesh_dataset(sFilename_mpas_mesh_in,
                      sFilename_uraster_in,
                      sFilename_flow_directions_in,
                      sFilename_out = sFilename_out,
                         iFlag_full_mesh = 1,
                         iFlag_elevation = 1,
                         iFlag_flow_directions = 0,
                      dLongitude_focus_in = -97,
                      dLatitude_focus_in=43,
                      dZoom_factor=0.7,
                      sPolyline_width_attribute = None, # 'drainage_area',
                      iFlag_enable_profiling=True,
                      window_size=(8000, 7000) )

    pass