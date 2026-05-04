# Antimeridian Buffer Implementation Plan

## Objective
Improve the `fix_raster_antimeridian_issue` function to select multiple buffer columns on both sides of the antimeridian using the `iRaster_buffer_pixel` parameter.

## Current Implementation Issues

1. Only selects exact dateline columns
2. Doesn't utilize the `iRaster_buffer_pixel` parameter
3. Has placeholder comments for left/right processing
4. Line 110 references undefined `updated_columns` variable
5. **CRITICAL**: For global rasters (-180° to 180°), the dateline is at BOTH column 0 AND the last column (wrap-around)

## Key Insight: Global Raster Wrap-Around

For global rasters with extent from -180° to 180°:
- **Column 0 (index 0)**: Represents longitude ≈ -180°
- **Column ncols-1 (last column)**: Represents longitude ≈ +180°
- These are the SAME longitude line (antimeridian) but at opposite edges of the raster

Therefore, we need to select:
- **Left buffer**: Last N columns of the raster (right edge, wrapping to -180°)
- **Right buffer**: First N columns of the raster (left edge, starting from -180°)

## Proposed Solution

### Algorithm Overview

```
1. Identify dateline column(s) where |longitude| ≈ 180°
2. Use leftmost dateline column as reference point
3. Expand selection to include buffer columns:
   - Left side: [reference_col - buffer_pixels : reference_col]
   - Right side: [reference_col : reference_col + buffer_pixels + 1]
4. Apply boundary checks to prevent index errors
5. Combine and return all buffer column indices
```

### Detailed Implementation Steps

#### Step 1: Identify Left-Side Buffer Columns (Last N columns of raster)
```python
# For global rasters, left side of dateline is the RIGHT EDGE of the raster
# These are the last iRaster_buffer_pixel columns (near +180° wrapping to -180°)
left_buffer_cols = np.arange(ncols - iRaster_buffer_pixel, ncols)
```

**Logic:**
- Left side of dateline (near -180°) is at the **right edge** of the raster
- Select the last `iRaster_buffer_pixel` columns
- Example: if ncols=360 and buffer=2, select columns [358, 359]

#### Step 2: Identify Right-Side Buffer Columns (First N columns of raster)
```python
# For global rasters, right side of dateline is the LEFT EDGE of the raster
# These are the first iRaster_buffer_pixel columns (starting from -180°)
right_buffer_cols = np.arange(0, iRaster_buffer_pixel)
```

**Logic:**
- Right side of dateline is at the **left edge** of the raster
- Select the first `iRaster_buffer_pixel` columns
- Example: if buffer=2, select columns [0, 1]

#### Step 3: Combine All Buffer Columns
```python
# Combine left and right buffer columns
all_buffer_cols = np.concatenate([left_buffer_cols, right_buffer_cols])

# Already sorted and unique by construction (no overlap possible)
# all_buffer_cols will be [ncols-buffer, ncols-buffer+1, ..., ncols-1, 0, 1, ..., buffer-1]
```

**Key Point:** No need for `np.unique()` since left and right buffers cannot overlap (they're at opposite ends)

#### Step 4: Replace Current Line 110
Replace:
```python
raster_data[:, dateline_columns] = updated_columns
```

With:
```python
raster_data[:, all_buffer_cols] = target_value
```

Or if applying a special function:
```python
# Extract the columns that need processing
columns_to_process = raster_data[:, all_buffer_cols]
# Apply special function (user-defined)
updated_columns = special_function(columns_to_process)
# Write back
raster_data[:, all_buffer_cols] = updated_columns
```

**Important Note:** The order of columns in `all_buffer_cols` is: [right-edge columns, left-edge columns]
This matches the logical order: [near +180°, near -180°]

## Edge Cases to Consider

### Case 1: Buffer exceeds raster width
**Scenario:** ncols = 360, buffer = 400
- Left buffer would try: columns [ncols-400, ncols] → negative start index
- Right buffer would try: columns [0, 400] → exceeds ncols
- **Solution:**
  ```python
  left_buffer_cols = np.arange(max(0, ncols - iRaster_buffer_pixel), ncols)
  right_buffer_cols = np.arange(0, min(ncols, iRaster_buffer_pixel))
  ```
- **Result:** Safely clips to available columns

### Case 2: Buffer = 0
**Scenario:** `iRaster_buffer_pixel = 0`
- Left buffer: `np.arange(ncols-0, ncols)` → empty array
- Right buffer: `np.arange(0, 0)` → empty array
- **Result:** No columns selected, raster unchanged (correct behavior)

### Case 3: Buffer = 1 (default)
**Scenario:** `iRaster_buffer_pixel = 1`, ncols = 360
- Left buffer: columns [359] (last column at +180°)
- Right buffer: columns [0] (first column at -180°)
- **Result:** Exactly the two dateline columns are selected

### Case 4: Large buffer overlapping in middle
**Scenario:** ncols = 100, buffer = 60
- Left buffer: columns [40, 41, ..., 99]
- Right buffer: columns [0, 1, ..., 59]
- **Total:** 120 column indices (more than exists!)
- **Impact:** All columns get modified, which may be unintended
- **Recommendation:** Add validation check if buffer > ncols/2

### Case 5: Regional (non-global) rasters
**Scenario:** Raster covers only North America (-170° to -50°), not global
- The current `dateline_columns` detection would find NO dateline columns
- **Result:** Function returns early (if dateline_columns.size == 0), no modification needed
- **Status:** Already handled correctly by existing code

## Function Documentation Update

Update the docstring to include:

```python
def fix_raster_antimeridian_issue(sFilename_in,
                                  sFilename_out,
                                  target_value,
                                  iRaster_buffer_pixel=1):
    """
    Replace raster values at the international dateline and surrounding buffer columns
    to avoid artifacts when creating buffer zones.

    This function identifies columns at the antimeridian (±180° longitude) and replaces
    values in those columns plus a configurable buffer on each side. This prevents
    buffer dilation operations from incorrectly propagating across the dateline.

    Parameters:
    -----------
    sFilename_in : str
        Input raster file path
    sFilename_out : str
        Output raster file path
    target_value : int or float
        Value to assign to antimeridian buffer columns
    iRaster_buffer_pixel : int, optional (default=1)
        Number of buffer columns to include on each side of the dateline

    Returns:
    --------
    None

    Example:
    --------
    With iRaster_buffer_pixel=2 and dateline at column 100:
    - Columns selected: [98, 99, 100, 101, 102]
    - All selected columns set to target_value
    """
```

## Testing Recommendations

1. **Test with buffer=0**: Should select no columns (empty arrays)
2. **Test with buffer=1**: Should select columns [ncols-1, 0] (the two dateline columns)
3. **Test with buffer=5**: Should select columns [ncols-5 to ncols-1, 0 to 4] (10 columns total)
4. **Test with large buffer**: buffer=50 on a 360-column raster
5. **Test with regional raster**: Raster that doesn't cross dateline (should exit early)
6. **Test with buffer > ncols/2**: Ensure graceful handling or warning

## Final Implementation Code

### Complete Code for lines 96-112 in raster_process.py

Replace the existing implementation with:

```python
# Read the raster data as a numpy array
raster_data = band.ReadAsArray()

# Identify columns corresponding to the antimeridian using normalized longitude
ncols = raster_data.shape[1]
lon = geo_transform[0] + np.arange(ncols) * geo_transform[1]
lon180 = ((lon + 180.0) % 360.0) - 180.0

# Check if this raster crosses the antimeridian
tolerance = abs(geo_transform[1]) * 0.5
dateline_columns = np.where(np.isclose(np.abs(lon180), 180.0, atol=tolerance))[0]

# Process only if the antimeridian exists in this raster extent
if dateline_columns.size > 0:
    # For global rasters (-180 to 180), the antimeridian wraps around:
    # - Left side (near -180°): last iRaster_buffer_pixel columns (right edge)
    # - Right side (near +180°): first iRaster_buffer_pixel columns (left edge)

    # Select buffer columns from the right edge (wrapping to left/west at -180°)
    left_start = max(0, ncols - iRaster_buffer_pixel)
    left_buffer_cols = np.arange(left_start, ncols)

    # Select buffer columns from the left edge (east side at -180° going toward 0°)
    right_end = min(ncols, iRaster_buffer_pixel)
    right_buffer_cols = np.arange(0, right_end)

    # Combine both buffer regions
    all_buffer_cols = np.concatenate([left_buffer_cols, right_buffer_cols])

    # Optional: Warn if buffer is too large
    if len(all_buffer_cols) > ncols:
        print(f"Warning: Buffer ({iRaster_buffer_pixel}) exceeds raster width ({ncols})")

    # Apply target value to all buffer columns
    raster_data[:, all_buffer_cols] = target_value

    print(f"Fixed {len(all_buffer_cols)} antimeridian buffer columns: "
          f"right-edge [{left_start}:{ncols}], left-edge [0:{right_end}]")
else:
    print("No antimeridian detected in this raster extent")
```

### Alternative: If using a special function

If you need to apply a custom function instead of directly assigning `target_value`:

```python
if dateline_columns.size > 0:
    # ... (same buffer column calculation as above) ...

    # Extract columns for processing
    columns_to_process = raster_data[:, all_buffer_cols]

    # Apply your special function here
    # For example: set to ocean value, or compute based on neighbors, etc.
    updated_columns = your_special_function(columns_to_process, all_buffer_cols, raster_data)

    # Write back the processed columns
    raster_data[:, all_buffer_cols] = updated_columns

    print(f"Processed {len(all_buffer_cols)} antimeridian buffer columns")
```

## Summary

The improved implementation will:
- ✅ Use the `iRaster_buffer_pixel` parameter effectively
- ✅ Select buffer columns on both sides of the antimeridian (wrap-around at edges)
- ✅ Handle boundary conditions safely with max/min checks
- ✅ Support user-defined processing of selected columns
- ✅ Maintain backward compatibility (default buffer=1)
- ✅ Handle regional rasters correctly (early exit if no dateline detected)
- ✅ Provide informative console output for debugging

### Key Behavior Summary

| Buffer Value | ncols=360 | Selected Columns | Total Count |
|--------------|-----------|------------------|-------------|
| 0 | 360 | None (empty) | 0 |
| 1 | 360 | [359, 0] | 2 |
| 5 | 360 | [355-359, 0-4] | 10 |
| 10 | 360 | [350-359, 0-9] | 20 |
| 180 | 360 | [180-359, 0-179] | 360 (all) |
