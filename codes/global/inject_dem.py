import os, sys, platform

sPlatform_os = platform.system()
from pyearth.toolbox.conversion.convert_vector_format import convert_vector_format
# Get the directory of the current script

sPath_library = '/qfs/people/liao313/workspace/python/uraster'
sys.path.append(sPath_library)
from uraster.classes.uraster import uraster
# Download input data using Pooch (downloads to system cache)

sFolder_output = "/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20260203001"

sFilename_source_mesh = os.path.join(sFolder_output, "mpas.geojson")
sFilename_target_mesh =  os.path.join(sFolder_output, "uraster.geojson")
sFilename_mesh_png = os.path.join(sFolder_output, "mesh.png")
sFilename_variable_png = os.path.join(sFolder_output, "variable.png")
sFilename_variable_animation = os.path.join(sFolder_output, "variable.mp4")

sFolder_gebco = "/compyfs/liao313/00raw/dem/global/"

def main():
    aConfig = dict()
    aConfig["sFilename_source_mesh"] = (
        sFilename_source_mesh  # use the L10-100 test mesh
    )
    aFilename_source_raster = []
    sFilename_dem0 = os.path.join(sFolder_gebco, "gebco_2025_n0.0_s-90.0_w0.0_e90.0.tif")
    sFilename_dem1 = os.path.join(sFolder_gebco, "gebco_2025_n0.0_s-90.0_w-180.0_e-90.0.tif")
    sFilename_dem2 = os.path.join(sFolder_gebco, "gebco_2025_n0.0_s-90.0_w-90.0_e0.0.tif")
    sFilename_dem3 = os.path.join(sFolder_gebco, "gebco_2025_n0.0_s-90.0_w90.0_e180.0.tif")
    sFilename_dem4 = os.path.join(sFolder_gebco, "gebco_2025_n90.0_s0.0_w-180.0_e-90.0.tif")
    sFilename_dem5 = os.path.join(sFolder_gebco, "gebco_2025_n90.0_s0.0_w-90.0_e0.0.tif")
    sFilename_dem6 = os.path.join(sFolder_gebco, "gebco_2025_n90.0_s0.0_w0.0_e90.0.tif")
    sFilename_dem7 = os.path.join(sFolder_gebco, "gebco_2025_n90.0_s0.0_w90.0_e180.0.tif")
    aFilename_source_raster.append(sFilename_dem0)  # global dem from gebco
    aFilename_source_raster.append(sFilename_dem1)  # global dem from gebco
    aFilename_source_raster.append(sFilename_dem2)  # global dem from gebco
    aFilename_source_raster.append(sFilename_dem3)  # global dem from gebco
    aFilename_source_raster.append(sFilename_dem4)  # global dem from gebco
    aFilename_source_raster.append(sFilename_dem5)  # global dem from gebco
    aFilename_source_raster.append(sFilename_dem6)  # global dem from gebco
    aFilename_source_raster.append(sFilename_dem7)  # global dem from gebco

    aConfig["aFilename_source_raster"] = aFilename_source_raster
    aConfig["sFilename_target_mesh"] = sFilename_target_mesh
    pRaster = uraster(aConfig)

    if not os.path.exists(pRaster.sFilename_target_mesh):
        print(f"Target mesh file does not exist: {pRaster.sFilename_target_mesh}")


    sFilename_source_mesh_parquet = os.path.join(sFolder_output, "mpas_source.parquet"
    )
    # convert_vector_format(pRaster.sFilename_source_mesh, sFilename_source_mesh_parquet)

    pRaster.setup(iFlag_verbose_in=True)
    # pRaster.report_inputs()
    dLongitude_focus_in = -112.033964
    dLatitude_focus_in = 43.491977
    #pRaster.visualize_source_mesh(
    #    sFilename_out=sFilename_mesh_png,
    #    dLongitude_focus_in=dLongitude_focus_in,
    #    dLatitude_focus_in=dLatitude_focus_in,
    #    dImage_scale_in=10,
    #    iFlag_show_graticule=False,
    #    iFlag_wireframe_only=True,
    #)

    pRaster.run_remap(iFlag_verbose_in=True)
    sFilename_mesh_parquet = os.path.join(sFolder_output, "mpas_uraster.parquet"
    )
    convert_vector_format(pRaster.sFilename_target_mesh, sFilename_mesh_parquet)

    # pRaster.report_outputs()
    sColormap = "terrain"

    #pRaster.visualize_target_mesh(
    #    sFilename_out=sFilename_variable_png,
    #    dLongitude_focus_in=dLongitude_focus_in,
    #    dLatitude_focus_in=dLatitude_focus_in,
    #    dImage_scale_in=10,
    #    sColormap=sColormap,
    #    iFlag_show_graticule=False,
    #)

    # pRaster.visualize_target_mesh(
    #    sFilename_out=sFilename_variable_animation,
    #    sColormap=sColormap,
    #    dLongitude_focus_in=dLongitude_focus_in,
    #    dLatitude_focus_in=dLatitude_focus_in,
    #    iFlag_create_animation=True,
    #    iAnimation_frames=360,       # 1° longitude per frame
    #    sAnimation_format='mp4')

    pRaster.cleanup()

    print("done")


if __name__ == "__main__":
    main()
