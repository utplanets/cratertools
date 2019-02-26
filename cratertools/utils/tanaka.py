import logging
import sys


def rasterize(filename, output_filename, hdf_output,resolution=1./128,scale=128):
    # setup for geopython
    import numpy as np
    import geopandas as gpd
    import pandas as pd
    import io
    import rasterio
    from rasterio import features
    from rasterio.transform import from_origin
    import pkg_resources
    import imageio
    from PIL import Image
    import h5py
    
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)#    logger.setLevel(logging.INFO)

    # load the shapefile
    logger.info("Reading Shapefile")
    shp_df = gpd.read_file(filename)

    # load the unit type data
    logger.info("Reading unit types")
    unit_file = pkg_resources.resource_filename('cratertools',
                                                'data/tanaka_unittype.csv')
    logger.info("Merging data")
    df = pd.read_csv(unit_file, index_col=0)
    #
    index = dict(zip(*np.unique(df.index, return_index=True)))

    df = pd.merge(shp_df[["Unit", "UnitDesc"]].groupby("Unit")
                  .agg(lambda x: x.iloc[0]),
                  df,
                  left_index=True,
                  right_index=True).reset_index()

    df["number"] = [index[c] for c in df['index']]  # map numbers
    df = (df.set_index(["unit", "index"])
            .sort_index()
            .sort_values("number"))

    # setup the transform function
    x = [-180, 180]  # -180+res/2, 180-res/2
    y = [-90, 90]  # -89.5, 89.5]
    transform = from_origin(x[0], y[-1], resolution, resolution)
    logger.info("Setting up rasterization {}".format(resolution))

    # setup the output file
    meta = dict(driver='PNG',
                height=-1+(y[1]-y[0])/resolution,
                width=-1+(x[1]-x[0])/resolution,
                count=1,
                dtype=np.uint16,
                crs='+proj=latlong',
                transform=transform)

    logger.info("Rasterizing")
    with rasterio.open(output_filename, 'w', **meta) as out:
        out_arr = out.read(1)
        unit_num = [index[x] for x in shp_df.Unit]

        # this is where we create a generator of geom
        # value pairs to use in rasterizing
        shapes = ((geom, value) for geom, value in zip(shp_df.geometry,
                                                       unit_num))
        burned = features.rasterize(shapes=shapes,
                                    fill=0,
                                    out=out_arr,
                                    transform=out.transform)
        out.write_band(1, burned)

    # Now read that file and save to HDF
    Image.MAX_IMAGE_PIXELS = None
    logger.info("Saving HDF5 file")
    with pd.HDFStore(hdf_output,'w') as h5f:
        h5f.append("table", df.reset_index())
        raster = pd.DataFrame(
                              np.array(
                                       imageio.imread(output_filename)
                                      , dtype=np.int8)
                              )

        h5f.put("map", raster, compression='zlib')
    logger.info("Reading HDF5 file")
    with pd.HDFStore(hdf_output, 'r') as h5f:
        tab = h5f["table"].sort_values('number')
        lmap = h5f["map"].values[::scale, ::scale]
        red, green, blue = (np.array(tab.red)/255,
                            np.array(tab.green)/255,
                            np.array(tab.blue)/255)
        im = (np.array([red[lmap.ravel()],
                        green[lmap.ravel()],
                        blue[lmap.ravel()]],
                       dtype=np.float32).
              reshape(3,
                      lmap.shape[0],
                      lmap.shape[1]).transpose([1, 2, 0]))
    logger.info("Drawing rastered image")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.imshow(im)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.savefig("img.png")

#@click.command()
#@click.argument("input_filename", nargs=1)
#@click.argument("png_filename",nargs=1)
#@click.argument("hdf_filename",nargs=1)
#@click.option("--resolution",default=1./128)
#def main(input_filename, png_filename, hdf_filename, resolution=1./128):
#    """Console script for cratertools."""
#    
#    rasterize_tanaka(input_filename,
#                     png_filename,
#                     hdf_filename,
#                     resolution=resolution)
#
#if __name__ == "__main__":
#    sys.exit(main(*sys.argv[1:]))  # pragma: no cover
