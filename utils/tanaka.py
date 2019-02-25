def rasterize_tanaka(filename, output_filename, hdf_output,resolution=128):
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
    import PIL
    import h5py
    
    # load the shapefile
    shp_df = gpd.read_file(filename)
    # load the unit type data
    unit_file = pkg_resources.resource_filename('cratertools',
                                                'data/tanaka_unittype.csv')
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

    # setup the output file
    meta = dict(driver='PNG',
                height=-1+(y[1]-y[0])/resolution,
                width=-1+(x[1]-x[0])/resolution,
                count=1,
                dtype=np.uint16,
                crs='+proj=latlong',
                transform=transform)

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
    PIL.Image.MAX_IMAGE_PIXELS = None
    with pd.HDFStore(hdf_output,'w') as h5f:
        h5f.append("table", df.reset_index())
        raster = pd.DataFrame(
                              np.array(
                                       imageio.imread('raster.png')
                                      )
                              )

        h5f.put("map", raster, compression='zlib')

    with pd.HDFStore(hdf_output, 'r') as h5f:
        tab = h5f["table"].sort_values('number')
        lmap = h5f["map"].values[::128, ::128]
        red, green, blue = (np.array(tab.red)/255,
                            np.array(t.green)/255,
                            np.array(t.blue)/255)
        im = (np.array([red[lmap.ravel()],
                       green[lmap.ravel()],
                       blue[lmap.ravel()]],
                       dtype=np.float32).
              reshape(3,
                      lmap.shape[0],
                      lmap.shape[1]).transpose([1, 2, 0]))
    plt.imshow(im)


import sys
import click

@click.command()
def main(input_filename, png_filename, hdf_filename, resolution=128):
    """Console script for cratertools."""
    click.echo("Rasterize Tanaka data")
    
    rasterize_tanaka(input_filename,
                     png_filename,
                     hdf_filename,
                     resolution=resolution):

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
