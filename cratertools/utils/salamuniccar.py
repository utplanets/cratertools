# extract the Salamunnicar data from the XLS file
import pandas as pd
import pkg_resources
import logging
import os

def extract_salamuniccar(filename, tables=None,
                         output_prefix=None,
                         output_filename=None):
    """Extract the lat,long, diameter from the Salamuniccar catalogs."""
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    logger.info("Reading Excel file")
    logger.info(output_filename)
    dfe = pd.ExcelFile(filename)
    names = [x for x in dfe.sheet_names if
             x != "YourCatalogue" and x != "Macros"]

    tables = tables or names
    if isinstance(tables, str):
        tables = [tables]

    output_prefix = output_prefix or "GS_"

    mapping_name = pkg_resources.resource_filename('cratertools',
                                                   'data/salamuniccar_mapping.csv',)

    mapping = pd.read_csv(mapping_name, index_col=0)

    for name in tables:
        logger.info("Processing table : {}".format(name))
        df = pd.read_excel(filename, name)
        outname = output_prefix+name
        df.to_hdf(outname, "/"+name)

        if output_filename is None:
            continue
        print(name, mapping.index)
        if name in mapping.index:
            d = mapping[mapping.index == name]
            v, k = d.columns.values, d.values[0]
            df = df.loc[:, k]
            df.rename(columns=dict(zip(k, v)),
                      inplace=True)
            # warp the longitude
            df["Long"][df["Long"] > 180] -= 360
            df = df.dropna()
            df.to_hdf(output_filename, name,
                      append=os.path.exists(output_filename), complevel=5)


def extract_robbins(filename, output_filename=None):
    """Extract the lat,long, diameter from the robbins catalog."""
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    logger.info("Reading Robbins data")
    robbins = pd.read_table(filename, engine="python", delimiter="\t")
    mapping_name = pkg_resources.resource_filename('cratertools',
                                                   'data/salamuniccar_mapping.csv',)

    mapping = pd.read_csv(mapping_name, index_col=0)

    d = mapping[mapping.index == "Robbins"]
    v, k = d.columns.values, d.values[0]
    robbins = robbins[k]
    robbins.rename(columns=dict(zip(k, v)), inplace=True)
    if output_filename is not None:
        robbins.to_hdf(output_filename, "/Robbins",
                       append=os.path.exists(output_filename), index=False)
