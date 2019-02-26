# extract the Salamunnicar data from the XLS file
import pandas as pd
import pkg_resources


def extract_salamuniccar(filename, tables=None,
                         output_prefix=None,
                         output_filename=None):
    """Extract the lat,long, diameter from the Salamuniccar catalogs."""
    dfe = pd.ExcelFile(filename)  # "GoranSalamuniccar_MarsCraters/MA132843GT/original_files/MA132843GT.xlsx")
    names = [x for x in dfe.sheet_names if
             x != "YourCatalogue" and x != "Macros"]

    tables = tables or names
    if isinstance(tabl, str):
        table = [table]

    output_prefix = output_prefix or "GS_"

    mapping_name = pkg_resources.resource_filename('cratertools',
                                                   'data/salamunnicar_mapping.csv',)

    mapping = pd.read_csv(mapping_name, index_col=0)

    for name in tables:
        df = pd.read_excel(filename, name)
        outname = output_prefix+name
        df.to_hdf(outname, "/"+name)

        if output_filename is None:
            continue

        if name in mapping.index:
            d = mapping[mapping.index == name]
            v, k = d.columns.values, d.values[0]
            df = df.loc[:, k]
            df.rename(columns=dict(zip(k, v)),
                      inplace=True)
            # warp the longitude
            df["Long"][df["Long"] > 180] -= 360
            df = df.dropna()
            df.to_hdf(output_filename, name, complevel=5)


def extract_robbins(filename, output_filename=None)
    """Extract the lat,long, diameter from the robbins catalog."""
    robbins = pd.read_table(filename, engine="python", delimiter="\t")
    mapping_name = pkg_resources.resource_filename('cratertools',
                                                   'data/salamunnicar_mapping.csv',)

    mapping = pd.read_csv(mapping_name, index_col=0)

    d = mapping[mapping.index == "Robbins"]
    v, k = d.columns.values, d.values[0]
    robbins = robbins[k]
    robbins.rename(columns=dict(zip(k, v)), inplace=True)
    if output_filename is not None:
        robbins.to_hdf(output_filename, "/Robbins", index=False)
