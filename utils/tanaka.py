# -*- coding: utf-8 -*-

"""Console script for cratertools."""
import sys
import click
import cratertools.utils.tanaka as tanaka
import cratertools.utils.salamuniccar as salamuniccar


@click.group()
@click.option('--debug/--no-debug', default=False)
def cli(debug):
    pass


@cli.command()
@click.argument("filename", type=str)
@click.argument("output_filename", type=str)
@click.argument("hdf_output", type=str)
@click.option("--scale",default=128)
def tanaka_rasterize(filename, output_filename, hdf_output, scale):
    tanaka.rasterize(filename, output_filename, hdf_output,resolution=1./scale, scale=scale)
    return 0


@cli.command()
@click.argument("filename", type=str)
@click.option("--tables", nargs=+, default=None)
@click.option("--output_prefix", type=str, default=None)
@click.option("--output_filename", type=str, default=None)
def salamuniccar(filename, tables, output_prefix, output_filename):
    salamuniccar.extract_salamuniccar(filename,
                                      tables=tables,
                                      output_prefix=output_prefix,
                                      output_filename=output_filename)
    return 0


@cli.command()
@click.argument("filename", type=str)
@click.option("--output_filename", type=str, default=None)
def robbins(filename, output_filename):
    salamuniccar.extract_robbins(filename, output_filename)
    return 0


if __name__ == "__main__":
    sys.exit(cli())  # pragma: no cover
