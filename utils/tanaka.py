# -*- coding: utf-8 -*-

"""Console script for cratertools."""
import sys
import click
import cratertools.utils.tanaka as tanaka

@click.group()
@click.option('--debug/--no-debug', default=False)
def cli(debug):
    pass

@cli.command()
@click.argument("filename", type=str)
@click.argument("output_filename", type=str)
@click.argument("hdf_output", type=str)
@click.option("--scale",default=128)
def rasterize(filename, output_filename, hdf_output, scale):
    tanaka.rasterize(filename, output_filename, hdf_output,resolution=1./scale, scale=scale)
    return 0


if __name__ == "__main__":
    sys.exit(cli())  # pragma: no cover
