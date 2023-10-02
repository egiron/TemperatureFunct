# __main__.py

# import sys
from pathlib import Path
import click
# import tfunct

class PathType(click.Path):
    """A Click path argument that returns a pathlib Path, not a string"""

    def convert(self, value, param, ctx):
        return Path(super().convert(value, param, ctx))


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    click.echo(click.style('Welcome to Temperature response funtions library', fg='green'))
    # click.echo("{}".format(tfunct.__version__))
    pass

"""Read the Temperature response function example dataset"""
@cli.command()
@click.option('-l', '--length', type=int, help='Length of password to be generated')
# @click.option('-s', '--length', type=str, default="", help="")
@click.option('-o', '--option', type=click.Choice(['1', '2', '3', '4']), default = '4',
    help='''Options\n
    1 - alphabetic lowercase\n
    2 - alphabetic both cases\n
    3 - alphanumeric\n
    4 - alphanumeric + special characters'''
)
@click.option(
    "-i",
    "--input",
    type=PathType(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path to input nz-building-outlines.gpkg",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=PathType(file_okay=True, dir_okay=False, writable=True),
    help="Path to output Parquet file.",
    required=True,
)
def main(length: int, option: str, input: Path, output: Path):
    click.echo("CLI had not implemented yet...")
    # click.echo(click.style('ATTENTION!', blink=True))
    # click.echo(click.style('Some things', reverse=True, fg='cyan'))

def main_old():
    # If an site ID is given, then show the location
    # if len(sys.argv) > 1:
    #     location = iwin.get_location(sys.argv[1])
    #     viewer.show(location)
    # If no ID is given, then show a list of all location
    # else:
    #     site = iwin.get_site()
    #     titles = iwin.get_titles()
    #     viewer.show_list(site, titles)
    click.echo("CLI had not implemented yet")

if __name__ == "__main__":
    cli()