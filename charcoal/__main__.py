"Enable python -m charcoal"
import sys
import os
import click
import subprocess

def get_snakefile_path():
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'Snakefile')
    return snakefile


def get_package_configfile(filename):
    thisdir = os.path.dirname(__file__)
    configfile = os.path.join(thisdir, 'conf', filename)
    return configfile


@click.group()
def cli():
    pass

# create a run subcommand that by default passes all of its arguments
# on to snakemake (after setting Snakefile and config)
@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('configfile')
@click.option('--no-use-conda', is_flag=True, default=False)
@click.option('--verbose', is_flag=True)
@click.argument('snakemake_args', nargs=-1)
def run(configfile, snakemake_args, no_use_conda, verbose):
    "execute charcoal workflow (using snakemake underneath)"
    # find the Snakefile relative to package path
    snakefile = get_snakefile_path()

    # basic command
    cmd = ["snakemake", "-s", snakefile]

    # add --use-conda
    if not no_use_conda:
        cmd += ["--use-conda"]

    # add defaults and system config files, in that order
    cmd += ["--configfile", get_package_configfile("defaults.conf")]
    cmd += ["--configfile", get_package_configfile("system.conf")]

    # add explicit configfile from command line
    cmd += ["--configfile", configfile]
    # snakemake sometimes seems to want a default -j; set it to 1 for now.
    # can overridden later on command line.
    cmd += ["-j", "1"]

    # add rest of snakemake arguments
    cmd += list(snakemake_args)

    if verbose:
        print('final command:', cmd)

    # runme
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        return e.returncode

cli.add_command(run)

def main():
    cli()

if __name__ == '__main__':
    main()
