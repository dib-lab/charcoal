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


def run_snakemake(configfile, no_use_conda=False, verbose=False,
                  extra_args=[]):
    # find the Snakefile relative to package path
    snakefile = get_snakefile_path()

    # basic command
    cmd = ["snakemake", "-s", snakefile]

    # add --use-conda
    if not no_use_conda:
        cmd += ["--use-conda"]

    # snakemake sometimes seems to want a default -j; set it to 1 for now.
    # can overridden later on command line.
    cmd += ["-j", "1"]

    # add rest of snakemake arguments
    cmd += list(extra_args)

    # add defaults and system config files, in that order
    configfiles = [get_package_configfile("defaults.conf"),
                   get_package_configfile("system.conf"),
                   configfile]
    cmd += ["--configfile"] + configfiles

    if verbose:
        print('final command:', cmd)

    # runme
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        return e.returncode


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
    run_snakemake(configfile, no_use_conda, verbose, snakemake_args)

@click.command()
@click.argument('configfile')
def check(configfile):
    "check configuration"
    run_snakemake(configfile, extra_args=['check'])

cli.add_command(run)
cli.add_command(check)


def main():
    cli()

if __name__ == '__main__':
    main()
