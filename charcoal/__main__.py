"Enable python -m charcoal"
import sys
import os
import subprocess
import glob

import click

def get_snakefile_path(name):
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, name)
    return snakefile


def get_package_configfile(filename):
    thisdir = os.path.dirname(__file__)
    configfile = os.path.join(thisdir, 'conf', filename)
    return configfile


def run_snakemake(configfile, no_use_conda=False, verbose=False,
                  snakefile_name='Snakefile', extra_args=[]):
    # find the Snakefile relative to package path
    snakefile = get_snakefile_path(snakefile_name)

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

    if configfile:
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

#
# actual command line functions
#

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
    run_snakemake(configfile, snakefile_name='Snakefile',
                  no_use_conda=no_use_conda, verbose=verbose,
                  extra_args=snakemake_args)

# download databases using a special Snakefile
@click.command()
def download_db():
    "download the necessary databases"
    run_snakemake(None, snakefile_name='Snakefile.download_db',
                  no_use_conda=True)

# 'check' command
@click.command()
@click.argument('configfile')
def check(configfile):
    "check configuration"
    run_snakemake(configfile, extra_args=['check'])

# 'showconf' command
@click.command()
@click.argument('configfile')
def showconf(configfile):
    "show full configuration across default, system and project config files"
    run_snakemake(configfile, extra_args=['showconf'])

# 'info' command
@click.command()
def info():
    "provide basic install/config file info"
    from .version import version
    print(f"""
This is charcoal version v{version}

Package install path: {os.path.dirname(__file__)}
Install-wide config file: {get_package_configfile('system.conf')}
snakemake Snakefile: {get_snakefile_path('Snakefile')}
""")

# 'init' command
@click.command()
@click.argument('configfile')
@click.option('--genome-dir', nargs=1)
@click.option('--lineages', nargs=1, default="")
@click.option('-f', '--force', is_flag=True)
def init(configfile, genome_dir, lineages, force):
    "create a new, empty config file."
    stubname = os.path.basename(configfile)
    if configfile.endswith('.conf'):
        stubname = stubname[:-5]
    else:
        configfile += '.conf'

    if os.path.exists(configfile) and not force:
        print(f"** ERROR: configfile '{configfile}' already exists.")
        return -1

    genome_list = f'{stubname}.genome-list.txt'
    if genome_dir:
        if os.path.exists(genome_list) and not force:
            print(f"** ERROR: genome list file '{genome_list}' already exists.")
            return -1
        genomes = glob.glob(f'{genome_dir}/*.fa')
        genomes += glob.glob(f'{genome_dir}/*.fna')
        genomes += glob.glob(f'{genome_dir}/*.fa.gz')
        genomes += glob.glob(f'{genome_dir}/*.fna.gz')
        print(f'found {len(genomes)} genomes in {genome_dir}/*.{{fa,fna}}{{,.gz}}')
        genomes = [ os.path.basename(g) for g in genomes ]
        with open(genome_list, 'wt') as fp:
            fp.write("\n".join(genomes))
        print(f"created '{genome_list}' with {len(genomes)} genomes in it.")

    if lineages:
        print(f"Using provided lineages from '{lineages}'")
    else:
        print("(No provided lineages file given.)")

    print(f"creating configfile '{configfile}' for project '{stubname}'")
    with open(configfile, 'wt') as fp:
        fp.write(\
f"""\
# location for all generated files
output_dir: output.{stubname}/

# list of genome filenames to decontaminate
genome_list: {stubname}.genome-list.txt

# directory in which genome filenames live
genome_dir: {genome_dir}

# (optional) list of lineages for input genomes. comment out or leave
# blank if none.
provided_lineages: {lineages}

# match_rank is the rank _above_ which cross-lineage matches are considered
# contamination. e.g. if set to 'superkingdom', then Archaeal matches in
# Bacterial genomes will be contamination, but nothing else.
#
# values can be superkingdom, phylum, class, order, family, or genus.
match_rank: order
""")

cli.add_command(run)
cli.add_command(check)
cli.add_command(showconf)
cli.add_command(info)
cli.add_command(init)
cli.add_command(download_db)

def main():
    cli()

if __name__ == '__main__':
    main()
