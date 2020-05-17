# charcoal documentation

## Introduction

## Getting help

Please file issues on
[the charcoal issue tracker](https://github.com/dib-lab/charcoal/issues).

## Installing charcoal

TODO: installable via pip or conda

## Running charcoal

Run charcoal like so:

```
charcoal run <config file>
```

A demo config file can be found in `demo/demo.conf`. Note that this uses
precomputed matches; see the Quickstart (@ctb doesn't exist yet!) for
downloading a complete bacterial/archaeal database.

An empty config file for a new project can be generated using
`charcoal init`.

Underneath, charcoal uses
[snakemake](https://snakemake.readthedocs.io/en/stable/) to run its
jobs, and `charcoal run` will pass on any parameters given after the
required config file.

### Running with multiple processes

We generally recommend running charcoal on a single computer using
multiple processes.  You can do this with `charcoal run <config file>
-j NUM`, where NUM is the number of processes to run.

### Needed resources

In general, each charcoal job needs less than 5 GB of memory, and
should take far less than 5 minutes per genome.  If you run with
multiple processes using `-j` (see above), the necessary memory
requirements will just multiply; e.g. request 80 GB of RAM if you are
running with `-j 16`.

charcoal's output files will use approximately the same amount of disk
space as the set of input genomes. charcoal compresses genomic output
(both cleaned and dirty) automatically using gzip.

### Cluster and cloud configuration

charcoal is lightweight lightweight and generally we do recommend running it on
a cluster; instead, run it with multiple processes using `-j` (see above).

However, if you want to run it across a cluster, you can do so!
charcoal uses snakemake underneath, so you can follow the
[snakemake cluster and cloud configuration instructions](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html).

## Configuring charcoal

@CTB: add gather-db and lineages configuration command.

You can generate a configuration file using `charcoal init <project>`:

```
charcoal init new-project
```
This will create a file `new-project.conf` with reasonable defaults.
You will need to fill in `genome_list` and `genome_dir` yourself.

If you have a directory with your genome bins in it already, you can
provide that to init with `--genome-dir`:
```
charcoal init new-project --genome-dir path/to/genomes/
```

You can check your configuration with `charcoal check <project>.conf`,
and show the aggregated configuration (defaults + system + project-specific
configuration) with `charcoal showconf <project>.conf`.

### Per-project configuration



### Installation-wide ccnfiguration

Each installation of charcoal (system-wide or in a conda
environment) has a configuration file that can provide
default settings for that installation. This is a good place to configure
the search database information.

Note that the installation-wide configuration is overridden by values
in the project configuration, so you can change any setting in the
project configuration. Any value not in the project config file will be
taken from the installation-wide configuration.

@ctb more - how to find it.

## Developer info

### Running charcoal from the development repo

You will need to run `pip install -e .` in the charcoal devleopment
repo in order to use the `charcoal` command from the development version.
