# charcoal documentation

May 2020

Contents:
* [Introduction](#introduction)
* [How charcoal finds contamination](#how-charcoal-finds-contamination)
* [Authorship and Acknowledgements](#authorship-and-acknowledgements)
* [Getting help](#getting-help)
* [Installing charcoal](#installing-charcoal)
* [Running charcoal](#running-charcoal)
* [Configuring charcoal](#configuring-charcoal)
* [Frequently Asked Questions](#frequently-asked-questions)
* [Developing charcoal](#developing-charcoal)

## Introduction

charcoal is a pipeline that removes contaminant contigs from genome
assemblies.

charcoal is designed for de-contaminating metagenome-assembled genomes (MAGs).
It is focused on removing bacterial and archaeal contaminants for now.

charcoal uses k-mer-based methods to identify potential
contaminants. Contigs that are taxonomically inconsistent with the
rest of the genome are identified and then removed.

charcoal relies on reference databases of genomes with a high quality
taxonomy to flag likely contamination.  We provide a database and set
of lineages from the [GTDB taxonomy](https://gtdb.ecogenomic.org/) for
this purpose.

charcoal uses relatively little memory (~2 GB per genome), takes less than
5 minutes per genome, and is fully parallelizable per genome. We've analyzed
1,000 genomes in well under an hour.

We are working on validating charcoal now (May 2020).

charcoal is open source under the BSD 3-clause license, and is free for
use, reuse, modification, and remixing. Please enjoy responsibly!

We provide a [Code of Conduct](CODE_OF_CONDUCT.rst) that we expect
developers, contributors, and users to follow when requesting help,
filing issues, writing documentation, and engaging in the project.

## How charcoal finds contamination

charcoal assigns a lineage to each input genome based on the search
database(s), as well as any lineages that are provided.

charcoal then uses a relatively straightforward heuristic to determine
if a contig is a contaminant: it searches 25,000 GTDB genomes using
`sourmash search --containment` to see where the contig matches best
against known genomes. If the contig places within a lineage that has
a different lineage than the containing genome, it is flagged as a
contaminant.

Importantly, charcoal defaults to assuming that a contig is "clean" -
if it has no information on a contig, it does not remove it.

So, charcoal will fail to detect contamination in very short contigs
for which no k-mers are chosen, as well as contigs that are completely
novel in their DNA content.

And, of course, charcoal is database dependent. So if genomes or
databases in the GTDB 25k collection are contaminated, charcoal will
not be able to detect those contaminants.

## Authorship and Acknowledgements

charcoal development is led by Titus Brown and Taylor Reiter, and is
heavily based on the [sourmash](http://sourmash.rtfd.io/) software. We
would especially like to thank Luiz Irber and Tessa Pierce for their
intellectual contributions to charcoal development!

Erich Schwarz demonstrated the utility of a k-mer based approach to
eukaryotic genome assembly decontamination a while back, and [contributed
significantly](https://github.com/dib-lab/sourmash/issues/940) to pushing
this forward in sourmash.

Donovan Parks suggested using sourmash to develop systematic
decontamination approaches.

Boris Vinatzer, Lenwood Heath, Leighton Pritchard, and Jonathan Eisen
also served as a useful sounding board in developing an understanding
of the taxonomic issues at hand.

We would especially like to thank the
[GTDB project](https://gtdb.ecogenomic.org/) for all their hard work
and their software and databases!

charcoal development is funded by the Moore Foundation through grant
GBMF4551 to C. Titus Brown. The codebase is Copyright 2020, Regents of
the University of California (as of May 2020).

This initial version of charcoal is a product of the [Lab for Data-
Intensive Biology](http://ivory.idyll.org/lab/).

## Getting help

We're happy to help you with any problems or questions you have!
Please file them as issues on
[the charcoal issue tracker](https://github.com/dib-lab/charcoal/issues).

## Installing charcoal

Please follow the installation instructions in
[the project README](https://github.com/dib-lab/charcoal/blob/master/README.md).

TODO: installable via pip or conda

## Running charcoal

Run charcoal like so:

```
charcoal run <config file>
```

A demo config file can be found in `demo/demo.conf`. Note that this
uses precomputed matches; see
[the Quickstart docs](https://github.com/dib-lab/charcoal/blob/master/README.md)
for downloading a complete bacterial/archaeal database.

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

### Output files

(UNDER DEVELOPMENT - FIX ME)

* discuss hit list
* discuss editing hit list
* output files etc

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
You can also pass in a provided lineages file with `--lineage`.

You can check your configuration with `charcoal check <project>.conf`,
and show the aggregated configuration (defaults + system + project-specific
configuration) with `charcoal showconf <project>.conf`.

### Per-project configuration

The project-specific settings are as follows:

* `output_dir`: the directory in which all of the output files will be placed. This will be created and populated by charcoal. This parameter is required.
* `genome_dir`: the directory in which all of the input genomes lie. Note, soft linked (`ln -s`) genome files are permitted. This parameter is required.
* `genome_list`: the basenames (e.g. `ls -1`) of the genome files to run decontamination on. They must all be located in `genome_dir`. FASTA, gzipped FASTA, or bzip2'ed FASTA are all supported, with no particular naming format required. This parameter is required.
* `provided_lineages`: an optional spreadsheet of `genome_filename,superkingdom,phylum,...` lines used to provide lineages for the input genomes. Here `genome_filename` must exactly match a filename in `genome_list`.  The provided lineage overrides the automatic lineage detection done by charcoal, and can be used to label genomes as lineages that do not belong to the sourmash databases specified in `gather_db`. For example, if you specify `genome.fa.gz,d__Eukaryota,` for a genome in this file while using the GTDB database for gather, then charcoal will remove bacterial and archaeal contaminants but not Eukaryotal.

Database parameters (for intermediate users):

* `gather_db`: a list of sourmash databases (SBT, LCA, or collections of signatures) against which to "collect" relevant genomic matches to the query genome. See [the sourmash documentation](http://sourmash.rtfd.io/) for more information. By default, we suggest using the GTDB .sbt.zip database here, as it is low memory and quick to search. Custom databases are completely supported as long as you supply an accompanying set of lineages.
* `lineages_csv`: a lineage spreadsheet (see `sourmash lca index` documentation in [the sourmash docs](http://sourmash.rtfd.io/)) specifying a mapping from identifiers to a fully resolved lineage. Any taxonomy can be used for this lineage, including NCBI or GTDB taxonomies; you probably shouldn't mix them though.
* `scaled`: the scaled resolution at which you want to detect contamination. This must be no smaller than the scaled parameter of the sourmash database(s) listed in `gather_db`.
* `ksize`: the k-mer size at which you want to detect contamination. This must be matched by the k-mer size of the sourmash database(s) listed in `gather_db`.

Other settings:
* `strict` (0 or 1, default 1) -- check and validate config settings & filenames strictly.
* `force` (0 or 1, default 0) -- continue past survive-able errors in decontamination.

### Rerunning charcoal with different parameters

Unless you change the database and lineage spreadsheet, the ksize, or the
scaled, you can rerun the filtering with different match ranks and provided
lineages.  To do this, remove `*.clean.fa.gz` in the output directory.

### Installation-wide configuration

Each installation of charcoal (system-wide or in a conda
environment) has a configuration file that can provide
default settings for that installation. This is a good place to configure
the search database information.

Note that the installation-wide configuration is overridden by values
in the project configuration, so you can change any setting in the
project configuration. Any value not in the project config file will be
taken from the installation-wide configuration.

You can find the install-wide config file location by running `charcoal info`.

### Downloading databases

The command
```
charcoal download-db
```
will download a sourmash database containing 25k genomes from GTDB, along
with the GTDB lineage spreadsheet for those genomes. These will take up about
1.5GB and will be placed in the `db/` subdirectory of the current working
directory.

## Frequently Asked Questions

#### Do I need to use the GTDB taxonomy with charcoal?

charcoal uses DNA similarity to identify potential contamination, and
relies on taxonomic annotations to decide if segments of shared DNA
are contaminants are not. We think GTDB is the best choice for this
purpose because their species-level clusters are based largely
entirely on whole-genome similarity (Average Nucleotide Identity).

However, outside of charcoal, there is no need to use GTDB for your genomes.
charcoal will automatically classify your genomes for you and do the
contamination analysis based on its own classification. The only
time you need to work with the GTDB taxonomy is when you override this
classification using the provided-lineages file.

#### How do I use my own classification databases?

You can absolutely use your own classification databases! You'll need
to provide one or more collections of signatures calculated by
sourmash (SBT or LCA databases, or multiple signatures), along with a
lineage spreadsheet that connects sequence identifiers to taxonomy.
The lineages spreadsheet should be in the form used as input by
`sourmash lca index`.

If you're adding to the default database used by charcoal, you'll need
to provide your lineages in a GTDB-compatible way
(e.g. `d__Eukaryota`).  But you can also easily use (e.g.) the NCBI
taxonomy if you want; there are scripts available through
[the sourmash project](https://github.com/dib-lab/sourmash) to help
you with this, please ask
[on the issue tracker](https://github.com/dib-lab/sourmash/issues).

One warning: we do not yet have a simple way to evaluate the impact of
confused taxonomies on charcoal's performance, and it's not
immediately clear what would happen if a "bad" set of genomes or
confounded set of taxonomies are provided to charcoal. See
[charcoal issue #77](https://github.com/dib-lab/charcoal/issues/77).

#### How do I decontaminate genomes that charcoal can't classify?

With the default GTDB database, charcoal cannot automatically classify
any eukaryote genomes, and may also miss unknown archaea and bacteria.
However, this doesn't mean that charcoal cannot decontaminate these
genomes.  You can provide your own lineages to charcoal and it will
happily remove sequences that _disagree_ with those lineages -- see the
`provided_lineages` option in the config file! charcoal will put
unidentified sequences in the clean genome bin.

## Developing charcoal

charcoal is developed collaboratively at https://github.com/dib-lab/charcoal/.

### Contributions and authorship

We welcome contributions - please feel free to open a Pull Request or [file an
issue](https://github.com/dib-lab/charcoal/issues) to get started!

Scientific, engineering and community contributors to charcoal are
eligible for authorship on publications about the charcoal
software. This can include (but is not limited to) documentation
updates, detailed bug reporting, performance and benchmark reporting,
validation and evaluation, and code contributions.

### Our philosophy

We adhere to the following guiding principles in developing charcoal:

* we default to high specificity when removing contamination.
* decisions made by charcoal will be clearly documented and supported by reporting.
* decisions made by the software and the developers will be made transparently and documented via issues and pull requests.
* the software will be documented and supported.
* as of a v1.0 release, we will use [semantic versioning](https://semver.org/) to support stable use of charcoal in pipelines.
* the developers of charcoal will strive to ["eat our own dogfood"](https://en.wikipedia.org/wiki/Eating_your_own_dog_food) - we will document how and why we use it ourselves.
* charcoal should work well with other programs and be a good member of the bioinformatics ecosystem.
* the charcoal software should be well tested and stable.

### Running charcoal from the development repo

You will need to run `pip install -e .` in the charcoal development
repo in order to use the `charcoal` command from the development version.
