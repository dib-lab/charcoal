# charcoal documentation

## Introduction

charcoal is a pipeline that cleans genome assemblies of contaminated contigs.

charcoal is designed for decontaminating metagenome-assembled genomes (MAGs).

charcoal uses k-mer-based methods to identify potential contaminants. It
identifies contigs that are taxonomically inconsistent with the rest of the
genome, and removes them. It relies on reference databases of genomes with
a high quality taxonomy.  We currently recommend using the GTDB taxonomy,
and we provide databases for that.

charcoal uses relatively little memory (~2 GB per genome), takes less than
5 minutes per genome, and is fully parallelizable per genome. We've analyzed
1,000 genomes in well under an hour.

We are working on validating charcoal now.

charcoal is open source under the BSD 3-clause license, and is free for
use, reuse, modification, and remixing. Please enjoy responsibly!

## Authorship and Acknowledgements

charcoal development is led by Titus Brown and Taylor Reiter, and is
based on the [sourmash](http://sourmash.rtfd.io/) software. We would
especially like to thank Luiz Irber and Tessa Pierce for their
intellectual contributions to charcoal development.

## Getting help

We're happy to help you with any problems or questions you have!
Please file them as issues on
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

### Output files

All output files will be placed in the output directory specified in
the config file.

A summary across all genomes will be in `combined_summary.csv`.
The columns are explained below.


For each genome, there are three key output files:
* `<filename>.report.txt` - a detailed report of the decontamination process.
* `<filename>.clean.fa.gz` - all of the kept ("clean") contigs.
* `<filename>.dirty.fa.gz` - all of the removed ("dirty") contigs.

The columns in the combined summary are:

1. `genomefile` - the genome file name
2. `brieftax` - a brief form of the lineage used (provided lineage if given; or guessed taxonomy, if not)
3. `f_major` - the fraction of k-mers belonging to the majority lineage in the genome
4. `f_ident` - the fraction of k-mers for which a match in the database was found
5. `f_removed` - the fraction of total base pairs removed as contaminated (in contigs) and put in the `.dirty.fa.gz` file
6. `n_reason_1` - see report.txt
7. `n_reason_2` - see report.txt
8. `n_reason_3` - see report.txt
9. `refsize` - the approximate size of the nearest match genome in the database, if any
10. `ratio` - the ratio between the size of this genome and the refsize
11. `clean_bp` - total bp in "clean" contigs in the `.clean.fa.gz` file
12. `clean_n` - total number of "clean" contigs
13. `dirty_n` - total number of "dirty" contigs
14. `dirty_bp` - total bp in "dirty" contigs in the `.clean.fa.gz` file
15. `missed_n` - total number of contigs for which no hash values were found (these are ignored by charcoal, and placed in the clean contigs file)
16. `missed_bp` - total number of bp in contigs for which no hash values were found
17. `taxguessed` - the lineage guessed by charcoal based on majority taxonomy (see `f_ident` and `f_match`)
18. `taxprovided` - the lineage given in the provided-lineages file, if any.
19. `comment` - a comment explaining why this genome was not processed.

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

The project-specific settings are as follows:

* `output_dir`: the directory in which all of the output files will be placed. This will be created and populated by charcoal. This parameter is required.
* `genome_dir`: the directory in which all of the input genomes lie. Note, soft linked (`ln -s`) genome files are permitted. This parameter is required.
* `genome_list`: the basenames (e.g. `ls -1`) of the genome files to run decontamination on. They must all be located in `genome_dir`. FASTA, gzipped FASTA, or bzip2'ed FASTA are all supported, with no particular naming format required. This parameter is required.
* `provided_lineages`: an optional spreadsheet of `genome_filename,superkingdom,phylum,...` lines used to provide lineages for the input genomes. Here `genome_filename` must exactly match a filename in `genome_list`.  The provided lineage overrides the automatic lineage detection done by charcoal, and can be used to label genomes as lineages that do not belong to the sourmash databases specified in `gather_db`. For example, if you specify `genome.fa.gz,d__Eukaryota,` for a genome in this file while using the GTDB database for gather, then charcoal will remove bacterial and archaeal contaminants but not Eukaryotal.
* `match_rank`: rank at or below which a contigs are **not** removed as contaminants. For example, if `match_rank` is genus, then contigs in a genome file belonging to different families than the genome will be removed as contaminants. Defaults to order.

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

## Developer info

charcoal is developed collaboratively at https://github.com/dib-lab/charcoal/.

### Contributions and authorship

We welcome contributions - please feel free to open a Pull Request!

Both scientific and engineering contributions shall be considered for
authorship on software publications. This can include feature ideas,
debugging, validation approaches, and documentation updates.

### Running charcoal from the development repo

You will need to run `pip install -e .` in the charcoal development
repo in order to use the `charcoal` command from the development version.
