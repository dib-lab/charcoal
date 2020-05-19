# charcoal

Remove contaminated bits of genomes using k-mer based taxonomic analysis
with sourmash.

**Still early in development.** Buyer beware! Here be dragons!!

## Installing!

In brief: clone this repository and change into the top-level repo
directory.  The file `environment.yml` contains the necessary conda
packages (python and snakemake) to run charcoal; see the Quickstart
section for explicit instructions.

### Quickstart:

Clone the repository, change into it, create the environment, and activate it:

```
git clone https://github.com/dib-lab/charcoal
cd ./charcoal/
conda env create -f environment.yml -n charcoal
conda activate charcoal
```

### Run the demo! (~1 minute)

To run, execute (in the top-level directory):

```
pip install -e .
charcoal run demo/demo.conf -j 8
```

You will end up with clean genomes in `output.demo/*.clean.fa.gz`, and
a summary of the demo genomes in `output.demo/combined_summary.csv`.

### Do a full configure & run! (~10 minutes)

Now that the demo runs, you've got all the software installed! Hooray!
Now let's see how to set up a real run, with real databases!

We'll use a set of 10 genomes taken from
[Nitrogen-fixing populations of Planctomycetes and Proteobacteria are abundant in surface ocean metagenomes, Delmont et al., 2018](https://www.nature.com/articles/s41564-018-0176-9). These
10 genomes have three eukaryotic and seven bacterial bins.

#### Install the database.

First, install the sourmash database for GTDB. 

```
charcoal download-db
```

This will put two files in the `db/` directory. You can run this
command multiple times and it should only download the databases once.

#### Download the example genomes

Next, download and unpack the example genomes:

```
curl -L https://osf.io/5pej8/download > example-genomes.tar.gz
tar xzf example-genomes.tar.gz
ls example-genomes/
```
The `example-genomes/` directory should have 10 genomes in it. It also
has a file `provided-lineages.csv` which labels the three eukaryotic
genomes as `d__Eukaryota`. This is needed because eukaryotes cannot be
automatically classified by charcoal, but they can be decontaminated.

#### Initiate a new project

Next, create a new project configuration:
```
charcoal init newproject --genome-dir example-genomes \
    --lineages example-genomes/provided-lineages.csv
```

This creates two files, `newproject.genome-list.txt` and
`newproject.conf`. The genome-list file contains the names of all
genome files, and the `newproject.conf` file contains the
configuration options for charcoal.

To do a "dry run" of charcoal, which lists out the jobs that will be
run, execute:
```
charcoal run newproject.conf -n
```

#### Decontaminate!

And, finally, run the decontamination routine. This will run four
processes in parallel (`-j 4`)
```
charcoal run newproject.conf -j 4
```

#### Examine the results

The results will be in `output.newproject/`; see the file
`combined_summary.csv`, as well as the `*.clean.fa.gz` files, which
contain the cleaned contigs. You might also take a look at the
`*.report.txt` files which contain individual genome cleaning reports.

For one example, the summary spreadsheet shows that approximately 10%
of `TARA_PSW_MAG_00136.fa` was removed (column `f_removed`), and the
report in `output.newproject/TARA_PSW_MAG_00136.fa.report.txt` shows
that contigs were removed for being members of a variety of different
bacterial lineages.

## Need help?

There's [more documentation](doc/README.md) under the `doc/` directory.

Please ask questions and file issues on [the GitHub issue tracker](https://github.com/dib-lab/charcoal/issues)!

[@ctb](https://github.com/ctb/) [@taylorreiter](https://github.com/taylorreiter)
May 2020
