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
pip install -e .
```

### Run the demo! (~2 minutes)

To run, execute (in the top-level directory):

```
python -m charcoal run demo/demo.conf -j 4
```

This will create two summary files in `output.demo/`,
`genome_summary.csv` and `hit_list_for_filtering.csv`. You can open
these in your favorite spreadsheet program.

For a friendlier summary, run:
```
python -m charcoal run demo/demo.conf -j 4 report
```

This will create a directory `output.demo/report/` that contains an
index page, `index.html`, that summarizes the charcoal run. This directory
also contain individual genome reports that you can reach through links
in the index.

Finally, you can run
```
python -m charcoal run demo/demo.conf -j 4 clean
```
and this will produce "cleaned" genomes based on the information in
`output.demo/hit_list_for_filtering.csv`.

### Do a full configure & run! (~10 minutes)

Now that the demo runs, you've got everything working! Hooray!
Now let's see how to set up a real run, with real databases!

This will take under 10 minutes and under 2 GB of disk space. You'll
need about 8 GB of RAM (change `-j 4` to `-j 1`, below, to run it in 2
GB of RAM, albeit 4x slower).

We'll use a set of 10 genomes taken from
[Nitrogen-fixing populations of Planctomycetes and Proteobacteria are abundant in surface ocean metagenomes, Delmont et al., 2018](https://www.nature.com/articles/s41564-018-0176-9). These
10 genomes have three eukaryotic and seven bacterial bins.

#### Install the database.

First, install the sourmash database for GTDB. 

```
charcoal download-db
```

This will put two files in the `db/` directory, totalling 1.5 GB. (You
can run this command multiple times and it should only download the
databases once.)

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
python -m charcoal run newproject.conf -n
```

#### Decontaminate!

And, finally, run the first round of analysis! This will run four
processes in parallel (`-j 4`)
```
python -m charcoal run newproject.conf -j 4
```

#### Examine the results

**out of date**

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
