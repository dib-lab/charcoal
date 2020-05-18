# charcoal

Remove contaminated bits of genomes using k-mer based taxonomic analysis
with sourmash.

**Still early in development.** Buyer beware! Here be dragons!!

## Installing!

Clone this repository and change into the top-level repo directory.
The file `environment.yml` contains the necessary conda packages
(python and snakemake) to run charcoal; see the Quickstart section
for explicit instructions.

### Quickstart:

Clone the repository, change into it, create the environment, and activate it:

```
git clone https://github.com/dib-lab/charcoal
cd ./charcoal/
conda env create -f environment.yml -n charcoal
conda activate charcoal
```

### Running the demo!

To run, execute (in the top-level directory):

```
pip install -e .
charcoal run demo/demo.conf -j 8
```

You will end up with clean genomes in `output.demo/*.clean.fa.gz`, and
a summary of the demo genomes in `output.demo/combined_summary.csv`.

Once that succeeds, you can configure it yourself by copying
`demo/demo.conf` to a new file and editing it, or creating a new
configuration file with `charcoal init <project name>`.

## Explanation of output files.

All output files are in the output directory specified in the config
file.

A summary across all genomes will be placed in
`combined_summary.csv`.

For each genome, there are three key output files:
* `<filename>.report.txt` - a detailed report of the decontamination process.
* `<filename>.clean.fa.gz` - all of the kept ("clean") contigs.
* `<filename>.dirty.fa.gz` - all of the removed ("dirty") contigs.

## Need help?

There's [more documentation](doc/README.md) under the `doc/` directory.

Ask questions and file issues on [the GitHub issue tracker](https://github.com/dib-lab/charcoal/issues).

[@ctb](https://github.com/ctb/) [@taylorreiter](https://github.com/taylorreiter)
May 2020
