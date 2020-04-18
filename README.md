# charcoal

Remove contaminated bits of genomes using "togetherness" correlations
across metagenomes, combined with k-mer based taxonomic
analysis.

**Still early in development.** Buyer beware! Here be dragons!!

## Installing!

Clone this repository and change into the top-level repo directory.
The file `environment.yml` contains the necessary conda packages
(python and snakemake) to run charcoal; see the Quickstart section
for explicit instructions.

### Quickstart:

Clone the repository, change into it, create the environment, and activate it:

```
git clone https://github.com/ctb/charcoal
cd ./charcoal/
conda env create -f environment.yml -n charcoal
conda activate charcoal
```

## Running!

To run, execute (in the top-level directory):

```
snakemake --use-conda -j 8 --configfile test-data/conf-test.yml
```

Once that works, you can configure it yourself by copying
`test-data/conf-test.yml` to a new file and editing it.

A few important notes --

* metagenome signatures need to be calculated with the following parameters:
   `--track-abundance --scaled=1000 -k 31`
   
## Explanation of output files.

TBW.
   
## Need help?

Ask questions and file issues on [the GitHub issue tracker](https://github.com/dib-lab/charcoal/issues).

[@ctb](https://github.com/ctb/) [@taylorreiter](https://github.com/taylorreiter)
March 2020
