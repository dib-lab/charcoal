# charcoal

Remove contaminated bits of genomes using "togetherness" correlations
across metagenomes, combined with k-mer based taxonomic
analysis.

**Still early in development.** Buyer beware! Here be dragons!!

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

Ask questions and file issues on [the GitHub issue tracker](https://github.com/dib-lab/charcoal/issues)

@ctb @taylorreiter
March 2020
