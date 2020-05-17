# charcoal demo

This directory (`demo/`) contains a quick and lightweight demonstration
of some of charcoal's features.

If you want to take a look at some output without running it yourself,
you can see the summary spreadsheet at
[../example-output/demo.combined_summary.csv](../example-output/demo.combined_summary.csv) and see two genome reports, one for [a LoombaR bin](../example-output/LoombaR_2017__SID1050_bax__bin.11.fa.gz.report.txt) and one for a [TARA bin from Delmont et al., 2019](../example-output/TARA_ANE_MAG_00014.fa.gz.report.txt).

To run the demo yourself, execute `charcoal run demo/demo.conf` from
the top-level directory, and then look at the output in `output.demo/`.

Briefly,

* the `LoombaR_2017__SID1050_bax__bin.11.fa.gz` could be identified
  against GTDB release 89 based on k-mers (76% were identifiable
  (`f_ident`), and of those, 95% (`f_major`) belonged to
  *s__Anaeromassilibacillus sp002159845*).
* for `TARA_ANE_MAG_00014.fa.gz`, we provided a lineage (*s__Salipiger
  thiooxidans*) in `provided-lineages.csv`, and contaminants were
  removed based on that lineage.
* no matches at all were found for `GCF_000005845-subset.fa.gz` in the
  GTDB database, so nothing is reported.
* too few of the k-mers in `TARA_ANE_MAG_00014.fa.gz` could be
  identified in the GTDB database (`f_ident < 0.10`), so nothing was
  done.
* of the identifiable k-mers in `TARA_PON_MAG_00084.fa.gz`, too few
  belonged to the major lineage (`f_major < 0.20`), so nothing was
  done.
