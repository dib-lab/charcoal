#
# run with --use-conda for maximal froodiness.
#
import csv, sys

# override this with --configfile on command line
configfile: 'test-data/conf-test.yml'

strict_val = config.get('strict', '1')
strict_mode = int(strict_val)
if not strict_mode:
    print('** WARNING: strict mode is OFF. Config errors will not force exit.')

### config stuff loaded from config file
genome_list_file = config['genome_list']
genome_list = [ line.strip() for line in open(genome_list_file, 'rt') ]
genome_list = [ line for line in genome_list if line ]   # remove empty lines

genome_dir = config['genome_dir'].rstrip('/')
output_dir = config['output_dir'].rstrip('/')

# read in provided lineages, if any.
provided_lineages_file = config.get('provided_lineages', '')
provided_lineages = {}
if provided_lineages_file:
    with open(provided_lineages_file, 'rt') as fp:
        r = csv.reader(fp)
        for row in r:
            genome_filename = row[0]
            if genome_filename not in genome_list:
                print(f'** WARNING: lineage was provided for unknown genome {genome_filename}')
                print(f'** in provided lineages file {provided_lineages_file}')
                print(f'** ({genome_filename} not in {genome_list_file})')
                if strict_mode:
                    sys.exit(-1)
            provided_lineages[genome_filename] = row[1:]

    print(f'** read {len(provided_lineages)} provided lineages')

lca_db = config['lca_db']

### utility functions
def output_files(filename_template, **kw):
    return expand(output_dir + filename_template, **kw)

def get_provided_lineage(w):
    "retrieve a lineage for this filename from provided_lineages dictionary"
    filename = w.f
    if filename in provided_lineages:
        lineage = provided_lineages[filename]
        lineage = ";".join(lineage)
        return lineage
    else:
        return "NA"

### rules!

wildcard_constraints:
    size="\d+"

rule all:
    input:
        expand(output_dir + '/{g}.clean.fa.gz', g=genome_list),
        output_dir + '/just_taxonomy.combined_summary.csv',

rule contigs_sig:
    input:
        genome_dir + '/{filename}'
    output:
        output_dir + '/{filename}.sig'
    conda: 'conf/env-sourmash.yml'
    params:
        scaled = config['sig_scaled'],
        ksize = config['sig_ksize']
    shell: """
        sourmash compute -k {params.ksize} --scaled {params.scaled} \
            {input} -o {output}
    """

rule gather_all:
    input:
        query = output_dir + '/{filename}.sig',
        databases = config['gather_db']
    output:
        csv = output_dir + '/{filename}.gather-matches.csv',
        matches = output_dir + '/{filename}.gather-matches.sig',
        txt = output_dir + '/{filename}.gather-matches.txt'
    conda: 'conf/env-sourmash.yml'
    shell: """
        sourmash gather {input.query} {input.databases} -o {output.csv} \
            --save-matches {output.matches} --threshold-bp=0 >& {output.txt}
        cat {output.txt}
        touch {output.csv} {output.matches}
    """

rule contigs_clean_just_taxonomy:
    input:
        genome = genome_dir + '/{f}',
        matches = output_dir + '/{f}.gather-matches.sig',
        lineages = config['lineages_csv']
    output:
        clean=output_dir + '/{f}.clean.fa.gz',
        dirty=output_dir + '/{f}.dirty.fa.gz',
        report=output_dir + '/{f}.report.txt',
        csv=output_dir + '/{f}.summary.csv'
    conda: 'conf/env-sourmash.yml'
    params:
        lineage = get_provided_lineage
    shell: """
        python -m charcoal.just_taxonomy \
            --genome {input.genome} --lineages_csv {input.lineages} \
            --matches_sig {input.matches} \
            --clean {output.clean} --dirty {output.dirty} \
            --report {output.report} --summary {output.csv} \
            --lineage {params.lineage:q}
    """
#            --lineage 'Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella'


rule combined_summary:
    input:
        expand(output_dir + '/{g}.summary.csv', g=genome_list),
    output:
        output_dir + '/just_taxonomy.combined_summary.csv',
    run:
        # combine all of the summary CSV files
        import csv
        with open(output[0], 'wt') as fp:
            w = csv.writer(fp)
            header = ["genomefile", "taxbrief", "taxfull", "refsize", "ratio", "clean_bp", "clean_n", "dirty_n", "dirty_bp", "missed_n", "missed_bp", "f_major","n_reason_1", "n_reason_2", "n_reason_3", "comment"]
            w.writerow(header)

            for i in input:
                with open(i, 'rt') as in_fp:
                   r = csv.reader(in_fp)
                   rows = list(r)
                   assert len(rows) == 1
                   row = rows[0]
                   assert len(row) == len(header), (len(row), len(header), i)

                   w.writerow(row)
