#
# run with --use-conda for maximal froodiness.
#
import csv, sys

# override this with --configfile on command line
configfile: 'test-data/00-test.conf'

strict_val = config.get('strict', '1')
strict_mode = int(strict_val)
if not strict_mode:
    print('** WARNING: strict mode is OFF. Config errors will not force exit.')

force = config.get('force', '0')
force = int(force)
force_param = ''
if force:
    force_param = '--force'

### config stuff loaded from config file
genome_list_file = config['genome_list']
genome_list = [ line.strip() for line in open(genome_list_file, 'rt') ]
genome_list = [ line for line in genome_list if line ]   # remove empty lines

genome_dir = config['genome_dir'].rstrip('/')
output_dir = config['output_dir'].rstrip('/')

### verification / strict mode

scaled = config['scaled']
try:
    scaled = int(scaled)
    if scaled < 1 or scaled > 100000:
        raise ValueError
except ValueError:
    print('** ERROR: scaled should be a number between 1 and 100000')
    print('** (it must also match the query database scaled value)')
    if strict_mode:
        sys.exit(-1)

ksize = config['ksize']
try:
    ksize = int(ksize)
    if ksize < 15 or ksize > 101:
        raise ValueError
except ValueError:
    print('** ERROR: ksize should be a nubmer between 15 and 101.')
    print('** (it must also match the query database ksize value)')
    if strict_mode:
        sys.exit(-1)

# verify that all genome files exist -
for filename in genome_list:
    fullpath = os.path.join(genome_dir, filename)
    if not os.path.exists(fullpath):
        print(f'** ERROR: genome file {filename} does not exist in {genome_dir}')
        if strict_mode:
            print('** exiting.')
            sys.exit(-1)

# verify that all query databases exist --
for filename in config['gather_db']:
    if not os.path.exists(filename):
        print(f'** ERROR: database {filename} does not exist.')
        if strict_mode:
            print('** exiting.')
            sys.exit(-1)

# does lineage csv exist?
filename = config['lineages_csv']
if not os.path.exists(filename):
    print(f'** ERROR: lineage CSV {filename} does not exist.')
    if strict_mode:
        print('** exiting.')
        sys.exit(-1)

# read in provided lineages, if any, and verify file.
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
            if len(row[1:]) <= 1:
                print(f'** cannot parse provided lineage for {genome_filename}')
                print(f'** ; is it comma separated?')
                sys.exit(-1)
            provided_lineages[genome_filename] = row[1:]

    print(f'** read {len(provided_lineages)} provided lineages')

print('** config file checks PASSED!')
print('** from here on out, it\'s all snakemake...')

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

###
### rules!
###

wildcard_constraints:
    size="\d+"

rule all:
    input:
        expand(output_dir + '/{g}.clean.fa.gz', g=genome_list),
        output_dir + '/combined_summary.csv',

rule contigs_sig:
    input:
        genome_dir + '/{filename}'
    output:
        output_dir + '/{filename}.sig'
    conda: 'conf/env-sourmash.yml'
    params:
        scaled = config['scaled'],
        ksize = config['ksize']
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
        lineage = get_provided_lineage,
        force = force_param
    shell: """
        python -m charcoal.just_taxonomy \
            --genome {input.genome} --lineages_csv {input.lineages} \
            --matches_sig {input.matches} \
            --clean {output.clean} --dirty {output.dirty} \
            --report {output.report} --summary {output.csv} \
            --lineage {params.lineage:q} {params.force}
    """

rule combined_summary:
    input:
        expand(output_dir + '/{g}.summary.csv', g=genome_list),
    output:
        output_dir + '/combined_summary.csv',
    run:
        # combine all of the summary CSV files
        import csv
        with open(output[0], 'wt') as fp:
            w = csv.writer(fp)

            header = ["genomefile", "brieftax",
                      "f_major", "f_ident", "f_removed",
                      "n_reason_1", "n_reason_2", "n_reason_3",
                      "refsize", "ratio",
                      "clean_bp", "clean_n", "dirty_n", "dirty_bp", "missed_n",
                      "missed_bp", "taxguessed", "taxprovided",
                      "comment"]

            w.writerow(header)

            for i in input:
                with open(i, 'rt') as in_fp:
                   r = csv.reader(in_fp)
                   rows = list(r)
                   assert len(rows) == 1
                   row = rows[0]
                   assert len(row) == len(header), (len(row), len(header), i)

                   w.writerow(row)
