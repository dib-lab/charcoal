#
# run with --use-conda for maximal froodiness.
#

# override this with --configfile on command line
configfile: 'test-data/conf-test.yml'

### config stuff loaded from config file
genome_list_file = config['genome_list']
genome_list = [ line.strip() for line in open(genome_list_file, 'rt') ]

genome_dir = config['genome_dir'].rstrip('/')
output_dir = config['output_dir'].rstrip('/')

metagenome_sig_list = config['metagenome_sig_list']
metagenome_sig_dir = config['metagenome_sig_dir'].rstrip('/')

lca_db = config['lca_db']

### rules!

rule all:
    input:
        expand(output_dir + '/{g}.hashes', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.10000', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.5000', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.matrix.csv', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.matrix.csv.mat.pdf', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.tax', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.tree', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.tree.newick', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.tree.png', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.tax.rm.clean.fa', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.tree.rm.clean.fa', g=genome_list),

rule make_hashes:
    input:
        genome_dir + '/{filename}'
    output:
        output_dir + '/{filename}.hashes'
    conda: 'conf/env-sourmash.yml'
    shell: """
        ./charcoal/process_genome.py {input} {output}
     """

rule make_hashes_fragment:
    input:
        genome_dir + '/{filename}'
    output:
        hashes=output_dir + '/{filename}.hashes.fragment.{size,\d+}',
        stats=output_dir + '/{filename}.hashes.fragment.{size,\d+}.stats'
    conda: 'conf/env-sourmash.yml'
    params:
        scaled=config['lca_scaled']
    shell: """
        ./charcoal/process_genome.py {input} {output.hashes} \
             --fragment {wildcards.size} --stats {output.stats} \
             --scaled={params.scaled}
     """

rule make_matrix_csv:
    input:
        hashes=output_dir + '/{filename}.hashes{postfix}',
        metag_list=metagenome_sig_list,
    output:
        output_dir + '/{filename}.hashes{postfix}.matrix.csv'
    params:
        metagenome_sig_dir=metagenome_sig_dir
    conda: 'conf/env-sourmash.yml'
    shell: """
        ./charcoal/match_metagenomes.py {input.hashes} {input.metag_list} \
            {output} -d {params.metagenome_sig_dir}
    """

rule make_matrix_pdf:
    input:
        output_dir + '/{g}.matrix.csv',
    output:
        matrix_pdf=output_dir + '/{g}.matrix.csv.mat.pdf',
        dendro_pdf=output_dir + '/{g}.matrix.csv.dendro.pdf',
        out=output_dir + '/{g}.matrix.csv.dendro.out'
    conda: 'conf/env-sourmash.yml'
    shell: """
        ./charcoal/cluster_and_plot.py {input} {output.matrix_pdf} \
            --dendro {output.dendro_pdf} > {output.out}
    """

rule make_taxhashes:
    input:
        genome_dir + '/{filename}'
    output:
        taxhashes=output_dir + '/{filename}.hashes.fragment.{size,\d+}.tax',
        taxcsv=   output_dir + '/{filename}.hashes.fragment.{size,\d+}.tax.csv'
    conda: 'conf/env-sourmash.yml'
    params:
        lca_db=lca_db,
    shell: """
        ./charcoal/genome_shred_to_tax.py {params.lca_db} {input} \
             {output.taxcsv} \
             --fragment {wildcards.size} --save-tax-hashes {output.taxhashes}
     """

rule make_tree:
    input:
        matrix=   output_dir + '/{f}.hashes.fragment.{size}.matrix.csv',
        hashes=   output_dir + '/{f}.hashes.fragment.{size}',
        taxhashes=output_dir + '/{f}.hashes.fragment.{size}.tax',
    output:
        output_dir + "/{f}.hashes.fragment.{size,\d+}.tree"
    conda: 'conf/env-sourmash.yml'
    params:
        lca_db=lca_db,
    shell: """
        ./charcoal/combine_tax_togetherness.py \
             {input.matrix} {input.hashes} {input.taxhashes} \
             --pickle-tree {output}
     """

rule make_tree_viz:
    input:
        output_dir + "/{f}.hashes.fragment.{size}.tree"
    output:
        output_dir + '/{f}.hashes.fragment.{size,\d+}.tree.newick',
    conda: 'conf/env-ete.yml'
    shell: """
        ./charcoal/ete_make_newick.py {input} {output}
     """

rule make_tree_viz_output:
    input:
        output_dir + '/{f}.hashes.fragment.{size}.tree.newick',
    output:
        png = output_dir + '/{f}.hashes.fragment.{size,\d+}.tree.png',
        svg = output_dir + '/{f}.hashes.fragment.{size,\d+}.tree.svg',
    conda: 'conf/env-ete.yml'
    shell: """
        export QT_QPA_PLATFORM=offscreen   # turn off interactive ete3 behavior
        ete3 view -t {input} --show_internal_names -i {output.png}
        ete3 view -t {input} --show_internal_names -i {output.svg}
     """

# generic rule for removing contigs/fragments by hash
rule separate_clean_dirty:
    input:
        genome  =genome_dir + '/{f}',
        rmhashes=output_dir + '/{f}.hashes.fragment.{size}.{suffix}',
    output:
        clean=output_dir + '/{f}.hashes.fragment.{size,\d+}.{suffix}.clean.fa',
        dirty=output_dir + '/{f}.hashes.fragment.{size,\d+}.{suffix}.dirty.fa'
    conda: 'conf/env-sourmash.yml'
    shell: """
        charcoal/remove_contigs_by_hash.py {input.genome} {input.rmhashes} \
            --fragment {wildcards.size} {output.clean} {output.dirty}
     """

# do dumb cleaning, based solely on majority order.
rule clean_tax_only:
    input:
        output_dir + '/{f}.hashes.fragment.{size}.tax',
    output:
        output_dir + '/{f}.hashes.fragment.{size,\d+}.tax.rm'
    conda: 'conf/env-sourmash.yml'
    shell: """
        charcoal/remove_tax_hashes.py {input} --rm-hashes {output}
     """

# do dumb cleaning, based solely on tip taxonomy
rule clean_tip_only:
    input:
        hashes=   output_dir + '/{f}.hashes.fragment.{size}',
        taxhashes=output_dir + '/{f}.hashes.fragment.{size}.tree',
    output:
        output_dir + '/{f}.hashes.fragment.{size,\d+}.tree.rm'
    conda: 'conf/env-sourmash.yml'
    shell: """
        charcoal/remove_tips.py {input.taxhashes} {input.hashes} \
              --rm-hashes {output}
     """

