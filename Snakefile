#
# run with --use-conda for maximal froodiness.
#

# override this with --configfile on command line
configfile: 'test-data/conf-test.yml'

### config stuff loaded from config file
genome_list_file = config['genome_list']
genome_list = [ line.strip() for line in open(genome_list_file, 'rt') ]
genome_list = [ line for line in genome_list if line ]   # remove empty lines

genome_dir = config['genome_dir'].rstrip('/')
output_dir = config['output_dir'].rstrip('/')

metagenome_sig_list = config['metagenome_sig_list']
metagenome_sig_dir = config['metagenome_sig_dir'].rstrip('/')

lca_db = config['lca_db']

### utility functions
def output_files(filename_template, **kw):
    return expand(output_dir + filename_template, **kw)

### rules!

wildcard_constraints:
    size="\d+"

rule all:
    input:
        expand(output_dir + '/{g}.clean.fa.gz', g=genome_list),
        output_dir + '/just_taxonomy.combined_summary.csv',

rule old:
    input:
        expand(output_dir + '/{g}.hash.100000', g=genome_list),
        expand(output_dir + '/{g}.hash.10000', g=genome_list),
        expand(output_dir + '/{g}.hash.5000', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.matrix.csv', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.matrix.csv.mat.pdf', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.tax', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.tree', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.tax.rm.clean.fa', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.tree.rm.clean.fa', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.tree.cut.clean.fa', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.tree.json', g=genome_list),

rule all_make_tree_viz:
    input:
        expand(output_dir + '/{g}.hash.100000.tree.newick', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.tree.png', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.tree.json.pdf', g=genome_list),
        expand(output_dir + '/{g}.hash.100000.tree.json.svg', g=genome_list),

rule reports:
    input:
        expand(output_dir + '/{g}.hash.100000.report.txt', g=genome_list)

####

rule make_hashes:
    input:
        genome_dir + '/{filename}'
    output:
        output_dir + '/{filename}.hashes'
    conda: 'conf/env-sourmash.yml'
    params:
        scaled=config['lca_scaled']
    shell: """
        python -m charcoal.process_genome \
             --genome {input} --save-hashes {output} \
             --scaled={params.scaled}
     """

rule make_hashes_fragment:
    input:
        genome_dir + '/{filename}'
    output:
        hashes=output_dir + '/{filename}.hash.{size}',
        stats=output_dir + '/{filename}.hash.{size}.stats'
    conda: 'conf/env-sourmash.yml'
    params:
        scaled=config['lca_scaled']
    shell: """
        python -m charcoal.process_genome --genome {input} \
             --save-hashes {output.hashes} \
             --fragment {wildcards.size} --stats {output.stats} \
             --scaled={params.scaled}
     """

rule make_matrix:
    input:
        hashes=output_dir + '/{filename}.hash{postfix}',
        metag_list=metagenome_sig_list,
    output:
        csv = output_dir + '/{filename}.hash{postfix}.matrix.csv',
        mat = output_dir + '/{filename}.hash{postfix}.matrix'
    params:
        metagenome_sig_dir=metagenome_sig_dir
    conda: 'conf/env-sourmash.yml'
    shell: """
        python -m charcoal.match_metagenomes --load-hashes {input.hashes} \
            --metagenome-sigs-list {input.metag_list} \
            --matrix-csv-out {output.csv} --matrix-pickle-out {output.mat} \
            --metagenome-sigs-dir {params.metagenome_sig_dir}
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
        python -m charcoal.cluster_and_plot --load-matrix-csv {input} \
            --output-fig {output.matrix_pdf} \
            --dendro-out {output.dendro_pdf} > {output.out}
    """

rule make_taxhashes_multi:
    input:
        expand(genome_dir + '/{f}', f=genome_list)
    output:
        taxhashes = output_files('/{f}.hash.{{size}}.tax', f=genome_list),
        taxcsv    = output_files('/{f}.hash.{{size}}.tax.csv', f=genome_list)
    conda: 'conf/env-sourmash.yml'
    resources:
        mem_mb=180000,
    params:
        lca_db=lca_db,
        output_dir=output_dir,
        tax_template = lambda wildcards: output_dir + '/{{genome}}.hash.{size}.tax'.format(size=wildcards.size),
        csv_template = lambda wildcards: output_dir + '/{{genome}}.hash.{size}.tax.csv'.format(size=wildcards.size),
    shell: """
        python -m charcoal.genome_shred_to_tax_multi \
             --lca-db {params.lca_db} --genomes {input} \
             --csv-output-template {params.csv_template} \
             --save-tax-hashes-template {params.tax_template} \
             --fragment {wildcards.size}
     """

rule make_tree:
    input:
        matrix    = output_dir + '/{f}.hash.{size}.matrix',
        taxhashes = output_dir + '/{f}.hash.{size}.tax',
    output:
        output_dir + "/{f}.hash.{size}.tree"
    conda: 'conf/env-sourmash.yml'
    params:
        lca_db=lca_db,
    shell: """
        python -m charcoal.combine_tax_togetherness \
             --load-matrix-pickle {input.matrix} \
             --load-tax-hashes {input.taxhashes} \
             --pickle-tree {output}
     """

rule make_tree_viz:
    input:
        output_dir + "/{f}.hash.{size}.tree"
    output:
        output_dir + '/{f}.hash.{size}.tree.newick',
    conda: 'conf/env-ete.yml'
    shell: """ ##
        python -m charcoal.ete_make_newick {input} {output}
     """

rule make_tree_viz_output:
    input:
        output_dir + '/{f}.hash.{size}.tree.newick',
    output:
        png = output_dir + '/{f}.hash.{size}.tree.png',
        svg = output_dir + '/{f}.hash.{size}.tree.svg',
    conda: 'conf/env-ete.yml'
    shell: """
        export QT_QPA_PLATFORM=offscreen   # turn off interactive ete3 behavior
        ete3 view -t {input} --show_internal_names -i {output.png}
        ete3 view -t {input} --show_internal_names -i {output.svg}
     """

rule make_tree_viz_output_pretty:
    input: 
        hash_to_frag = output_dir + '/{f}.hash.{size}.tree.json.hashes_to_fragment.csv',
        hash_to_tax = output_dir + '/{f}.hash.{size}.tree.json.hashes_to_tax.csv',
        leaves_to_hash = output_dir + '/{f}.hash.{size}.tree.json.leaves_to_hashval.csv',
        node_to_children = output_dir + '/{f}.hash.{size}.tree.json.node_to_children.csv',
    output:
        pdf = output_dir + '/{f}.hash.{size}.tree.json.pdf',
        svg = output_dir + '/{f}.hash.{size}.tree.json.svg',
    params: 
        output_dir = output_dir  
    conda:  'conf/env-ggtree.yml' 
    script: "charcoal/plot_json_pretty.R"

# generic rule for removing contigs/fragments by hash
rule separate_clean_dirty:
    input:
        genome   = genome_dir + '/{f}',
        rmhashes = output_dir + '/{f}.hash.{size}.{suffix}',
    output:
        clean=output_dir + '/{f}.hash.{size}.{suffix}.clean.fa',
        dirty=output_dir + '/{f}.hash.{size}.{suffix}.dirty.fa'
    conda: 'conf/env-sourmash.yml'
    shell: """
        python -m charcoal.remove_contigs_by_hash --genome {input.genome} \
            --hashlist {input.rmhashes} \
            --fragment {wildcards.size} \
            --clean-output {output.clean} --dirty-output {output.dirty}
     """

# do dumb cleaning, based solely on majority order.
rule clean_tax_only:
    input:
        output_dir + '/{f}.hash.{size}.tax',
    output:
        output_dir + '/{f}.hash.{size}.tax.rm'
    conda: 'conf/env-sourmash.yml'
    shell: """
        python -m charcoal.remove_tax_hashes {input} --rm-hashes {output}
     """

# do dumb cleaning, based solely on tip taxonomy
rule clean_tip_only:
    input:
        hashes= output_dir + '/{f}.hash.{size}',
        tree  = output_dir + '/{f}.hash.{size}.tree',
    output:
        output_dir + '/{f}.hash.{size}.tree.rm'
    conda: 'conf/env-sourmash.yml'
    shell: """ ##
        python -m charcoal.remove_tips {input.tree} {input.hashes} \
              --rm-hashes {output}
     """

# slightly cleverer cleaning, based on tree cutting
rule cut_tree:
    input:
        hashes= output_dir + '/{f}.hash.{size}',
        tree  = output_dir + '/{f}.hash.{size}.tree',
    output:
        output_dir + '/{f}.hash.{size}.tree.cut'
    conda: 'conf/env-sourmash.yml'
    shell: """ ##
        python -m charcoal.cut_tree_1 {input.tree} {input.hashes} \
              --rm-hashes {output}
     """

# JSON output
rule together_json:
    input:
        genome=genome_dir + "/{f}",
        taxhashes=output_dir + "/{f}.hash.{size}.tax",
        tree=output_dir + "/{f}.hash.{size}.tree"
    output:
        json = output_dir + '/{f}.hash.{size}.tree.json',
        csv1 = output_dir + '/{f}.hash.{size}.tree.json.hashes_to_fragment.csv',
        csv2 = output_dir + '/{f}.hash.{size}.tree.json.hashes_to_tax.csv',
        csv3 = output_dir + '/{f}.hash.{size}.tree.json.leaves_to_hashval.csv',
        csv4 = output_dir + '/{f}.hash.{size}.tree.json.node_id_to_tax.csv',
        csv5 = output_dir + '/{f}.hash.{size}.tree.json.node_to_children.csv',
    conda: 'conf/env-sourmash.yml'
    shell: """ ##
        python -m charcoal.together_tree_to_json \
               {input.genome} {input.taxhashes} {input.tree} {output.json} \
               --fragment {wildcards.size}
    """

rule make_report:
    input:
        genome=genome_dir + "/{f}",
        taxhashes=output_dir + "/{f}.hash.{size}.tax",
        tree=output_dir + "/{f}.hash.{size}.tree",
        matrix=output_dir + '/{f}.hash.{size}.matrix',
        tips_rm=output_dir + '/{f}.hash.{size}.tree.rm',
        tax_rm=output_dir + '/{f}.hash.{size}.tax.rm',
        cut1_rm=output_dir + '/{f}.hash.{size}.tree.cut',
    output:
        output_dir + '/{f}.hash.{size}.report.txt'
    shell: """ ##
        python -m charcoal.report_and_summarize {input.genome} \
            --tax-hashes {input.taxhashes} --matrix {input.matrix} \
            --tips-rm {input.tax_rm} \
            --tax-rm {input.tips_rm} \
            --cut1-rm {input.cut1_rm} \
            -o {output}
    """

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
        database = config['gather_db']
    output:
        csv = output_dir + '/{filename}.gather-matches.csv',
        matches = output_dir + '/{filename}.gather-matches.sig',
        txt = output_dir + '/{filename}.gather-matches.txt'
    conda: 'conf/env-sourmash.yml'
    shell: """
        sourmash gather {input.query} {input.database} -o {output.csv} \
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
    shell: """
        python -m charcoal.just_taxonomy \
            --genome {input.genome} --lineages_csv {input.lineages} \
            --matches_sig {input.matches} \
            --clean {output.clean} --dirty {output.dirty} \
            --report {output.report} --summary {output.csv}

    """


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
            header = ["genomefile","clean_n", "clean_bp", "dirty_n", "dirty_bp","f_major","n_reason_1", "n_reason_2", "n_reason_3"]
            w.writerow(header)

            for i in input:
                with open(i, 'rt') as in_fp:
                   r = csv.reader(in_fp)
                   rows = list(r)
                   assert len(rows) == 1
                   row = rows[0]
                   assert len(row) == len(header)

                   w.writerow(row)
