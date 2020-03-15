configfile: 'test-data/conf-test.yml'

genome_list_file = config['genome_list']
genome_list = [ line.strip() for line in open(genome_list_file, 'rt') ]

genome_dir = config['genome_dir'].rstrip('/')
output_dir = config['output_dir'].rstrip('/')

metagenome_sig_list = config['metagenome_sig_list']
metagenome_sig_dir = config['metagenome_sig_dir'].rstrip('/')

lca_db = config['lca_db']

rule all:
    input:
        expand(output_dir + '/{g}.hashes', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.10000', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.5000', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.matrix.csv', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.matrix.csv.mat.pdf', g=genome_list),
        expand(output_dir + '/{g}.hashes.fragment.100000.tax', g=genome_list)

rule make_hashes:
    input:
        genome_dir + '/{filename}'
    output:
        output_dir + '/{filename}.hashes'
    conda: 'env-sourmash.yml'
    shell: """
        ./process-genome.py {input} {output}
     """

rule make_hashes_fragment:
    input:
        genome_dir + '/{filename}'
    output:
        hashes=output_dir + '/{filename}.hashes.fragment.{size,\d+}',
        stats=output_dir + '/{filename}.hashes.fragment.{size,\d+}.stats'
    conda: 'env-sourmash.yml'
    params:
        scaled=config['lca_scaled']
    shell: """
        ./process-genome.py {input} {output.hashes} \
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
    conda: 'env-sourmash.yml'
    shell: """
        ./match-metagenomes.py {input.hashes} {input.metag_list} {output} \
            -d {params.metagenome_sig_dir}
    """

rule make_matrix_pdf:
    input:
        output_dir + '/{g}.matrix.csv',
    output:
        matrix_pdf=output_dir + '/{g}.matrix.csv.mat.pdf',
        dendro_pdf=output_dir + '/{g}.matrix.csv.dendro.pdf',
        out=output_dir + '/{g}.matrix.csv.dendro.out'
    conda: 'env-sourmash.yml'
    shell: """
        ./cluster-and-plot.py {input} {output.matrix_pdf} \
            --dendro {output.dendro_pdf} > {output.out}
    """

rule make_taxhashes:
    input:
        genome_dir + '/{filename}'
    output:
        taxhashes=output_dir + '/{filename}.hashes.fragment.{size,\d+}.tax',
        taxcsv=   output_dir + '/{filename}.hashes.fragment.{size,\d+}.tax.csv'
    conda: 'env-sourmash.yml'
    params:
        lca_db=lca_db,
    shell: """
        ./genome-shred-to-tax.py {params.lca_db} {input} {output.taxcsv} \
             --fragment {wildcards.size} --save-tax-hashes {output.taxhashes}
     """
