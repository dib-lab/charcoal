genome_list = [ line.strip() for line in open('genome-list.txt', 'rt') ]

gtdb_db = '/home/ctbrown/sourmash_databases/gtdb/gtdb-release89-k31.lca.json.gz'

rule all:
    input:
        expand('genomes/{g}.hashes', g=genome_list),
        expand('genomes/{g}.hashes.fragment.100000', g=genome_list),
        expand('genomes/{g}.hashes.fragment.10000', g=genome_list),
        expand('genomes/{g}.hashes.fragment.5000', g=genome_list),
        expand('genomes/{g}.hashes.fragment.100000.matrix.csv', g=genome_list),
        expand('genomes/{g}.hashes.fragment.100000.matrix.csv.mat.pdf', g=genome_list),
        expand('genomes/{g}.hashes.fragment.100000.tax', g=genome_list)

rule make_hashes:
    input:
        'genomes/{filename}.fna'
    output:
        'genomes/{filename}.fna.hashes'
    conda: 'env-sourmash.yml'
    shell: """
        ./process-genome.py {input} {output}
     """

rule make_hashes_fragment:
    input:
        'genomes/{filename}.fna'
    output:
        hashes='genomes/{filename}.fna.hashes.fragment.{size,\d+}',
        stats='genomes/{filename}.fna.hashes.fragment.{size,\d+}.stats'
    conda: 'env-sourmash.yml'
    shell: """
        ./process-genome.py {input} {output.hashes} \
             --fragment {wildcards.size} --stats {output.stats}
     """

rule make_matrix_csv:
    input:
        hashes='genomes/{filename}.hashes{postfix}',
        metag_list='ibd_metagenome_prefixes.txt'
    output:
        'genomes/{filename}.hashes{postfix}.matrix.csv'
    conda: 'env-sourmash.yml'
    shell: """
        ./match-metagenomes.py {input.hashes} {input.metag_list} {output}
    """

rule make_matrix_pdf:
    input:
        'genomes/{g}.matrix.csv',
    output:
        matrix_pdf='genomes/{g}.matrix.csv.mat.pdf',
        dendro_pdf='genomes/{g}.matrix.csv.dendro.pdf',
        out='genomes/{g}.matrix.csv.dendro.out'
    conda: 'env-sourmash.yml'
    shell: """
        ./cluster-and-plot.py {input} {output.matrix_pdf} \
            --dendro {output.dendro_pdf} > {output.out}
    """

rule make_taxhashes:
    input:
        'genomes/{filename}.fna'
    output:
        taxhashes='genomes/{filename}.fna.hashes.fragment.{size,\d+}.tax',
        taxcsv=   'genomes/{filename}.fna.hashes.fragment.{size,\d+}.tax.csv'
    conda: 'env-sourmash.yml'
    params:
        lca_db=gtdb_db,
    shell: """
        ./genome-shred-to-tax.py {params.lca_db} {input} {output.taxcsv} \
             --fragment {wildcards.size} --save-tax-hashes {output.taxhashes}
     """

