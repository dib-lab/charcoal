genome_list = [ line.strip() for line in open('genome-list.txt', 'rt') ]

rule all:
     input:
        expand('genomes/{g}.hashes', g=genome_list),
        expand('genomes/{g}.hashes.fragment.100000', g=genome_list),
        expand('genomes/{g}.hashes.fragment.10000', g=genome_list)

rule make_hashes:
     input:
        'genomes/{filename}.fna'
     output:
        'genomes/{filename}.fna.hashes'
     shell: """
        ./process-genome.py {input} {output}
     """

rule make_hashes_fragment:
     input:
        'genomes/{filename}.fna'
     output:
        'genomes/{filename}.fna.hashes.fragment.{size}'
     shell: """
        ./process-genome.py {input} {output} --fragment {wildcards.size}
     """
