##
## database download stuff
##

rule all:
    input:
        database = 'db/gtdb-release89-k31.sbt.zip',
        lineages = 'db/gtdb-release89-lineages.csv'

rule download_gtdb_lineages:
    output:
        lineages = 'db/gtdb-release89-lineages.csv'
    shell: """
        curl -L https://osf.io/q645v/download | gunzip > {output.lineages}
    """

rule download_gtdb_sbt_zip:
    output:
        database = 'db/gtdb-release89-k31.sbt.zip',
    shell: """
        curl -L https://osf.io/5mb9k/download > {output.database}
    """
