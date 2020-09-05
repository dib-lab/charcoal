import os.path
import collections

import plotly.graph_objects as go

import sourmash
from charcoal import utils
from .sourmash_sankey import GenomeSankeyFlow


def plot_genome_area_chart(genome, summary):
    row = summary[genome]

    x=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus']

    y_ignored_bp = []
    y_noident_bp = []
    y_clean_bp = []
    y_dirty_bp = []

    ignored_bp = int(row['ignored_contigs_bp'])
    noident_bp = int(row['noident_contigs_bp'])

    for rank in x:
        tax_good_bp = int(row[f'good_{rank}_bp'])
        tax_bad_bp = int(row[f'bad_{rank}_bp'])

        y_ignored_bp.append(ignored_bp)
        y_noident_bp.append(noident_bp)
        y_clean_bp.append(tax_good_bp)
        y_dirty_bp.append(tax_bad_bp)


    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=x, y=y_clean_bp,
        mode='lines',
        line=dict(width=0.5, color='rgb(184, 247, 212)'),
        stackgroup='one',
        name='clean bp',
        groupnorm='percent' # sets the normalization for the sum of the stackgroup
    ))
    fig.add_trace(go.Scatter(
        x=x, y=y_noident_bp,
        mode='lines',
        line=dict(width=0.5, color='rgb(111, 231, 219)'),
        stackgroup='one',
        name='no ident'
    ))
    fig.add_trace(go.Scatter(
        x=x, y=y_ignored_bp,
        mode='lines',
        line=dict(width=0.5, color='rgb(127, 166, 238)'),
        stackgroup='one',
        name='no hashes (ignored)'
    ))
    fig.add_trace(go.Scatter(
        x=x, y=y_dirty_bp,
        mode='lines',
        line=dict(width=0.5, color='rgb(131, 90, 241)'),
        stackgroup='one',
        name='dirty bp'
    ))

    fig.update_layout(
        showlegend=True,
        xaxis_type='category',
        yaxis=dict(
            type='linear',
            range=[1, 100],
            ticksuffix='%'))
    
    title = genome
    fig.update_layout(template='plotly',
                      title=title)

    return fig


def make_flow_fig(title, genome_lineage, contigs_info):
    obj = GenomeSankeyFlow()
    
    counts = collections.Counter()
    for contig_name, gather_info in contigs_info.items():
        contig_taxlist = gather_info.gather_tax

        # iterate over each contig match and summarize counts.              
        # note - here we can stop at first one, or track them all.          
        # note - b/c gather counts each hash only once, these               
        #     are non-overlapping                                           
        total_hashcount = 0
        for lin, hashcount in contig_taxlist:
            counts[lin] += hashcount
            total_hashcount += hashcount

        # track missing => unassigned lineage
        unident = gather_info.num_hashes - total_hashcount
        counts[obj.unassigned_lin] += unident
    
    # set the color of the main lineage
    genome_lineage = tuple(genome_lineage)
    obj.colors[genome_lineage] = "lightseagreen"
    
    # for phylum level disagreements, let's go with palevioletred
    for lin in counts:
        if not utils.is_lineage_match(lin, genome_lineage, 'phylum'):
            obj.colors[lin] = 'palevioletred'
            
    # assign unassigned to good lineage, maybe?
    counts[genome_lineage] += counts[obj.unassigned_lin]
    del counts[obj.unassigned_lin]
    
    obj.make_links(genome_lineage, counts)
    fig = obj.make_plotly_fig(title)
    
    return fig


def load_and_make_flow(dirname, genome_name):
    "Load charcoal output for one specific genome."
    contigs_filename = f'{dirname}/{genome_name}.contigs-tax.json'
    contigs_info = utils.load_contigs_gather_json(contigs_filename)

    hit_list_filename = f'{dirname}/hit_list_for_filtering.csv'
    hit_list = utils.CSV_DictHelper(hit_list_filename, 'genome')

    row = hit_list[genome_name]
    genome_lineage = utils.make_lineage(row['lineage'])
    
    title = sourmash.lca.display_lineage(genome_lineage)
    fig = make_flow_fig(title, genome_lineage, contigs_info)
    return fig


def load_all(dirname):
    "Load all charcoal genomes & plot."
    hit_list_filename = f'{dirname}/hit_list_for_filtering.csv'
    hit_list = utils.CSV_DictHelper(hit_list_filename, 'genome')
    
    for genome_name in hit_list.rows:
        if int(hit_list[genome_name]['total_bad_bp']):
            fig = load_and_make_flow(dirname, genome_name)
            contigs_filename = f'{dirname}/{genome_name}.contigs-tax.json'
            contigs_info = utils.load_contigs_gather_json(contigs_filename)

            row = hit_list[genome_name]
            genome_lineage = utils.make_lineage(row['lineage'])
            
            fig = make_sankey_fig(genome_name, genome_lineage, contigs_info)

            fig.show()
            try:
                fig.write_image(f"images/{genome_name}.png", width=1500)
            except ValueError: # orca may not be configured properly
                pass
