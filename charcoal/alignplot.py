#! /usr/bin/env python
"""
Class etc to produce a stacked dotplot and other genome overlap/contamination
stats.

TODO:
* argparse the thang
"""
import sys
import argparse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import csv
import tempfile
import shutil
import subprocess
import os
import glob
from collections import defaultdict, namedtuple
import numpy

import screed
from interval import interval


AlignedRegion = namedtuple(
    "AlignedRegion", "query, target, qstart, qend, tstart, tend, pident, qsize, tsize"
)


def glob_all(pattern, endings):
    g = []
    for end in endings:
        p = pattern + end
        g.extend(glob.glob(p))

    return g


def group_regions_by(regions, by_type):
    "gather regions by 'target' or 'query'"

    assert by_type in ("target", "query")

    regions_by = defaultdict(list)
    for region in regions:
        name = getattr(region, by_type)
        x = regions_by[name]
        x.append(region)

    return regions_by


def calc_regions_aligned_bp(regions_by, by_type, filter_by=None):
    assert by_type in ("target", "query")

    regions_sum_kb = {}
    for name, rr in regions_by.items():
        if filter_by:
            rr = [r for r in rr if filter_by(r)]
        if not rr:
            continue

        if by_type == "target":
            ivals = [interval[r.tstart, r.tend] for r in rr]
        elif by_type == "query":
            ivals = [interval[r.qstart, r.qend] for r in rr]
        else:
            raise Exception(by_type)

        ii = interval()
        for iv in ivals:
            ii |= iv

        sum_kb = 0
        for component in ii:
            (start, end) = component
            size = end - start
            sum_kb += size

        regions_sum_kb[name] = sum_kb
    return regions_sum_kb


def region_size(region, by_type):
    "Calculate the size of an alignment based on target or query."
    assert by_type in ("target", "query")

    if by_type == "target":
        return abs(region.tend - region.tstart)
    elif by_type == "query":
        return abs(region.qend - region.qstart)
    raise Exception(f"unhandled by_type {by_type}")


def load_contig_sizes(genomefile, init_d=None):
    "load in all the actual contig sizes for this genome (in kb)"
    if init_d is None:
        all_sizes = {}
    else:
        all_sizes = init_d

    for record in screed.open(genomefile):
        name = record.name.split()[0]
        assert name not in all_sizes
        all_sizes[name] = len(record.sequence) / 1e3

    return all_sizes


class AlignmentContainer:
    """
    Build or load, then store a set of alignments between a
    query and a bunch of targets genomes.

    Takes:
    * query accession,
    * multiple target accessions,
    * an optional info file containing mappings from accession to names
    """

    endings = ".gz", ".fa", ".fna"

    def __init__(
        self, q_acc, queryfile, targetfiles, info_file=None,
    ):
        self.q_acc = q_acc
        self.queryfile = queryfile

        self.queryfile = queryfile
        self.q_acc = q_acc

        self.t_acc_list = []
        self.targetfiles = []
        for (t_acc, targetfile) in targetfiles:
            self.t_acc_list.append(t_acc)
            self.targetfiles.append(targetfile)

        self.query_name = q_acc
        self.target_names = {}
        for acc in self.t_acc_list:
            self.target_names[acc] = acc

        if info_file:
            for row in csv.DictReader(open(info_file, "rt")):
                if self.q_acc == row["acc"]:
                    self.query_name = row["ncbi_tax_name"]
                if row["acc"] in self.target_names:
                    self.target_names[row["acc"]] = row["ncbi_tax_name"]

    def get_targetfile(self, t_acc):
        targetfile = None
        for find_t_acc, find_targetfile in zip(self.t_acc_list, self.targetfiles):
            if find_t_acc == t_acc:
                targetfile = find_targetfile
                break

        assert targetfile
        return targetfile

    def run_mashmap(self):
        "Run all the things, save the results."
        results = {}

        for t_acc, targetfile in zip(self.t_acc_list, self.targetfiles):
            name = self.target_names[t_acc]
            regions = self._run_mashmap(targetfile)
            results[t_acc] = regions

        self.results = results

    def run_nucmer(self):
        "Run all the things, save the results."
        results = {}

        for t_acc, targetfile in zip(self.t_acc_list, self.targetfiles):
            name = self.target_names[t_acc]
            regions = self._run_nucmer(targetfile)
            results[t_acc] = regions

        self.results = results

    def _run_mashmap(self, targetfile):
        "Run mashmap instead of nucmer."
        print("running mashmap...")
        tempdir = tempfile.mkdtemp()
        outfile = os.path.join(tempdir, "mashmap.out")
        cmd = f"mashmap -q {self.queryfile} -r {targetfile} -o {outfile} --pi 95"  # -f none -s 1000
        print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        print(f"...done! reading output from {outfile}.")

        results = self._read_mashmap(outfile)
        shutil.rmtree(tempdir)
        return results

    def _read_mashmap(self, filename):
        "Parse the mashmap output."
        fp = open(filename, "rt")

        regions = []
        for line in fp:
            line = line.strip().split()
            (
                query,
                qsize,
                qstart,
                qend,
                strand,
                target,
                tsize,
                tstart,
                tend,
                pident,
            ) = line
            if strand == '-':
                tstart, tend = tend, tstart

            region = AlignedRegion(
                qsize=int(qsize) / 1e3,
                qstart=int(qstart) / 1e3,
                qend=int(qend) / 1e3,
                tsize=int(tsize) / 1e3,
                tstart=int(tstart) / 1e3,
                tend=int(tend) / 1e3,
                pident=float(pident),
                query=query,
                target=target,
            )

            assert region.qend > region.qstart
            regions.append(region)

        return regions

    def _run_nucmer(self, targetfile):
        "Run nucmer and show coords."
        print(f"running nucmer & show-coords for {targetfile}...")
        tempdir = tempfile.mkdtemp()

        queryfile = self.queryfile
        if self.queryfile.endswith(".gz"):
            queryfile = os.path.join(tempdir, "query.fa")
            subprocess.check_call(
                f"gunzip -c {self.queryfile} > {queryfile}", shell=True
            )

        if targetfile.endswith(".gz"):
            newfile = os.path.join(tempdir, "target.fa")
            subprocess.check_call(f"gunzip -c {targetfile} > {newfile}", shell=True)
            targetfile = newfile

        cmd = f"nucmer -p {tempdir}/cmp {queryfile} {targetfile} 2> /dev/null"
        # print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        deltafile = f"{tempdir}/cmp.delta"
        coordsfile = f"{tempdir}/cmp.coords"

        cmd = f"show-coords -T {deltafile} > {coordsfile} 2> /dev/null"
        # print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        print(f"...done! reading output from {tempdir}.")

        results = self._read_nucmer(coordsfile)
        # shutil.rmtree(tempdir)
        return results

    def _read_nucmer(self, filename):
        "Parse the nucmer output."
        fp = open(filename, "rt")
        lines = fp.readlines()
        assert lines[1].startswith("NUCMER"), (filename, lines[0])
        assert not lines[2].strip()

        regions = []
        for line in lines[4:]:
            line = line.strip().split("\t")
            qstart, qend, tstart, tend, qsize, tsize, pident, query, target = line
            region = AlignedRegion(
                qsize=int(qsize) / 1e3,
                qstart=int(qstart) / 1e3,
                qend=int(qend) / 1e3,
                tsize=int(tsize) / 1e3,
                tstart=int(tstart) / 1e3,
                tend=int(tend) / 1e3,
                pident=float(pident),
                query=query,
                target=target,
            )

            # identity and length filter - @CTB move outside!
            #            if region.pident < 95 or abs(region.qend - region.qstart) < 0.5:
            #                continue

            regions.append(region)

        return regions

    def filter(self, pident=None, query_size=None):
        "Filter alignments at given pident and query size."
        if pident is None and query_size is None:
            return

        new_results = {}
        for t_acc, t_results in self.results.items():
            filtered = []
            for region in t_results:
                keep = True
                if pident and region.pident < pident:
                    keep = False
                if query_size and region_size(region, 'query') < query_size:
                    keep = False

                if keep:
                    filtered.append(region)
            new_results[t_acc] = filtered

    def calc_shared(self, t_acc=None):
        "Calculate the number of bases shared by query and ..."
        if t_acc:
            regions = self.results[t_acc]
        else:
            regions = []
            for t_acc in self.results:
                regions.extend(self.results[t_acc])

        regions_by_query = group_regions_by(regions, 'query')
        sum_kb = calc_regions_aligned_bp(regions_by_query, 'query')
        return sum_kb


class StackedDotPlot:
    def __init__(self, alignment):
        self.alignment = alignment

    def plot(self):
        "Do the actual stacked dotplot plotting."
        alignment = self.alignment

        if alignment.q_acc == alignment.query_name:
            ylabel_text = f"{alignment.q_acc} (coords in kb)"
        else:
            ylabel_text = f"{alignment.q_acc}: {alignment.query_name} (kb)"
        plt.ylabel(ylabel_text)
        plt.xlabel("coordinates of matches (scaled to kb)")

        colors = ("r-", "b-", "g-")

        q_starts = {}
        q_sofar = 0

        # the use of max_x is what makes it a stacked dotplot!! :)
        max_x = 0  # track where to start each target

        # iterate over each set of features, plotting lines.
        for t_acc, color in zip(alignment.t_acc_list, colors):
            name = alignment.target_names[t_acc]
            # CTB: if we move this out of the loop and plot self-x-self
            # there is an interesting effect of showing distribution.
            # exploreme!
            t_starts = {}
            t_sofar = 0

            line = None
            this_max_x = 0
            for region in alignment.results[t_acc]:
                # calculate the base y position for this query contig --
                q_base = q_starts.get(region.query)
                if q_base is None:
                    q_starts[region.query] = q_sofar
                    q_base = q_sofar
                    q_sofar += region.qsize

                # calculate the base x position for this target contig --
                t_base = t_starts.get(region.target)
                if t_base is None:
                    t_starts[region.target] = t_sofar
                    t_base = t_sofar
                    t_sofar += region.tsize

                x_0 = t_base + region.tstart
                y_0 = q_base + region.qstart

                x_1 = t_base + region.tend
                y_1 = q_base + region.qend

                # stack 'em horizontally with max_x
                line = plt.plot((x_0 + max_x, x_1 + max_x), (y_0, y_1), color)
                this_max_x = max(this_max_x, x_0, x_1)

            # label the last plotted line w/the right name to make legend
            if line:
                line[0].set_label(name)

            # "stack" the dotplots horizontally.
            max_x = this_max_x

        plt.legend(loc="lower right")

        return plt.gcf()

    def target_response_curve(self, t_acc):
        alignment = self.alignment
        regions = alignment.results[t_acc]

        # first, find the targetfile (genome) for this accession
        targetfile = alignment.get_targetfile(t_acc)

        # calculate and sort region summed kb in alignments over 95%
        regions_by_target = group_regions_by(regions, "target")
        regions_aligned_kb = calc_regions_aligned_bp(
            regions_by_target, "target", filter_by=lambda r: r.pident >= 95
        )
        region_items = list(regions_aligned_kb.items())
        region_items.sort(key=lambda x: -x[1])

        # load in all the actual contig sizes for this genome
        all_sizes = load_contig_sizes(targetfile)
        sum_bp = sum(all_sizes.values())

        # construct points for plot --
        x = [0]  # kb in target contigs
        y = [0]  # alignments in those target contigs
        sofar = 0
        aligned_sofar = 0

        # start with contigs with most aligned bases first - the sorting order matters here!
        for name, ani95_kb in region_items:
            # ok, track total kb and aligned kb added by this contig
            sofar += all_sizes[name]
            aligned_sofar += ani95_kb
            assert all_sizes[name] > 0

            x.append(sofar)
            y.append(aligned_sofar)

        saturation_point = sofar

        # add in the rest of the contigs that have no alignments in 'em'
        remaining_names = set(all_sizes) - set(region_items)
        for contig in remaining_names:
            sofar += all_sizes[contig]
            x.append(sofar)
            y.append(aligned_sofar)

        return numpy.array(x), numpy.array(y), saturation_point

    def query_response_curve(self):
        alignment = self.alignment

        # aggregate regions over _all_ results
        regions = []
        for k, v in alignment.results.items():
            regions.extend(v)

        queryfile = alignment.queryfile

        # calculate and sort region summed kb in alignments over 95%
        regions_by_query = group_regions_by(regions, "query")
        regions_aligned_kb = calc_regions_aligned_bp(
            regions_by_query, "query", filter_by=lambda r: r.pident >= 95
        )
        region_items = list(regions_aligned_kb.items())
        region_items.sort(key=lambda x: -x[1])

        # load in all the actual contig sizes for this genome
        all_sizes = load_contig_sizes(queryfile)
        sum_bp = sum(all_sizes.values())

        # construct points for plot --
        x = [0]  # kb in query contigs
        y = [0]  # alignments in those query contigs
        sofar = 0
        aligned_sofar = 0

        # start with contigs with most aligned bases first - the sorting order matters here!
        for name, ani95_kb in region_items:
            # ok, track total kb and aligned kb added by this contig
            sofar += all_sizes[name]
            aligned_sofar += ani95_kb
            assert all_sizes[name] > 0

            x.append(sofar)
            y.append(aligned_sofar)

        saturation_point = sofar

        # add in the rest of the contigs that have no alignments in 'em'
        remaining_names = set(all_sizes) - set(region_items)
        for contig in remaining_names:
            sofar += all_sizes[contig]
            x.append(sofar)
            y.append(aligned_sofar)

        return numpy.array(x), numpy.array(y), saturation_point


class AlignmentSlopeDiagram:
    """
    A FamilyRelations-style diagram showing proportional sloped lines
    for alignments.
    """
    def __init__(self, alignment):
        self.alignment = alignment
        
    def calculate(self, select_n=None, plot_all_contigs=False):
        alignment = self.alignment

        regions = []
        for k, v in alignment.results.items():
            regions.extend(v)

        queryfile = alignment.queryfile

        # calculate and sort region summed kb in alignments over 95%
        regions_by_query = group_regions_by(regions, "query")
        regions_aligned_kb = calc_regions_aligned_bp(
            regions_by_query, "query", filter_by=lambda r: r.pident >= 95
        )
        region_items = list(regions_aligned_kb.items())
        region_items.sort(key=lambda x: -x[1])
        
        if select_n:
            region_items = region_items[:select_n]

        query_sizes = load_contig_sizes(alignment.queryfile)

        target_sizes = {}
        for targetfile in alignment.targetfiles:
            load_contig_sizes(targetfile, target_sizes)


        query_list = []
        target_list = []
        align = []

        query_idx_d = {}
        target_idx_d = {}

        for name, aligned_kb in region_items:
            for alignment in regions_by_query[name]:
                a = alignment

                query_idx = query_idx_d.get(a.query)
                if query_idx is None:
                    query_idx_d[a.query] = len(query_list)
                    query_idx = len(query_list)
                    query_list.append(query_sizes[a.query])

                target_idx = target_idx_d.get(a.target)
                if target_idx is None:
                    target_idx_d[a.target] = len(target_list)
                    target_idx = len(target_list)
                    target_list.append(target_sizes[a.target])

                assert a.qstart <= query_list[query_idx]
                assert a.qstart >= 0
                assert a.qend <= query_list[query_idx]
                assert a.qend >= 0

                assert a.tstart <= target_list[target_idx], (a.tstart, target_idx, target_list[target_idx])
                assert a.tstart >= 0
                assert a.tend <= target_list[target_idx]
                assert a.tend >= 0

                align.append((query_idx, a.qstart, a.qend,
                              target_idx, a.tstart, a.tend))

        if plot_all_contigs:
            remaining_query = set(query_sizes) - set(query_idx_d)
            remaining_target = set(target_sizes) - set(target_idx_d)

            for query in remaining_query:
                query_list.append(query_sizes[query])
            for target in remaining_target:
                target_list.append(target_sizes[target])
                
        self.from_contig_sizes = query_list
        self.to_contig_sizes = target_list
        self.alignments = align
        
    def plot(self, select_n=None, plot_all_contigs=False, use_labels=True):
        self.calculate(select_n=select_n, plot_all_contigs=plot_all_contigs)

        from_contigs = self.from_contig_sizes
        to_contigs = self.to_contig_sizes
        alignments = self.alignments

        if not from_contigs or not to_contigs or not self.alignments:
            return
        
        fig, ax = plt.subplots()

        from_sum = sum(from_contigs)
        to_sum = sum(to_contigs)
        
        from_bases = [0]
        sofar = 0
        for i in from_contigs:
            sofar += i
            from_bases.append(sofar)
        
        to_bases = [0]
        sofar = 0
        for i in to_contigs:
            sofar += i
            to_bases.append(sofar)

        patches = []
        for (from_i, qstart, qend, to_i, tstart, tend) in alignments:
            assert qstart <= from_contigs[from_i]
            assert qstart >= 0
            assert qend <= from_contigs[from_i]
            assert qend >= 0

            assert tstart <= to_contigs[to_i], (tstart, to_i, to_contigs[to_i])
            assert tstart >= 0
            assert tend <= to_contigs[to_i]
            assert tend >= 0

            upper_left = (0, 1 - (from_bases[from_i] + qstart) / from_sum)
            upper_right = (0, 1 - (from_bases[from_i] + qend) / from_sum)
            bottom_right = (1, 1 - (to_bases[to_i] + tend) / to_sum)
            bottom_left = (1, 1 - (to_bases[to_i] + tstart) / to_sum)
 
            polygon = Polygon((upper_left, upper_right, bottom_right, bottom_left), True)
            patches.append(polygon)

        p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)

        colors = 100*numpy.random.rand(len(patches))
        p.set_array(numpy.array(colors))

        ax.add_collection(p)

        ax.set_ylabel('contaminated genome contigs')
        from_yticks = [ 1 - i / from_sum for i in from_bases]
        ax.set_yticks(from_yticks)
        if use_labels:
            ax.set_yticklabels([ f"{i:.0f}kb" for i in from_bases ])
        else:
            ax.set_yticklabels([])
        
        secax = ax.secondary_yaxis('right')       
        secax.set_ylabel('source genome contigs')
        to_yticks = [ 1 - i/to_sum for i in to_bases ]
        secax.set_yticks(to_yticks)
        if use_labels:
            secax.set_yticklabels([ f"{i:.0f}kb" for i in to_bases ])
        else:
            secax.set_yticklabels([])

        ax.set_xticks([])
        ax.set_xticklabels([])
        
        return plt.gcf()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("query_acc")
    p.add_argument("target_accs", nargs="+")
    p.add_argument(
        "-g",
        "--genomes-directory",
        default="./genomes",
        help="directory with genome files in it",
    )
    p.add_argument("-i", "--info-file", help="CSV file with nicer names for accessions")
    p.add_argument("-o", "--output-prefix", default="alignplot")
    args = p.parse_args()

    alignment = AlignmentContainer(
        args.query_acc, args.target_accs, args.info_file, args.genomes_directory,
    )
    alignment.run_nucmer()

    dotplot = StackedDotPlot(alignment)
    dotplot.plot()

    print(f"saving {args.output_prefix}-nucmer.png")
    plt.savefig(f"{args.output_prefix}-nucmer.png")
    plt.cla()

    alignment.run_mashmap()

    dotplot = StackedDotPlot(alignment)
    dotplot.plot()

    print(f"saving {args.output_prefix}-mashmap.png")
    plt.savefig(f"{args.output_prefix}-mashmap.png")
    plt.cla()

    t_acc = alignment.t_acc_list[0]
    x, y, sat1 = dotplot.target_response_curve(t_acc)
    x2, y2, sat2 = dotplot.query_response_curve()

    plt.plot(x, y / max(y), "b-", label=f"target loss ({t_acc})")
    plt.plot(x2, y2 / max(y2), "g-", label="query loss")

    plt.xlabel("kb in genome contigs removed")
    plt.ylabel("fraction of alignments removed")
    plt.legend(loc="lower right")

    print(f"saving {args.output_prefix}-response.png")
    plt.savefig(f"{args.output_prefix}-response.png")
    plt.cla()

    slope = AlignmentSlopeDiagram(alignment)
    fig = slope.plot()

    print(f"saving {args.output_prefix}-alignplot.png")
    plt.savefig(f"{args.output_prefix}-alignplot.png")
    plt.cla()

    return 0


if __name__ == "__main__":
    sys.exit(main())
