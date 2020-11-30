#! /usr/bin/env python
"""
TODO:
* only get filenames for accs that are not already in output
"""
import sys
import argparse
import urllib.request
import csv

from lxml import etree


def url_for_accession(accession):
    db, acc = accession.strip().split("_")
    if '.' in acc:
        number, version = acc.split(".")
    else:
        number, version = acc, '1'
    number = "/".join([number[p : p + 3] for p in range(0, len(number), 3)])
    url = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/{number}"

    with urllib.request.urlopen(url) as response:
        all_names = response.read()

    all_names = all_names.decode("utf-8")

    full_name = None
    for line in all_names.splitlines():
        name = line.split()[-1]
        db_, acc_, *_ = name.split("_")
        if db_ == db and acc_.startswith(acc):
            full_name = name
            break

    if full_name is None:
        return None
    else:
        url = "https" + url[3:]
        return (
            f"{url}/{full_name}/{full_name}_genomic.fna.gz",
            f"{url}/{full_name}/{full_name}_assembly_report.txt",
        )


def get_taxid_from_assembly_report(url):
    with urllib.request.urlopen(url) as response:
        content = response.read()

    content = content.decode("utf-8").splitlines()
    for line in content:
        if "Taxid:" in line:
            line = line.strip()
            pos = line.find("Taxid:")
            assert pos >= 0
            pos += len("Taxid:")
            taxid = line[pos:]
            taxid = taxid.strip()
            return taxid

    assert 0


def get_tax_name_for_taxid(taxid):
    tax_url = (
        f"https://www.ncbi.nlm.nih.gov/taxonomy/?term={taxid}&report=taxon&format=text"
    )
    with urllib.request.urlopen(tax_url) as response:
        content = response.read()

    root = etree.fromstring(content)
    notags = etree.tostring(root).decode("utf-8")
    if notags.startswith("<pre>"):
        notags = notags[5:]
    if notags.endswith("</pre>"):
        notags = notags[:-6]
    notags = notags.strip()

    return notags


def main():
    p = argparse.ArgumentParser()
    p.add_argument("accession")
    p.add_argument("-o", "--output")
    args = p.parse_args()

    fieldnames = ["acc", "genome_url", "assembly_report_url", "ncbi_tax_name"]
    if args.output:
        fp = open(args.output, "wt")
        w = csv.DictWriter(fp, fieldnames=fieldnames)
    else:
        w = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    w.writeheader()

    acc = args.accession

    genome_url, assembly_report_url = url_for_accession(acc)
    taxid = get_taxid_from_assembly_report(assembly_report_url)
    tax_name = get_tax_name_for_taxid(taxid)

    d = dict(
        acc=acc,
        genome_url=genome_url,
        assembly_report_url=assembly_report_url,
        ncbi_tax_name=tax_name,
    )

    w.writerow(d)
    print(f"retrieved for {acc} - {tax_name}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
