#! /usr/bin/env python
"""
Combine CSVs and sort by given field.
"""
import sys
import argparse
import csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sort-by', default=None)
    p.add_argument('--reverse', action='store_true')
    p.add_argument('csvs', nargs='+')
    args = p.parse_args()

    first_csv = args.csvs[0]
    csvs = args.csvs[1:]

    rows = []
    with open(first_csv, 'rt') as fp:
        r = csv.DictReader(fp)
        rows.extend(list(r))
        fieldnames = r.fieldnames

    if args.sort_by:
        sort_by = args.sort_by
    else:
        sort_by = fieldnames[0]

    try:
        float(rows[0][sort_by])
        key_fn = lambda x: float(x[sort_by])
    except ValueError:
        key_fn = lambda x: x[sort_by]

    for csvfile in csvs:
        with open(csvfile, 'rt') as fp:
            r = csv.DictReader(fp)
            if r.fieldnames != fieldnames:
                diff = set(r.fieldnames) ^ set(fieldnames)
                print(f"error! disjoint fieldnames b/t {first_csv} and {csvfile}: {str(diff)}", file=sys.stderr)
                sys.exit(-1)

            new_rows = list(r)
            rows.extend(new_rows)

    print(f'loaded {len(rows)} total. now sorting!', file=sys.stderr)

    rows.sort(key=key_fn, reverse=args.reverse)

    o = csv.DictWriter(sys.stdout, fieldnames)
    o.writeheader()
    for row in rows:
        o.writerow(row)

    return 0


if __name__ == '__main__':
    sys.exit(main())
