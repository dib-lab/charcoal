#! /usr/bin/env python
import sys
import csv
import argparse
import copy
import collections


HEADERS="genomefile,brieftax,f_major,f_ident,f_removed,n_reason_1,n_reason_2,n_reason_3,refsize,ratio,clean_bp,clean_n,dirty_n,dirty_bp,missed_n,missed_bp,noident_n,noident_bp,taxguessed,taxprovided,comment".split(',')


def main(args):
    # collect CSVs
    all_rows = {}
    n_dups = 0
    for filename in args.summary_csvs:
        with open(filename, 'rt') as fp:
            r = csv.DictReader(fp)

            for row in r:
                g = row['genomefile']
                if g in all_rows:
                    n_dups += 1
                else:
                    all_rows[g] = row

    print(f'loaded {len(all_rows)} records from {len(args.summary_csvs)} files; ignored {n_dups} duplicates.')

    # summarize: reasons, etc.

    cleaned = {}
    reasons = collections.Counter()
    for g, row in all_rows.items():
        f_removed = row['f_removed']
        if not f_removed: f_removed = 0.0
        f_removed = float(f_removed)
        if f_removed > 0.0:
            row = copy.copy(row)
            row['f_removed'] = f_removed
            cleaned[g] = row

        if row['comment']:
            reasons[row['comment']] += 1

    print('Summary of comments:')
    for reason, count in reasons.most_common():
        print(f'   {count} rows - {reason}')
    print('')
        
    # output a sorted version, with only places where f_removed is set
    cleaned_list = sorted(cleaned.values(), key=lambda x: -x['f_removed'])
    print(f'{len(cleaned_list)} rows had some contamination.')

    if args.only_cleaned:
        with open(args.only_cleaned, 'wt') as fp:
            w = csv.DictWriter(fp, fieldnames=HEADERS)
            w.writeheader()
            for row in cleaned_list:
                w.writerow(row)

        print(f"wrote only these rows to '{args.only_cleaned}'")



def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('summary_csvs', nargs='+')
    p.add_argument('-C', '--only-cleaned', help='only rows that were cleaned')
    args = p.parse_args()

    main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(0)
