import sys
import argparse
import Bio.SeqIO

import cps

def main():
    parser = argparse.ArgumentParser(
        description='Calculate codon pair score/bias.')

    parser.add_argument('fna_filename',
        metavar='fna_file',
        type=str,
        default=sys.stdin,
        help='FASTA-formatted nucleotide file of the Coding Sequences'
             '(stdin by default)')

    parser.add_argument('--ref',
        metavar='reference_file',
        type=str,
        default=None,
        help='Calculate reference pairs for provided file')

    parser.add_argument('--output',
        metavar='output_file',
        type=str,
        default=sys.stdout,
        help='Specifies name of file to write to (stdout by default)')

    parser.add_argument('--calc_reference',
        action='store_true',
        help='Calculate reference pairs for provided file')

    args = parser.parse_args()

    if args.ref is not None:
        cps.codon_pair_ref(args.ref)

    if isinstance(args.fna_filename, str):
        with open(args.fna_filename, 'r+') as fh:
            seq_records = list(Bio.SeqIO.parse(fh, "fasta"))
    else:
        print("Reading sequences from stdin", file=sys.stderr)
        seq_records = Bio.SeqIO.parse(args.fna_file, "fasta")

    if args.calc_reference:
        cps.write_reference(seq_records, args.output)
    else:
        for record in seq_records:
            print(cps.calc_cpb(str(record.seq)))

if __name__ == '__main__':
    main()
