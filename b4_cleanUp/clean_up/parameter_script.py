import argparse
import os

# check if file exists
def is_file_valid(f):
    if not os.path.exists(f):
        raise argparse.ArgumentTypeError("File provided '%s' does not exist!" % f)
    else:
        return os.path.abspath(f)  # return an open file handle


def fn_get_args():
    # arg parser
    parser = argparse.ArgumentParser(description='Mismatch analyzer: Analyzes residue conservation using multiple sequence alignment')

    # required main arguments
    parser.add_argument('-r', metavar='ref', dest='ref_f', type=is_file_valid, help='Reference fasta file with only one sequence', required=True)
    parser.add_argument('-q', metavar='query', dest='query_f', type=is_file_valid, help='Query sequence fasta file')
    parser.add_argument('-a', metavar='aligned', dest='aligned_f', type=is_file_valid, help='Pre-aligned MSA fasta file')

    # parse the arguments
    args = parser.parse_args()

    reference_f = args.ref_f
    query_f = args.query_f
    aligned_f = args.aligned_f

    return reference_f, query_f, aligned_f


if __name__ == '__main__':
    reference_f, query_f, aligned_f = fn_get_args()
    # print reference_f, query_f, aligned_f