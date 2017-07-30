'''
Script to get arguments from user.

Author: 'Mana'valan Gajapathy
'''

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
    parser.add_argument('-r', metavar='--ref_seq', dest='ref_f', type=is_file_valid, help='Reference fasta file with only one sequence', required=True)
    parser.add_argument('-qs', metavar='--query_seq', dest='query_f', type=is_file_valid, help='Query sequence fasta file')
    parser.add_argument('-a', metavar='--aligned', dest='aligned_f', type=is_file_valid, help='Pre-aligned MSA fasta file')
    parser.add_argument('-qp', metavar='--query_pos', dest='query_pos_list',
                                help='Query residue positions (based on reference sequence). Use comma as delimiter. Example: "20,30".'
                                     'Use "all" if all residues are query positions', required=True)

    # parse the arguments
    args = parser.parse_args()

    reference_f = args.ref_f
    query_f = args.query_f
    alignment_f = args.aligned_f

    try:
        if args.query_pos_list.strip() == 'all':
            query_pos_list = 'all'
        else:
            query_pos_list =[int(x)  for x in args.query_pos_list.strip().split(',')]
    except Exception as E:
        print 'Error in query positions submitted. Use comma as delimiter. Eg: "20,30"'
        exit()

    # checks arg requirements
    if query_f and alignment_f:
        print 'Error. Arguments -qs and -a cannot be used together. Only one can be used at once.'
        exit()

    print 'Done reading arguments from user'
    return reference_f, query_f, alignment_f, query_pos_list


if __name__ == '__main__':
    reference_f, query_f, alignment_f, query_pos_list = fn_get_args()
    print reference_f, query_f, alignment_f, query_pos_list