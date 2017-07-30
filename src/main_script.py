import os
import errno
from Bio import SeqIO
from parameter_script import fn_get_args
from run_clustalo import run_clustal_omega


# verifies if a directory exists in the path
def mkdir_if_not_exist(dir_path):
    try:
        os.makedirs(dir_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            pass


# Does following checks when reference and query files are provided:
# 1. Proper data type and files not being empty
# 2. Adds reference seq to query seq file if not already present
# 3. If already present, test reference sequence matches that in query file
def check_input_files(reference_f, query_f, out_dir):
    # read reference sequence file
    try:
        reference_seq = SeqIO.read(reference_f, "fasta")
    except Exception as e:
        print "Ran into problem reading reference sequence file:\n  '%s'" \
              'Error: %s' % (reference_f, e)
        exit()

    # read query sequence file
    try:
        query_seqs = list(SeqIO.parse(query_f, "fasta"))
    except Exception as e:
        print "Ran into problem reading query sequence file:\n  '%s'" \
              'Error: %s' % (query_f,e)
        exit()

    # testing reference sequence file
    if not len(reference_seq):
        print 'Reference sequence file seems to be empty. Make sure it is in FASTA format.' \
              'Try again after fixing this problem'
        exit()

    # testing query sequence file
    if not len(query_seqs):
        print 'Query sequences file seems to be empty. Make sure it is in FASTA format.' \
              'Try again after fixing this problem'
        exit()

    query_seqs_id_list = [item.id for item in query_seqs]


    # Verifies if Reference sequence is present in query sequences and adds to it, if non-existent.
    if reference_seq.id not in query_seqs_id_list:
        query_seqs.append(reference_seq)
        query_RefSeq_added_f = os.path.join(out_dir, os.path.splitext(os.path.basename(query_f))[0] + "_Reference_SeqAdded.fasta")
        with open(query_RefSeq_added_f, 'w') as add_ref_seq_handle:
            SeqIO.write(query_seqs, add_ref_seq_handle, 'fasta')

        print ("Reference seq is not present in file: '%s'. "
               "A new file '%s' is created with Reference Seq appended to it." % (query_f, query_RefSeq_added_f))

        # new reference seq added file becomes the query file
        query_f = query_RefSeq_added_f
    else:
        # print 'Reference seq is present in file: %s' % query_f

        # verify sequence in Reference file matches to that in Sequences file
        refseq_in_query = query_seqs[ query_seqs_id_list.index(reference_seq.id) ].seq
        if str(reference_seq.seq) != str(refseq_in_query):
            print "Reference sequence does not match to corresponding sequence in 'Sequences file'. " \
                    "Fix it and try again!"
            exit()
        # else:
        #     print "Reference sequence matches to corresponding sequence in Sequences file."

    return query_f


# def sddbd():
if __name__ == '__main__':
    # get input filenames
    # reference_f, query_f, aligned_f = fn_get_args()

    reference_f = 'test_data/Reference.fasta'
    query_f = 'test_data/query.fasta'


    # run clustal omega alignment
    if reference_f and query_f:
        out_dir = os.path.join(os.path.dirname(query_f), 'output')
        mkdir_if_not_exist(out_dir)

        # put reference and query file through some basic tests
        query_f = check_input_files(reference_f, query_f, out_dir)

        # align the sequences using clustal omega
        aligned_f = os.path.join(out_dir, 'aligned.fasta')
        run_clustal_omega(query_f, aligned_f)
