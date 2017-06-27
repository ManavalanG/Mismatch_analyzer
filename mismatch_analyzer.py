from Bio import AlignIO
from Bio import SeqIO
import logging


# Function to index each character (residue) in the Reference sequence
def indexing(string, res):
    index_list = []
    for index, char in enumerate(string):
        if char == res:
            index_list.append(index)
    return index_list


def get_ref_index_in_alignment(aligned_data, reference):
    # Checks if Reference sequence ID is present in alignment file.
    ref_index_in_alignment = None
    for index, _  in enumerate(aligned_data):
        if reference.id == aligned_data[index].id:
            ref_index_in_alignment = index
            break


    # Checks if Reference seq is present in MSA and if present, checks if they match or not
    if not ref_index_in_alignment:
        temp = 'Error: Alignment file provided does not have Sequence ID corresponding to your reference sequence. ' \
                'Fix it and try again!'
        print temp
        raise
    else:
        reference_seq_in_alignment = aligned_data[ref_index_in_alignment].seq
        reference_seq_in_alignment = str(reference_seq_in_alignment).replace('-', '')
        if str(reference.seq) == reference_seq_in_alignment:
            # print "Reference sequence matches to corresponding sequence in Alignment file provided."
            pass
        else:
            temp = "Reference sequence does not match to corresponding sequence in Alignment file provided. " \
                    "Fix it and try again!"
            print temp
            raise

    return ref_index_in_alignment



def get_query_sites_in_alignment(aligned_data, reference, ref_index_in_alignment, query_pos_actual, id_delimiter):
    # Indexes residues (characters) in reference sequence and also for reference sequence in alignment.
    # print 'Begins indexing residues in reference seq and also of reference sequence in alignment'
    reference_seq_aligned = str(aligned_data[ref_index_in_alignment].seq)

    # following line creates list of alphabets (both upper and lower).
    # Includes non-amino acid chars as well to support nucleotides other than A, T, G and C
    # works for both upper and lower case characters and considers non-amino acid and non-nucleotide residues as mismatch
    aa_list = map(chr, range(65, 91)) + map(chr, range(97,123))
    Dict_aa_Index_Reference_Seq = {k: [] for k in aa_list}
    Dict_aa_Index_Reference_Seq_inAlignment = {k: [] for k in aa_list}
    for aa in aa_list:
        Dict_aa_Index_Reference_Seq_inAlignment[aa] = indexing(reference_seq_aligned, aa)
        Dict_aa_Index_Reference_Seq[aa] = indexing(reference.seq, aa)
    # print 'Done indexing residues in Reference seq and also of Reference sequence in alignment'

    # This section is to get appropriate number of 'name components' in the title row.
    # Also this accounts even if the clustal record ids have excessive number of pipes without any info in them.
    aligned_ids_list = []
    pipes_total_list = []
    for record in aligned_data:
        aligned_ids_list.append(record.id)
        record_id_split = (record.id).split(id_delimiter)
        pipes_total_list.append(len(record_id_split))
    pipes_max = max(pipes_total_list)
    # print 'Done determining number of title components in the sequence record IDs.'

    # Extracts residues present at query sites in the reference sequence and also reference sequence in alignment
    # print 'Begins the part to extract residues at the requested query sites.'
    query_pos_list = []        # query_site_actual refers to actual positioning whereas query_pos_list to python 0-indices.
    query_pos_residues = []
    for n, _ in enumerate(query_pos_actual):
        query_pos_list.append(int(query_pos_actual[n]) - 1)
        query_pos_residues.append(reference.seq[query_pos_list[n]])

    # print 'Query sites requested: %s' % (str(query_site_actual))
    # print 'Residues at those query sites: %s' % str(query_pos_residues)


    # Finds the position of query site residues in reference's alignment sequence
    query_pos_in_alignment = []
    for aa in aa_list:
        for query_site in query_pos_list:
            if query_site in Dict_aa_Index_Reference_Seq[aa]:
                pos = Dict_aa_Index_Reference_Seq[aa].index(query_site)
                query_pos_in_alignment.append(Dict_aa_Index_Reference_Seq_inAlignment[aa][pos])

    query_pos_in_alignment = sorted(query_pos_in_alignment)

    # To count for the fact that python indexes from 0 instead of 1
    query_pos_in_alignment_actual = [(item + 1) for item in query_pos_in_alignment]

    # print "Corresponding sites in the alignment's Reference sequence : %s" % str(query_in_Alignment_actual)
    
    return query_pos_list, query_pos_residues, query_pos_in_alignment, query_pos_in_alignment_actual, aligned_ids_list



def get_mismatch_info(aligned_ids_list, aligned_data, query_pos_in_alignment, query_pos_residues):
    aa_set = ['AVFPMILW', 'DE', 'RK', 'STYHCNGQ', 'avfpmilw', 'de', 'rk', 'styhcngq']
    dict_alignment_residues = {k: [] for k in aligned_ids_list}
    dict_identifier_status = {}
    for record in aligned_data:
        dict_identifier_status[record.id] = []  # to assist in color coding IDs in html output
        # for no in range(0, len(query_pos_in_alignment)):
        for no, _ in enumerate(query_pos_in_alignment):

            # this obtains similarity count and percentage
            aa_similar = ''
            for group in aa_set:  # aa_set is defined in settings text file and allowed to be changed via GUI
                if query_pos_residues[no] in group:
                    aa_similar = group
                    break

            if record.seq[query_pos_in_alignment[no]] == query_pos_residues[no]:
                dict_alignment_residues[record.id].append('')
                dict_identifier_status[record.id].append('match')
            elif record.seq[query_pos_in_alignment[no]] in aa_similar:
                dict_identifier_status[record.id].append('similar')
                dict_alignment_residues[record.id].append(record.seq[query_pos_in_alignment[no]])
            else:
                dict_alignment_residues[record.id].append(record.seq[query_pos_in_alignment[no]])
                dict_identifier_status[record.id].append('mismatch')

    for key in dict_identifier_status:  # to assist in color coding IDs in html output
        dict_identifier_status[key] = list(set(dict_identifier_status[key]))

    return dict_alignment_residues, dict_identifier_status






# def xx():
if __name__ == '__main__':
    id_delimiter = '|'

    aligned_f = 'aligned.fasta'
    ref_f = 'Reference.fasta'
    query_pos_actual = [2,4]
    query_pos_actual = sorted(query_pos_actual)

    aligned_data = AlignIO.read(aligned_f, "fasta")
    reference = SeqIO.read(ref_f, "fasta")

    ref_index_in_alignment = get_ref_index_in_alignment(aligned_data, reference)

    query_pos_list, query_pos_residues, query_pos_in_alignment, query_pos_in_alignment_actual, aligned_ids_list = get_query_sites_in_alignment(aligned_data, reference, ref_index_in_alignment, query_pos_actual, id_delimiter)

    dict_alignment_residues, dict_identifier_status = get_mismatch_info(aligned_ids_list, aligned_data, query_pos_in_alignment, query_pos_residues)
