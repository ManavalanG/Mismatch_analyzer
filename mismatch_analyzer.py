from Bio import AlignIO
from Bio import SeqIO
import errno


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
    query_pos_list = []        # query_pos_actual refers to actual positioning whereas query_pos_list to python 0-indices.
    query_pos_residues = []
    for n, _ in enumerate(query_pos_actual):
        query_pos_list.append(int(query_pos_actual[n]) - 1)
        query_pos_residues.append(reference.seq[query_pos_list[n]])

    # print 'Query sites requested: %s' % (str(query_pos_actual))
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
    
    return query_pos_list, query_pos_residues, query_pos_in_alignment, query_pos_in_alignment_actual, aligned_ids_list, pipes_max



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



def write_csv_out(query_pos_actual, query_pos_residues, aligned_ids_list, dict_alignment_residues, pipes_max,
                  id_delimiter, query_pos_in_alignment, aligned_data, input_file, Reference_file, method_name):
    Output_filename_csv = 'Mismatches_Tabulated.csv'
    Output_Path = ''
    Output_file_csv = Output_Path + Output_filename_csv

    try:
        Output_handle_mismatch_csv = open(Output_file_csv, 'w')
    except IOError as exception:
        print 'Error that resulted when opening csv file for writing: %s' % exception   # this happens in Windows OS
        if exception.errno != errno.EEXIST:
            temp = ("A file titled \n\t'%s' \nseems to be open in MS Excel. Close that file and try again!" % Output_filename_csv)
            print temp
            raise

    # Writes the title line for csv file and txt file
    Title_components_csv = ['S.No']
    for no in range(0, pipes_max):
        Title_components_csv += ['Title_' + str(no + 1)]
    Title_components_csv += ',Sequence Length'
    Title_components_csv += ',No. of Mismatches'
    # print 'Determined the title to be used in csv and tect files.'


    for no in range(0, len(query_pos_actual)):
        Title_components_csv += (',' + str(query_pos_actual[no]) + ' "' + query_pos_residues[no] + '"')
    Title_components_csv += '\n'
    # print 'Created title line for using in csv file'

    # Gethers residues at query sites and decides if they are mismatch or not.
    # print 'Starting to extract and write mismatch details into csv file'

    # sno_1 = 0
    # sno_2 = 0
    Mismatches_details = ''
    No_Mismatches_details = ''
    mismatched_seqs_total = 0
    match_seqs_total = 0
    mismatches_none = True  # This var is to get info if there is no mismatches noted at any of the sites tested
    for identifier in aligned_ids_list:
        mismatch_count = 0  # to count number of mismatches
        for char in dict_alignment_residues[identifier]:
            if char != '':
                mismatch_count += 1

        pipe_split = identifier.split(id_delimiter)
        if pipe_split[-1] == '':  # Accounts for if pipe is in the end of the id name
            pipe_split = pipe_split[:-1]

        if len(pipe_split) < pipes_max:  # makes sure all IDs have same number of naming components
            for no in range(0, pipes_max - len(pipe_split)):
                pipe_split.append('-na-')

        id_components = ''  # adds name components
        for no in range(0, pipes_max):
            if ',' in pipe_split[no] or '"' in pipe_split[no]:  # standardizes as per csv file format since spreadsheet software will misbehave when a cell value has comma or double quote in it
                pipe_split[no] = pipe_split[no].replace('"', '""')
                pipe_split[no] = ('"%s"' % pipe_split[no])

            id_components += (',' + pipe_split[no])

        residues_data = ''  # obtains residue in that position, mismatched or not.
        for no in range(0, len(query_pos_in_alignment)):
            res = dict_alignment_residues[identifier][no]
            if dict_alignment_residues[identifier][no] == '':  # If residue matches, equal sign is used to improved visibility
                res = '='
            if dict_alignment_residues[identifier][no] == '-':  # Star symbol is used instead of gap symbols to improved visibility
                res = '*'
            residues_data += (',' + res)
        residues_data += '\n'

        len_seq = None  # this part obtains the length of sequence
        for no in range(0, len(aligned_data)):
            if aligned_data._records[no].id == identifier:
                seq = str(aligned_data._records[no].seq)
                seq = seq.replace('-', '')
                len_seq = len(seq)
                break

        if len_seq is None:
            print "Error obtaining sequence length of '%s'." % identifier

        if mismatch_count > 0:  # Writes into output file if at least one mismatch found in  query positions
            # serial_number_1= str (sno_1 + 1)
            # Output_handle_mismatch_csv.write(str (serial_number_1) +id_components + ',' + str(len_seq) + ',' +
            # 						str (mismatch_count) + residues_data)
            Mismatches_details += (id_components + ',' + str(len_seq) + ',' + str(mismatch_count) + residues_data)
            mismatched_seqs_total += 1
            mismatches_none = False

        else:  # Gathers details of records without any mismatches in query positions
            # serial_number_2= str (sno_2 + 1)
            No_Mismatches_details += (id_components + ',' + str(len_seq) + ',' + str(mismatch_count) + residues_data)
            match_seqs_total += 1

    Output_handle_mismatch_csv.write('Sequences file used:    "%s"\n' % input_file)
    Output_handle_mismatch_csv.write('Reference file used:    "%s"\n' % Reference_file)
    Output_handle_mismatch_csv.write('Residue conservation method used:    %s\n\n' % method_name)
    # Output_handle_mismatch_csv.write('Reference sequence included in calculations:    %s\n\n' %ref_included)

    if ref_included:  # reference sequence included in calculations
        total_no_seqs = len(aligned_data)
        ref_str = 'including Reference sequence'
    else:
        total_no_seqs = len(aligned_data) - 1
        match_seqs_total -= 1
        ref_str = 'excluding Reference sequence'

    perc_match = float(match_seqs_total) * 100 / total_no_seqs
    perc_mismatch = float(mismatched_seqs_total) * 100 / total_no_seqs

    Output_handle_mismatch_csv.write('Number of sequences\n' +
                                     '       in alignment (%s):      %i\n' % (ref_str, total_no_seqs) +
                                     '       that have all residues matching (%s):   %i (%2.1f %%)\n'
                                     % (ref_str, match_seqs_total, perc_match) +
                                     '       that have at least one mismatching residue:      %i (%2.1f %%)\n\n\n'
                                     % (mismatched_seqs_total, perc_mismatch))

    Output_handle_mismatch_csv.write("*** Records that have mismatches in at least one of the query sites ***\n\n")
    Output_handle_mismatch_csv.write(Title_components_csv)

    if mismatches_none:
        Output_handle_mismatch_csv.write(",*** No mismatches at any of the sites requested ***\n")
    else:
        # sorts the data that have at least one mismatch
        csv_data = csv.reader(StringIO.StringIO(Mismatches_details))
        low = len(id_components.split(',')) + 1
        high = low + len(query_site) + 1

        if natural_sorting:
            sortedlist1 = natsorted(csv_data, key=operator.itemgetter(*range(low, high)))
        else:
            sortedlist1 = sorted(csv_data, key=operator.itemgetter(*range(low, high)))
        sno = 1
        for line in sortedlist1:
            Output_handle_mismatch_csv.write(str(sno) + ','.join(line) + '\n')
            sno += 1

    # sorts the data that have no mismatches at all
    csv_data = csv.reader(StringIO.StringIO(No_Mismatches_details))
    low = 2  # will have to be changed in the end to 1
    high = len(id_components.split(',')) + 1
    if natural_sorting:
        sortedlist2 = natsorted(csv_data, key=operator.itemgetter(*range(low, high)))
    else:
        sortedlist2 = sorted(csv_data, key=operator.itemgetter(*range(low, high)))

    # to find unique residue in query sites
    mismatch_log.debug('Unique residues at each requested site will be extracted')

    Dict_Unique_Residues = {k: [] for k in query_pos_actual}
    for no in range(0, len(query_pos_actual)):
        for res_list in dict_alignment_residues.values():
            Dict_Unique_Residues[query_pos_actual[no]].append(res_list[no])

        # this part will enable calculating each unique residue's count and fraction w/o including reference seq's residue.
        # This is so as to offer true calculation as reference seq' residue is always matching anyway.
        if '' in Dict_Unique_Residues[query_pos_actual[no]] and not ref_included:
            Dict_Unique_Residues[query_pos_actual[no]].remove('')

    # if in case user sets 'similar aa sets' empty in settings (either through gui or settings file), the script will
    # treat rest of analysis as if it is in dna mode (where similarity is not considered)
    if aa_set_str == '':
        protein_mode = False

    Output_handle_mismatch_csv.write("\n\n*** Unique residues seen at the query sites and their count. ***\n")
    if ref_included:  # tells whether reference sequence is involved or not in calculations
        Output_handle_mismatch_csv.write(
            "** (Note: Calculation includes Reference sequences's residue at that position) **\n\n")
    else:
        Output_handle_mismatch_csv.write(
            "** (Note: Calculation DOESN'T include Reference sequences's residue at that position) **\n\n")

    if protein_mode:
        Output_handle_mismatch_csv.write((
                                         ",Expected Residue, ,Identity_count,%% Identity,%% Conservation (%s),Unique residues' count and fraction\n") % method_name)
    else:
        Output_handle_mismatch_csv.write(",Expected Residue, ,Identity_count,% Identity,"
                                         "Unique residues' count and fraction\n")

    # aa_set = ['AVFPMILW', 'DE', 'RK', 'STYHCNGQ', 'avfpmilw', 'de', 'rk', 'styhcngq']
    unique_residues_line = ''
    unique_residues_list = []
    identity_count_list = []
    identity_perc_list = []
    # similarity_count_list = []
    similarity_perc_list = []

    # List 'conserve_score_list' will be used to write data in output files based on residue conservation method needed
    if conserve_method == 'amino_acid_grouping':
        conserve_score_list = similarity_perc_list  # both lists will change when one changes
    elif conserve_method == 'liu08_non_seq_weighted':
        conserve_score_list = liu08_simple_score_list
    elif conserve_method == 'liu08_seq_weighted':
        conserve_score_list = liu08_weighted_score_list

    # get stats about %Identity (and optionally %similarity) at column positions requested
    for no in range(0, len(query_pos_actual)):
        for item in range(0, len(Dict_Unique_Residues[query_pos_actual[no]])):
            if Dict_Unique_Residues[query_pos_actual[no]][item] == '':
                Dict_Unique_Residues[query_pos_actual[no]][item] = query_pos_residues[no]

        unique_res_detail = ','
        unique_res_detail += (str(query_pos_actual[no]) + ' "' + query_pos_residues[no] + '",')

        # this writes unique residues for each site in descending order of their count
        Dict_sorting_temp = {}
        for elem in sorted(set(Dict_Unique_Residues[query_pos_actual[no]])):
            Dict_sorting_temp[elem] = Dict_Unique_Residues[query_pos_actual[no]].count(elem)

        # this obtains identity count and percentage
        if query_pos_residues[no] in Dict_sorting_temp:
            identity_count = Dict_sorting_temp[query_pos_residues[no]]
            identity_count_list.append(identity_count)
            if not ref_included:
                identity_perc = float(identity_count) / (len(aligned_ids_list) - 1) * 100
            else:
                identity_perc = float(identity_count) / (len(aligned_ids_list)) * 100
            # this is to avoid Python rounding up, for eg., 99.999 to 100.0. This is a bit of hard coding but this
            # list will not be used anywhere else for calculation purposes. It'll be just used for writing in output
            if identity_perc > 99.94 and identity_perc < 100:
                identity_perc = 99.9
            identity_perc_list.append(identity_perc)
        else:
            identity_count = 0
            identity_count_list.append(identity_count)
            identity_perc = 0
            identity_perc_list.append(identity_perc)

        # this obtains residue conservation score (similarity count and percentage) by amino acid grouping method
        aa_similar = ''
        for group in aa_set:  # Note: aa_set is defined in settings text file and allowed to be changed via GUI
            if query_pos_residues[no] in group:
                aa_similar = group
                break

        # if residue conservation score needs to be determined based on amino acid grouping set
        if conserve_method == 'amino_acid_grouping':
            similarity_count = 0
            for key in Dict_sorting_temp:
                if key in aa_similar:
                    similarity_count += Dict_sorting_temp[key]
            # similarity_count_list.append(similarity_count)
            if not ref_included:
                similarity_perc = float(similarity_count) / (len(aligned_ids_list) - 1) * 100
            else:
                similarity_perc = float(similarity_count) / (len(aligned_ids_list)) * 100
            # this is to avoid Python rounding up, for eg., 99.999 to 100.0. This is a bit of hard coding but this
            # list will not be used anywhere else for calculation purposes. It'll be just used for writing in output
            if similarity_perc > 99.94 and similarity_perc < 100:
                similarity_perc = 99.9
            similarity_perc_list.append(similarity_perc)
        # method_name = "Amino acid Grouping" 			    # for csv output


        # calculates count of unique amino acidss in column and sorts them in ascending order
        unique_each = ''
        for key, value in sorted(Dict_sorting_temp.iteritems(), key=lambda (k, v): (v, k), reverse=True):
            perc = float(value) / (len(aligned_ids_list) - 1) * 100
            perc = '%3.1f' % perc
            unique_each += (',' + key + ': ' + str(value) + ' (' + str(perc) + '%)')

        # Collects data, which needs to be written into a output file, in a variable
        if protein_mode:  # residue conservation details added only if protein mode is used
            unique_res_detail += (',%i,%3.1f,%0.1f' % (identity_count, identity_perc, conserve_score_list[no]))
        else:
            unique_res_detail += (',%i,%3.1f' % (identity_count, identity_perc))
        unique_res_detail += unique_each
        unique_residues_line += (unique_res_detail + '\n')
        unique_residues_list.append(unique_each)  # for text file output

    mismatch_log.debug('Done extracting unique residues')

    Output_handle_mismatch_csv.write(unique_residues_line)
    Output_handle_mismatch_csv.write("\n\n*** Records that Do Not have mismatches at any of the query sites ***")
    # Output_handle_mismatch_csv.write('\n\n' + Title_components_csv + No_Mismatches_details)
    Output_handle_mismatch_csv.write('\n\n' + Title_components_csv)
    sno = 1
    for ele in sortedlist2:
        Output_handle_mismatch_csv.write(str(sno) + ','.join(ele) + '\n')
        sno += 1
    Output_handle_mismatch_csv.close()
    mismatch_log.info('Done writing mismatch details in to csv file mentioned above.')


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

    query_pos_list, query_pos_residues, query_pos_in_alignment, query_pos_in_alignment_actual, aligned_ids_list, pipes_max = get_query_sites_in_alignment(aligned_data, reference, ref_index_in_alignment, query_pos_actual, id_delimiter)

    dict_alignment_residues, dict_identifier_status = get_mismatch_info(aligned_ids_list, aligned_data, query_pos_in_alignment, query_pos_residues)

    write_csv_out(query_pos_actual, query_pos_residues, aligned_ids_list, dict_alignment_residues, pipes_max, id_delimiter, query_pos_in_alignment)