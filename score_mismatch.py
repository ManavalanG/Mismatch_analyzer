from Bio import AlignIO
from Bio import SeqIO
import pandas as pd
import pandascharm as pc
import numpy as np

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 100000)
pd.set_option('max_colwidth', 100000)

pd.options.mode.chained_assignment = None  # default='warn'

def example_df():
    # aa = pd.DataFrame({'a': [1,2], 'b':[3,4], 'c':[5,2]})
    # aa = pd.DataFrame({'a': ['1','2', '7'], 'b':['3','4', '8'], 'c':['2','2', '8']})
    aa = pd.DataFrame({'a': ['a','c', 'c'], 'b':['d','e', 'f'], 'c':['g','h', 'c']})
    print aa
    # compared = aa.apply(lambda x:x.isin(aa['c']))
    # print compared

    print aa.apply(lambda x:x==aa['a'])

    for i, x in aa.iteritems():
        print x


# converts fasta files to pandas dataframes
def fasta_to_pandas(f):
    data_fasta = AlignIO.read(f, "fasta")
    data_pd = pc.from_bioalignment(data_fasta)       # converted to pandas dataframe

    # convert categories datatype to objects
    for col_name in data_pd.columns:
        data_pd[col_name] = data_pd[col_name].astype(object)

    return data_pd


# Makes following checks with respect to reference sequence:
# 1. Does reference file has only one sequence? Is its ID part of Sequences file?
# 2. Does seq in Reference file and corresponding seq in Sequence file match?
def check_ref_seq_in_alignment(reference_pd, alignment_pd):
    ref_id = list(reference_pd.columns)
    alignment_IDs = list(alignment_pd.columns)

    if len(ref_id) > 1:
        print 'Reference file has more than one sequence. Only one sequence is allowed'
        raise
    else:
        ref_id = ref_id[0]

    if ref_id not in alignment_IDs:
        print "Reference sequence's ID not found in Query sequences file. Fix it and try again."
        raise


    # ref_seq_in_alignment = alignment_pd.loc[(alignment_pd[ref_id] != '-'), ref_id].reset_index(drop=True)
    # ref_seq_in_alignment = alignment_pd.loc[(alignment_pd[ref_id] != '-'), ref_id].reset_index()
    ref_seq_in_alignment = alignment_pd.loc[(alignment_pd[ref_id] != '-'), [ref_id]]
    ref_seq_in_alignment['index_in_alignment'] = ref_seq_in_alignment.index
    ref_seq_in_alignment = ref_seq_in_alignment.reset_index(drop=True)
    # print ref_seq_in_alignment


    if len(ref_seq_in_alignment) != len(reference_pd[ref_id]):
        print "Reference sequence in Query sequences file does not match with sequence in Reference File. Also, their length don't match. Fix it and try again."
        raise
    elif not ref_seq_in_alignment[ref_id].equals(reference_pd[ref_id]):
        print 'Reference sequence in Query sequences file does not match with sequence in Reference File. Fix it and try again.'
        raise

    return ref_id, ref_seq_in_alignment


# Based on query positions in input, make pandas dataframes - reading reference file and alignment file
def get_query_residue_info(query_pos_actual, reference_pd, ref_id, ref_seq_in_alignment):
    query_pos_0indexed = [x-1  for x in query_pos_actual]

    # get query residues
    query_ref_pd = reference_pd.iloc[query_pos_0indexed, :]
    query_ref_pd['actual_pos'] = query_ref_pd.index + 1

    query_info_in_alignment_pd = ref_seq_in_alignment.iloc[query_ref_pd.index, :]
    query_info_in_alignment_pd.rename(columns={ref_id : 'query_residues', 'index_in_alignment': 'pos_in_alignment'}, inplace=True)
    query_info_in_alignment_pd['query_pos_in_ref'] = query_info_in_alignment_pd.index
    query_info_in_alignment_pd['query_pos_actual'] = query_info_in_alignment_pd.index + 1
    query_info_in_alignment_pd['query_label_csv'] = query_info_in_alignment_pd['query_pos_actual'].astype(str) + ' "' + query_info_in_alignment_pd['query_residues'] + '"'
    query_info_in_alignment_pd['query_label_html'] = query_info_in_alignment_pd['query_pos_actual'].astype(str) + '<BR>' + query_info_in_alignment_pd['query_residues']
    query_info_in_alignment_pd['query_line_html'] = '<a name="' + query_info_in_alignment_pd['query_pos_actual'].astype(str) \
                                                    + '"><vert>' + query_info_in_alignment_pd['query_pos_actual'].astype(str) + '</vert></a>'

    return query_ref_pd, query_info_in_alignment_pd


# Get DF of alignment at query positions. Match query seqs to reference seq, and return a boolean DF.
def booleante(alignment_pd, query_pos_aligned, ref_id):
    query_region_pd = alignment_pd.loc[query_pos_aligned, :]

    bool_query_region_pd = query_region_pd.apply(lambda x:x==query_region_pd[ref_id])
    inverse_bool_query_region_pd = query_region_pd.apply(lambda x:x!=query_region_pd[ref_id])    # True and False are switched

    return query_region_pd, bool_query_region_pd, inverse_bool_query_region_pd


def split_id_by_delimiter(results_csv_pd):
    title_df = pd.DataFrame(index=results_csv_pd.index)
    title_df['title'] = results_csv_pd.index
    title_df = title_df['title'].str.split('|', expand=True)  # split using a delimiter symbol
    title_df.columns = pd.Series(['ID'] * 2).map(str) + '_' + pd.Series(title_df.columns + 1).map(str)

    return title_df


# adds serial number to DF in the front
def add_serial_no(df):
    df.insert(0, 'S.No', [x for x in range(1, len(df) + 1)])


def summarize_unique_residues(query_info_in_alignment_pd, query_region_pd, ref_id):
    # print query_info_in_alignment_pd
    query_region_pd = query_region_pd.T     # transpose and save
    query_region_pd.columns = query_info_in_alignment_pd['query_pos_in_ref']    # rename columns from pos-in-alignment to pos-in-reference-seq
    # print query_region_pd
    query_region_pd.drop(ref_id, axis=0, inplace=True)
    # print query_region_pd

    unique_pd = pd.DataFrame(index=query_info_in_alignment_pd.index)
    unique_pd['Expected Residue'] = query_info_in_alignment_pd['query_label_csv']
    for query_pos in query_region_pd.columns:
        a = query_region_pd[query_pos].value_counts()
        a_perc = (query_region_pd[query_pos].value_counts(normalize=True)*100).round(1)

        value_count_pos_pd = pd.DataFrame()
        value_count_pos_pd['count'] = a
        value_count_pos_pd['perc'] = a_perc
        value_count_pos_pd['residue'] = a_perc.index

        temp = ''
        for pos in value_count_pos_pd.index:
            temp += '%s: %i (%0.1f%%), ' % (pos, value_count_pos_pd.loc[pos, 'count'], value_count_pos_pd.loc[pos, 'perc'])
        temp = temp.strip(', ')

        query_residue = query_info_in_alignment_pd.loc[query_pos, 'query_residues']

        if query_residue in a:
            unique_pd.loc[query_pos, 'Identity_count'] = a[query_residue]
            unique_pd.loc[query_pos, '% Identity'] = a_perc[query_residue]
        else:
            unique_pd.loc[query_pos, 'Identity_count'] = 0
            unique_pd.loc[query_pos, '% Identity'] = 0

        unique_pd.loc[query_pos, "Unique residues' count and fraction"] = temp

    # unique_pd.to_csv('valuecount.tsv', sep='\t', index=False)
    return unique_pd


def prepare_for_csv_output(query_region_pd, bool_query_region_pd, inverse_bool_query_region_pd, seq_length_alignment_pd, query_info_in_alignment_pd):
    results_csv_pd = query_region_pd.T * inverse_bool_query_region_pd.T

    for col_name in results_csv_pd.columns:
        results_csv_pd.ix[(results_csv_pd[col_name]==''), col_name] = '='

    dict_name = query_info_in_alignment_pd.set_index('pos_in_alignment')['query_label_csv'].to_dict()
    results_csv_pd.rename(columns=dict_name, inplace=True)

    results_csv_pd.insert(0, 'mismatch_count', bool_query_region_pd[bool_query_region_pd==0].count())
    results_csv_pd.insert(0, 'seq_length_alignment_pd', seq_length_alignment_pd)

    title_df = split_id_by_delimiter(results_csv_pd)
    results_csv_pd = pd.concat([title_df, results_csv_pd], axis=1)

    return results_csv_pd


# write results to a csv file
def write_csv_out(results_csv_pd, unique_pd):
    match_only_pd = results_csv_pd[results_csv_pd['mismatch_count'] == 0]
    add_serial_no(match_only_pd)

    mismatch_only_pd = results_csv_pd[results_csv_pd['mismatch_count'] != 0]
    mismatch_only_pd.sort_values('mismatch_count', ascending=False)
    add_serial_no(mismatch_only_pd)

    # unique_pd = summarize_unique_residues(query_info_in_alignment_pd, query_region_pd, ref_id)

    csv_outfile = 'csv_out.tsv'
    with open(csv_outfile, 'w') as csv_handle:
        csv_handle.write('\n*** Records that have mismatches in at least one of the query sites ***\n')
        mismatch_only_pd.to_csv(csv_handle, sep='\t', index=False)

        csv_handle.write('\n\n*** Unique residues seen at the query sites and their count. ***\n')
        unique_pd.to_csv(csv_handle, sep='\t', index=False)

        csv_handle.write('\n\n*** Records that Do Not have mismatches at any of the query sites ***\n')
        match_only_pd.to_csv(csv_handle, sep='\t', index=False)


def write_html(alignment_pd, id_maxLength, query_info_in_alignment_pd):
    alignment_len = len(alignment_pd.columns)
    block_size = 60
    if alignment_len % block_size:
        no_of_blocks = (alignment_len/block_size + 1)
    else:
        no_of_blocks = alignment_len/block_size

    with open('aa.txt', 'w') as handle:

        query_line_pd = pd.DataFrame(index=alignment_pd.columns)
        query_line_pd['trial'] = ' '
        # print query_info_in_alignment_pd
        query_line_pd.loc[(query_info_in_alignment_pd['pos_in_alignment']), 'trial'] = list(query_info_in_alignment_pd['query_line_html'])
        query_line_pd['trial'].fillna(' ')


        alignment_pd = pd.concat([query_line_pd.T, alignment_pd], axis=0)
        html_blocks_pd = pd.DataFrame(index=alignment_pd.index)
        for block_no in range(0, no_of_blocks):
            min = block_no * block_size
            max = (block_no +1) * block_size
            if alignment_len <= max:
                max = alignment_len

            html_blocks_pd[block_no] = alignment_pd[range(min, max)].apply(lambda x: ''.join(x), axis=1)
            # print html_blocks_pd

            for i, item in html_blocks_pd[block_no].iteritems():
                # print i
                handle.write('%s    <spaced>%s</spaced>\n' % (i.ljust(id_maxLength), item))

            handle.write('\n\n')

    # query_pos_index(html_blocks_pd[block_no])


def process_for_html(bool_query_region_pd, alignment_pd, query_pos_aligned, query_info_in_alignment_pd):
    # alignment_len = len(alignment_pd)

    bool_query_region_pd = bool_query_region_pd.T
    alignment_pd = alignment_pd.T
    id_maxLength = alignment_pd.index.str.len().max()

    match_prefix = '<strong><span style="background-color: #42DB33">'
    match_suffix = '</span></strong>'
    mismatch_prefix = '<strong><span style="background-color: #F73D94">'

    for pos in query_pos_aligned:
        prefix = bool_query_region_pd[pos].apply(lambda x:match_prefix if x else mismatch_prefix)
        suffix = bool_query_region_pd[pos].apply(lambda x:match_suffix)

        alignment_pd[pos] = prefix + alignment_pd[pos] + suffix

    write_html(alignment_pd, id_maxLength, query_info_in_alignment_pd)



# def xxxx():
if __name__ == '__main__':
    query_pos_actual = [1,3,13]

    reference_f = 'test_data/Reference.fasta'
    alignment_f = 'test_data/aligned.fasta'

    reference_pd = fasta_to_pandas(reference_f)
    alignment_pd = fasta_to_pandas(alignment_f)

    # check reference sequence criteria pass prelimiary test
    ref_id, ref_seq_in_alignment = check_ref_seq_in_alignment(reference_pd, alignment_pd)

    # process query pos and residue info in reference seq and in alignment
    query_ref_pd, query_info_in_alignment_pd = get_query_residue_info(query_pos_actual, reference_pd, ref_id, ref_seq_in_alignment)

    # get query positions in alignment
    query_pos_aligned = list(query_info_in_alignment_pd['pos_in_alignment'])

    # get query region of alignment
    query_region_pd, bool_query_region_pd, inverse_bool_query_region_pd = booleante(alignment_pd, query_pos_aligned, ref_id)

    # get length of sequences in alignment
    seq_length_alignment_pd = alignment_pd[alignment_pd != '-'].count()

    # summarize residues at query positions
    unique_pd = summarize_unique_residues(query_info_in_alignment_pd, query_region_pd, ref_id)

    # process data for csv output writing
    results_csv_pd = prepare_for_csv_output(query_region_pd, bool_query_region_pd, inverse_bool_query_region_pd, seq_length_alignment_pd, query_info_in_alignment_pd)

    # write csv output
    write_csv_out(results_csv_pd, unique_pd)

    process_for_html(bool_query_region_pd, alignment_pd, query_pos_aligned, query_info_in_alignment_pd)
