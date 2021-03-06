'''
Reads MSA and analyzes mismatches at query positions which are in reference to residue positions in reference sequence.
Writes results in CSV and HTML.

Author: 'Mana'valan Gajapathy
'''

from Bio import AlignIO
import pandas as pd
import pandascharm as pc
from shutil import copyfile
import os

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 100000)
pd.set_option('max_colwidth', 100000)

pd.options.mode.chained_assignment = None  # default='warn'


# converts fasta files to pandas dataframes
def fasta_to_pandas(f):
    data_fasta = AlignIO.read(f, "fasta")
    data_pd = pc.from_bioalignment(data_fasta)       # converted to pandas dataframe

    # convert categories datatype to objects
    for col_name in data_pd.columns:
        data_pd[col_name] = data_pd[col_name].astype(object)

    return data_pd


# Makes following checks with respect to reference sequence:
# 1. Does reference file have only one sequence? Is its ID part of Sequences file?
# 2. Does seq in Reference file and corresponding seq in Sequence file match?
def check_ref_seq_in_alignment_pd(reference_pd, alignment_pd):
    ref_id = list(reference_pd.columns)
    alignment_IDs = list(alignment_pd.columns)

    if len(ref_id) > 1:
        print 'Reference file has more than one sequence. Only one sequence is allowed'
        exit()
    else:
        ref_id = ref_id[0]

    if ref_id not in alignment_IDs:
        print "Reference sequence's ID not found in Query sequences file. Fix it and try again."
        exit()


    # ref_seq_in_alignment_pd = alignment_pd.loc[(alignment_pd[ref_id] != '-'), ref_id].reset_index(drop=True)
    # ref_seq_in_alignment_pd = alignment_pd.loc[(alignment_pd[ref_id] != '-'), ref_id].reset_index()
    ref_seq_in_alignment_pd = alignment_pd.loc[(alignment_pd[ref_id] != '-'), [ref_id]]
    ref_seq_in_alignment_pd['index_in_alignment'] = ref_seq_in_alignment_pd.index
    ref_seq_in_alignment_pd = ref_seq_in_alignment_pd.reset_index(drop=True)
    # print ref_seq_in_alignment_pd


    if len(ref_seq_in_alignment_pd) != len(reference_pd[ref_id]):
        print "Reference sequence in Query sequences file does not match with sequence in Reference File. Also, their length don't match. Fix it and try again."
        exit()
    elif not ref_seq_in_alignment_pd[ref_id].equals(reference_pd[ref_id]):
        print 'Reference sequence in Query sequences file does not match with sequence in Reference File. Fix it and try again.'
        exit()

    return ref_id, ref_seq_in_alignment_pd


# Based on query positions in input, make pandas dataframes - reading reference file and alignment file
def get_query_residue_info(query_pos_actual, reference_pd, ref_id, ref_seq_in_alignment_pd):
    query_pos_0indexed = [x-1  for x in query_pos_actual]

    # get query residues
    query_ref_pd = reference_pd.iloc[query_pos_0indexed, :]
    query_ref_pd['actual_pos'] = query_ref_pd.index + 1

    query_info_in_alignment_pd = ref_seq_in_alignment_pd.iloc[query_ref_pd.index, :]
    query_info_in_alignment_pd.rename(columns={ref_id : 'query_residues', 'index_in_alignment': 'pos_in_alignment'}, inplace=True)
    query_info_in_alignment_pd['query_pos_in_ref'] = query_info_in_alignment_pd.index
    query_info_in_alignment_pd['query_pos_actual'] = query_info_in_alignment_pd.index + 1
    query_info_in_alignment_pd['query_label_csv'] = query_info_in_alignment_pd['query_pos_actual'].astype(str) + ' "' + query_info_in_alignment_pd['query_residues'] + '"'
    query_info_in_alignment_pd['table_label_html'] = query_info_in_alignment_pd['query_pos_actual'].astype(str) + '<BR>' + query_info_in_alignment_pd['query_residues']
    query_info_in_alignment_pd['query_label_html'] = query_info_in_alignment_pd['query_pos_actual'].apply(lambda x:'<a name="%s"><vert>%s</vert></a>' % (x,'<br />'.join(list(str(x)))))

    return query_ref_pd, query_info_in_alignment_pd


# Get DF of alignment at query positions. Match query seqs to reference seq, and return a boolean DF.
def booleante(alignment_pd, query_pos_aligned, ref_id):
    query_region_pd = alignment_pd.loc[query_pos_aligned, :]

    bool_query_region_pd = query_region_pd.apply(lambda x:x==query_region_pd[ref_id])
    # inverse_bool_query_region_pd = query_region_pd.apply(lambda x:x!=query_region_pd[ref_id])    # True and False are switched
    inverse_bool_query_region_pd = ~bool_query_region_pd    # True and False are switched

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
        # print value_count_pos_pd

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


def prepare_for_csv_output(query_region_pd, bool_query_region_pd, inverse_bool_query_region_pd, seq_length_alignment_s, query_info_in_alignment_pd):
    results_csv_pd = query_region_pd.T * inverse_bool_query_region_pd.T

    for col_name in results_csv_pd.columns:
        results_csv_pd.ix[(results_csv_pd[col_name]==''), col_name] = '='

    dict_name = query_info_in_alignment_pd.set_index('pos_in_alignment')['query_label_csv'].to_dict()
    results_csv_pd.rename(columns=dict_name, inplace=True)

    results_csv_pd.insert(0, 'mismatch_count', bool_query_region_pd[bool_query_region_pd==0].count())
    results_csv_pd.insert(0, 'seq_length_alignment_s', seq_length_alignment_s)

    title_df = split_id_by_delimiter(results_csv_pd)
    results_csv_pd = pd.concat([title_df, results_csv_pd], axis=1)

    return results_csv_pd


# write results to a csv file
def write_csv_out(results_csv_pd, unique_pd, bool_query_region_pd, reference_f, query_seq_f, alignment_f, out_dir):
    match_only_pd = results_csv_pd[results_csv_pd['mismatch_count'] == 0]
    add_serial_no(match_only_pd)

    mismatch_only_pd = results_csv_pd[results_csv_pd['mismatch_count'] != 0]
    mismatch_only_pd.sort_values('mismatch_count', ascending=False)
    add_serial_no(mismatch_only_pd)

    csv_outfile = os.path.join(out_dir, 'csv_out.tsv')
    with open(csv_outfile, 'w') as csv_handle:
        csv_handle.write('Reference sequences file used:    "%s"\n'
                         'Alignment file:    "%s"\n' % (reference_f, alignment_f))
        if query_seq_f:
            csv_handle.write('Query sequences file used:    "%s"\n' % query_seq_f)
        csv_handle.write('\n')


        # summarize mismatch analysis for output writing
        compare_length_pd = bool_query_region_pd.sum() == len(unique_pd)
        temp_1 = compare_length_pd.sum() - 1
        temp_2 = len(compare_length_pd) - 1 - temp_1
        temp_1_perc = float(temp_1) / (len(compare_length_pd)-1) * 100
        temp_2_perc = float(temp_2) / (len(compare_length_pd)-1) * 100

        mismsatch_summary_info = ('Number of sequences (excluding Reference sequence)\n'
                         '       in alignment:                               %i\n'
                         '       that have all residues matching:            %i (%0.1f %%)\n'
                         '       that have at least one mismatching residue: %i (%0.1f %%)\n\n'
                         % (len(results_csv_pd)-1, temp_1, round(temp_1_perc, 1), temp_2, temp_2_perc))
        csv_handle.write(mismsatch_summary_info)


        csv_handle.write('\n*** Records that have mismatches in at least one of the query sites ***\n')
        mismatch_only_pd.to_csv(csv_handle, sep='\t', index=False)

        csv_handle.write('\n\n*** Unique residues seen at the query sites and their count. ***\n')
        unique_pd.to_csv(csv_handle, sep='\t', index=False)

        csv_handle.write('\n\n*** Records that Do Not have mismatches at any of the query sites ***\n')
        match_only_pd.to_csv(csv_handle, sep='\t', index=False)

    return mismsatch_summary_info


def html_chart_text(info_lines, width_canvas, query_site_actual, identity_perc_list, out_html_handle):
    text_for_chart = ('\n</HEAD>\n<BODY LANG="en-US" DIR="LTR">\n'
                      '<PRE CLASS="western">\n%s'
                      '<div>\n  <canvas id="canvas_bar" height="350", width = "%i" ></canvas>\n'
                      '<div id="legend"></div>\n</div>\n\n<script>\n'
                      '\tvar barChartData = \n\t{\n'
                      '\t\tlabels : %s,\n'
                      '\t\tdatasets : [\n'
                      '\t\t{	fillColor : "rgba(51,51,255,0.8)",\n\t\t\ttitle: "%% Identity",\n'
                      '\t\t\tdata : %s  },\n\n'
                      % (info_lines, width_canvas, query_site_actual, identity_perc_list))

    out_html_handle.write(text_for_chart)
    out_html_handle.write('\t\t]\n\t}\n')


def html_table_text(query_site_residues, query_site_actual, identity_perc_list, out_html_handle):
    limit_col = 20  # this variable determines how many columns are allowed in each table
    no_of_tables = (len(query_site_residues) - 1 + limit_col) / limit_col
    limit_col_temp = limit_col

    row1 = []
    for n in range(0, len(query_site_actual)):
        temp = '<a href="#%i">%i<BR>%s</a>' % (
            query_site_actual[n], query_site_actual[n], query_site_residues[n])  # hyperlink included
        row1.append(temp)
    row2 = []
    for n in range(0, len(identity_perc_list)):
        row2.append("%3.1f" % identity_perc_list[n])

    rows = [row1, row2]
    heading = ['Position', '% Identity', '% Conservation']
    for table_no in range(0, no_of_tables):
        table_string = '<TABLE border="1" cellpadding="0" cellspacing="0" HEIGHT="90px" style="table-layout:fixed">\n'
        if table_no < (no_of_tables - 1):
            end_col = limit_col_temp
        else:
            end_col = len(query_site_residues)

        for no in range(0, len(rows)):
            if not no:
                table_string += ('\t<TR ALIGN="CENTER" BGCOLOR="#E3E3E0">\n')
            else:
                table_string += ('\t<TR ALIGN="CENTER">\n')

            table_string += ('\t\t<TH scope="row" width=130>%s</TH>\n' % heading[no])
            for ele in range((limit_col_temp - limit_col), end_col):
                table_string += ('\t\t<TD width=60>%s</TD>\n' % rows[no][ele])

            table_string += '\t</TR>\n'
        table_string += ('</TABLE>\n')
        out_html_handle.write(table_string)
        limit_col_temp += limit_col



def html_color_legend(match_color, similar_color, mismatch_color, out_html_handle):
    temp = '-' * 100 + '\n\n'
    temp += '<strong> Color           For Residues                      For Sequence IDs  </strong>\n\n'
    temp += (
    '<span style="background-color: %s">        </span>      Matching                        Seq. has no mismatches at all\n\n' % match_color)

    protein_mode = False
    if protein_mode:
        temp += (
        '<span style="background-color: %s">        </span>      Mismatching but <ins>similar</ins>         Seq. has at least one mismatch but <ins>similar</ins> residue\n\n' % similar_color)
        temp += (
        '<span style="background-color: %s">        </span>      Mismatching and <ins>dissimilar</ins>      Seq. has at least one mismatch and <ins>dissimilar</ins> residue\n\n' % mismatch_color)

    else:
        temp += (
        '<span style="background-color: %s">        </span>      Mismatching                     Seq. has at least one mismatching residue\n\n' % mismatch_color)

    temp += "<span style='background-color: #00FFFF'>        </span>      Not applicable                  Reference Sequence's ID\n\n"

    out_html_handle.write(temp + '-' * 100 + '\n\n')


# get length of df excluding gaps
def get_length_so_far(df):
    return df[df!='-'].count(axis=1)


def html_alignment_text(alignment_highlighted_pd, alignment_pd, id_maxLength, colored_id_pd, query_info_in_alignment_pd, out_html_handle):
    alignment_len = len(alignment_highlighted_pd.columns)
    block_size = 60
    if alignment_len % block_size:
        no_of_blocks = (alignment_len/block_size + 1)
    else:
        no_of_blocks = alignment_len/block_size

    # with open('aa.txt', 'w') as handle:
    if True:
        # for labelling query residue positions above alignment in html output
        query_line_pd = pd.DataFrame(index=alignment_highlighted_pd.columns)
        query_line_pd['pos_label'] = ' '
        query_line_pd.loc[(query_info_in_alignment_pd['pos_in_alignment']), 'pos_label'] = list(query_info_in_alignment_pd['query_label_html'])
        query_line_pd['pos_label'].fillna(' ')


        html_blocks_pd = pd.DataFrame(index=alignment_highlighted_pd.index)
        for block_no in range(0, no_of_blocks):
            min = block_no * block_size
            max = (block_no +1) * block_size
            if alignment_len <= max:
                max = alignment_len

            # finds the length of sequence so far excluding gaps
            length_s = get_length_so_far(alignment_pd.T[range(0, max)])

            # residue positions label line
            label_line = '%s    <spaced>%s</spaced>\n' % (''.ljust(id_maxLength), ''.join(query_line_pd.loc[min:max, 'pos_label']))
            out_html_handle.write(label_line)

            # get required block of alignment
            html_blocks_pd[block_no] = alignment_highlighted_pd[range(min, max)].apply(lambda x: ''.join(x), axis=1)

            for i, item in html_blocks_pd[block_no].iteritems():
                out_html_handle.write('%s    <spaced>%s</spaced>   %i\n' % (colored_id_pd['ID'][i].ljust(id_maxLength), item, length_s[i]))

            out_html_handle.write('\n\n')



def html_text_out(alignment_highlighted_pd, alignment_pd, id_maxLength, query_info_in_alignment_pd, unique_pd, colored_id_pd, match_color,
                  similar_color, mismatch_color, mismsatch_summary_info, out_dir):
    html_template_head_f = 'Template_html5_ChartNewjs_head.txt'
    html_template_tail_f = 'Template_html5_ChartNewjs_tail.txt'
    html_outfile = os.path.join(out_dir, 'html_out.html')

    # write starting block of html file to outfile
    copyfile(html_template_head_f, html_outfile)

    with open(html_outfile, 'a') as out_html_handle, open(html_template_tail_f, 'r') as html_tail_handle:
        info_lines = mismsatch_summary_info
        info_lines = info_lines.replace(' matching:', ' <ins>matching</ins>:')
        info_lines = info_lines.replace(' mismatching ', ' <ins>mismatching</ins> ')
        info_lines += '-' * 100 + '\n\n'

        if (63 * len(unique_pd)) < 325:
            width_canvas = 325
        else:
            width_canvas = 63 * len(unique_pd)

        query_site_actual = list(query_info_in_alignment_pd['query_pos_actual'])
        query_site_residues = list(query_info_in_alignment_pd['query_residues'])
        identity_perc_list = list(unique_pd['% Identity'])


        # makes chart
        html_chart_text(info_lines, width_canvas, query_site_actual, identity_perc_list, out_html_handle)

        # write tailing text
        out_html_handle.write(html_tail_handle.read())

        #write table data
        html_table_text(query_site_residues, query_site_actual, identity_perc_list, out_html_handle)

        # write color legend data
        html_color_legend(match_color, similar_color, mismatch_color, out_html_handle)

        # write alignment text
        html_alignment_text(alignment_highlighted_pd, alignment_pd, id_maxLength, colored_id_pd, query_info_in_alignment_pd, out_html_handle)

        # write closing lines to html
        out_html_handle.write('</PRE>\n</BODY>\n</HTML>')



def process_for_html(ref_id, bool_query_region_pd, alignment_pd, query_pos_aligned, query_info_in_alignment_pd, unique_pd, mismsatch_summary_info, out_dir):
    # alignment_len = len(alignment_pd)
    bool_query_region_pd = bool_query_region_pd.T
    alignment_highlighted_pd = alignment_pd.T.copy()    # deep copy
    id_maxLength = alignment_highlighted_pd.index.str.len().max()

    match_color = '#42DB33'
    similar_color = '#FFFF00'
    mismatch_color = '#F73D94'
    ref_color = '#00FFFF'

    match_prefix = '<strong><span style="background-color: %s">' % match_color
    suffix_string = '</span></strong>'
    mismatch_prefix = '<strong><span style="background-color: %s">' % mismatch_color
    ref_prefix = '<strong><span style="background-color: %s">' % ref_color

    for pos in query_pos_aligned:
        prefix = bool_query_region_pd[pos].apply(lambda x:match_prefix if x else mismatch_prefix)
        suffix = bool_query_region_pd[pos].apply(lambda x:suffix_string)
        alignment_highlighted_pd[pos] = prefix + alignment_highlighted_pd[pos] + suffix


    # To color IDs in html output
    colored_id_pd = pd.DataFrame(index=bool_query_region_pd.index)
    colored_id_pd['match_status'] = bool_query_region_pd.sum(axis=1) == len(unique_pd)
    colored_id_pd['ID'] = bool_query_region_pd.index

    colored_id_pd.loc[~(colored_id_pd['match_status']), 'ID'] = colored_id_pd['ID']\
        .apply(lambda x:'%s%s%s%s' % (mismatch_prefix.replace('<strong>', ''), x[:2], suffix_string.replace('</strong>', ''), x[2:].ljust(id_maxLength-2)))
    colored_id_pd.loc[(colored_id_pd['match_status']), 'ID'] = colored_id_pd['ID']\
        .apply(lambda x: '%s%s%s%s' % (match_prefix.replace('<strong>', ''), x[:2], suffix_string.replace('</strong>', ''), x[2:].ljust(id_maxLength-2)))
    colored_id_pd.loc[ref_id, 'ID'] = '%s%s%s' % (ref_prefix, ref_id.ljust(id_maxLength), suffix_string)


    html_text_out(alignment_highlighted_pd, alignment_pd, id_maxLength, query_info_in_alignment_pd, unique_pd, colored_id_pd, match_color, similar_color, mismatch_color, mismsatch_summary_info, out_dir)



def get_query_pos_list(query_pos_list, reference_pd):
    if query_pos_list == 'all':
        query_pos_actual = range(0, len(reference_pd))
    else:
        query_pos_actual = query_pos_list
        if max(query_pos_actual) > len(reference_pd):
            print 'Error. The highest of query residue positions provided exceeds length of reference sequence (%i).' % (len(reference_pd))
            exit()

    return query_pos_actual


def main_script(reference_f, query_seq_f, alignment_f, query_pos_list, out_dir):
    reference_pd = fasta_to_pandas(reference_f)
    alignment_pd = fasta_to_pandas(alignment_f)

    # get query residue positions
    query_pos_actual = get_query_pos_list(query_pos_list, reference_pd)

    # check reference sequence criteria pass prelimiary test
    ref_id, ref_seq_in_alignment_pd = check_ref_seq_in_alignment_pd(reference_pd, alignment_pd)

    # process query pos and residue info in reference seq and in alignment
    query_ref_pd, query_info_in_alignment_pd = get_query_residue_info(query_pos_actual, reference_pd, ref_id, ref_seq_in_alignment_pd)

    # get query positions in alignment
    query_pos_aligned = list(query_info_in_alignment_pd['pos_in_alignment'])

    # get query region of alignment
    query_region_pd, bool_query_region_pd, inverse_bool_query_region_pd = booleante(alignment_pd, query_pos_aligned, ref_id)

    # get length of sequences in alignment
    seq_length_alignment_s = alignment_pd[alignment_pd != '-'].count()

    # summarize residues at query positions
    unique_pd = summarize_unique_residues(query_info_in_alignment_pd, query_region_pd, ref_id)

    # process data for csv output writing
    results_csv_pd = prepare_for_csv_output(query_region_pd, bool_query_region_pd, inverse_bool_query_region_pd, seq_length_alignment_s, query_info_in_alignment_pd)

    # write csv output
    mismsatch_summary_info = write_csv_out(results_csv_pd, unique_pd, bool_query_region_pd, reference_f, query_seq_f, alignment_f, out_dir)
    print 'Done writing CSV output'

    # write html output
    process_for_html(ref_id, bool_query_region_pd, alignment_pd, query_pos_aligned, query_info_in_alignment_pd, unique_pd, mismsatch_summary_info, out_dir)
    print 'Done writing HTML output'


if __name__ == '__main__':
    query_pos_list = [1,58]
    reference_f = '../data/Reference.fasta'
    query_seq_f = None
    alignment_f = '../data/aligned.fasta'
    out_dir = '../data/output/'

    main_script(reference_f, query_seq_f, alignment_f, query_pos_list, out_dir)