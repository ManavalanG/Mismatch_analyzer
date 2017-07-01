from Bio import AlignIO
import pandas as pd
import pandascharm as pc
import numpy as np

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 100000)
pd.set_option('max_colwidth', 100000)

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


def split_id_by_delimiter(results_csv_pd):
    title_df = pd.DataFrame(index=results_csv_pd.index)
    title_df['title'] = results_csv_pd.index
    title_df = title_df['title'].str.split('|', expand=True)  # split using a delimiter symbol
    title_df.columns = pd.Series(['Title'] * 2).map(str) + '_' + pd.Series(title_df.columns).map(str)

    return title_df


def csv_data(query_pd, bool_query_pd, inverse_bool_query_pd, seq_length):
    results_csv_pd = query_pd.T * inverse_bool_query_pd.T
    for col_name in results_csv_pd.columns:
        results_csv_pd.ix[(results_csv_pd[col_name]==''), col_name] = '='

    results_csv_pd['mismatch_count'] = bool_query_pd[bool_query_pd==0].count()
    results_csv_pd['seq_length'] = seq_length

    title_df = split_id_by_delimiter(results_csv_pd)

    results_csv_pd = pd.concat([title_df, results_csv_pd], axis=1)

    return results_csv_pd


def booleante(alignment_pd, query_pos_aligned, ref_id):
    query_pd = alignment_pd.loc[query_pos_aligned, :]

    bool_query_pd = query_pd.apply(lambda x:x==query_pd[ref_id])
    inverse_bool_query_pd = query_pd.apply(lambda x:x!=query_pd[ref_id])    # True and False are switched

    return query_pd, bool_query_pd, inverse_bool_query_pd


def write_html(alignment_pd, id_maxLength):
    alignment_len = len(alignment_pd.columns)
    block_size = 60
    if alignment_len % block_size:
        no_of_blocks = (alignment_len/block_size + 1)
    else:
        no_of_blocks = alignment_len/block_size

    with open('aa.txt', 'w') as handle:

        html_blocks_pd = pd.DataFrame(index=alignment_pd.index)
        for block_no in range(0, no_of_blocks):
            min = block_no * block_size
            max = (block_no +1) * block_size
            if alignment_len <= max:
                max = alignment_len
            # print range(min, max)
            # print

            html_blocks_pd[block_no] = alignment_pd[range(min, max)].apply(lambda x: ''.join(x), axis=1)
            # html_block = alignment_pd[range(min, max)].apply(lambda x: ''.join(x), axis=1)

        # print html_blocks_pd

                # print type(html_blocks_pd[[0]])
                # print html_blocks_pd.index.str.len().max()
            for i, item in html_blocks_pd[block_no].iteritems():
                handle.write('%s    <spaced>%s</spaced>\n' % (i.ljust(id_maxLength), item))

            handle.write('\n\n')
                # print i, '\t\t2222\t', row
            # xx = html_blocks_pd[[0]].to_string(index=False, justify='right', header=None)
            # handle.write(xx)


def process_for_html(bool_query_pd, alignment_pd, query_pos_aligned):
    # alignment_len = len(alignment_pd)

    print 'adsf'
    bool_query_pd = bool_query_pd.T
    alignment_pd = alignment_pd.T
    id_maxLength = alignment_pd.index.str.len().max()

    match_prefix = '<strong><span style="background-color: green">'
    match_suffix = '</span></strong>'
    mismatch_prefix = '<strong><span style="background-color: red">'

    for pos in query_pos_aligned:
        prefix = bool_query_pd[pos].apply(lambda x:match_prefix if x else mismatch_prefix)
        suffix = bool_query_pd[pos].apply(lambda x:match_suffix)

        alignment_pd[pos] = prefix + alignment_pd[pos] + suffix

    write_html(alignment_pd, id_maxLength)




# def xxxx():
if __name__ == '__main__':
    ref_id = 'id_1|a'
    query_pos_aligned = [1,2]
    aligned_data = AlignIO.read('test_data/aligned.fasta', "fasta")
    alignment_pd = pc.from_bioalignment(aligned_data)       # converted to pandas dataframe



    # changes datatype of all columns from category to object
    for col_name in alignment_pd.columns:
        alignment_pd[col_name] = alignment_pd[col_name].astype(object)

    query_pd, bool_query_pd, inverse_bool_query_pd = booleante(alignment_pd, query_pos_aligned, ref_id)

    seq_length = alignment_pd[alignment_pd != '-'].count()
    results_csv_pd = csv_data(query_pd, bool_query_pd, inverse_bool_query_pd, seq_length)

    process_for_html(bool_query_pd, alignment_pd, query_pos_aligned)



    # if True:    # for html output
    if False:
        bool_query_pd = bool_query_pd.T
        alignment_pd = alignment_pd.T
        id_maxLength = alignment_pd.index.str.len().max()

        match_prefix = '<strong><span style="background-color: green">'
        match_suffix = '</span></strong>'
        mismatch_prefix = '<spaced><strong><span style="background-color: red">'

        for pos in query_pos_aligned:
            prefix = bool_query_pd[pos].apply(lambda x:match_prefix if x else mismatch_prefix)
            suffix = bool_query_pd[pos].apply(lambda x:match_suffix)

            alignment_pd[pos] = prefix + alignment_pd[pos] + suffix

        # print alignment_pd[[0,1,2,3]]

        # alignment_pd['merged'] = alignment_pd[range(5)].apply(lambda x: ''.join(x), axis=1)
        # print alignment_pd['merged']

        # alignment_len = 61
        block_size = 30
        if alignment_len % block_size:
            no_of_blocks = (alignment_len/block_size + 1)
        else:
            no_of_blocks = alignment_len/block_size

        with open('aa.txt', 'w') as handle:

            html_blocks_pd = pd.DataFrame(index=alignment_pd.index)
            for block_no in range(0, no_of_blocks):
                min = block_no * block_size
                max = (block_no +1) * block_size
                if alignment_len <= max:
                    max = alignment_len
                # print range(min, max)
                # print

                html_blocks_pd[block_no] = alignment_pd[range(min, max)].apply(lambda x: ''.join(x), axis=1)
                # html_block = alignment_pd[range(min, max)].apply(lambda x: ''.join(x), axis=1)

            # print html_blocks_pd

                    # print type(html_blocks_pd[[0]])
                    # print html_blocks_pd.index.str.len().max()
                for i, item in html_blocks_pd[block_no].iteritems():
                    handle.write('%s    %s\n' % (i.ljust(id_maxLength), item))

                handle.write('\n\n')
                    # print i, '\t\t2222\t', row
                # xx = html_blocks_pd[[0]].to_string(index=False, justify='right', header=None)
                # handle.write(xx)