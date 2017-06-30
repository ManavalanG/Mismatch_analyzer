from Bio import AlignIO
import pandas as pd
import pandascharm as pc
import numpy as np

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 100000)

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

# def xxxx():
if __name__ == '__main__':
    ref_id = 'id_1|a'
    query_pos_aligned = [0,1,2]
    aligned_data = AlignIO.read('test_data/aligned.fasta', "fasta")
    alignment_pd = pc.from_bioalignment(aligned_data)       # converted to pandas dataframe

    # changes datatype of all columns from category to object
    for col_name in alignment_pd.columns:
        alignment_pd[col_name] = alignment_pd[col_name].astype(object)

    query_pd, bool_query_pd, inverse_bool_query_pd = booleante(alignment_pd, query_pos_aligned, ref_id)

    seq_length = alignment_pd[alignment_pd != '-'].count()
    results_csv_pd = csv_data(query_pd, bool_query_pd, inverse_bool_query_pd, seq_length)



def xxwef():
    seq_length = alignment_pd[alignment_pd!='-'].count()


    query_df = alignment_pd.loc[query_pos_aligned, :]
    print query_df


    bool_query_df = query_df.apply(lambda x:x==query_df[ref_id])
    bool_query_inverse_df = query_df.apply(lambda x:x!=query_df[ref_id])
    print bool_query_df
    # bool_query_df = bool_query_df.astype(int)

    mismatch_count = bool_query_df[bool_query_df<1].count()

    # print query_df.values * bool_query_df.values
    # print pd.DataFrame(query_df.values*bool_query_df, columns=query_df.columns, index=query_df.index)
    mult = query_df * bool_query_inverse_df

    # for col_name in alignment_pd.columns:
    #     mult.ix[mult[col_name] == '', col_name] = '='
    # out_pd = mult.T
    # out_pd['mismatch_count'] = mismatch_count
    # out_pd['seq_length'] = seq_length

    # out_pd['title'] = out_pd.index
    title_df = out_pd['title'].str.split('|', expand=True)
    # print pd.Series(title_df.columns, dtype=str)

    title_df.columns = pd.Series(['a']*2).map(str) + '-' + pd.Series(title_df.columns).map(str)
    print title_df.columns
    # out_pd.drop(['title'])
    out_pd.drop('title', axis=1, inplace=True)

    final = pd.concat([title_df, out_pd], axis=1)
    # print final.columns

