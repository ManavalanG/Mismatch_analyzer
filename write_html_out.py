from shutil import copyfile


def html_chart_text(info_lines, width_canvas, query_site_actual, identity_perc_list, out_html_handle):
    text_for_chart = ('</HEAD>\n<BODY LANG="en-US" DIR="LTR">\n'
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


# def asdf():
if __name__ == '__main__':
    html_template_head_f = '/Users/mana/Documents/GitHub_downloaded/Mismatch_analyzer/test_data/Template_html5_ChartNewjs_head.txt'
    html_template_tail_f = '/Users/mana/Documents/GitHub_downloaded/Mismatch_analyzer/test_data/Template_html5_ChartNewjs_tail.txt'
    html_outfile = 'html_out.html'

    copyfile(html_template_head_f, html_outfile)

    with open(html_outfile, 'a') as out_html_handle, open(html_template_tail_f, 'r') as html_tail_handle:
        out_html_handle.write('sadfghj\sdfgh\dsfghj\n')

        info_lines ='***info_lines***'
        width_canvas = 325
        query_site_actual = [1,2,3]
        query_site_residues = ['A', 'D', 'R']
        identity_perc_list = [100, 33, 50]

        # makes chart
        html_chart_text(info_lines, width_canvas, query_site_actual, identity_perc_list, out_html_handle)

        # write tailing text
        out_html_handle.write(html_tail_handle.read())

        #write table data
        html_table_text(query_site_residues, query_site_actual, identity_perc_list, out_html_handle)

        # write color legend data
        match_color = '#42DB33'
        similar_color = '#FFFF00'
        mismatch_color = '#F73D94'
        html_color_legend(match_color, similar_color, mismatch_color, out_html_handle)
