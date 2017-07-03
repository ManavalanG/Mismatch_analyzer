from shutil import copyfile

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
    out_html_handle.write('</HEAD>\n<BODY LANG="en-US" DIR="LTR">\n<PRE CLASS="western">\n%s'
                          '<div>\n  <canvas id="canvas_bar" height="350", width = "%i" ></canvas>\n'
                          '<div id="legend"></div>\n</div>\n\n<script>\n'
                          '\tvar barChartData = \n\t{\n'
                          '\t\tlabels : %s,\n\t\tdatasets : [\n'
                          '\t\t{	fillColor : "rgba(51,51,255,0.8)",\n\t\t\ttitle: "%% Identity",\n'
                          '\t\t\tdata : %s  },\n\n'
                          % (info_lines, width_canvas, query_site_actual, identity_perc_list))
    out_html_handle.write('\t\t]\n\t}\n')

    out_html_handle.write(html_tail_handle.read())


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


