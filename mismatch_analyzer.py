def fetch_mismatch():
	global unique_residues_line
	mismatch_log.debug("Module 'Parsing for Aligned residues' initiated.")
	global resi_positions
	global Aligned_Filename
	global Output_Path
	global query_in_Alignment
	global query_site_residues
	global query_site_actual
	global identity_perc_list
	# global similarity_perc_list
	global conserve_score_list
	global protein_mode
	global dict_identifier_status
	global match_seqs_total
	global mismatched_seqs_total
	global method_name

	mismatches_text_output = False		# Default: False. If True, writes a text output file along with csv output file.

	mismatch_log.debug("Alignment file '%s' will be read." % Aligned_Filename)
	# data_alignment = AlignIO.read(Aligned_Filename, "clustal")
	# clustal_aligned = False
	format_supported = False

	try:
		# this block identifies if alignment is in any of the supported formats
		format_type_list = ['clustal','emboss','fasta','fasta-m10','ig','maf','nexus','phylip','phylip-sequential','phylip-relaxed','stockholm']
		for format_type in format_type_list:
			try:
				data_alignment = AlignIO.read(Aligned_Filename, format_type)
				# if format_type == "clustal":
					# clustal_aligned = True
				mismatch_log.info('Alignment format is detected as "%s" format' % format_type)
				format_supported = True
				break
			except ValueError:
				pass
	except Exception as e:
		mismatch_log.error("Error that resulted: %s" % e)
		temp = ("Could not open file: \n\t'%s'" % Aligned_Filename)
		mismatch_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	if not format_supported:
		temp = "Your Alignment file does not have alignment format supported by ResCons"
		mismatch_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	mismatch_log.debug("Alignment file '%s' was read successfully." % Aligned_Filename)

	# Checks if Reference sequence ID is present in alignment file.
	Reference_index = None
	for serial_no in range(len(data_alignment)):
		if Reference.id == data_alignment[serial_no].id:
			Reference_index = serial_no

	# Checks if Reference seq is present in MSA and if present, checks if they match or not
	if Reference_index is None:
		temp = 'Error: Alignment file provided does not have Sequence ID corresponding to your reference sequence. ' \
				'Fix it and try again!'
		mismatch_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		refseq_in_alignment = data_alignment[Reference_index].seq
		refseq_in_alignment = str(refseq_in_alignment).replace('-', '')
		if str(Reference.seq) == refseq_in_alignment:
			clustal_log.info("Reference sequence matches to corresponding sequence in Alignment file provided.")
		else:
			temp = "Reference sequence does not match to corresponding sequence in Alignment file provided. " \
					"Fix it and try again!"
			mismatch_log.error(temp)
			popup_error(temp)
			raise_enabler('stop')

	mismatch_log.info("Reference sequence's ID is present in alignment file and their sequences match.")

	# Indexes residues (characters) in reference sequence and also for reference sequence in alignment.
	mismatch_log.debug('Begins indexing residues in reference seq and also of reference sequence in alignment')
	reference_seq_aligned = str(data_alignment[Reference_index].seq)

	# following line creates list of upper case alphabets.
	# Includes non-amino acid chars as well to support nucleotides other than A, T, G and C
	# works for both upper and lower case characters and considers non-amino acid and non-nucleotide residues as mismatch
	aa_list = map(chr, range(65, 91)) + map(chr, range(97,123))
	Dict_aa_Index_Reference_Seq = {k: [] for k in aa_list}
	Dict_aa_Index_Reference_Seq_inAlignment = {k: [] for k in aa_list}
	for aa in aa_list:
		Dict_aa_Index_Reference_Seq_inAlignment[aa] = indexing(reference_seq_aligned, aa)
		Dict_aa_Index_Reference_Seq[aa] = indexing(Reference.seq, aa)
	mismatch_log.debug('Done indexing residues in Reference seq and also of Reference sequence in alignment')


	# This section is to get appropriate number of 'name components' in the title row.
	# Also this accounts even if the clustal record ids have excessive number of pipes without any info in them.
	aligned_ids_list = []
	pipes_total_list = []
	justified_list = []
	for record in data_alignment:
		aligned_ids_list.append(record.id)
		pipes = (record.id).split(id_delimiter)

		if len(pipes) == 1:
			justified_list.append(len(pipes[0]))

		pipes_temp = pipes
		for no in range(1, len(pipes)+1):
			if pipes[-no] == '':
				pipes_temp = pipes[:(-no)]
			else:
				break
		pipes_total_list.append(len(pipes_temp))
	pipes_most = max(pipes_total_list)
	mismatch_log.debug('Done determining number of title components in the sequence record IDs.')

	# If sequence IDs are made only of few components, this will enable reasonable alignment in the text file output
	if pipes_most == 1:
		justified = max(justified_list)
	elif pipes_most == 2:
		justified = 45
	elif pipes_most ==3:
		justified = 33
	else:
		justified = 23


	# Writes the title line for csv file and txt file
	Title_components_csv = 'S.No '
	Title_components_txt = '\t\t' + 'S.No'.ljust(7) + '\t'
	for no in range(0, pipes_most):
		Title_components_csv += (',' + 'Title_' + str(no+1))
		Title_components_txt += (('Title_' + str(no + 1)).ljust(justified) + '\t')
	Title_components_csv += ',Sequence Length'
	Title_components_csv += ',No. of Mismatches'
	Title_components_txt += ('Mismatch'.ljust(10)+ '\t' + 'Expected'.ljust(10) + '\n')
	mismatch_log.debug('Determined the title to be used in csv and tect files.')


	if checkboxval_Clustal.get():	# To account for input file used depending on clustal alignment was selected or not
		input_file = SeqQuery_file
	else:
		input_file = Aligned_Filename


	# Extracts residues present at query sites in the reference sequence and also reference sequence in alignment
	mismatch_log.debug('Begins the part to extract residues at the requested query sites.')
	query_site_actual = resi_positions		# query_site_actual refers to actual positioning whereas query_site to python indexes.
	query_site_actual = sorted(query_site_actual)
	query_site = []
	query_site_residues = []
	for n in range(0, len(query_site_actual)):
		query_site.append(int(query_site_actual[n]) - 1)
		query_site_residues.append(Reference.seq[query_site[n]])

	mismatch_log.info('Query sites requested: %s' % (str(query_site_actual)))
	mismatch_log.info('Residues at those query sites: %s' % str(query_site_residues))


	# Finds the position of query site residues in reference's alignment sequence
	query_in_Alignment = []
	for aa in aa_list:
		for query_site_element in query_site:
			if query_site_element in Dict_aa_Index_Reference_Seq[aa]:
				pos = Dict_aa_Index_Reference_Seq[aa].index(query_site_element)
				query_in_Alignment.append(Dict_aa_Index_Reference_Seq_inAlignment[aa][pos])

	query_in_Alignment = sorted(query_in_Alignment)

	#To count for the fact that python indexes from 0 instead of 1
	query_in_Alignment_actual = [(item + 1) for item in query_in_Alignment]

	mismatch_log.info("Corresponding sites in the alignment's Reference sequence : %s" %str(query_in_Alignment_actual))

	# data_alignment = AlignIO.read(Aligned_Filename, "clustal")		# seek function doesn't work here

	# Extracts residues across alignment at the requested query sites
	mismatch_log.info('Begin extracting residues across alignment at the requested query sites')

	dict_clustal_residues = {k: [] for k in aligned_ids_list}
	dict_identifier_status = {}
	for record in data_alignment:
		dict_identifier_status[record.id] = []			# to assist in color coding IDs in html output
		for no in range(0, len(query_in_Alignment)):

			# this obtains similarity count and percentage
			aa_similar = ''
			for group in aa_set:		# aa_set is defined in settings text file and allowed to be changed via GUI
				if query_site_residues[no] in group:
					aa_similar = group
					break

			if record.seq[query_in_Alignment[no]] == query_site_residues[no]:
				dict_clustal_residues[record.id].append('')
				dict_identifier_status[record.id].append('match')
			elif record.seq[query_in_Alignment[no]] in aa_similar:
				dict_identifier_status[record.id].append('similar')
				dict_clustal_residues[record.id].append(record.seq[query_in_Alignment[no]])
			else:
				dict_clustal_residues[record.id].append(record.seq[query_in_Alignment[no]])
				dict_identifier_status[record.id].append('mismatch')

	for key in dict_identifier_status:		# to assist in color coding IDs in html output
		dict_identifier_status[key] = list(set(dict_identifier_status[key]))


	# Calculate "Residue Conservation score" by Liu08 method using Similarity matrix S obtained from BLOSUM62 matrix
	if protein_mode and conserve_method != 'amino_acid_grouping':
		# Since Liu08 similarity (conservation) score needs to be calculated from MSA w/o reference seq, a new MSA is
		# created that will not have reference sequence in it.
		if not ref_included:
			msa_data = []
			for item_no, item in enumerate(data_alignment._records):
				if item_no != Reference_index:
					msa_data.append(item)
			msa_data = MultipleSeqAlignment(msa_data)		# List of SeqRecords is now converted to be a MSA
		else:
			msa_data = data_alignment

		alignment_len = len(msa_data[0].seq)

		# Choose which Liu08 scoring method (Sequence weighted or not)
		if conserve_method  == 'liu08_non_seq_weighted':
			liu08_simple_score_list = []
			positions_concerned = query_in_Alignment
			method_name = "Liu08 Non-Sequence-weighted"		# for csv and html output
		elif conserve_method == 'liu08_seq_weighted':
			liu08_weighted_score_list = []
			positions_concerned = range(alignment_len)
			method_name = "Liu08 Sequence-weighted"			# for csv and html output


		# get frequency and unique amino acid count in a column
		# To calculate sequence weight, all columns need to be processed here
		# if conserve_method in ['liu08_non_seq_weighted', 'liu08_seq_weighted']:
		dict_pos_freq = {}
		dict_pos_unique_aa = {}
		for column_no in positions_concerned:
			column_aa = msa_data[:,column_no]
			column_aa = column_aa.upper()		# to help if seqs are in lower case or mix of lower and upper cases

			unique_aa = list( set( column_aa ) )
			dict_pos_unique_aa[column_no] = [ len(unique_aa) ]
			dict_pos_unique_aa[column_no].append(unique_aa)

			dict_pos_freq[column_no] = {}
			for aa in aa_label:
				column_aa = column_aa.upper()		# to help if seqs are in lower case or mix of lower and upper cases
				dict_pos_freq[column_no][aa] = column_aa.count(aa)


		# calculate seq weight for each seq when Liu08 seq weighted score needs to be determined
		if conserve_method == 'liu08_seq_weighted':
			dict_seq_weight = {}
			for seq_no in range(0, len(msa_data)):
				sequence = msa_data[seq_no].seq
				sequence = sequence.upper()				# to help if seqs are in lower case or mix of lower and upper cases
				seq_weight = 0
				for aa_no in range(0, len(sequence)):
					freq = float( dict_pos_freq[aa_no][ sequence[aa_no] ] )
					if freq:
						seq_weight += 1 / (freq * dict_pos_unique_aa[aa_no][0])
				dict_seq_weight[seq_no] = seq_weight / alignment_len        # Normalized sequence weight

		# determine residue conservation score by Liu08 method
		# for column_no in range(0, len(query_in_Alignment)):
		for column_no in query_in_Alignment:
			column_aa = msa_data[:,column_no]
			column_aa = str(column_aa).upper()		# to help if seqs are in lower case or mix of lower and upper cases

			# Identifies most frequent residue. If it is gap, next most common residue is chosen.
			# At the moment, doesn't account for the issue where more than one residue are sharing second most common status.
			sorted_aa_count= sorted(dict_pos_freq[column_no].iteritems(), key = lambda (k,v): (v,k), reverse = True)
			aa_most = sorted_aa_count[0][0]
			if sorted_aa_count[0][0] == '-' and sorted_aa_count[1][0] != 0:
				aa_most = sorted_aa_count[1][0]

			# Calculates Liu08 Sequence Weighted residue conservation score
			# While Liu08 score ranges from 0 to 10, here we present it as percentage from 0 to 100
			if conserve_method == 'liu08_seq_weighted':
				liu08_weighted_score = 0
				for row_no in range(0, len(column_aa)):
					matrix_score = float( dict_similarity_matrix[aa_most][column_aa[row_no]] )
					liu08_weighted_score += (dict_seq_weight[row_no] * matrix_score)
				temp = liu08_weighted_score * 10
				# this is to avoid Python rounding up, for eg., 99.999 to 100.0. This is a bit of hard coding but this
				# list will not be used anywhere else for calculation purposes. It'll be just used for writing in output
				if temp > 99.94 and temp < 100:
					temp = 99.9
				liu08_weighted_score_list.append(temp)	# multiplied by 10 to present in percentage
				# print '%i Liu08  %s  %0.2f' %(column_no,aa_most, liu08_weighted_score)

			# Calculates simple Liu08 NON-Sequence Weighted residue conservation score
			# Score is presented in percentage
			if conserve_method == 'liu08_non_seq_weighted':
				liu08_simple_score = 0
				for aa in aa_label:
					matrix_score = float(dict_similarity_matrix[aa_most][aa])
					liu08_simple_score += ( matrix_score * dict_pos_freq[column_no][aa] )
				temp = liu08_simple_score/len(column_aa) * 10
				# this is to avoid Python rounding up, for eg., 99.999 to 100.0. This is a bit of hard coding but this
				# list will not be used anywhere else for calculation purposes. It'll be just used for writing in output
				if temp > 99.94 and temp < 100:
					temp = 99.9
				liu08_simple_score_list.append(temp)


		# # calculates shannon entropy - Turned off for now as it needs to be tested further
		# # Non-standard amino acids are treated as 'gaps'; Similar to Scorecons server
		# If this need to be used import this ->  from math import log as math_log
		# score_shannon = 0
		# gaps_rel_freq = 0
		# for symbol in aa_label:
		# 	if symbol not in ['-', 'B', 'J', 'O', 'U', 'X', 'Z']:	# non-std amino acids
		# 		rel_freq_aa = float( dict_pos_freq[column_no][symbol] ) / len(data_alignment)
		# 		if rel_freq_aa > 0:
		# 			score_shannon -= ( rel_freq_aa * math_log(rel_freq_aa, len(data_alignment)) )
		# 	else:						# non-std amino acids treated as gaps
		# 		gaps_rel_freq += float( dict_pos_freq[column_no][symbol] ) / len(data_alignment)
		#
		# if gaps_rel_freq > 0:
		# 	score_shannon -= ( gaps_rel_freq * math_log(gaps_rel_freq, len(data_alignment)) )
		#
		# score_shannon = 1 - score_shannon
		# # print "%0.3f" %score_shannon, dict_pos_unique_aa[column_no][1]

	if conserve_method == 'amino_acid_grouping':
		method_name = "Amino acid Grouping"             # for csv and html output


	# Creates and writes  mismatch details in to a csv file
	mismatch_log.debug('Creates and writes  mismatch details in to a csv file')

	Output_filename_csv = 'Mismatches_Tabulated.csv'
	Output_file_csv = Output_Path + Output_filename_csv

	mismatch_log.info('Mismatch details will be written in to file:\n  "\t%s"' %Output_file_csv )
	try:
		Output_handle_mismatch_csv = open(Output_file_csv, 'w')
	except IOError as exception:
		mismatch_log.error('Error that resulted when opening csv file for writing: %s' % exception)
		if exception.errno != errno.EEXIST:
			temp = ("A file titled \n\t'%s' \nseems to be open in MS Excel. Close that file and try again!" %Output_filename_csv)
			mismatch_log.error(temp)
			popup_error(temp)
			raise_enabler('stop')

	for no in range(0, len(query_site_actual)):
		Title_components_csv += (',' + str(query_site_actual[no]) + ' "' + query_site_residues[no] + '"')
	Title_components_csv += '\n'
	mismatch_log.debug('Created title line for using in csv file')

	# Gethers residues at query sites and decides if they are mismatch or not.
	mismatch_log.debug('Starting to extract and write mismatch details into csv file')

	# sno_1 = 0
	# sno_2 = 0
	Mismatches_details = ''
	No_Mismatches_details = ''
	mismatched_seqs_total = 0
	match_seqs_total = 0
	mismatches_none = True		# This var is to get info if there is no mismatches noted at any of the sites tested
	for identifier in aligned_ids_list:
		mismatch_count = 0					# to count number of mismatches
		for char in dict_clustal_residues[identifier]:
			if char != '':
				mismatch_count += 1

		pipe_split = identifier.split(id_delimiter)
		if pipe_split[-1] == '':			# Accounts for if pipe is in the end of the id name
			pipe_split = pipe_split[:-1]

		if len(pipe_split) < pipes_most:	# makes sure all IDs have same number of naming components
			for no in range(0, pipes_most - len(pipe_split)):
				pipe_split.append('-na-')

		id_components = ''					# adds name components
		for no in range(0, pipes_most):
			if ',' in pipe_split[no] or '"' in pipe_split[no]:	# standardizes as per csv file format since spreadsheet software will misbehave when a cell value has comma or double quote in it
				pipe_split[no] = pipe_split[no].replace('"','""')
				pipe_split[no] = ('"%s"' % pipe_split[no])

			id_components += (',' + pipe_split[no])

		residues_data = ''					# obtains residue in that position, mismatched or not.
		for no in range(0, len(query_in_Alignment)):
			res = dict_clustal_residues[identifier][no]
			if dict_clustal_residues[identifier][no] == '':		# If residue matches, equal sign is used to improved visibility
				res = '='
			if dict_clustal_residues[identifier][no] == '-':	# Star symbol is used instead of gap symbols to improved visibility
				res = '*'
			residues_data += (',' + res)
		residues_data += '\n'

		len_seq = None			# this part obtains the length of sequence
		for no in range(0, len(data_alignment)):
				if data_alignment._records[no].id == identifier:
					seq = str(data_alignment._records[no].seq)
					seq = seq.replace('-', '')
					len_seq = len(seq)
					break

		if len_seq is None:
			mismatch_log.error("Error obtaining sequence length of '%s'." % identifier)

		if mismatch_count > 0:				# Writes into output file if at least one mismatch found in  query positions
			# serial_number_1= str (sno_1 + 1)
			# Output_handle_mismatch_csv.write(str (serial_number_1) +id_components + ',' + str(len_seq) + ',' +
			# 						str (mismatch_count) + residues_data)
			Mismatches_details += (id_components + ',' + str(len_seq) + ',' + str (mismatch_count) + residues_data)
			mismatched_seqs_total += 1
			mismatches_none = False

		else:								# Gathers details of records without any mismatches in query positions
			# serial_number_2= str (sno_2 + 1)
			No_Mismatches_details += (id_components  + ',' + str(len_seq)  + ',' + str(mismatch_count) + residues_data)
			match_seqs_total += 1

	Output_handle_mismatch_csv.write('Sequences file used:    "%s"\n' % input_file)
	Output_handle_mismatch_csv.write('Reference file used:    "%s"\n' % Reference_file)
	Output_handle_mismatch_csv.write('Residue conservation method used:    %s\n\n' %method_name)
	# Output_handle_mismatch_csv.write('Reference sequence included in calculations:    %s\n\n' %ref_included)

	if ref_included:		# reference sequence included in calculations
		total_no_seqs = len(data_alignment)
		ref_str = 'including Reference sequence'
	else:
		total_no_seqs = len(data_alignment)-1
		match_seqs_total -= 1
		ref_str = 'excluding Reference sequence'

	perc_match = float(match_seqs_total) * 100 / total_no_seqs
	perc_mismatch = float(mismatched_seqs_total) * 100 / total_no_seqs

	Output_handle_mismatch_csv.write('Number of sequences\n' +
									 '       in alignment (%s):      %i\n' %(ref_str, total_no_seqs)+
									 '       that have all residues matching (%s):   %i (%2.1f %%)\n'
									 %( ref_str, match_seqs_total, perc_match ) +
									 '       that have at least one mismatching residue:      %i (%2.1f %%)\n\n\n'
									 %( mismatched_seqs_total, perc_mismatch ) )

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
			Output_handle_mismatch_csv.write( str(sno) + ','.join(line) + '\n')
			sno += 1

	# sorts the data that have no mismatches at all
	csv_data = csv.reader(StringIO.StringIO(No_Mismatches_details))
	low = 2						# will have to be changed in the end to 1
	high = len(id_components.split(',')) + 1
	if natural_sorting:
		sortedlist2 = natsorted(csv_data, key=operator.itemgetter(*range(low, high)))
	else:
		sortedlist2 = sorted(csv_data, key=operator.itemgetter(*range(low, high)))


	# to find unique residue in query sites
	mismatch_log.debug('Unique residues at each requested site will be extracted')

	Dict_Unique_Residues = {k: [] for k in query_site_actual}
	for no in range(0, len(query_site_actual)):
		for res_list in dict_clustal_residues.values():
			Dict_Unique_Residues[query_site_actual[no]].append(res_list[no])

		# this part will enable calculating each unique residue's count and fraction w/o including reference seq's residue.
		# This is so as to offer true calculation as reference seq' residue is always matching anyway.
		if '' in Dict_Unique_Residues[query_site_actual[no]] and not ref_included:
			Dict_Unique_Residues[query_site_actual[no]].remove('')

	# if in case user sets 'similar aa sets' empty in settings (either through gui or settings file), the script will
	# treat rest of analysis as if it is in dna mode (where similarity is not considered)
	if aa_set_str == '':
		protein_mode = False

	Output_handle_mismatch_csv.write("\n\n*** Unique residues seen at the query sites and their count. ***\n")
	if ref_included:		# tells whether reference sequence is involved or not in calculations
		Output_handle_mismatch_csv.write("** (Note: Calculation includes Reference sequences's residue at that position) **\n\n")
	else:
		Output_handle_mismatch_csv.write("** (Note: Calculation DOESN'T include Reference sequences's residue at that position) **\n\n")

	if protein_mode:
		Output_handle_mismatch_csv.write( (",Expected Residue, ,Identity_count,%% Identity,%% Conservation (%s),Unique residues' count and fraction\n") % method_name)
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
		conserve_score_list = similarity_perc_list      # both lists will change when one changes
	elif conserve_method == 'liu08_non_seq_weighted':
		conserve_score_list = liu08_simple_score_list
	elif conserve_method == 'liu08_seq_weighted':
		conserve_score_list = liu08_weighted_score_list

	# get stats about %Identity (and optionally %similarity) at column positions requested
	for no in range(0, len(query_site_actual)):
		for item in range(0, len(Dict_Unique_Residues[query_site_actual[no]])):
			if Dict_Unique_Residues[query_site_actual[no]][item] == '':
				Dict_Unique_Residues[query_site_actual[no]][item] = query_site_residues[no]

		unique_res_detail = ','
		unique_res_detail += (str(query_site_actual[no]) + ' "' + query_site_residues[no] + '",')

		# this writes unique residues for each site in descending order of their count
		Dict_sorting_temp = {}
		for elem in sorted(set(Dict_Unique_Residues[query_site_actual[no]])):
			Dict_sorting_temp[elem] = Dict_Unique_Residues[query_site_actual[no]].count(elem)

		# this obtains identity count and percentage
		if query_site_residues[no] in Dict_sorting_temp:
			identity_count = Dict_sorting_temp[query_site_residues[no]]
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
		for group in aa_set:		# Note: aa_set is defined in settings text file and allowed to be changed via GUI
			if query_site_residues[no] in group:
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
		for key, value in sorted(Dict_sorting_temp.iteritems(), key = lambda (k,v): (v,k), reverse = True):
			perc = float (value) / (len(aligned_ids_list) - 1) * 100
			perc = '%3.1f' %perc
			unique_each += (',' + key + ': ' + str(value) + ' (' + str(perc) + '%)')

		# Collects data, which needs to be written into a output file, in a variable
		if protein_mode:			# residue conservation details added only if protein mode is used
			unique_res_detail += (',%i,%3.1f,%0.1f' % (identity_count, identity_perc, conserve_score_list[no]))
		else:
			unique_res_detail += (',%i,%3.1f' % (identity_count, identity_perc))
		unique_res_detail += unique_each
		unique_residues_line += ( unique_res_detail + '\n')
		unique_residues_list.append(unique_each)		# for text file output

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

	# This part gathers mismatch details to write into a text file.
	if mismatches_text_output:
		mismatch_log.info('Begin process to write mismatch details in to a "text file". People find this useful ah?')

		Output_file_Mismatch_txt = Output_Path + 'Mismatches_Detailed.txt'
		Output_handle_mismatch_txt = open(Output_file_Mismatch_txt, 'w')
		Output_handle_mismatch_txt.write('Sequences file used'.ljust(38) + ':  "' + input_file + '"\n' +
					'Reference file used'.ljust(38) + ':  "' + Reference_file + '"\n\n' +
					'Query site(s) for which mismatching are requested'.ljust(55) + ':\t' + str(query_site_actual) + '\n' +
					'Corresponding Query site residue(s) in Reference seq'.ljust(55) + ':\t' +  str(query_site_residues) + '\n' +
					'Their corresponding position(s) in alignment:'.ljust(55) + ':\t' +  str(query_in_Alignment_actual) + '\n\n' +
					"Number of sequences in the provided alignment:".ljust(55) + ':\t' +  str(len(data_alignment)) + '\n\n\n')

		mismatch_log.info('Mismatch details will be written in to text file: \n\t%s' %Output_file_Mismatch_txt)

		# If residues are not matching, this fetches the details about corresponding sequence records
		# for writing in to text file
		mismatch_log.debug('Begins extracting mismatch details to be written in to text file')

		siteno = 0
		for num in range(0, len(query_in_Alignment)):
			aligned_residues = data_alignment[:, query_in_Alignment[num]]		# Obtains residues across alignment in requested position

			element_index = 0
			Mismatch = ''
			increment = 1
			for element in aligned_residues:
				if element != query_site_residues[siteno]:		# Obtains details if residue is a mismatch
					ID_target = data_alignment[element_index].id
					ID_target = str(ID_target)

					pipe_split = ID_target.split(id_delimiter)
					if pipe_split[-1] == '':			# Accounts for if pipe is in the end of the id name
						pipe_split = pipe_split[:-1]

					if len(pipe_split) < pipes_most:	# makes sure all IDs have same number of naming components
						for no in range(0, pipes_most - len(pipe_split)):
							pipe_split.append('-na-')

					id_components = ''					# adds name components
					for no in range(0, pipes_most):
						id_components += (pipe_split[no].ljust(justified) + '\t')

					Details = ('\t\t' + str(increment).ljust(5) + '\t' +
							   id_components + element.ljust(10) + '\t' + query_site_residues[siteno].ljust(10) + '\n')

					Mismatch += Details
					increment += 1

				element_index += 1

			unique_res = unique_residues_list[num].split(',')
			unique_res = unique_res[1:]

			Output_handle_mismatch_txt.write('Site#	' + str(siteno + 1) + '/' + str(len(query_site)) + '\n')
			Output_handle_mismatch_txt.write('Expected Residue'.ljust(45) + ':\t"' +  query_site_residues[siteno]+ '"\n' +
						'Position of this residue in reference seq'.ljust(45) + ':\t' + str(query_site_actual[siteno]) + '\n' +
						'Position of this residue in alignment'.ljust(45) + ':\t' +  str(query_in_Alignment[num] + 1) + '\n' +
						'Residues across alignment at this site'.ljust(45) + ':\t' +  aligned_residues + '\n' +
						'Unique residues present and their count'.ljust(45) + ':\t' +  str(unique_res) + '\n' +
						'% Identity'.ljust(45) + ':\t' +  str(identity_perc_list[num]) +
											 '%%\t(%i residues are identical)' % identity_count_list[num] + '\n')

			if protein_mode:		# residue conservation details only if proteins sequences are used
				Output_handle_mismatch_txt.write('% Conservation_Score'.ljust(45) + ':\t' +
													str(conserve_score_list[num]) + '%%\n')

			if len(Mismatch) != 0:					# This is to write into output file if mismatches were found
				Output_handle_mismatch_txt.write('\n*** Following are the mismatches noted ***\n\n')
				Output_handle_mismatch_txt.write(Title_components_txt + Mismatch + '\n')
			Output_handle_mismatch_txt.write('\n')
			siteno += 1

			if len(Mismatch) == 0:					# This is if there are no mismatches for a residue
				Output_handle_mismatch_txt.write('*** No mismatches at this site ***\n\n\n')

		Output_handle_mismatch_txt.close()
		mismatch_log.info('Done extracting and writing mismatch details in to text file mentioned above')
