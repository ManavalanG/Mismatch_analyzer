from Bio.Align.Applications import ClustalOmegaCommandline

seq_f = '/Users/mana/Documents/GitHub_downloaded/Mismatch_analyzer/test_data/clustal/Sequences.fasta'
aligned_f = '/Users/mana/Documents/GitHub_downloaded/Mismatch_analyzer/test_data/clustal/aligned_clustalo.fasta'

def run_clustal_omega(seq_f, aligned_f):
    clustal_commandline = {'infile': seq_f,
                           'outfile': aligned_f,
                           'outfmt': 'fasta',
                           'iterations': 3,
                           'force': True}

    try:
        clustalomega_cline = ClustalOmegaCommandline(**clustal_commandline)
    except Exception as e:  # This will catch all the major errors
        print 'Error in command directed to ClustalO. Check the command entered!'
        print 'Error produed:\n%s' % e
        raise

    try:
        clustalomega_cline()
    except Exception as e:  # This will catch all the major errors
        print 'Error when running Clustal Omega. Make sure Clustal Omega is installed and working properly'
        print 'Error produed:\n%s' % e
        raise

run_clustal_omega(seq_f, aligned_f)
