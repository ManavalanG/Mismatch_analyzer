# Mismatch Analyzer
Analyzes residue conservation in multiple sequence
alignment.

### How to use?
1. Clone/[download](https://github.com/ManavalanG/Mismatch_analyzer/archive/master.zip) the repository
2. Run tool '[src/mismatch analyzer](src/mismatch_analyzer.py)' with appropriate arguments. For example:
```bash
cd src
python mismatch_analyzer.py -r ../data/reference.fasta  -qs ../data/query.fasta -qp "5,10,15"
```
3. Required arguments:
    * A reference sequence file
    * A query sequences file OR a pre-aligned MSA file
    * Query residue positions
    
4. To see arguments available, refer to help
```bash'
python mismatch_analyzer.py -h
```

### Requirements
* Python 2.7
* Required Python modules are listed in [requirements.txt](requirements.txt)
* [Clustal Omega](http://www.clustal.org/omega/) if performing MSA is desired.  


### Example output
See following files for output produced by this tool.
1. [CSV output](data/output/csv_out.tsv)
2. [HTML output](data/output/html_out.html)

### License
Data and code in this repository is licensed under [MIT license](LICENSE.md). 

