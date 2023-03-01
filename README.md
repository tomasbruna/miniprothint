# miniprothint

### Usage

    miniprot_boundary_scorer -o miniprot.gff -s blosum62.csv < miniprot.aln
    miniprothint.py miniprot.gff --workdir miniprothint
    
### Outputs

* `miniprothint.gff` A set of hints passing a relaxed set of thresholds.
* `hc.gff`           A set of hints passing stringent thresholds.


### Notes

Note that `miniprot.gff` itself (the output of `miniprot_boundary_scorer`) is not identical to running miniprot and outputting to a gff file because a set of basic filters is already applied during parsing. This is done to keep the output file size small by removing the most obvious noise. This `miniprot_boundary_scorer` filtering can be turned off by passing the following parameters:

    miniprot_boundary_scorer -e -9999 -x -9999 -i -9999 -o miniprot.gff -s blosum62.csv < miniprot.aln
