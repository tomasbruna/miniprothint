# Miniprothint

Miniprothint selects a set of reliable gene prediction hints from [miniprot](https://github.com/lh3/miniprot) alignments scored by [miniprot boundary scorer](https://github.com/tomasbruna/miniprot-boundary-scorer).

## Usage

We recommend running [miniprot](https://github.com/lh3/miniprot) and [miniprot boundary scorer](https://github.com/tomasbruna/miniprot-boundary-scorer) separately for easier debugging:

    miniprot genome.fasta proteins.fasta --aln > miniprot.aln
    miniprot_boundary_scorer -o miniprot_parsed.gff -s blosum62.csv < miniprot.aln
    miniprothint.py miniprot_parsed.gff --workdir miniprothint

**Coming soon**

However with a protein fasta file on input, miniprothint can take care of running [miniprot](https://github.com/lh3/miniprot) and [miniprot boundary scorer](https://github.com/tomasbruna/miniprot-boundary-scorer):

    miniprothint.py --genome genome.fasta --proteins proteins.fasta --workdir miniprothint
    
or just [miniprot boundary scorer](https://github.com/tomasbruna/miniprot-boundary-scorer):

    miniprot genome.fasta proteins.fasta --aln > miniprot.aln
    miniprothint.py --alignment miniprot.aln --workdir miniprothint

If the latter two ways of running miniprothint are used, [miniprot](https://github.com/lh3/miniprot) and [miniprot boundary scorer](https://github.com/tomasbruna/miniprot-boundary-scorer) must be in the `$PATH`.
    
## Outputs

* `miniprothint.gff` A set of hints passing a relaxed set of thresholds.
* `hc.gff`           A set of hints passing stringent thresholds.

If [miniprot](https://github.com/lh3/miniprot) and/or [miniprot boundary scorer](https://github.com/tomasbruna/miniprot-boundary-scorer) are run by miniprothint, their outputs are saved to:

* `miniprot.aln`
* `miniprot_parsed.gff`
