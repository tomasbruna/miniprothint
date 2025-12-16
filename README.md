# Miniprothint

Miniprothint selects a set of reliable gene prediction hints from [miniprot](https://github.com/lh3/miniprot) alignments scored by [miniprot boundary scorer](https://github.com/tomasbruna/miniprot-boundary-scorer).

## Usage

We recommend running [miniprot](https://github.com/lh3/miniprot) and [miniprot boundary scorer](https://github.com/tomasbruna/miniprot-boundary-scorer) separately for easier debugging:

    miniprot genome.fasta proteins.fasta --aln > miniprot.aln
    miniprot_boundary_scorer -o miniprot_parsed.gff -s blosum62.csv < miniprot.aln
    miniprothint.py miniprot_parsed.gff --workdir miniprothint
    
### Running with Apptainer/Singularity

The required tools are available in the [GALBA](https://github.com/Gaius-Augustus/GALBA) container:
    
    apptainer exec docker://katharinahoff/galba-notebook:latest miniprot genome.fasta proteins.fasta --aln > miniprot.aln
    apptainer exec docker://katharinahoff/galba-notebook:latest miniprot_boundary_scorer -o miniprot_parsed.gff -s /opt/miniprot-boundary-scorer/blosum62.csv < miniprot.aln
    apptainer exec docker://katharinahoff/galba-notebook:latest miniprothint.py miniprot_parsed.gff --workdir miniprothint

## Outputs

* `miniprothint.gff` A set of hints passing a relaxed set of thresholds.
* `hc.gff`           A set of hints passing stringent thresholds.

If [miniprot](https://github.com/lh3/miniprot) and/or [miniprot boundary scorer](https://github.com/tomasbruna/miniprot-boundary-scorer) are run by miniprothint, their outputs are saved to:

* `miniprot.aln`
* `miniprot_parsed.gff`
