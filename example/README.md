# Miniprothint Usage Example

To run ProtHint, use the following command

    ../bin/prothint.py inputs/miniprot_parsed.gff  --workdir miniprothint

The results in the `miniprothint` folder should match the contents of the
`reference` folder. Note that the order of hints in gff files might differ.

To check the results, use for example:

    diff <(sort reference/miniprothint.gff) <(sort miniprothint/miniprothint.gff)
    diff <(sort reference/hc.gff) <(sort reference/hc.gff)

The files `inputs/genome.fasta` and `inputs/proteins.fasta` can be used to test the [alternative modes](https://github.com/tomasbruna/miniprothint#usage) of miniprothint execution.
