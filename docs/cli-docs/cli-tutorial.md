# CLI Tutorial

## Overview

The `tsplit` command is used to split transposable elements into their internal and external segments. It can be used for both long terminal repeats (LTRs) and terminal inverted repeats (TIRs).

The command line tools use sequence alignment to identify the boundaries of the segments based on user-defined parameters.

## Example usage

tSplit can be run in two modes: `tsplit LTR` and `tsplit TIR`, for extracting long terminal repeats or terminal inverted repeats, respectively.

Options are the same for each.

## tsplit TIR

### Basic TIR processing

For each element in _TIR_element.fa_ split into internal and external (TIR) segments.

Split segments will be written to _TIR_split_tsplit_output.fasta_ with suffix "\_I" for internal or "\_TIR" for external segments.

If the `--both` flag is set, then both left and right termini will be reported with suffixes "\_L\_TIR" and "\_R\_TIR".

In this example, TIRs must be at least 10bp in length and share 80%
identity and occur within 10bp of each end of the input element.

```bash
tsplit TIR -i tests/data/TIR_element.fa -d results -p TIR_split

# Equivalet to defaults
tsplit TIR -i tests/data/TIR_element.fa -d results -p TIR_split --method blastn --maxdist 10 --minid 80.0 --minterm 10 --blast_evalue 0.001 --method blastn --splitmode split
```

Output: `results/TIR_split_tsplit_output.fasta`

### Report both TIRs

With `--splitmode external` only the outermost TIR is returned. Setting `--both` returns both TIRs - useful if not identical.

```bash
tsplit TIR -i tests/data/TIR_element.fa -d results -p TIR_external_both --splitmode external --both
```

Output: `results/TIR_external_both_tsplit_output.fasta`

### Output PAF alignments and GFF TIR annotations

```bash
tsplit TIR -i tests/data/TIR_element.fa -d results --splitmode split --blast_evalue 0.001 --minid 60.0 --paf --gff
```

Output:

- `results/TIR_element_tsplit_output.fasta`
- `results/TIR_element.gff3`
- `results/TIR_element.paf`

### Generate TIR annotated dotplot

```bash
flexidot -i tests/data/TIR_element.fa -a results/TIR_element.paf -m 2 -o results/blast_dotplot --gff results/TIR_element.gff3
```

Output:

<p align="center">

<img src="https://github.com/Adamtaranto/tSplit/blob/main/docs/images/TIR_blast_dotplot-Polydotplot.png?raw=true" alt="TIR element blastn self-align dotplot with detected TIRs highlighted in red on the diagonal." style="width:80%;height:auto;">

</p>

## tsplit LTR

### Basic LTR processing

For each element in _LTR_retrotransposon.fa_ split into internal and external segments.

Split segments will be written to _LTR_split_tsplit_output.fasta_ with suffix "\_I" for internal or "\_LTR" for external segments.

If the `--both` flag is set, then both left and right termini will be reported with suffixes "\_L\_LTR" and "\_R\_LTR".

By default, LTRs must be at least 10bp in length and share 80% identity and occur within 10bp of each end of the input element.

```bash
tsplit LTR -i tests/data/LTR_retrotransposon.fa -d results
```

Output:

- `results/LTR_retrotransposon_tsplit_output.fasta`

### Output PAF alignments and GFF LTR annotations

```bash
tsplit LTR -i tests/data/LTR_retrotransposon.fa -d results --splitmode split --blast_evalue 0.001 --minid 60.0 --paf --gff
```

Output:

- `results/LTR_retrotransposon_tsplit_output.fasta`
- `results/LTR_retrotransposon.gff3`
- `results/LTR_retrotransposon.paf`

### Generate LTR annotated dotplot

```bash
flexidot -i tests/data/LTR_retrotransposon.fa -a results/LTR_retrotransposon.paf -m 2 -o results/blast_dotplot --gff results/LTR_retrotransposon.gff3
```

<p align="center">

<img src="https://github.com/Adamtaranto/tSplit/blob/main/docs/images/LTR_blast_dotplot-Polydotplot.png?raw=true" alt="LTR element blastn self-align dotplot with detected LTRs highlighted in red on the diagonal." style="width:80%;height:auto;">

</p>
