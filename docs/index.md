[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/tSplit.svg)](https://badge.fury.io/py/tSplit)
[![codecov](https://codecov.io/gh/Adamtaranto/tSplit/graph/badge.svg?token=24AGM1OWS5)](https://codecov.io/gh/Adamtaranto/tSplit)
[![Downloads](https://pepy.tech/badge/tsplit)](https://pepy.tech/project/tsplit)

# tSplit the TE-splitter

Extract terminal repeats from retrotransposons (LTRs) or DNA transposons (TIRs). Returns compontent segments of the element for use with transposon mapping tools.

Optionally, `tsplit TIR` can also compose synthetic MITES from complete DNA transposons.

## Table of contents

- [Algorithm overview](#algorithm-overview)
- [Options and usage](#options-and-usage)
  - [Installing tSplit](#installing-tsplit)

## Algorithm overview

tSplit attempts to identify terminal repeats in transposable elements by
first aligning each element to itself using `blastn` or `nucmer`, and then applying a set of
tuneable heuristics to select an alignment pair most likely to represent an LTR or TIR, as follows:

1. Exclude all diagonal/self-matches
2. If `tsplit LTR`: Retain only alignment pairs on the same strand (tandem repeats)
3. If `tsplit TIR`: Retain only alignment pairs on opposite strands (inverse repeats)
4. Retain pairs for which the 5' match begins within x bases of element start
   and whose 3' match ends within x bases of element end
5. If multiple candidates remain select alignment pair with largest internal segment
   (i.e. closest to element ends)

## Installing tSplit

Requirements:

- [pymummer](https://pypi.python.org/pypi/pymummer) version >= 0.10.3 with wrapper for nucmer option _--diagfactor_.
- [MUMmer4](https://github.com/mummer4/mummer)
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Optional)

You can create a Conda environment with these dependencies using the YAML file in this repo.

```bash
conda env create -f environment.yml

conda activate tsplit
```

After activating the tsplit environment you can use pip to install the latest version of tsplit.

Installation options:

1. Install from PyPi.
   This will get you the latest stable release.

```bash
pip install tsplit
```

2. Pip install directly from this git repository.

This is the best way to ensure you have the latest development version.

```bash
pip install git+https://github.com/Adamtaranto/tSplit.git
```
