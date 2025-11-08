# BioSeqToolkit — GitHub repo skeleton
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![Issues](https://img.shields.io/github/issues/HimanshuNeve/BioSeqToolkit)](https://github.com/HimanshuNeve/BioSeqToolkit/issues)

This document contains everything you need to create a clean GitHub repository for your sequence parsing pipeline. It includes:

* `README.md` (detailed usage and examples)
* `LICENSE` (MIT)
* `.gitignore`
* `scripts/` (primary script `genbank_to_fasta_no_biopy_v2.py` — single-file, Biopython-free)
* `examples/` (small example `seq.gb` and `seq.fasta` you can test with)
* `analysis/` (small `analyze_results.py` helper)
* `.vscode/launch.json` (optional VS Code run config)

Copy each file into your repo folder, then run the git commands in the README to push to GitHub.

---

## README.md

````markdown
# BioSeqToolkit

This repository provides a compact, Biopython-free toolkit to extract, parse, and analyze biological sequence data from GenBank and FASTA files. It generates comprehensive FASTA outputs including full sequences, CDS nucleotide sequences, translated proteins, and peptide feature sets.

## Features

- Works for all biological sequences (plants, microbes, animals, viruses, etc.).
- Accepts single files, directories, or glob patterns.
- Auto-detects `ORIGIN` and `SEQUENCE` blocks in GenBank files.
- Parses `FEATURES` including `CDS`, `mat_peptide`, and `sig_peptide` and extracts their sequences.
- Supports `join(...)` and `complement(...)` location expressions.
- Translates CDS with built-in codon table (honors `/codon_start`).
- Produces timestamped output folders to avoid overwrites.
- No Biopython dependency (standard library only). Optionally uses `pandas` for CSV summaries.

## Quick start

1. Clone this repo (or create a new one and add files):

```bash
git clone <your-repo-url>
cd BioSeqToolkit
````

2. Put your input files into `data/` (create that folder):

```
data/
  example.gb
  example.fasta
```

3. Run the script:

```bash
python genbank_to_fasta_no_biopy_v2.py data/example.gb --translate --to-stop
# or process whole folder
python genbank_to_fasta_no_biopy_v2.py data --translate
```

4. Check outputs in `results/` (timestamped subfolders) and `results/aggregate_summary.csv`.

## Requirements

* Python 3.8+ (no external packages required).
* Optional: `pandas` for CSV summary output: `pip install pandas`.

## Files

* `genbank_to_fasta_no_biopy_v2.py` — main parser script (see `scripts/`)
* `analyze_results.py` — analysis and QC helper
* `examples/` — small example GenBank and FASTA files for testing

## Development & Contribution

* Keep code simple and modular.
* Add tests under `tests/` for parser functions.
* Follow consistent naming for new CLI arguments.

## License

MIT

````

---

## LICENSE (MIT)

```text
MIT License

Copyright (c) 2025 <Your Name>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
````

---

## .gitignore

```text
__pycache__/
*.py[cod]
*$py.class
.env/
.venv/
venv/
.vscode/
results/
.DS_Store
Thumbs.db
```

---

## scripts/genbank_to_fasta_no_biopy_v2.py

(Place the full script you use here.)

---

## examples/seq.gb (small test file)

```text
LOCUS       NM_001204686             968 bp    mRNA    linear   INV 27-APR-2025
DEFINITION  Aplysia californica insulin precursor (PIN), mRNA.
...
ORIGIN
        1 cctgaatata gccaactaaa ttctaggaac tctaagagga ctacgcttgt ctccaacatc
       61 ttatcgtcaa catcttctgc aagcgataac tatatttctg gtccgccaaa gtagtatacg
      ...
//
```

---

## examples/seq.fasta

```fasta
>example_seq
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
```

---

## analysis/analyze_results.py

```python
import csv
from pathlib import Path

agg = Path('results/aggregate_summary.csv')
if not agg.exists():
    print('No aggregate summary found. Run the main script first.')
    raise SystemExit(1)

with agg.open() as fh:
    reader = csv.DictReader(fh)
    rows = list(reader)

print('Processed files:')
for r in rows:
    print('-', r['file'], '->', r.get('outdir'), 'qc:', r.get('qc'))

try:
    import pandas as pd
    df = pd.read_csv(agg)
    print('\nSummary (pandas):')
    print(df.describe(include='all'))
except Exception:
    pass
```

---

## .vscode/launch.json

```json
{
  "version": "0.2.0",
  "configurations": [
    {
      "name": "Run parser",
      "type": "python",
      "request": "launch",
      "program": "${workspaceFolder}/genbank_to_fasta_no_biopy_v2.py",
      "args": ["data/seq.gb", "--translate", "--to-stop"],
      "console": "integratedTerminal"
    }
  ]
}
```

---

## Git initialization

```bash
git init
git add .
git commit -m "Initial commit: BioSeqToolkit skeleton"
git branch -M main
```

To add and push to GitHub:

```bash
git remote add origin https://github.com/<your-username>/BioSeqToolkit.git
git push -u origin main
```

---

## Notes

* `results/` is ignored by Git.
* Use semantic commit messages.
* Keep example files small (<10 KB) for fast tests.

---
