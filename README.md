# Faster Algorithms for Subset Sum using Filtered Representation Trees

Complete source package, version 25.

## Contents

- `main.tex`: the complete manuscript, including all sections, tables and the AI statement.
- `references.bib` and `main.bbl`: the bibliography source and compiled bibliography.
- `verify_classical.py`, `verify_qracm.py`, `verify_qraqm.py`: three complete verification programs with all exact parameters embedded.
- `Makefile`: commands for compiling the paper and running the checks.

## Compile the paper

Upload this ZIP directly to Overleaf and select `main.tex` as the main document, with pdfLaTeX as compiler. All manuscript sources are in the ZIP root.

Locally, with a standard TeX Live installation, run:

```sh
make paper
```

This runs pdfLaTeX, BibTeX and two further pdfLaTeX passes.

## Verify the numerical constructions

Use Python 3.10 or newer:

```sh
make verify
```

The verification programs check exact probability identities and resource inequalities using 100-digit directed interval arithmetic. The conclusions retain the algorithmic assumptions stated in the paper.
