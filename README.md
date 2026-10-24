# @deckyard/bioinformatics

[![test](https://github.com/deckyard/deck-bioinformatics/actions/workflows/test.yml/badge.svg)](https://github.com/deckyard/deck-bioinformatics/actions/workflows/test.yml)
[![Documentation](https://img.shields.io/badge/docs-github%20pages-blue)](https://deckyard.github.io/deck-bioinformatics/)
[![release](https://img.shields.io/github/v/release/deckyard/deck-bioinformatics)](https://github.com/deckyard/deck-bioinformatics/releases/latest)

Sequence manipulation for DNA, RNA and protein in
[Jennifer](https://jennifer-lang.dev/), modelled on the Sequence Manipulation
Suite: transforms, composition and physical properties, pattern, ORF and
restriction-site search, FASTA and FASTQ, random sequences, and pairwise
alignment. Sequences are plain strings. Nothing under `src/` touches the
network, a subprocess or the filesystem, so the deck needs no capability and
runs on `jennifer-tiny`.

**[Documentation](https://deckyard.github.io/deck-bioinformatics/)**

## Install

```sh
jvc add @deckyard/bioinformatics
```

This vendors the deck to `vendor/deckyard/bioinformatics/`. Without `jvc`, unpack
the `src/` subtree of a release archive
(`deckyard-bioinformatics-<version>.tar.gz`) there. Requires Jennifer 0.25.0 or
newer; no deck dependencies.

## Usage

```jennifer
import "@deckyard/bioinformatics/" as bioseq;

def gene as string init bioseq.filterDna("1 atg acg tta cat gga taa");

bioseq.reverseComplement($gene);            # TTATCCATGTAACGTCAT
bioseq.translate($gene);                    # MTLHG*
bioseq.gcContent($gene);                    # 33.33
bioseq.molecularWeight($gene, bioseq.DNA);  # 5546.66
```

## Development

```sh
tools/check-all.sh                          # the full gate, as CI runs it
jennifer test src/bioinformatics_test.j     # the test suite
jennifer run examples/transforms.j          # one example per topic in examples/
jennifer lint src/bioinformatics.j src/bioinformatics_test.j
```

A topic file cannot be linted alone: it is `include`d, so lint the entry module
and the test overlay, which splice every topic. `src/topics/tabledata.j` and
`src/topics/enzymedata.j` are generated; see [`tools/README.md`](tools/README.md).

The book's source is [`docs/`](docs/), built with
[Grimoire](https://github.com/jennifer-language/grimoire). To read it locally:

```sh
git clone https://github.com/jennifer-language/grimoire .grimoire
jennifer run .grimoire/bin/grimoire serve
```

## Licence

LGPL-3.0-only, see [LICENSE](LICENSE).
