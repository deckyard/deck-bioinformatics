# @deckyard/bioinformatics

Sequence manipulation for DNA, RNA and protein in
[Jennifer](https://jennifer-lang.dev/), modelled on the
[Sequence Manipulation Suite](https://www.bioinformatics.org/sms2/).

This deck is part of [Deckyard](https://github.com/deckyard), community decks for [Jennifer](https://jennifer-lang.dev/).

```jennifer
import "@deckyard/bioinformatics/" as bioseq;

def gene as string init bioseq.filterDna("1 atg acg tta cat gga taa");

bioseq.reverseComplement($gene);            # TTATCCATGTAACGTCAT
bioseq.translate($gene);                    # MTLHG*
bioseq.gcContent($gene);                    # 33.33
bioseq.molecularWeight($gene, bioseq.DNA);  # 5546.66
```

Sequences are plain strings, offsets are 0-based half-open, and every error is
`Error{kind: "bioinformatics"}`. Read [Sequences and alphabets](sequences.md)
first: it states the conventions every other page relies on. Every snippet in
this book imports the deck as `bioseq`; nothing depends on the alias.

| | |
| - | - |
| [Transforms](transforms.md) | complement, transcribe, translate, six frames, residue spellings |
| [Composition](composition.md) | GC content, molecular weight, melting temperature, pI, GRAVY, codon usage |
| [Search](search.md) | exact, IUPAC and approximate patterns; ORFs; CpG islands; restriction mapping and digestion |
| [Formats](formats.md) | FASTA and FASTQ, quality scores |
| [Manipulation](manipulate.md) | random sequences, shuffle, mutation, windows |
| [Comparison](compare.md) | Needleman-Wunsch, Smith-Waterman, Hamming and edit distance |
| [Genetic codes and residues](tables.md) | nine NCBI codes, IUPAC codes, residue properties, parameter constants |
| [Restriction enzymes](enzymes.md) | 76 hand-curated type II enzymes |
| [Cheatsheet](cheatsheet.md) | every export on one page |
| [Design notes](design.md) | measurements behind the implementation, correctness anchors |

## No I/O, so it runs on `jennifer-tiny`

Nothing under `src/` touches the network, a subprocess or the filesystem. FASTA
and FASTQ are parsed from and written to strings; reading the file is the
caller's job ([Formats](formats.md#the-io-boundary)). There is no
`pragma-jennifer-capability` in the deck. `[engines]` in `deck.toml` is an
allowlist, not a version floor, so it lists `jennifer-tiny` beside `jennifer`;
CI runs every example under both.

Melting temperatures, isoelectric points and masses depend on published
parameter sets that disagree with each other. Each function names its set in
its docblock and the parameters are exported constants
([Parameter constants](tables.md#parameter-constants)).
