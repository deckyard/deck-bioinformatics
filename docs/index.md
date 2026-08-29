# @jennifer/bioinformatics

**Sequence manipulation for DNA, RNA and protein**, in the
[Jennifer](https://jennifer-lang.dev/) language.

Transforms, composition and physical properties, pattern and feature search,
FASTA and FASTQ, sequence generation, and pairwise alignment.

```jennifer
import "@jennifer/bioinformatics/" as bioseq;

def gene as string init bioseq.filterDna("1 atg acg tta cat gga taa");

bioseq.reverseComplement($gene);            # TTATCCATGTAACGTCAT
bioseq.translate($gene);                    # MTLHG*
bioseq.gcContent($gene);                    # 33.33
bioseq.molecularWeight($gene, bioseq.DNA);  # 5546.66
```

## Sequences are strings

There is no `Seq` type to construct and no wrapper to unwrap.
`bioseq.reverseComplement("ATGC")` is `"GCAT"`, and every verb in the deck takes
and returns plain strings.

The package is named for the field; the alias is named for what the calls act
on. `bioseq.reverseComplement($s)` reads better in a pipeline than
`bioinformatics.reverseComplement($s)`, and every snippet in the documentation
binds `bioseq`, so they copy-paste directly. Nothing in the deck depends on the
name the caller gives it.

The rules that follow from strings-in, strings-out:

- **Input is accepted in any case**, and residues come back uppercase.
- **A molecule kind is an argument only where it changes the answer** -
  `molecularWeight($s, bioseq.DNA)` takes one, `reverseComplement($s)` does not.
- **Offsets are 0-based and half-open**, matching Jennifer's own slicing, so
  `$s[$m.start..$m.end]` is the match. (`subrange` is the one deliberate
  exception: it is 1-based inclusive.)
- **Two tiers of tolerance.** `filterDna` / `filterRna` / `filterProtein`
  silently strip anything unrecognised; every other verb tolerates whitespace,
  digits and case but raises on a genuine foreign character, naming it and its
  position.
- **One error shape** - `Error{kind: "bioinformatics"}`, with a message written
  for a human.

[Sequences and alphabets](sequences.md) has the full treatment, and is the
page to read first.

## What is here

| | |
| - | - |
| [Transforms](transforms.md) | complement, transcribe, translate, six-frame reading, residue spellings |
| [Composition](composition.md) | GC content, molecular weight, melting temperature, isoelectric point, GRAVY, codon usage |
| [Search](search.md) | exact, ambiguous and approximate patterns; open reading frames; CpG islands; restriction mapping and digestion |
| [Formats](formats.md) | FASTA and FASTQ, parsing and writing, quality scores |
| [Manipulation](manipulate.md) | random sequences, composition-preserving shuffle, mutation, windows |
| [Comparison](compare.md) | Needleman-Wunsch, Smith-Waterman, Hamming and edit distance |
| [Genetic codes and residues](tables.md) | nine NCBI codes, IUPAC ambiguity, per-residue properties |
| [Restriction enzymes](enzymes.md) | 76 hand-curated enzymes and their cut geometry |

The [cheatsheet](cheatsheet.md) is the whole surface on one page. The
[design notes](design.md) explain the decisions behind it, including the
ones that were measured rather than assumed.

## It runs on `jennifer-tiny`

Nothing under `src/` touches the network, a subprocess, or the filesystem.
FASTA and FASTQ are parsed from and written to **strings**; opening the file is
the caller's job, which `examples/` shows. There is no
`pragma-jennifer-capability` anywhere in the deck, by design - and the
consequence is that every function is a pure transformation that can be tested
without a fixture on disk.

`deck.toml` says so out loud: `[engines]` lists `jennifer-tiny` beside
`jennifer`. That table is an allowlist of engines rather than a version floor,
so a manifest naming only the full build would have a tiny run *refused* - and
CI runs the whole example set under both on every push.

## Every number names its source

This deck reports melting temperatures, isoelectric points and molecular
weights. Those depend on published parameter sets that disagree with each other
by small but real amounts, so each function names the set it uses in its own
documentation, and the parameters themselves are exported constants. Where a
value can be checked against an independent source it is, and the tests say
which checks are external anchors and which are not - see
[Design notes](design.md#correctness).
