# Comparison

Runnable: `jennifer run examples/compare.j`.

## Global or local

- **Global**, [Needleman-Wunsch](https://doi.org/10.1016/0022-2836(70)90057-4)
  (`alignGlobal`): both sequences end to end, every residue in the result. For
  sequences homologous over their whole length.
- **Local**, [Smith-Waterman](https://doi.org/10.1016/0022-2836(81)90087-5)
  (`alignLocal`): the single best-matching stretch. For a short sequence against
  a long one, or a shared domain. A global alignment of a probe against a long
  sequence spends its score on terminal gaps.

```jennifer
def opts as bioseq.AlignOptions init bioseq.dnaOptions();
def al as bioseq.Alignment init bioseq.alignGlobal($one, $other, $opts);
def best as bioseq.Alignment init bioseq.alignLocal($one, $other, $opts);
```

## Scoring

```jennifer
bioseq.dnaOptions();       # +2 match, -1 mismatch, -2 gap
bioseq.proteinOptions();   # BLOSUM62, -8 gap
```

A non-empty `matrix` replaces `matchScore` and `mismatchScore`. Gaps are linear:
every gap position costs `gapScore`. Affine gaps are not implemented.

A substitution matrix scores conservative substitutions above disruptive ones,
which a match/mismatch pair cannot:

| pair | BLOSUM62 |
| - | - |
| L-L | +4 |
| L-I | +2 |
| K-R | +2 |
| F-Y | +3 |
| L-D | -4 |
| W-P | -4 |

```jennifer
bioseq.substitutionScore(bioseq.BLOSUM62, "W", "W");   # 11
bioseq.substitutionScore(bioseq.BLOSUM62, "L", "I");   # 2
```

`BLOSUM62` stores one order of each pair (300 entries over the twenty residues,
`B`, `Z`, `X` and `*`), so read it through `substitutionScore`, which tries both.
An uncovered pair scores the matrix's `XX` entry (-1), or 0 for a matrix without
one.

## Reading an alignment

`a` and `b` are the aligned sequences with `-` for gaps. For a global alignment,
removing the gaps gives the normalised inputs back. `markup` has one character per column: `|` identical,
`:` different but scoring positively, space otherwise.

```
MALWMRLLPLLALLALWGPDPAAA
||||:::|||:|||||||||||||
MALWLKILPLVALLALWGPDPAAA

83.3% identical, 100.0% similar
```

```jennifer
bioseq.percentIdentity($al);
bioseq.percentSimilarity($al);
bioseq.percentGaps($al);
bioseq.formatAlignment($al, 60);   # EMBOSS-style numbered blocks
```

The percentages are over the alignment's length, gap columns included, so a
global alignment of sequences of different lengths reports low identity even
when the shorter matches perfectly. Without a matrix, similarity equals
identity. `formatAlignment` numbers rows by 1-based residue position in the
original sequences, gaps not counted. These four raise on an empty `Alignment`
or one whose row lengths disagree with `length` (possible when built by hand).

### Local alignment coordinates

`startA` and `startB` are 0-based offsets of the aligned region (0 for a global
alignment):

```jennifer
def matched as string init strings.replace($al.a, "-", "");
$one[$al.startA..$al.startA + len($matched)] == $matched;   # true
```

Local boundaries are soft. With +2/-1 scoring, random sequence still scores
above zero on average, so the reported span can reach past the region of
interest.

## Distances

```jennifer
bioseq.hammingDistance("ATGC", "ATCC");   # 1
bioseq.editDistance("ATGC", "AGC");       # 1
bioseq.editDistance("kitten", "sitting"); # 3
```

**Hamming** counts differing positions and raises on unequal lengths.
**Edit distance** (Levenshtein) counts insertions, deletions and substitutions,
each costing one. A sequence shifted by one position (12 bases in
`examples/compare.j`):

| | |
| - | - |
| Hamming | 12 |
| edit | 2 (one deletion, one insertion) |

Edit distance keeps two rows of `len(b) + 1` entries: O(m) memory in the second
sequence's length, O(n*m) time, and no `maxCells` limit.

## Performance and the memory limit

Alignment is the deck's one O(n*m) surface: about 380 000 cells a second.

| | cells | time |
| - | - | - |
| 50 x 50 | 2 500 | 8 ms |
| 100 x 100 | 10 000 | 31 ms |
| 200 x 200 | 40 000 | 120 ms |
| 400 x 400 | 160 000 | 420 ms |

Memory runs out first: every matrix entry is a boxed interpreter value, about a
kilobyte per cell.

| | cells | resident |
| - | - | - |
| 100 x 100 | 10 000 | 35 MB |
| 350 x 350 | 122 500 | 151 MB |
| 700 x 700 | 490 000 | 570 MB |

1000 x 1000 is about a gigabyte; 2000 x 2000 is an out-of-memory kill. So
`AlignOptions.maxCells` defaults to `DEFAULT_MAX_ALIGNMENT_CELLS` (1 000 000).
The matrix has `(len(a) + 1) * (len(b) + 1)` cells, and an alignment over the
limit raises before allocating. Without the limit, aligning two untrusted
sequences can take the process down.

```jennifer
def opts as bioseq.AlignOptions init bioseq.dnaOptions();
$opts.maxCells = 4000000;   # about 4 GB
```

`0` removes the limit. `jennifer run examples/bench.j` reproduces the numbers;
see [Design notes](design.md#alignment-is-the-one-quadratic-surface).
