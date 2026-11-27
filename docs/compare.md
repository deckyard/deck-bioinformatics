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

When no residue pair scores above zero the local alignment is empty: `length`
and `score` 0, empty rows, `startA` and `startB` equal to the two lengths (an
empty slice at the end). The readouts and `formatAlignment` raise on it.

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
sequence's length, O(n*m) time. `maxCells` does not apply; instead
`len(a) * len(b)` above `MAX_EDIT_DISTANCE_CELLS` (2 500 000; 1500 x 1500 fits)
raises kind `"limit"` before any work.

| `editDistance`, 1500 x 1500 (2 250 000 cells) | time | per cell |
| - | - | - |
| Jennifer 0.27.0 | 5.5 s | about 2.4 microseconds |
| Jennifer 0.28.0 | 3.0 s | about 1.3 microseconds |

## Performance and the memory limit

Alignment is the deck's one O(n*m) surface: about 380 000 cells a second on
Jennifer 0.27.0, about 500 000 on 0.28.0.

| | cells | 0.27.0 | 0.28.0 |
| - | - | - | - |
| 50 x 50 | 2 500 | 8 ms | 4 ms |
| 100 x 100 | 10 000 | 31 ms | 18 ms |
| 200 x 200 | 40 000 | 120 ms | 79 ms |
| 400 x 400 | 160 000 | 420 ms | 300 ms |

Memory runs out first: every matrix entry is a boxed interpreter value, about
1.1 KB a cell on 0.27.0 and about 0.2 KB on 0.28.0.

| | cells | 0.27.0 resident | 0.28.0 resident | 0.28.0 time |
| - | - | - | - | - |
| 100 x 100 | 10 000 | 35 MB | 40 MB | 0.07 s |
| 350 x 350 | 122 500 | 151 MB | 50 MB | 0.3 s |
| 700 x 700 | 490 000 | 570 MB | 108 MB | 1.0 s |
| 1000 x 1000 | 1 000 000 | 1135 MB | 221 MB | 2.4 s |
| 1414 x 1414 | 2 000 000 | not measured | 339 MB | 4.2 s |
| 2000 x 2000 | 4 000 000 | out of memory | 636 MB | 7.1 s |

So `AlignOptions.maxCells` defaults to `DEFAULT_MAX_ALIGNMENT_CELLS`
(1 000 000). The matrix has `(len(a) + 1) * (len(b) + 1)` cells, and an
alignment over the limit raises kind `"limit"` before allocating. Without the
limit, aligning two untrusted sequences can take the process down.

```jennifer
def opts as bioseq.AlignOptions init bioseq.dnaOptions();
$opts.maxCells = 4000000;   # about 4 GB on 0.27.0, 0.6 GB on 0.28.0
```

On 0.28.0 a limit of 2 000 000 cells is reasonable; the default stays at
1 000 000 while the deck's floor is 0.27.0.

`0` removes the limit. `jennifer run examples/bench.j` reproduces the numbers;
see [Design notes](design.md#alignment-is-the-one-quadratic-surface).
