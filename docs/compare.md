# Comparison

Pairwise alignment and distance.

Runnable version: `jennifer run examples/compare.j`.

## Which alignment

**Global** - [Needleman-Wunsch](https://doi.org/10.1016/0022-2836(70)90057-4) -
aligns both sequences end to end, padding with gaps. Every residue of both
appears in the result. This is the right question for two sequences believed
homologous over their whole length: two alleles of a gene, or the same protein
from two species.

**Local** - [Smith-Waterman](https://doi.org/10.1016/0022-2836(81)90087-5) -
finds the single best-matching stretch and discards everything outside it. This
is the right question for a short sequence against a long one, or for a shared
domain inside two otherwise unrelated proteins.

Using the wrong one is the most common mistake: a global alignment of a probe
against a chromosome spends its entire score on terminal gaps.

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

`AlignOptions` carries `matchScore`, `mismatchScore`, `gapScore` and `matrix`.
When `matrix` is non-empty it replaces the match/mismatch pair entirely.

Gap scoring is **linear**: every gap position costs the same. Affine gaps - a
larger penalty to open than to extend, which models indels better - are the
obvious next step and are deliberately not in this release.

### Why protein needs a matrix

Leucine against isoleucine is a conservative substitution that barely changes a
protein. Leucine against aspartate is not. A match/mismatch score cannot express
that difference; a substitution matrix can:

| pair | BLOSUM62 | |
| - | - | - |
| L → L | +4 | identical |
| L → I | +2 | conservative |
| K → R | +2 | conservative |
| F → Y | +3 | conservative |
| L → D | -4 | disruptive |
| W → P | -4 | disruptive |

```jennifer
bioseq.substitutionScore(bioseq.BLOSUM62, "W", "W");   # 11
bioseq.substitutionScore(bioseq.BLOSUM62, "L", "I");   # 2
```

`BLOSUM62` is stored as a symmetric upper triangle, so read it through
`substitutionScore`, which tries both orders, rather than indexing it directly.

## Reading an alignment

An `Alignment` carries `a`, `b`, `markup`, `score`, `length`, `identities`,
`similarities`, `gaps`, `startA` and `startB`.

`a` and `b` are the sequences with gaps inserted; strip the gaps and you have the
originals back. `markup` is one character per column:

| | |
| - | - |
| `\|` | identical residues |
| `:` | different, but scoring positively - a conservative substitution |
| space | anything else, including a gap |

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
bioseq.formatAlignment($al, 60);   # numbered blocks, EMBOSS style
```

Identity is measured over the alignment's own length, gaps included - so a
global alignment of two sequences of very different lengths reports a low
identity even when the shorter matches perfectly. That is the honest answer to
the question asked; ask `alignLocal` a different one.

Similarity only differs from identity under a substitution matrix. With plain
match/mismatch scoring nothing but an identity can score positively, so the two
are equal.

### Local alignment coordinates

`startA` and `startB` locate the alignment in the sequences you passed in:

```jennifer
def matched as string init strings.replace($al.a, "-", "");
$one[$al.startA..$al.startA + len($matched)] == $matched;   # true
```

Note that a local alignment's **boundaries are soft** unless the flanks are
genuinely dissimilar. Smith-Waterman keeps extending outward while the running
score rises, and with +2 for a match against -1 for a mismatch, random sequence
still scores above zero on average. Expect the reported span to reach a little
past the region you were looking for.

## Distances

```jennifer
bioseq.hammingDistance("ATGC", "ATCC");   # 1
bioseq.editDistance("ATGC", "AGC");       # 1
bioseq.editDistance("kitten", "sitting"); # 3
```

**Hamming** counts differing positions and therefore requires equal lengths - it
has no way to express an insertion. It refuses unequal lengths rather than
guessing.

**Edit distance** (Levenshtein) allows insertions and deletions as well, each
costing one.

Where the two part company is instructive. Shift a sequence by one position:

| | |
| - | - |
| Hamming | **12** - every position looks wrong |
| edit | **2** - one deletion and one insertion |

Computed with two rows rather than a full matrix, so edit distance costs
`O(min(n, m))` memory. The time is still `O(n·m)`.

## Performance

Alignment is this deck's one `O(n·m)` surface. Measured at about **70 000 cells
a second**, with cost quadratic in sequence length:

| | cells | time |
| - | - | - |
| 50 × 50 | 2 500 | 86 ms |
| 100 × 100 | 10 000 | 150 ms |
| 200 × 200 | 40 000 | 458 ms |
| 400 × 400 | 160 000 | 2.28 s |

But **memory binds first, and binds harder**. Every matrix entry is a boxed
interpreter value, measured at roughly **a kilobyte per cell**:

| | cells | resident |
| - | - | - |
| 100 × 100 | 10 000 | 35 MB |
| 350 × 350 | 122 500 | 151 MB |
| 700 × 700 | 490 000 | 570 MB |

So 1000 × 1000 is about a gigabyte, and 2000 × 2000 is not a wait but an
out-of-memory kill. That is why `AlignOptions.maxCells` exists and is set by
default to `DEFAULT_MAX_ALIGNMENT_CELLS` - a million cells, or a thousand
residues square. Aligning two untrusted sequences would otherwise be a reliable
way to take the process down.

```jennifer
def opts as bioseq.AlignOptions init bioseq.dnaOptions();
$opts.maxCells = 4000000;   # you have the memory and you mean it
```

Setting it to `0` removes the limit. Reproduce the numbers with
`jennifer run examples/bench.j`; the reasoning is in
[Design notes](design.md#alignment-is-the-one-quadratic-surface).

## The surface

| | |
| - | - |
| `AlignOptions` | `matchScore`, `mismatchScore`, `gapScore`, `matrix` |
| `dnaOptions()` / `proteinOptions()` | defaults |
| `substitutionScore(matrix, a, b)` | order-independent lookup |
| `alignGlobal(a, b, opts)` | Needleman-Wunsch |
| `alignLocal(a, b, opts)` | Smith-Waterman |
| `percentIdentity(al)` / `percentSimilarity(al)` / `percentGaps(al)` | |
| `formatAlignment(al, width)` | numbered blocks |
| `hammingDistance(a, b)` / `editDistance(a, b)` | |
| `BLOSUM62` | the substitution matrix |
