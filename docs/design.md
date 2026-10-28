# Design notes

The measurements come from `examples/bench.j`, best of several runs on one
developer machine. Compare ratios, not milliseconds: on a machine about twice
as fast every conclusion still holds, so a different absolute number is not a
regression. `jennifer run examples/bench.j` shows yours.

## Byte loops, not string indexing

Jennifer strings are rune-indexed, so `$s[$i..$i + 1]` rescans from the start on
every access and a per-residue loop goes quadratic. Sequences are ASCII, so every
hot loop converts once to `bytes`, indexes that, and converts back.

| over 20 000 residues | |
| - | - |
| `$s[$i..$i + 1]` | 1545 ms |
| `$b[$i]` on a `bytes` buffer | 38 ms |

## Accretion into a bytes buffer

`lists.push` returns a new list each time, so building a result with it copies
the whole result per step:

| building a 20 000-residue result | |
| - | - |
| `lists.push` in a loop | 54 700 ms |
| `$out[] = code;` append | 19 ms |
| `$s + "A"` concatenation | 77 ms |

## Plain maps beat clever encodings

Three implementations of `translate`, over 10 000 codons:

| | |
| - | - |
| string-slice codon key | 1094 ms |
| key rebuilt from bytes, map lookup | **63 ms** |
| base-4 packed index into a 64-entry string | 78 ms |

The packed index lost because its three helper calls per codon cost more than
the lookup they replaced. So the genetic codes stay plain
`map of string to string`, and a helper call inside a per-residue loop is
avoided: case folding is written inline as `if ($c >= 97 and $c <= 122)`.

## Alignment is the one quadratic surface

Everything else is linear in sequence length.

| | cells | |
| - | - | - |
| 50 x 50 | 2 500 | 8 ms |
| 100 x 100 | 10 000 | 31 ms |
| 200 x 200 | 40 000 | 120 ms |
| 400 x 400 | 160 000 | 420 ms |

About 380 000 cells a second: linear in cells, so quadratic in sequence length.
A few hundred residues a side is comfortable; a chromosome is out of scope. The
50 x 50 row carries some fixed overhead; the linear shape shows from 100 x 100
on. Memory and the
`maxCells` guard are in [Comparison](compare.md#performance-and-the-memory-limit).

The DP matrix is a flat `list of int` indexed `i * width + j`, and substitution
scores for the two sequences are precomputed into a map keyed by the two bytes,
so the inner loop is indexed reads and one lookup. The traceback recomputes
which move produced each cell instead of keeping a direction matrix: same
answer, half the memory.

## Helpers return, they do not fill

Jennifer passes arguments by value. A helper cannot fill in a list handed to
it, because it appends to a copy; helpers return what they found and the
caller concatenates.

## No filesystem, no network

Nothing under `src/` opens a file, makes a request or runs a subprocess. The
deck therefore runs on `jennifer-tiny` with no capability, and every test
fixture is embedded in the test files. The examples that read `testdata/`
(`formats.j`, `composition.j`, `search.j`, `manipulate.j`, `compare.j`) do it
with `fs` as the caller; `formats.j` is the one written to show that boundary.

## Generated data, hand-written API

The genetic codes, residue properties and enzyme table are wide, repetitive data
where a transcription slip is invisible in review. `tools/gentables.j` and
`tools/genenzymes.j` generate them and hold the data in reviewable form: each
alternative code as a short diff from the standard one. They assert invariants
before writing: BLOSUM62 symmetric across all 576 pairs, every enzyme site
palindromic with its cut offsets inside the site.

Generated data and hand-written accessors live in separate files (`tabledata.j`
and `tables.j`, `enzymedata.j` and `enzymes.j`), so a generator only overwrites
files holding nothing but `def const`. The first version had no split, and
regenerating deleted accessors appended to the output.

## Correctness

The tests distinguish three kinds of check:

- **External anchors**, values published independently of the deck.
  Bradykinin's average mass, a mass-spectrometry calibration standard: the deck
  computes 1060.22 against a published ~1060.2, which validates the residue
  masses and the residues-plus-water convention together. Translation is
  anchored on the human preproinsulin CDS and its published protein, so a slip
  in the codon table shows up as a wrong protein.
- **Cross-implementation checks**, values computed separately from the same
  published parameters: the nearest-neighbour Tm and the isoelectric points.
  They catch coding errors, not parameter errors, and are labelled as such.
- **Property checks**, invariants independent of parameters: charge falls
  monotonically with pH; Tm rises with GC, length and salt; each residue's
  synonymous codon fractions sum to 1.0; a digest's fragments reassemble to the
  original; an alignment's rows are the inputs with gaps added.

## Not included

| | |
| - | - |
| Affine gap penalties | linear gaps only; Gotoh would be next |
| Type IIS enzymes | their cut geometry does not fit two in-site offsets |
| PDB and 3-D structure | a much larger effort |
| Multiple sequence alignment | outside the catalogue |
| Streaming and high-throughput | a different deck |
