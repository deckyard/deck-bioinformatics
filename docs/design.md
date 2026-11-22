# Design

## Scope

Sequence manipulation for DNA, RNA and protein, modelled on the tool catalogue
of the [Sequence Manipulation Suite](https://stothardresearch.ca/sequence-manipulation-suite/)
(SMS): transforms, composition and physical properties, pattern / ORF /
restriction search, FASTA and FASTQ over strings, random sequences and
mutation, pairwise alignment and distances, plus the reference tables these
need (genetic codes, IUPAC codes, residue properties, BLOSUM62, a restriction
enzyme table).

It is pure computation over plain strings. Nothing under `src/` opens a file,
makes a request or runs a subprocess, so the deck needs no capability and runs
on `jennifer-tiny`.

Not included:

| | |
| - | - |
| Affine gap penalties | linear gaps only; Gotoh would be next |
| Type IIS enzymes | their cut geometry does not fit two in-site offsets |
| PDB and 3-D structure | a much larger effort |
| Multiple sequence alignment | outside the catalogue |
| Streaming and high-throughput | a different deck |
| Whole-genome alignment | alignment is quadratic; see below |

## Decisions

- **Pure `.j`, a deck and not core.** Almost all of the catalogue is string,
  list and map work over `strings`, `regex`, `lists`, `maps`, `math`,
  `decimal` and `binary`. Writing it in Jennifer keeps the tables and
  algorithms extendable by the community without Go.
- **Sequences are `string`.** A molecule kind (`DNA`, `RNA`, `PROTEIN`) is an
  argument only where it changes the answer, as in `molecularWeight`.
- **Validate first, strictly.** Every public verb normalises its input
  (upper case, whitespace and digits removed) and raises on any other
  character, naming the character and its 0-based position in the input as
  given. Gap characters (`-`, `.`) raise rather than being dropped, so
  coordinates never shift silently. `filterDna` and `filterProtein` are the
  explicit way to drop foreign characters.
- **Ambiguity codes are rejected** by mass, pI and Tm style verbs rather than
  averaged.
- **Coordinates** in result structs (`Match`, `Orf`, `Site`, `Region`,
  `Fragment`) are 0-based, half-open and on the forward strand, reverse-strand
  hits included. `subrange` is the one 1-based inclusive verb, matching how
  biologists write ranges.
- **Matches overlap.** `find("AAA", "AA")` gives two matches. IUPAC matching
  therefore does not use `regex`, whose scan is non-overlapping; it expands
  degenerate patterns into concrete sequences (up to 1024) scanned natively
  with `binary.indexOf`, and only beyond that falls back to an interpreted
  scanner. IUPAC matching reads `U` as `T` on both sides.
- **Exact CpG thresholds.** `cpgIslands` turns its `minGc` and `minObsExp`
  thresholds into integer fractions through `decimal` and cross-multiplies
  counts, so a window exactly at a threshold is decided exactly, not by float
  rounding.
- **Formats keep what they read.** FASTA and FASTQ parsing preserves case
  (lower case is soft-masking) and does not check the alphabet. FASTQ quality
  is read by length, not by spotting `@`, which is a valid quality character.
  Writers raise on any field that would not parse back (whitespace in an id,
  line breaks, `>` or `;` in a FASTA sequence, a leading `+` in a FASTQ one).
- **No RNG of its own.** Random verbs draw from `math.rand`, so
  `math.randSeed` makes a run reproducible.
- **`BLOSUM62` stores one order per pair** (300 entries); an uncovered pair
  scores as the `XX` entry.
- **`translate` drops** a trailing one or two bases rather than raising.

### Generated data, hand-written API

The genetic codes, residue properties and enzyme table are wide, repetitive data
where a transcription slip is invisible in review. `tools/gentables.j` and
`tools/genenzymes.j` generate them and hold the data in reviewable form: each
alternative code as a short diff from the standard one. They assert invariants
before writing: BLOSUM62 symmetric across all 576 pairs, every enzyme site
palindromic with its cut offsets inside the site and `cutBottom == len(site) -
cutTop`.

Generated data and hand-written accessors live in separate files (`tabledata.j`
and `tables.j`, `enzymedata.j` and `enzymes.j`), so a generator only overwrites
files holding nothing but `def const`. The first version had no split, and
regenerating deleted accessors appended to the output. CI regenerates both
files and fails when they differ from the committed ones.

### Helpers return, they do not fill

Jennifer passes arguments by value. A helper cannot fill in a list handed to
it, because it appends to a copy; helpers return what they found and the
caller concatenates. Passing the whole sequence to a helper once per
iteration copies it each time, so hot loops inline small tests instead.

## Performance

The measurements come from `examples/bench.j`, best of several runs on one
developer machine. Compare ratios, not milliseconds: on a machine about twice
as fast every conclusion still holds, so a different absolute number is not a
regression. `jennifer run examples/bench.j` shows yours.

### Byte loops, not string indexing

Jennifer strings are rune-indexed, so `$s[$i..$i + 1]` rescans from the start on
every access and a per-residue loop goes quadratic. Sequences are ASCII, so every
hot loop converts once to `bytes`, indexes that, and converts back.

| over 20 000 residues | |
| - | - |
| `$s[$i..$i + 1]` | 1545 ms |
| `$b[$i]` on a `bytes` buffer | 38 ms |

### Accretion into a bytes buffer

`lists.push` returns a new list each time, so building a result with it copies
the whole result per step:

| building a 20 000-residue result | |
| - | - |
| `lists.push` in a loop | 54 700 ms |
| `$out[] = code;` append | 19 ms |
| `$s + "A"` concatenation | 77 ms |

Parsers and writers collect pieces in a `list of string` and join once:
string `+` onto a growing result made an 8000-read FASTQ write take 28 s and an
800 kb FASTA parse 4.5 s before the change.

### Plain maps beat clever encodings

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

### Clean input takes a fast path

Validation first tries one anchored `regex.matches` against a pattern of the
admitted characters and returns the input untouched when it is already clean;
only dirty input goes through the byte loop. Building the pattern per call cost
more than the match, so the patterns are generated constants, tested against
the admitted-character lists.

### Alignment is the one quadratic surface

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
on. Memory and the `maxCells` guard are in
[Comparison](compare.md#performance-and-the-memory-limit).

The DP matrix is a flat `list of int` indexed `i * width + j`, and substitution
scores for the two sequences are precomputed into a map keyed by the two bytes,
so the inner loop is indexed reads and one lookup. The traceback recomputes
which move produced each cell instead of keeping a direction matrix: same
answer, half the memory.

## Errors

Every failure is `Error{kind: "bioinformatics"}`, raised before any work is
done, with a message naming the offending character and position where there
is one, except a reached cap below, which is the generic
`Error{kind: "limit"}` every library uses for it (since 0.4.0; until then the
caps were `"bioinformatics"` too). Invalid caller input, such as a window
longer than the sequence, keeps the deck's kind. The deck does no I/O, so no
transport or decoder kind reaches the caller.

## Security and limits

The input is treated as untrusted text; every bound is on memory or time.

- `MAX_RECORDS` (1 000 000) caps FASTA and FASTQ parsing: a file of bare
  `>` lines expands about a hundredfold into records.
- `AlignOptions.maxCells`, by default `DEFAULT_MAX_ALIGNMENT_CELLS`
  (1 000 000), bounds an alignment's DP matrix. At about a kilobyte per cell
  memory binds before time; 2000 x 2000 was an out-of-memory kill.
- `MAX_EDIT_DISTANCE_CELLS` (2 500 000, about 5 s) bounds `editDistance`,
  whose memory is linear, so time is what binds.
- IUPAC expansion is bounded (1024 concrete sequences) before the
  interpreted fallback takes over.

## Testing and validation

The overlays run in-process with fixtures embedded in the test files, so the
suite runs on `jennifer-tiny` and from a bare copy of `src/`. The tests
distinguish three kinds of check:

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

The generated tables are also checked independently of their generators: the
tests restate the genetic-code differences and re-derive enzyme overhangs and
palindromes. The examples read sample files from `testdata/` through `fs`, as
a caller would; `bench.j` runs in CI.

## Known limitations and open questions

- The protein alphabet is a superset of both nucleic-acid alphabets, so
  `isProtein("ACGT")` is true.
- Alignment uses linear gaps only, and its quadratic cost limits it to
  sequences of a few hundred residues to about a thousand.

## Dependencies on the interpreter

Floor `>=0.27.0` for `jennifer` and `jennifer-tiny`, set by the UTF-8 codecs
of `encoding.encode` / `decode` and `strings.replicate`. The exact CpG
thresholds use `decimal` (0.26). No capability.
