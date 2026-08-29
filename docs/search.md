# Search

Finding things inside a sequence: patterns, reading frames, islands, and where
an enzyme cuts.

Runnable version: `jennifer run examples/search.j`.

## Patterns

```jennifer
bioseq.find($dna, "GGATCC");            # exact
bioseq.findIupac($dna, "CCWGG");        # ambiguity codes in the pattern
bioseq.findFuzzy($dna, "GAATTC", 1);    # up to n substitutions
```

Each returns `list of Match`, with `start`, `end`, `text`, `strand` and
`mismatches`. Offsets are 0-based and half-open, so the match slices back out:

```jennifer
$dna[$m.start..$m.end] == $m.text;   # always true
```

### Matches overlap

```jennifer
len(bioseq.find("AAAA", "AA"));   # 3, not 1
```

A tandem repeat and a site nested inside another site are both real findings, and
a non-overlapping scan would silently drop the second. This is also why the
IUPAC matcher walks the sequence itself rather than delegating to `regex`, whose
scan is non-overlapping; [`iupacClass`](tables.md#iupac-ambiguity-codes) stays
available for composing your own patterns.

### Ambiguity is one-way

A pattern's `N` admits a subject's `A`. A pattern's `A` does not admit a
subject's `N`:

```jennifer
len(bioseq.findIupac("AAAA", "N"));   # 4
len(bioseq.findIupac("NNNN", "A"));   # 0
```

`findIupac` reports the **concrete** match it found, not the pattern - so
searching `CCWGG` tells you whether each hit was `CCAGG` or `CCTGG`.

### Approximate matching

`findFuzzy` counts substitutions only - no insertions or deletions - which is
what a primer-binding or probe-hybridisation question usually means. Each match
carries its `mismatches`, so exact hits are distinguishable from near ones. When
indels matter, use [`editDistance`](compare.md#distances) instead.

## Open reading frames

```jennifer
for (def o in bioseq.orfs($dna, 50)) {
    io.printf("frame %d  %s\n", $o.frame, $o.protein);
}
```

All six frames are scanned. Within a frame, reading begins at the first start
codon and runs to the first in-frame stop; the next ORF begins at the next start
after that stop. So ORFs never overlap within a frame, and each is the longest
reading from its region rather than one report per nested start codon.

An `Orf` carries `frame`, `strand`, `start`, `end`, `length`, `nucleotides`,
`protein` and `stopped`.

### Frames that run off the end

`stopped` is false when the frame reached the end of the sequence without a stop
codon. That is reported rather than discarded, because on a fragment it usually
means the gene simply continues past what you have:

```jennifer
def o as bioseq.Orf init bioseq.orfs($insulinFragment, 20)[0];
$o.stopped;   # false
$o.protein;   # MALWMRLLPLLALLALWGPDPAAA - the CDS continues
```

### Reverse-strand coordinates

An ORF found on the minus strand reports `start` and `end` in the coordinates of
the sequence **you passed in**, not of its reverse complement. `nucleotides`
holds the frame as read, so:

```jennifer
bioseq.reverseComplement($dna[$o.start..$o.end]) == $o.nucleotides;   # strand -1
$dna[$o.start..$o.end] == $o.nucleotides;                            # strand +1
```

This is the detail most easily got wrong in an ORF finder, so it is asserted
directly against the sequence in the test suite rather than against a remembered
number.

### Start codons follow the table

```jennifer
bioseq.orfsWith($dna, 50, 11);   # bacterial
```

`orfsWith` takes an NCBI table **id**, not a codon map - because an ORF search
needs the table's start codons as well as its assignments, and that is precisely
what separates some tables. `BACTERIAL` assigns every codon identically to
`STANDARD` but recognises seven starts rather than three, so the same sequence
yields different ORFs under each.

Worth knowing: `CTG` and `TTG` are start codons under the *standard* code too,
not only `ATG`. That is the sort of thing a table lookup gets right and an
assumption gets wrong.

## CpG islands

```jennifer
bioseq.cpgIslands($promoter, 200, 50.0, 0.6);   # the classic thresholds
```

CpG is scarce in vertebrate genomes: methylated cytosine deaminates to thymine,
so the dinucleotide has been eroded almost everywhere except where methylation is
kept off - which is where promoters sit. An island is a stretch that has escaped
that erosion.

Three thresholds must be met at once, following Gardiner-Garden and Frommer
(1987): the window must be long enough, GC-rich enough, and carry at least the
given ratio of observed to expected CpG, where expected is
`(C count × G count) / length`.

That third test is what separates *CpG-rich* from merely *GC-rich*:

```jennifer
def gcOnly as string init strings.repeat("G", 60) + strings.repeat("C", 60);
bioseq.gcContent($gcOnly);                             # 100.0
len(bioseq.cpgIslands($gcOnly, 20, 50.0, 0.6));        # 0
```

Qualifying windows that touch or overlap are merged, and each `Region`'s `gc`
and `obsExp` are recomputed across the merged span rather than carried over from
one window. Island edges are approximate by nature: a window straddling a
boundary is still half-dense and still qualifies, so the reported span runs a
little wider than the dense core.

## Restriction sites

```jennifer
bioseq.restrictionSites($plasmid, ["EcoRI", "BamHI", "HindIII"]);
```

Each `Site` gives `enzyme`, `start`, `end`, `cut`, `site` and `strand`, ordered
by cut position. `site` is the sequence actually found, which for a degenerate
recognition sequence is the concrete match:

```jennifer
bioseq.restrictionSites("AAACCAGGTTTCCTGGAAA", ["EcoRII"]);
# EcoRII reads CCWGG; the two hits report "CCAGG" and "CCTGG"
```

Every recognition sequence in the deck's table is palindromic, so a forward scan
finds every site. A non-palindromic site added later is searched on both strands
and reported with `strand: -1` where it is found only on the reverse.

### Which enzymes to reach for

```jennifer
bioseq.enzymesCutting($insert);   # every enzyme in the table with a site
bioseq.nonCutters($insert);       # every enzyme with none
```

The second is usually the more useful: the job is often to linearise a construct
exactly once, or to leave an insert intact.

## Digestion

```jennifer
bioseq.restrictionDigest($plasmid, ["EcoRI", "BamHI"], false);   # linear
bioseq.restrictionDigest($plasmid, ["EcoRI", "BamHI"], true);    # circular
```

Cuts are taken on the top strand at each enzyme's `cutTop` offset, which is what
determines fragment sizes on a gel; the staggered cut on the bottom strand
decides the ends' chemistry rather than the lengths, and is reported through
[`Enzyme.overhang`](enzymes.md) instead.

The invariant: the fragments, in order, are the original sequence.

| | linear | circular |
| - | - | - |
| no sites | one fragment, the whole sequence | one fragment, still circular |
| one site | two fragments | one fragment of full length - a linearised plasmid |
| *n* sites | *n* + 1 fragments | *n* fragments |

On a circle, the fragment spanning the origin wraps: it is the one whose `end`
is less than its `start`.

## Performance

Searching is handed to the interpreter's native substring scan wherever it can
be. An exact pattern is one `binary.indexOf`; a degenerate one is expanded into
the concrete sequences it stands for - `CCWGG` is two searches, `GANTC` four -
and each of those is native too. Sixty-seven of the seventy-six recognition
sequences are exact, so a restriction map takes the fast path throughout.

Measured on a 200 kb sequence:

| | |
| - | - |
| `enzymesCutting` over the whole table | ~1.2 s |
| the same, before this was done | ~13.6 s |

What is left is **validation**, not searching: every public verb normalises its
input once, which is an interpreted pass costing roughly half a second per
200 kb. That is now the floor for the topic.

Two costs still scale with more than length. `restrictionSites` repeats a scan
per enzyme, so asking it for all 76 over a long sequence is 76 passes - prefer
`enzymesCutting` for the *which enzymes* question, which is what it is for. And
the returned sites are sorted by insertion, which is quadratic in the number of
sites; that only shows up in the thousands.

## The surface

| | |
| - | - |
| `find(s, pattern)` | exact, overlapping |
| `findIupac(s, pattern)` | ambiguity codes in the pattern |
| `findFuzzy(s, pattern, maxMismatches)` | substitutions only |
| `orfs(nt, minAa)` / `orfsWith(nt, minAa, tableId)` | → `list of Orf` |
| `cpgIslands(s, window, minGc, minObsExp)` | → `list of Region` |
| `restrictionSites(s, names)` | → `list of Site` |
| `enzymesCutting(s)` / `nonCutters(s)` | → `list of string` |
| `restrictionDigest(s, names, circular)` | → `list of Fragment` |
