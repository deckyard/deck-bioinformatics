# Search

Runnable: `jennifer run examples/search.j`.

## Patterns

```jennifer
bioseq.find($dna, "GGATCC");            # exact
bioseq.findIupac($dna, "CCWGG");        # IUPAC codes in the pattern
bioseq.findFuzzy($dna, "GAATTC", 1);    # up to 1 substitution
```

Each returns `list of Match`, in ascending order of `start`, with
`$dna[$m.start..$m.end] == $m.text`. Only the given strand is searched;
reverse-complement the pattern for the other. `find` and `findFuzzy` compare
residues literally and work on protein too; `findIupac` takes nucleic acids.

Matches overlap:

```jennifer
len(bioseq.find("AAAA", "AA"));   # 3
```

That is why `findIupac` scans the sequence itself instead of using `regex`.
[`iupacClass`](tables.md#iupac-ambiguity-codes) is there for callers composing
their own patterns.

IUPAC matching is one-way: a pattern's `N` admits a subject's `A`, a pattern's
`A` does not admit a subject's `N`. `text` is the concrete subsequence found, so
a `CCWGG` hit says whether it was `CCAGG` or `CCTGG`.

```jennifer
len(bioseq.findIupac("AAAA", "N"));   # 4
len(bioseq.findIupac("NNNN", "A"));   # 0
```

`findFuzzy` counts substitutions only; each match carries `mismatches`. Use
[`editDistance`](compare.md#distances) when indels matter. A pattern longer
than the sequence, or a negative tolerance, raises.

## Open reading frames

```jennifer
for (def o in bioseq.orfs($dna, 50)) {
    io.printf("frame %d  %s\n", $o.frame, $o.protein);
}
bioseq.orfsWith($dna, 50, 11);   # NCBI table 11, bacterial starts
```

All six frames are scanned, results in frame order `1, 2, 3, -1, -2, -3`.
Within a frame an ORF starts at the first start codon and runs to the first
in-frame stop; the next starts at the next start codon after that stop. ORFs do
not overlap within a frame and nested starts are not reported separately.
`minAa` is the shortest protein reported, in residues, stop excluded.

`orfsWith` takes an NCBI table id, not a codon map, because it needs the table's
start codons. `BACTERIAL` (11) assigns codons as `STANDARD` (1) does but has
seven starts instead of three. `CTG` and `TTG` are starts under table 1 too.

`stopped` is false when a frame reaches the end of the sequence without a stop.
Such ORFs are reported; on a fragment the gene usually continues:

```jennifer
def o as bioseq.Orf init bioseq.orfs($insulinFragment, 20)[0];
$o.stopped;   # false
$o.protein;   # MALWMRLLPLLALLALWGPDPAAA
```

### Reverse-strand coordinates

`start` and `end` are offsets into the sequence you passed in, for both strands.
`nucleotides` is the frame as read (DNA-spelled, `U` read as `T`):

```jennifer
bioseq.reverseComplement($dna[$o.start..$o.end]) == $o.nucleotides;   # strand -1
$dna[$o.start..$o.end] == $o.nucleotides;                            # strand +1
```

The test suite asserts this against the sequence itself.

## CpG islands

```jennifer
bioseq.cpgIslands($promoter, 200, 50.0, 0.6);   # Gardiner-Garden and Frommer (1987)
```

A window qualifies when its GC percentage is at least `minGc` and its observed
CpG count over the expected `(C count * G count) / length` is at least
`minObsExp`. The obs/exp test separates CpG-rich from merely GC-rich:

```jennifer
def gcOnly as string init strings.repeat("G", 60) + strings.repeat("C", 60);
bioseq.gcContent($gcOnly);                             # 100.0
len(bioseq.cpgIslands($gcOnly, 20, 50.0, 0.6));        # 0
```

Qualifying windows that touch or overlap are merged into one `Region`, whose
`gc` and `obsExp` are recomputed over the merged span. A window straddling an
island edge can still qualify, so reported spans run a little wider than the
dense core. A sequence shorter than `window` gives an empty list.

## Restriction sites

```jennifer
bioseq.restrictionSites($plasmid, ["EcoRI", "BamHI", "HindIII"]);
bioseq.restrictionSites("AAACCAGGTTTCCTGGAAA", ["EcoRII"]);
# two sites: "CCAGG" at 3 and "CCTGG" at 11
```

Sites are ordered by `cut`, then enzyme name. `site` is the concrete match for a
degenerate recognition sequence. Names are case-insensitive; an unknown name
raises. Every site in the table is palindromic, so one forward scan finds them
all (a non-palindromic site would be searched on both strands and reported with
`strand: -1`). The sequence is scanned as linear: a site spanning the origin of
a circular sequence is not found.

```jennifer
bioseq.enzymesCutting($insert);   # enzymes with at least one site, alphabetical
bioseq.nonCutters($insert);       # enzymes with none
```

`nonCutters` is the list to pick from when a construct must be cut exactly once
or an insert left intact.

## Digestion

```jennifer
bioseq.restrictionDigest($plasmid, ["EcoRI", "BamHI"], false);   # linear
bioseq.restrictionDigest($plasmid, ["EcoRI", "BamHI"], true);    # circular
```

Cuts are on the top strand at each enzyme's `cutTop`, which sets the fragment
lengths; the bottom-strand stagger is described by
[`Enzyme.overhang`](enzymes.md#cut-geometry). Coinciding cuts count once; a cut
at offset 0 or at the end is ignored. Fragments are ordered by `start` and,
concatenated in order, give back the original sequence.

| | linear | circular |
| - | - | - |
| no sites | one fragment, the whole sequence | one fragment, `0` to length |
| one site | two fragments | one full-length fragment starting at the cut, `end == start` |
| *n* sites | *n* + 1 fragments | *n* fragments |

On a circle the last fragment runs through the origin: its `end` is not greater
than its `start`.

## Performance

An exact pattern is one native `binary.indexOf` scan. A degenerate one is
expanded into its concrete sequences (`CCWGG` is two, `GANTC` four), each
scanned natively; only one expanding past 1024 variants falls back to an
interpreted scan of about 3 microseconds per position. 67 of the 76 recognition
sequences are exact.

| on a 200 kb sequence | |
| - | - |
| `enzymesCutting` over the whole table | ~1.2 s |
| the same, before native scanning | ~13.6 s |

The rest is validation: every public verb normalises its input once, about half
a second per 200 kb, which is the floor for this topic. `restrictionSites`
scans once per enzyme, so for "which enzymes cut" use `enzymesCutting`, which
stops at each enzyme's first site. Sites are merge-sorted, O(n log n) in the
number of sites.
