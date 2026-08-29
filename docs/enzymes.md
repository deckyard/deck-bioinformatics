# Restriction enzymes

76 hand-curated type II enzymes, their recognition sequences, and where they
cut.

For finding and cutting, see [Search](search.md#restriction-sites). This page is
the catalogue.

## Looking one up

```jennifer
def e as bioseq.Enzyme init bioseq.enzyme("EcoRI");
$e.site;             # "GAATTC"
$e.cutTop;           # 1
$e.cutBottom;        # 5
$e.overhang;         # "5'"
$e.overhangLength;   # 4
```

Lookup is case-insensitive - `enzyme("ecori")` works - but the record carries the
conventional spelling. `hasEnzyme` is the guard, `enzymes()` and `enzymeNames()`
give the whole table in alphabetical order.

## Cut geometry

Both offsets are counted **along the top strand**, from the first base of the
recognition site:

```
EcoRI  G^AATTC     site "GAATTC", cutTop 1, cutBottom 5   →  5′ overhang, 4 nt
PstI   CTGCA^G     site "CTGCAG", cutTop 5, cutBottom 1   →  3′ overhang, 4 nt
SmaI   CCC^GGG     site "CCCGGG", cutTop 3, cutBottom 3   →  blunt
```

`overhang` is `"5'"`, `"3'"` or `"blunt"`, and is derived from the two offsets
rather than stored independently - the tests check that it agrees.

`cutTop` is what determines fragment sizes on a gel, which is what
[`restrictionDigest`](search.md#digestion) cuts on. `cutBottom` decides the
ends' chemistry: whether a fragment will ligate to another, and to which.

## Degenerate sites

A recognition sequence is stored in IUPAC codes, exactly as it is specified,
rather than expanded into several rows:

```jennifer
bioseq.enzyme("EcoRII").site;   # "CCWGG" - W is A or T
bioseq.enzyme("HinfI").site;    # "GANTC"
bioseq.enzyme("XmnI").site;     # "GAANNNNTTC"
```

The search treats them as patterns and reports the concrete match it found.

## Every site is palindromic

A type II enzyme binds as a homodimer, so its recognition site reads the same on
both strands. This holds for the degenerate sites too - `N` complements to `N`
and `W` to `W`, so the specificity in the outer bases pairs up:

```
GANTC  →  reverse complement  →  GANTC
```

There is no exception list. A site that fails this check is a typo, and both the
generator and the test suite assert it for all 76 entries.

## Isoschizomers

```jennifer
bioseq.isoschizomers("SmaI");   # ["XmaI"]
```

Enzymes sharing a recognition sequence. **True isoschizomers** recognise and cut
identically (`HpaII` / `MspI`, `Sau3AI` / `MboI`); **neoschizomers** recognise
the same site but cut it at different positions, which changes the ends you get:

| | site | cut | ends |
| - | - | - | - |
| `SmaI` | `CCCGGG` | `CCC^GGG` | blunt |
| `XmaI` | `CCCGGG` | `C^CCGGG` | 5′, 4 nt |

Both are returned; read `cutTop` to tell them apart.

## Eight-cutters

The eight-base recognition sites are the ones that rarely cut at all, which is
what makes them useful for linearising a construct exactly once:

`AscI`, `FseI`, `NotI`, `PacI`, `PmeI`, `SbfI`, `SrfI`, `SwaI`

## What is not here

**Type IIS enzymes** - `BsaI`, `BsmBI` and the rest, which cut at a distance
*outside* their recognition site - are deliberately absent. Their geometry does
not fit two offsets measured inside the site, and modelling them properly is a
separate job. `hasEnzyme("BsaI")` is false rather than quietly wrong.

## Provenance

Recognition sequences and cut positions are **hand-curated from published enzyme
specifications**, not imported in bulk from REBASE, whose data carries its own
usage terms. The table is meant to be extended by whoever needs an enzyme it
lacks rather than to be exhaustive on day one: a new entry is a four-field row in
`tools/genenzymes.j`, and the generator's assertions catch a malformed one.

Those assertions check that every site is spelled in the DNA IUPAC alphabet, that
both cut offsets fall inside the site, and that the site is palindromic. The test
suite independently re-derives the overhang label from the offsets and re-checks
the palindromes, so the data is verified on both sides of the generator.

## The surface

| | |
| - | - |
| `Enzyme` | `name`, `site`, `cutTop`, `cutBottom`, `overhang`, `overhangLength` |
| `enzyme(name)` | case-insensitive lookup |
| `hasEnzyme(name)` | the guard |
| `enzymes()` / `enzymeNames()` | the whole table |
| `isoschizomers(name)` | others sharing the site |
| `isBlunt(e)` | both strands cut at the same offset |
