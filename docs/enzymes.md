# Restriction enzymes

76 hand-curated type II enzymes. Finding and cutting sites is in
[Search](search.md#restriction-sites).

## Lookup

```jennifer
def e as bioseq.Enzyme init bioseq.enzyme("EcoRI");
$e.site;             # "GAATTC"
$e.cutTop;           # 1
$e.cutBottom;        # 5
$e.overhang;         # "5'"
$e.overhangLength;   # 4
```

Names are case-insensitive (`enzyme("ecori")`); the record carries the
conventional spelling. An unknown name raises; `hasEnzyme` is the guard.
`enzymes()` and `enzymeNames()` return the whole table in alphabetical order.

## Cut geometry

Both offsets count along the top strand from the first base of the site:

```
EcoRI  G^AATTC     site "GAATTC", cutTop 1, cutBottom 5   ->  5' overhang, 4 nt
PstI   CTGCA^G     site "CTGCAG", cutTop 5, cutBottom 1   ->  3' overhang, 4 nt
SmaI   CCC^GGG     site "CCCGGG", cutTop 3, cutBottom 3   ->  blunt
```

`overhang` (`"5'"`, `"3'"` or `"blunt"`) and `overhangLength` follow from the
two offsets; the tests check they agree. `isBlunt(e)` is
`cutTop == cutBottom`. [`restrictionDigest`](search.md#digestion) cuts at
`cutTop`.

## Degenerate sites

Sites are stored in IUPAC codes as specified:

```jennifer
bioseq.enzyme("EcoRII").site;   # "CCWGG"
bioseq.enzyme("HinfI").site;    # "GANTC"
bioseq.enzyme("XmnI").site;     # "GAANNNNTTC"
```

The search reports the concrete match. 67 of the 76 sites are exact.

## Every site is palindromic

Each site equals its own reverse complement, degenerate ones included (`N`
complements to `N`, `W` to `W`):

```
GANTC  ->  reverse complement  ->  GANTC
```

There is no exception list: the generator and the test suite both assert it for
all 76 entries.

## Isoschizomers

```jennifer
bioseq.isoschizomers("SmaI");   # ["XmaI"]
```

Other enzymes with the same site string, alphabetical. True isoschizomers cut
identically (`HpaII` / `MspI`, `Sau3AI` / `MboI`); neoschizomers cut the same
site at a different position:

| | site | cut | ends |
| - | - | - | - |
| `SmaI` | `CCCGGG` | `CCC^GGG` | blunt |
| `XmaI` | `CCCGGG` | `C^CCGGG` | 5', 4 nt |

Compare `cutTop` to tell them apart.

## Eight-base cutters

`AscI`, `FseI`, `NotI`, `PacI`, `PmeI`, `SbfI`, `SrfI`, `SwaI`. They rarely cut,
which makes them useful for linearising a construct once.

## Not included

Type IIS enzymes (`BsaI`, `BsmBI`, ...) cut outside their recognition site,
which two in-site offsets cannot describe. They are absent:
`hasEnzyme("BsaI")` is false.

## Provenance

Sites and cut positions are hand-curated from published enzyme specifications,
not imported from REBASE, whose data carries its own terms. A new enzyme is a
four-field row in `tools/genenzymes.j`. The generator asserts that names are
unique, every site is in the DNA IUPAC alphabet and palindromic, and both cut
offsets fall inside the site; the test suite re-derives the overhang labels and
re-checks the palindromes.
