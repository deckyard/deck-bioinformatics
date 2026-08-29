# Genetic codes and residues

The reference data the rest of the deck reads. All of it is exported, so a
caller who disagrees with a parameter can see exactly what was used.

Runnable tour: `jennifer run examples/tables.j`.

## Genetic codes

Nine NCBI genetic codes, each a `map of string to string` from codon to
one-letter residue, with a stop as `*`:

| id | constant | name |
| - | - | - |
| 1 | `STANDARD` | Standard |
| 2 | `VERTEBRATE_MITOCHONDRIAL` | Vertebrate Mitochondrial |
| 3 | `YEAST_MITOCHONDRIAL` | Yeast Mitochondrial |
| 4 | `MOLD_MITOCHONDRIAL` | Mold, Protozoan and Coelenterate Mitochondrial |
| 5 | `INVERTEBRATE_MITOCHONDRIAL` | Invertebrate Mitochondrial |
| 6 | `CILIATE_NUCLEAR` | Ciliate, Dasycladacean and Hexamita Nuclear |
| 9 | `ECHINODERM_MITOCHONDRIAL` | Echinoderm and Flatworm Mitochondrial |
| 10 | `EUPLOTID_NUCLEAR` | Euplotid Nuclear |
| 11 | `BACTERIAL` | Bacterial, Archaeal and Plant Plastid |

The numbering has gaps - 7 and 8 were retired - so walk `codonTableIds()` rather
than a range.

```jennifer
bioseq.STANDARD["ATG"];             # "M"
bioseq.codonTable(2)["AGA"];        # "*"
bioseq.codonTableName(2);           # "Vertebrate Mitochondrial"
bioseq.codonTableIds();             # [1, 2, 3, 4, 5, 6, 9, 10, 11]
```

### What actually differs

Each alternative code reassigns a handful of codons. That short list is the
reviewable form of the data, and the test suite states it independently and
fails if a table stops matching:

| code | reassigns |
| - | - |
| Vertebrate mito. | `AGA`→`*`, `AGG`→`*`, `ATA`→`M`, `TGA`→`W` |
| Yeast mito. | `ATA`→`M`, the whole `CTN` block→`T`, `TGA`→`W` |
| Mold mito. | `TGA`→`W` |
| Invertebrate mito. | `AGA`→`S`, `AGG`→`S`, `ATA`→`M`, `TGA`→`W` |
| Ciliate nuclear | `TAA`→`Q`, `TAG`→`Q` |
| Echinoderm mito. | `AAA`→`N`, `AGA`→`S`, `AGG`→`S`, `TGA`→`W` |
| Euplotid nuclear | `TGA`→`C` |
| Bacterial | nothing - see below |

The ciliate code is the famous one: two of the three stops read glutamine, so a
naive translation of a *Paramecium* gene stops far too early.

> NCBI additionally marks `CGA` and `CGC` as *absent* in the yeast mitochondrial
> code rather than as arginine. A map cannot express "absent", so they are kept
> as `R`; a sequence relying on that distinction needs a hand check.

### Start codons

```jennifer
bioseq.startCodons(1);    # TTG, CTG, ATG
bioseq.startCodons(11);   # TTG, CTG, ATT, ATC, ATA, ATG, GTG
```

This is what distinguishes `BACTERIAL` from `STANDARD`: identical codon
assignments, seven start codons rather than three. It is also why
[`orfsWith`](search.md#start-codons-follow-the-table) takes a table id rather
than a codon map.

Worth noting that `TTG` and `CTG` initiate under the standard code too - `ATG`
is not the only start codon, even in a eukaryotic nucleus.

## IUPAC ambiguity codes

```jennifer
bioseq.expandIupac("R");            # "AG"
bioseq.iupacClass("R");             # "[AG]"
bioseq.matchesIupac("R", "G");      # true
bioseq.matchesIupac("A", "N");      # false - the relation is one-way
```

`IUPAC_DNA` and `IUPAC_RNA` map each code to the bases it admits. The four
unambiguous bases map to themselves, so a lookup never needs a special case.

| | | | |
| - | - | - | - |
| `R` = A G | `Y` = C T | `S` = C G | `W` = A T |
| `K` = G T | `M` = A C | `B` = C G T | `D` = A G T |
| `H` = A C T | `V` = A C G | `N` = A C G T | |

`iupacClass` renders a code as a regular-expression character class, for callers
composing their own patterns. The deck's own
[`findIupac`](search.md#patterns) does not use it - it needs overlapping matches,
which a regex scan does not give.

## Complements

`COMPLEMENT_DNA` and `COMPLEMENT_RNA` cover the whole alphabet, ambiguity codes
included. They complement as sets - `R` (A or G) becomes `Y` (C or T), `B`
(not A) becomes `V` (not T) - and `S`, `W` and `N` are self-complementary.

Each is an involution: complementing twice returns the original.

## Residue properties

```jennifer
def k as bioseq.AminoAcid init bioseq.aminoAcid("K");
$k.three;        # "Lys"
$k.name;         # "Lysine"
$k.monoMass;     # 128.09496
$k.avgMass;      # 128.1741
$k.pKa;          # 10.8
$k.charge;       # 1
$k.hydropathy;   # -3.9
```

`aminoAcid` takes either spelling, in any case - `aminoAcid("Lys")` and
`aminoAcid("k")` are the same lookup. `aminoAcids()` returns all 22 records: the
standard twenty, plus selenocysteine and pyrrolysine.

**The ambiguity codes `B`, `Z`, `X` and `J` are deliberately absent.** They have
no single mass or pKa, so a caller who means one of them has to decide what it
means rather than receive a plausible wrong number. `hasAminoAcid` is the cheap
guard.

A few properties worth knowing, all asserted in the test suite:

- Leucine and isoleucine are **isobaric** - identical mass - which is why mass
  spectrometry reports "Leu/Ile".
- Glycine is the lightest residue.
- Monoisotopic mass is lighter than average mass for every residue **except
  selenocysteine**: selenium's most abundant isotope is heavier than the
  element's average, so the ordering inverts.
- A side chain ionises exactly when it carries a charge: three basic residues
  (H, K, R), four acidic (C, D, E, Y).

Three-letter spellings for the whole protein alphabet, ambiguity codes and the
stop included, live in `THREE_LETTER` and are what
[`oneToThree`](transforms.md#residue-spellings) reads.

## Parameter constants

| | |
| - | - |
| `PKA_N_TERMINUS`, `PKA_C_TERMINUS` | 8.6 and 3.6 - the EMBOSS `iep` set |
| `WATER_MONO_MASS`, `WATER_AVG_MASS` | added once per chain |
| `DNA_RESIDUE_MASS`, `RNA_RESIDUE_MASS` | per-base residue masses |
| `PHOSPHATE_MASS` | subtracted for a 5′-OH oligo |
| `NN_ENTHALPY`, `NN_ENTROPY` | SantaLucia (1998) nearest-neighbour parameters |
| `GAS_CONSTANT` | 1.987 cal/(mol·K) |
| `BLOSUM62` | the [substitution matrix](compare.md#why-protein-needs-a-matrix) |

## How this data is maintained

`src/topics/tabledata.j` is **generated** by `tools/gentables.j`, and the
generator is where the review happens. It expresses each alternative genetic code
as the standard code plus its documented differences - a four-line diff rather
than a 64-line blob - and asserts the anchors before emitting anything: that
BLOSUM62 is symmetric across all 576 pairs, that ten well-known matrix values are
right, and that every residue's identity score beats all nineteen of its
substitutions.

The tests then verify from the other side, restating the difference lists
independently. See [Design notes](design.md#generated-data-hand-written-api).
