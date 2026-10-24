# Genetic codes and residues

The reference data the deck reads, all exported. Runnable tour:
`jennifer run examples/tables.j`.

## Genetic codes

Nine NCBI genetic codes, each a `map of string to string` from DNA codon to
one-letter residue, stop as `*`:

| id | constant | NCBI name |
| - | - | - |
| 1 | `STANDARD` | Standard |
| 2 | `VERTEBRATE_MITOCHONDRIAL` | Vertebrate Mitochondrial |
| 3 | `YEAST_MITOCHONDRIAL` | Yeast Mitochondrial |
| 4 | `MOLD_MITOCHONDRIAL` | Mold, Protozoan and Coelenterate Mitochondrial; Mycoplasma/Spiroplasma |
| 5 | `INVERTEBRATE_MITOCHONDRIAL` | Invertebrate Mitochondrial |
| 6 | `CILIATE_NUCLEAR` | Ciliate, Dasycladacean and Hexamita Nuclear |
| 9 | `ECHINODERM_MITOCHONDRIAL` | Echinoderm and Flatworm Mitochondrial |
| 10 | `EUPLOTID_NUCLEAR` | Euplotid Nuclear |
| 11 | `BACTERIAL` | Bacterial, Archaeal and Plant Plastid |

Ids 7 and 8 were retired, so walk `codonTableIds()` rather than a range. An
unknown id raises.

```jennifer
bioseq.STANDARD["ATG"];             # "M"
bioseq.codonTable(2)["AGA"];        # "*"
bioseq.codonTableName(2);           # "Vertebrate Mitochondrial"
bioseq.codonTableIds();             # [1, 2, 3, 4, 5, 6, 9, 10, 11]
```

### Differences from the standard code

| code | reassigns |
| - | - |
| 2 Vertebrate mito. | `AGA`, `AGG` to `*`; `ATA` to `M`; `TGA` to `W` |
| 3 Yeast mito. | `ATA` to `M`; `CTN` to `T`; `TGA` to `W` |
| 4 Mold mito. | `TGA` to `W` |
| 5 Invertebrate mito. | `AGA`, `AGG` to `S`; `ATA` to `M`; `TGA` to `W` |
| 6 Ciliate nuclear | `TAA`, `TAG` to `Q` |
| 9 Echinoderm mito. | `AAA` to `N`; `AGA`, `AGG` to `S`; `TGA` to `W` |
| 10 Euplotid nuclear | `TGA` to `C` |
| 11 Bacterial | none; only the start codons differ |

The test suite restates these lists independently of the generator.

NCBI marks `CGA` and `CGC` as absent in code 3 rather than arginine. A map
cannot express "absent", so they are kept as `R`; a sequence that depends on the
distinction needs a hand check.

### Start codons

```jennifer
bioseq.startCodons(1);    # ["TTG", "CTG", "ATG"]
bioseq.startCodons(11);   # ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"]
```

This is why [`orfsWith`](search.md#open-reading-frames) takes a table id rather
than a codon map.

## IUPAC ambiguity codes

```jennifer
bioseq.expandIupac("R");            # "AG"
bioseq.iupacClass("R");             # "[AG]"
bioseq.matchesIupac("R", "G");      # true
bioseq.matchesIupac("A", "N");      # false - one-way
```

`IUPAC_DNA` and `IUPAC_RNA` map each code to the bases it admits; `A`, `C`, `G`,
`T`/`U` map to themselves. Lookups take any case; a non-IUPAC character raises.

| | | | |
| - | - | - | - |
| `R` = A G | `Y` = C T | `S` = C G | `W` = A T |
| `K` = G T | `M` = A C | `B` = C G T | `D` = A G T |
| `H` = A C T | `V` = A C G | `N` = A C G T | |

`iupacClass` renders a code as a regex character class for your own patterns;
[`findIupac`](search.md#patterns) does not use it.

`COMPLEMENT_DNA` and `COMPLEMENT_RNA` cover the whole alphabet. Codes complement
as sets (`R` to `Y`, `B` to `V`); `S`, `W` and `N` map to themselves.
Complementing twice returns the original.

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

`aminoAcid` takes a one- or three-letter code in any case (`"k"`, `"Lys"`).
`aminoAcids()` returns 22 records: the twenty standard residues in
three-letter alphabetical order (`ARNDCEQGHILKMFPSTWYV`), then selenocysteine
(`U`) and pyrrolysine (`O`). Masses are residue masses (free amino acid less one
water). `pKa` is the EMBOSS `iep` side-chain value, 0.0 for a side chain that
does not ionise.

`B`, `Z`, `X`, `J` and `*` have no record, since they have no single mass or
pKa; `aminoAcid` raises on them and `hasAminoAcid` is the guard.

Asserted in the test suite:

- Leucine and isoleucine have identical mass (mass spectrometry reports
  "Leu/Ile").
- Glycine is the lightest residue.
- Monoisotopic mass is below average mass for every residue except
  selenocysteine, because selenium's most abundant isotope is heavier than its
  average.
- A side chain has a charge exactly when it has a pKa: H, K, R basic; C, D, E, Y
  acidic.

`THREE_LETTER` maps the whole protein alphabet, ambiguity codes (`Asx`, `Glx`,
`Xaa`, `Xle`) and stop (`Ter`) included, to three-letter codes;
[`oneToThree`](transforms.md#residue-spellings) reads it.

## Parameter constants

| | |
| - | - |
| `PKA_N_TERMINUS`, `PKA_C_TERMINUS` | 8.6, 3.6 (EMBOSS `iep`) |
| `WATER_MONO_MASS`, `WATER_AVG_MASS` | 18.010565, 18.01528 Da, added once per chain |
| `DNA_RESIDUE_MASS`, `RNA_RESIDUE_MASS` | average residue mass per base (nucleoside monophosphate less water) |
| `PHOSPHATE_MASS` | 79.98 Da (HPO3), subtracted for a 5'-OH oligo |
| `NN_ENTHALPY`, `NN_ENTROPY` | SantaLucia (1998) nearest-neighbour dH (kcal/mol) and dS (cal/(mol K)), all 16 dinucleotides |
| `GAS_CONSTANT` | 1.987 cal/(mol K) |
| `BLOSUM62` | [substitution matrix](compare.md#scoring), Henikoff and Henikoff (1992) |

## How the data is maintained

`src/topics/tabledata.j` is generated by `tools/gentables.j`. The generator
writes each alternative code as the standard code plus its documented
reassignments, and before emitting asserts that BLOSUM62 is symmetric across
all 576 pairs, that ten known matrix values are right, and that every residue's
identity score beats all nineteen of its substitutions. See
[Design notes](design.md#generated-data-hand-written-api).
