# Composition and properties

Runnable: `jennifer run examples/composition.j`. The parameter sets behind these
numbers are exported; see [Parameter constants](tables.md#parameter-constants).

## Composition

```jennifer
bioseq.baseCount("AATTGGCC");       # {"A": 2, "T": 2, "G": 2, "C": 2}
bioseq.baseFrequency($gene);        # the same, as percentages
bioseq.residueCount($protein);
bioseq.residueFrequency($protein);
```

Only residues present get an entry; guard a lookup with `maps.has`.

## GC content

```jennifer
bioseq.gcContent("ATGCGC");   # 66.67
bioseq.atContent("ATGCGC");   # 33.33
bioseq.gcWindows($genome, 100, 50);   # list of float, one per window
```

Ambiguity codes count toward the length but toward neither G/C nor A/T/U, so on
a sequence with `N` (or `S`) the two sum to less than 100.

`gcWindows(s, size, step)` starts a window every `step` bases and drops the last
partial window, so every value covers exactly `size` bases;
[`windows`](manipulate.md#windows) cuts identically. A sequence shorter than
one window raises.

## Molecular weight

```jennifer
bioseq.molecularWeight("ATGC", bioseq.DNA);       # 1173.83
bioseq.molecularWeight($rna, bioseq.RNA);
bioseq.molecularWeight($protein, bioseq.PROTEIN);
```

Average masses, in daltons.

- `DNA`, `RNA`: one strand with a free 5' hydroxyl (an unphosphorylated
  synthetic oligo): residue masses + `WATER_AVG_MASS` - `PHOSPHATE_MASS`. Roughly
  double it for a duplex.
- `PROTEIN`: residue masses (`avgMass` of each [`AminoAcid`](tables.md#residue-properties))
  + one water, the ExPASy convention.

Ambiguity codes raise ([why](sequences.md#ambiguity-codes-are-rejected-not-averaged)).

## Melting temperature

```jennifer
bioseq.meltingTemp("ATGCATGCATGCATGC");                     # 48.0, Wallace rule
bioseq.meltingTempNn("ATGCATGCATGCATGC");                   # 50.7, 250 nM primer, 50 mM Na+
bioseq.meltingTempNnWith("ATGCATGCATGCATGC", 250.0, 100.0); # 53.9, 100 mM Na+
```

Both take DNA only, in degrees Celsius, and raise on ambiguity codes.

**Wallace rule**: `2(A+T) + 4(G+C)`. Calibrated for oligos of about 14 to 20
bases; on longer sequences it returns a number with no physical meaning.

**Nearest-neighbour**: SantaLucia (1998) unified parameters (`NN_ENTHALPY`,
`NN_ENTROPY`, `GAS_CONSTANT`). Enthalpy and entropy are summed over the
dinucleotides, with helix initiation per terminal base, a salt correction and a
self-complementarity correction, then `Tm = dH / (dS + R ln(Ct/4)) - 273.15`.
Needs at least two bases; concentrations must be positive (`primerNm` in nM,
`saltMm` in mM). Other parameter sets differ by a degree or two.

Order matters to the stacking model, not to Wallace:

| oligo | Wallace | nearest-neighbour |
| - | - | - |
| `GCGCGC` | 24.0 | 22.8 |
| `GGGCCC` | 24.0 | 10.4 |

## Protein chemistry

```jennifer
bioseq.netCharge($protein, 7.0);
bioseq.isoelectricPoint($protein);
bioseq.gravy($protein);
bioseq.aromaticity($protein);
```

All four raise on an empty or non-protein sequence; the first three also on
`B`, `Z`, `X`, `J` or `*`.

**Net charge** sums the fractional charge of the two termini and the side chains
of D, E, C, Y (acidic) and H, K, R (basic). It decreases monotonically with pH.

**Isoelectric point**: bisection on the net charge over pH 0 to 14, stopping
within 1e-4 of zero or after 100 steps. The pKa values are the EMBOSS `iep` set
(`PKA_N_TERMINUS`, `PKA_C_TERMINUS`, `AminoAcid.pKa`); EMBOSS, ExPASy and
Sillero differ by a few tenths of a unit, so compare against `iep`, not
ProtParam.

**GRAVY**: mean Kyte-Doolittle hydropathy. Positive is hydrophobic; soluble
globular proteins typically score -0.5 to -1.5. The preproinsulin signal
fragment in `testdata/` scores +1.25.

**Aromaticity**: fraction of F, W and Y, 0 to 1.

## Summaries

```jennifer
def s as bioseq.DnaStats init bioseq.dnaStats("ATGCNNNN");
$s.gc;               # 25.0
$s.ambiguous;        # true
$s.molecularWeight;  # 0.0
def p as bioseq.ProteinStats init bioseq.proteinStats($protein);
```

One call per report row. With ambiguity codes the summary still has length,
counts and GC/AT (DNA) or aromaticity (protein); mass, Tm, pI and GRAVY are
`0.0` and `ambiguous` is true. `dnaStats` rejects RNA. `DnaStats.meltingTemp`
is the Wallace figure. Fields: [cheatsheet](cheatsheet.md#dnastats).

## Codon usage

```jennifer
for (def u in bioseq.codonUsage($cds)) {
    io.printf("%s %s %d %f|prec=2\n", $u.codon, $u.aminoAcid, $u.count, $u.fraction);
}
```

One `CodonUse` per codon present, in first-seen order. `fraction` is the share
among codons for the same residue (each residue's fractions sum to 1.0);
`perThousand` sums to 1000. In the insulin fragment `CTG` is 6 of the 8
leucines (0.75). Codons are reported DNA-spelled; a codon the table does not
hold (one with an ambiguity code) counts under residue `X`. `codonUsageWith`
takes a codon table.
