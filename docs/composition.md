# Composition and properties

What a sequence is made of, and what that implies about the molecule.

This is the topic where the numbers come from published parameter sets rather
than from arithmetic on the sequence alone. Two sets can disagree about the same
sequence, so **every formula here names its source**, and the parameters are
exported constants you can read.

Runnable version: `jennifer run examples/composition.j`.

## Composition

```jennifer
bioseq.baseCount("AATTGGCC");       # {"A": 2, "T": 2, "G": 2, "C": 2}
bioseq.baseFrequency($gene);        # the same, as percentages
bioseq.residueCount($protein);      # the protein counterpart
bioseq.residueFrequency($protein);
```

A residue that does not occur has no entry, rather than an entry of zero - guard
a lookup with `maps.has`.

## GC content

```jennifer
bioseq.gcContent("ATGCGC");   # 66.67
bioseq.atContent("ATGCGC");   # 33.33
```

Ambiguity codes count toward the length but toward neither figure, which is what
keeps the two from summing above 100 on a sequence carrying `N`. On an
unambiguous sequence they partition it exactly.

### Along the sequence

A genome's average GC says little; its local variation says a great deal -
isochores, horizontally transferred regions and promoters all show up as
departures from the local baseline.

```jennifer
bioseq.gcWindows($genome, 100, 50);   # list of float, one per window
```

The last partial window is dropped, so every value covers exactly `size` bases
and the series is comparable along its length. [`windows`](manipulate.md#windows)
cuts the sequence the same way when you want the subsequences rather than their
GC.

## Molecular weight

```jennifer
bioseq.molecularWeight("ATGC", bioseq.DNA);       # 1173.83
bioseq.molecularWeight($rna, bioseq.RNA);
bioseq.molecularWeight($protein, bioseq.PROTEIN);
```

For `DNA` and `RNA` this is a **single strand with a free 5′ hydroxyl** - a
synthesised oligo, not a duplex and not a phosphorylated one. Residue masses are
summed, one water is added for the termini, and one metaphosphate is removed.
Double it, roughly, for a duplex. For `PROTEIN` it is the residue masses plus
one water, the ExPASy convention.

Ambiguity codes raise. See [Sequences](sequences.md#ambiguity-codes-are-rejected-not-averaged).

**Parameters:** `DNA_RESIDUE_MASS`, `RNA_RESIDUE_MASS`, `PHOSPHATE_MASS`,
`WATER_AVG_MASS`, and `avgMass` on each [`AminoAcid`](tables.md#residue-properties)
record.

## Melting temperature

Two models, and the difference between them is the point.

```jennifer
bioseq.meltingTemp("ATGCATGCATGCATGC");     # 48.0  - the Wallace rule
bioseq.meltingTempNn("ATGCATGCATGCATGC");   # the nearest-neighbour model
```

**The Wallace rule** is `2(A+T) + 4(G+C)`: crude, instant, and still the number
quoted on a primer tube. It is calibrated for oligos of roughly 14 to 20 bases
and is merely extrapolated beyond that - applied to a gene it returns a large
number with no physical meaning.

**The nearest-neighbour model** reads the bases in order, because stacking energy
depends on which base follows which. Enthalpy and entropy are summed over the
sequence's dinucleotides, corrected for helix initiation at each end, for salt,
and for self-complementarity.

The consequence, and the reason the model exists:

| oligo | Wallace | nearest-neighbour |
| - | - | - |
| `GCGCGC` | 24.0 | 22.8 |
| `GGGCCC` | 24.0 | 10.4 |

Identical composition. Wallace cannot tell them apart; the real model puts them
**12.4 °C** apart.

Conditions are parameters rather than constants, because they change the answer:

```jennifer
bioseq.meltingTempNn($primer);                       # 250 nM primer, 50 mM Na+
bioseq.meltingTempNnWith($primer, 250.0, 100.0);     # ~3.4 C higher
```

**Parameters:** SantaLucia (1998) unified set, exported as `NN_ENTHALPY` and
`NN_ENTROPY`, with `GAS_CONSTANT`. Other published sets differ by a degree or
two; this one is named so the difference is attributable.

## Protein chemistry

```jennifer
bioseq.netCharge($protein, 7.0);
bioseq.isoelectricPoint($protein);
bioseq.gravy($protein);
bioseq.aromaticity($protein);
```

**Net charge** sums the fractional protonation of every ionisable group: the two
termini and the side chains of D, E, C, Y (acidic) and H, K, R (basic). It falls
monotonically with pH.

**The isoelectric point** is where that curve crosses zero, found by bisection.
Because the curve is monotone the search always converges.

> The value depends on the pKa set. EMBOSS, ExPASy and Sillero disagree by a few
> tenths of a pH unit on the same sequence. This deck uses the **EMBOSS `iep`**
> set, exported as `PKA_N_TERMINUS`, `PKA_C_TERMINUS` and the `pKa` field on each
> residue record - so compare against `iep`, not against ProtParam.

**GRAVY** is the mean Kyte-Doolittle hydropathy. Positive means hydrophobic
overall; a soluble globular protein usually lands between -0.5 and -1.5, and a
membrane-spanning one runs positive. The insulin signal peptide in the test data
scores **+1.25**, which is what a sequence whose job is to cross a membrane
should look like.

**Aromaticity** is the fraction of F, W and Y.

## Summaries

```jennifer
def s as bioseq.DnaStats init bioseq.dnaStats($gene);
def p as bioseq.ProteinStats init bioseq.proteinStats($protein);
```

`DnaStats` carries `length`, `gc`, `at`, `counts`, `molecularWeight`,
`meltingTemp` and `ambiguous`. `ProteinStats` carries `length`,
`molecularWeight`, `isoelectricPoint`, `gravy`, `aromaticity`, `counts` and
`ambiguous`.

These are the convenience form, for a report or a table. When the sequence
carries a code with no defined mass or pKa, the summary still gives you length,
composition and GC - it sets `ambiguous` and reports the undefined figures as
`0.0` rather than raising or inventing them:

```jennifer
def s as bioseq.DnaStats init bioseq.dnaStats("ATGCNNNN");
$s.gc;               # 25.0
$s.ambiguous;        # true
$s.molecularWeight;  # 0.0 - not a guess
```

Note that `DnaStats.meltingTemp` is the Wallace figure, so read it together with
`length`; on anything gene-sized it is extrapolated well past the rule's range.

## Codon usage

```jennifer
for (def u in bioseq.codonUsage($cds)) {
    io.printf("%s %s %d %f|prec=2\n", $u.codon, $u.aminoAcid, $u.count, $u.fraction);
}
```

Each `CodonUse` gives `codon`, `aminoAcid`, `count`, `fraction` and
`perThousand`. `fraction` is the share of that residue's *synonymous* codons,
which is the figure a codon-usage table reports and the one that matters: a
highly expressed gene leans hard on a few preferred codons, and that bias is
what codon-optimisation software rewrites. In the insulin fragment, `CTG` takes
6 of the 8 leucines.

Each residue's fractions sum to 1.0, and the `perThousand` column sums to 1000.
`codonUsageWith` takes an alternative table.

## The surface

| | |
| - | - |
| `baseCount(s)` / `baseFrequency(s)` | nucleotide composition |
| `residueCount(p)` / `residueFrequency(p)` | protein composition |
| `gcContent(s)` / `atContent(s)` | percentages |
| `gcWindows(s, size, step)` | the GC-plot series |
| `molecularWeight(s, kind)` | daltons |
| `meltingTemp(s)` | Wallace rule |
| `meltingTempNn(s)` / `meltingTempNnWith(s, primerNm, saltMm)` | nearest-neighbour |
| `netCharge(p, ph)` / `isoelectricPoint(p)` | |
| `gravy(p)` / `aromaticity(p)` | |
| `dnaStats(s)` / `proteinStats(p)` | one-pass summaries |
| `codonUsage(nt)` / `codonUsageWith(nt, table)` | → `list of CodonUse` |
