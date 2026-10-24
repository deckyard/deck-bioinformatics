# Manipulation

Runnable: `jennifer run examples/manipulate.j`, which uses these verbs to build
null models.

## Reproducibility

The random verbs draw from `math.rand`; the deck keeps no generator of its own.
One `math.randSeed` makes the whole program reproducible:

```jennifer
math.randSeed(20260828);
bioseq.randomDna(40);   # the same sequence every run
```

## Generating

```jennifer
bioseq.randomDna(40);       # uniform over DNA_BASES, "ACGT"
bioseq.randomRna(40);       # RNA_BASES, "ACGU"
bioseq.randomProtein(40);   # STANDARD_RESIDUES, the twenty
```

Ambiguity codes are never generated. The length must be positive; there is no
upper bound, so `randomDna(2000000000)` tries to allocate two gigabytes.

### With a composition

```jennifer
bioseq.randomDnaWith(100, {"A": 35.0, "T": 35.0, "C": 15.0, "G": 15.0});
bioseq.randomDnaWith(40, {"A": 1.0, "G": 1.0});   # purines only
```

Weights are relative and normalised: `{"A": 3.0, "T": 1.0}` and
`{"A": 75.0, "T": 25.0}` are the same distribution. A residue left out (or
weighted 0) never appears. A key outside the generator's alphabet, a negative
weight, or all-zero weights raise. `randomRnaWith` and `randomProteinWith` work
the same way.

Real proteins are far from uniform. For a realistic null model, take the
composition of a real one:

```jennifer
def composition as map of string to float init bioseq.residueFrequency($real);
bioseq.randomProteinWith(len($real), $composition);
```

## Shuffling

```jennifer
bioseq.shuffle($sequence);
```

Fisher-Yates: every permutation equally likely, residue counts unchanged, so
`baseCount`, GC content, molecular weight and GRAVY are unchanged too. A
shuffle differs from the original only in order, which makes it the null model
for "is this motif surprising?":

```jennifer
def trials as int init 500;
def withSite as int init 0;
for (def i as int init 0; $i < $trials; $i = $i + 1) {
    if (len(bioseq.restrictionSites(bioseq.shuffle($construct), ["EcoRI"])) > 0) {
        $withSite = $withSite + 1;
    }
}
# examples/manipulate.j: 15 of 500 shuffles carried a site, about 3%
```

## Mutating

```jennifer
bioseq.mutate($gene, 0.01, bioseq.DNA);
```

Substitutions only, each position with probability `rate` (0 to 1), so the
number of changes varies between calls: at 0.1 a 100-base sequence gets about
ten. A mutated position always changes: the replacement is drawn from the
alphabet minus the current residue, so the realised rate matches `rate`. `kind`
selects that alphabet: `DNA_BASES`, `RNA_BASES` or `STANDARD_RESIDUES`.

## Extracting

### `subrange`

```jennifer
bioseq.subrange($gene, 13, 24);   # bases 13 to 24
$gene[12..24];                    # the same
```

The deck's one 1-based, inclusive verb, so positions quoted in a paper
transcribe directly. Out-of-range or backwards ranges raise.

### `windows`

```jennifer
bioseq.windows($gene, 24, 12);   # list of string
```

A window of `size` residues every `step`; they overlap when `step < size`. The
last partial window is dropped, as in
[`gcWindows`](composition.md#gc-content), so the two cut a sequence
identically. A sequence shorter than one window raises.
