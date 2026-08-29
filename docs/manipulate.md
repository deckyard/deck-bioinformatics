# Manipulation

Generating, shuffling, mutating and slicing.

Runnable version: `jennifer run examples/manipulate.j`.

## What these are for

Null models. A motif that looks striking in a real sequence has to be compared
against something, and the honest comparison is sequence with the same
composition and no structure. That is what `shuffle` produces and what
`randomDnaWith` approximates.

## Reproducibility

The random verbs draw from `math.rand`, so a single `math.randSeed` governs the
whole run:

```jennifer
math.randSeed(20260828);
bioseq.randomDna(40);   # the same sequence every time
```

The deck keeps no generator of its own, deliberately: one seed for the program is
easier to reason about than a hidden second stream.

## Generating

```jennifer
bioseq.randomDna(40);
bioseq.randomRna(40);
bioseq.randomProtein(40);
```

Uniform over the alphabet. The alphabets are `DNA_BASES` (`ACGT`), `RNA_BASES`
(`ACGU`) and `STANDARD_RESIDUES` (the twenty). Ambiguity codes are excluded: a
generator emitting `N` would be producing missing data, not sequence.

### With a composition

```jennifer
bioseq.randomDnaWith(100, {"A": 35.0, "T": 35.0, "C": 15.0, "G": 15.0});
```

Weights are **relative, not percentages** - they are normalised - so
`{"A": 3.0, "T": 1.0}` and `{"A": 75.0, "T": 25.0}` describe the same
distribution. A residue left out of the map never appears, which is how to
generate a purine-only sequence:

```jennifer
bioseq.randomDnaWith(40, {"A": 1.0, "G": 1.0});
```

A uniform random protein is not what a real one looks like - residue frequencies
in nature are far from flat. Pair `randomProteinWith` with
[`residueFrequency`](composition.md#composition) on a real protein when the null
model needs to be realistic:

```jennifer
def composition as map of string to float init bioseq.residueFrequency($real);
bioseq.randomProteinWith(len($real), $composition);
```

## Shuffling

```jennifer
bioseq.shuffle($sequence);
```

A Fisher-Yates shuffle: every permutation equally likely, and the residue counts
unchanged by construction. `baseCount(shuffle(s))` always equals
`baseCount(s)`, and so does anything derived from composition alone - GC
content, molecular weight, GRAVY.

This is the honest null model for *is this motif surprising?* A random sequence
of the same length answers a weaker question, because it differs from the real
one in composition as well as order. A shuffle differs only in order, which is
usually the hypothesis being tested:

```jennifer
def trials as int init 500;
def withSite as int init 0;
for (def i as int init 0; $i < $trials; $i = $i + 1) {
    if (len(bioseq.restrictionSites(bioseq.shuffle($construct), ["EcoRI"])) > 0) {
        $withSite = $withSite + 1;
    }
}
# of 500 shuffles with identical composition, 15 carried a site:
# a site like this arises by chance about 3% of the time
```

## Mutating

```jennifer
bioseq.mutate($gene, 0.01, bioseq.DNA);
```

Substitutions only - no insertions or deletions - at a per-position probability.
The number of changes varies between calls: at rate 0.1 a hundred-base sequence
takes about ten substitutions, not exactly ten.

**A mutated position always changes.** If a substitution could redraw the
original residue, the realised rate would come out lower than the one asked for,
so the draw is taken from the alphabet less the current residue.

The `kind` argument decides which alphabet substitutions come from: the four
bases for `DNA` and `RNA`, the twenty standard residues for `PROTEIN`.

## Extracting

### `subrange`

```jennifer
bioseq.subrange($gene, 13, 24);   # "bases 13 to 24"
```

**This is the deck's one 1-based inclusive verb**, and it is named for it.
Everything else is 0-based and half-open, matching Jennifer's own slicing - but
biology quotes positions 1-based inclusive everywhere, so a paper saying
"residues 25 to 40" transcribes directly. The equivalent native slice is
`$gene[12..24]`.

### `windows`

```jennifer
bioseq.windows($gene, 24, 12);   # list of string
```

Overlapping windows, advancing by `step`. The last partial window is dropped, so
every window is exactly `size` long and the series stays comparable - the same
rule [`gcWindows`](composition.md#along-the-sequence) follows, and the two cut a
sequence identically.

## The surface

| | |
| - | - |
| `randomDna(n)` / `randomRna(n)` / `randomProtein(n)` | uniform |
| `randomDnaWith(n, weights)` / `randomRnaWith` / `randomProteinWith` | weighted |
| `shuffle(s)` | composition-preserving |
| `mutate(s, rate, kind)` | point substitutions |
| `subrange(s, from, through)` | **1-based inclusive** |
| `windows(s, size, step)` | → `list of string` |
| `DNA_BASES`, `RNA_BASES`, `STANDARD_RESIDUES` | the generator alphabets |
