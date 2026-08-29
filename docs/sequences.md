# Sequences and alphabets

Read this page first. It describes the five conventions every other topic
follows, and most surprises with the deck come from not knowing one of them.

## A sequence is a string

```jennifer
bioseq.reverseComplement("ATGC");   # "GCAT"
```

No constructor, no wrapper type. A sequence read out of a FASTA record, typed
into a literal, or built by concatenation is the same thing to every verb here.
The cost of that choice is that the deck cannot tell DNA from protein by looking
at the type - so it looks at the residues instead, and says so when it cannot
tell.

## Case in, uppercase out

Input is accepted in any case. Residues come back uppercase.

```jennifer
bioseq.complement("atgc");   # "TACG"
```

The one deliberate exception is [FASTA parsing](formats.md#case-is-preserved),
where lowercase carries meaning - it marks soft-masked repeats - and is
preserved.

## Two tiers of tolerance

The deck is strict about what a sequence is, but it distinguishes *junk you
meant to remove* from *junk that means you have the wrong data*.

**The `filter*` verbs strip silently.** They are for text that arrived from
somewhere untidy - a paper, a web page, a numbered listing:

```jennifer
bioseq.filterDna("1 atg acg > x");   # "ATGACG"
bioseq.filterDna("...???...");       # "" - nothing recognisable, no error
```

**Everything else validates.** Whitespace, digits and case are tolerated,
because a pasted, wrapped, numbered sequence is still a sequence. A genuine
non-residue is an error naming the character and its position:

```jennifer
bioseq.reverseComplement("  atg acg\n  60 tta  ");   # fine: "TAACGTCAT"
bioseq.reverseComplement("ATG?C");
# Error{kind: "bioinformatics"}: not a DNA sequence: '?' at position 3
#   (use filterDna / filterRna / filterProtein to strip it)
```

Gap characters raise rather than being dropped, because silently removing a `-`
would shift every coordinate downstream of it. Filter an aligned sequence
deliberately if that is what you want.

## The alphabets

| | |
| - | - |
| `DNA_ALPHABET` | `ACGT` plus the eleven IUPAC ambiguity codes `RYSWKMBDHVN` |
| `RNA_ALPHABET` | `ACGU` plus the same eleven |
| `PROTEIN_ALPHABET` | the twenty standard residues, the ambiguity codes `BZXJ`, the rare `U` and `O`, and `*` for a stop |

Note that the protein alphabet is a superset of both nucleic-acid alphabets, so
`isProtein("ACGT")` and `isDna("ACGT")` are both true. A four-residue peptide
and a four-base oligo are spelled the same way; nothing can tell them apart, and
the deck does not pretend otherwise.

```jennifer
bioseq.isDna("ATGC");     # true
bioseq.isDna("AUGC");     # false - U is RNA
bioseq.isDna("");         # false - an empty string is not a sequence
bioseq.isRna("AUGC");     # true
bioseq.isProtein("MKV*"); # true - a stop is a protein-alphabet character
```

`filterDna` drops `U` rather than guessing that it meant `T`. A sequence spelled
with `U` is RNA: use `filterRna`, or [`backTranscribe`](transforms.md) to move
it to DNA deliberately.

## Coordinates are 0-based and half-open

Every offset the deck reports - `Match.start`, `Orf.start`, `Fragment.start`,
`Alignment.startA` - indexes the sequence the way Jennifer's own slicing does:

```jennifer
def m as bioseq.Match init bioseq.find($dna, "GAATTC")[0];
$dna[$m.start..$m.end];   # the match itself
```

That holds for features found on the reverse strand too: their offsets index the
*forward* sequence you passed in, not the reverse complement. See
[Search](search.md#reverse-strand-coordinates).

**One verb breaks this rule, and is named for it.**
[`subrange`](manipulate.md#subrange) is 1-based and inclusive, because biology
quotes positions that way:

```jennifer
bioseq.subrange($gene, 13, 24);   # "bases 13 to 24", as a paper would say
$gene[12..24];                    # the same slice, deck-style
```

## Molecule kind, only where it matters

A kind argument appears where the answer depends on which molecule it is, and
nowhere else:

```jennifer
bioseq.molecularWeight($s, bioseq.DNA);   # takes a kind
bioseq.reverseComplement($s);             # does not
```

The kinds are `bioseq.DNA`, `bioseq.RNA` and `bioseq.PROTEIN`.

Where a verb accepts either nucleic acid and does not care which - `gcContent`,
`find`, `windows` - it simply accepts both. Where the alphabet must be settled
and the sequence does not settle it, the deck raises rather than guessing:

```jennifer
bioseq.complement("ATGCU");
# Error: sequence mixes T and U, so it is neither DNA nor RNA
#   (use transcribe or backTranscribe to settle on one)
```

## One error shape

Everything raises `Error{kind: "bioinformatics"}` with a message written for a person:

```jennifer
try {
    bioseq.hammingDistance("ATGC", "ATG");
} catch (e) {
    io.printf("%s\n", $e.message);
    # Hamming distance needs equal lengths, got 4 and 3
    #   (use editDistance when the lengths differ)
}
```

Arguments are validated before any work is done, so a verb either produces a
correct answer or raises - it never returns a partial one.

## Ambiguity codes are rejected, not averaged

`N` means "any base". There is no defensible molecular weight for a residue that
might be any of four, and returning the mean of the possibilities would produce
a number that looks authoritative and is not. So the verbs that produce physical
measurements raise:

```jennifer
bioseq.molecularWeight("ATGN", bioseq.DNA);
# Error: no molecular weight for the ambiguity code 'N'
#   - resolve it, or drop it, before weighing the sequence
```

The summary structs take a softer route, because a summary of a real sequence is
still useful: [`dnaStats`](composition.md#summaries) reports length, composition
and GC content as usual, sets an `ambiguous` flag, and reports the undefined
numbers as `0.0` rather than inventing them.

## Value semantics

Jennifer copies on assignment and on argument passing, and this deck relies on
it: every verb returns a new string or list, and nothing you pass in is ever
modified.

```jennifer
def original as string init "ATGCATGC";
def shuffled as string init bioseq.shuffle($original);
# $original is untouched
```
