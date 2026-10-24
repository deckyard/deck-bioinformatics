# Sequences and alphabets

The conventions every other page relies on.

## A sequence is a string

```jennifer
bioseq.reverseComplement("ATGC");   # "GCAT"
```

There is no sequence type. Because the type cannot say whether a string is DNA
or protein, the deck looks at the residues, and raises when they do not settle
it.

## Case in, uppercase out

```jennifer
bioseq.complement("atgc");   # "TACG"
```

The exception is [FASTA and FASTQ parsing](formats.md#case-is-preserved), which
preserves case because lowercase marks soft-masked regions.

## Lenient filters, strict everything else

The `filter*` verbs drop anything outside their alphabet without error:

```jennifer
bioseq.filterDna("1 atg acg > x");   # "ATGACG"
bioseq.filterDna("...???...");       # ""
```

Every other verb validates. Whitespace (space, tab, CR, LF), digits and case are
tolerated, so a pasted, numbered, wrapped listing works. Anything else raises,
naming the character and its 0-based byte offset in the input as given:

```jennifer
bioseq.reverseComplement("  atg acg\n  60 tta  ");   # "TAACGTCAT"
bioseq.reverseComplement("ATG?C");
# Error{kind: "bioinformatics"}: not a nucleotide sequence: '?' at position 3
#   (use filterDna / filterRna / filterProtein to strip it)
```

Gap characters (`-`, `.`) raise too: dropping them would shift every coordinate
after them. Filter an aligned sequence deliberately.

## The alphabets

| | |
| - | - |
| `DNA_ALPHABET` | `ACGT` plus the eleven IUPAC ambiguity codes `RYSWKMBDHVN` |
| `RNA_ALPHABET` | `ACGU` plus the same eleven |
| `PROTEIN_ALPHABET` | the twenty standard residues, ambiguity codes `BZXJ`, `U` (selenocysteine), `O` (pyrrolysine) and `*` (stop) |

The protein alphabet contains both nucleic-acid alphabets, so `ACGT` passes
`isProtein` as well as `isDna`.

```jennifer
bioseq.isDna("ATGC");     # true
bioseq.isDna("AUGC");     # false - U is RNA
bioseq.isDna("");         # false
bioseq.isRna("AUGC");     # true
bioseq.isProtein("MKV*"); # true
```

`filterDna` drops `U`. Use `filterRna` for an RNA-spelled sequence, or
[`backTranscribe`](transforms.md#transcription) to convert it.

## Coordinates are 0-based and half-open

Every offset in a struct (`Match.start`, `Orf.start`, `Site.cut`,
`Fragment.start`, `Alignment.startA`, ...) indexes the sequence the way
Jennifer slicing does:

```jennifer
def m as bioseq.Match init bioseq.find($dna, "GAATTC")[0];
$dna[$m.start..$m.end];   # the match
```

Offsets count residues after whitespace and digits are removed. Features found
on the reverse strand report offsets into the forward sequence you passed in
([Search](search.md#reverse-strand-coordinates)).

[`subrange`](manipulate.md#subrange) is the one 1-based, inclusive verb:

```jennifer
bioseq.subrange($gene, 13, 24);   # bases 13 to 24
$gene[12..24];                    # the same slice
```

## Molecule kind only where it matters

A kind argument (`bioseq.DNA`, `bioseq.RNA`, `bioseq.PROTEIN`) appears only where
it changes the answer:

```jennifer
bioseq.molecularWeight($s, bioseq.DNA);   # takes a kind
bioseq.reverseComplement($s);             # does not
```

Verbs that accept either nucleic acid (`gcContent`, `windows`, `complement`)
accept both. Where the alphabet must be known and the sequence holds both `T`
and `U`, the verb raises:

```jennifer
bioseq.complement("ATGCU");
# Error: sequence mixes T and U, so it is neither DNA nor RNA
#   (use transcribe or backTranscribe to settle on one)
```

A sequence with neither `T` nor `U` is treated as DNA.

## One error shape

Every error is `Error{kind: "bioinformatics"}` with a readable message.
Arguments are validated before any work, so a call either returns a complete
answer or raises.

```jennifer
try {
    bioseq.hammingDistance("ATGC", "ATG");
} catch (e) {
    io.printf("%s\n", $e.message);
    # Hamming distance needs equal lengths, got 4 and 3
    #   (use editDistance when the lengths differ)
}
```

## Ambiguity codes are rejected, not averaged

Verbs that compute a physical quantity (mass, Tm, pI, GRAVY) raise on an
ambiguity code:

```jennifer
bioseq.molecularWeight("ATGN", bioseq.DNA);
# Error: no molecular weight for the ambiguity code 'N'
#   - resolve it, or drop it, before weighing the sequence
```

The summaries do not raise: [`dnaStats`](composition.md#summaries) and
`proteinStats` set `ambiguous` and report the undefined figures as `0.0`.

## Value semantics

Jennifer copies on assignment and argument passing. Every verb returns a new
string or list; nothing passed in is modified.
