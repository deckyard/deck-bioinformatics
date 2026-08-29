# Transforms

Turning one sequence into another: the other strand, the other alphabet, the
protein it encodes.

Runnable version of everything here: `jennifer run examples/transforms.j`.

## The two strands

```jennifer
bioseq.complement("ATGC");          # "TACG"
bioseq.reverseComplement("ATGC");   # "GCAT"
bioseq.reverse("ATGC");             # "CGTA"
```

`complement` pairs the bases and leaves them where they are, which means the
result reads 3′-to-5′ - almost never what you want on its own.
`reverseComplement` is the strand you would actually write down, and applying it
twice returns the original:

```jennifer
bioseq.reverseComplement(bioseq.reverseComplement($gene)) == $gene;   # true
```

`reverse` is a residue operation rather than a strand operation, so it works on
protein too. It exists mostly because `reverseComplement` decomposes into it:

```jennifer
bioseq.reverseComplement($s) == bioseq.reverse(bioseq.complement($s));   # true
```

### Alphabet and ambiguity

The alphabet is preserved. A sequence spelled with `U` complements with `U`; one
spelled with neither `T` nor `U` is read as DNA, which is what every tool does:

```jennifer
bioseq.complement("AUGC");     # "UACG"
bioseq.complement("ACGACG");   # "TGCTGC"
```

Ambiguity codes complement as sets - `R` is A-or-G, so it becomes `Y`, which is
C-or-T - and the self-complementary codes are fixed points:

```jennifer
bioseq.complement("RYSWKM");   # "YRSWMK"
bioseq.complement("BDHV");     # "VHDB"
bioseq.complement("SWN");      # "SWN"
```

A sequence carrying both `T` and `U` has no alphabet, and rather than pick one,
the deck raises.

### Palindromes

A sequence that is its own reverse complement is what a type II restriction
enzyme recognises, because the enzyme binds as a homodimer:

```jennifer
bioseq.reverseComplement("GAATTC") == "GAATTC";   # true - EcoRI
```

## Transcription

```jennifer
bioseq.transcribe("ATGC");       # "AUGC"
bioseq.backTranscribe("AUGC");   # "ATGC"
```

This is the coding-strand shorthand, not the biological event: nothing is
complemented, because an mRNA matches the coding strand and only the thymines
are spelled differently. Each verb insists on its own input alphabet, so
`transcribe` on RNA raises rather than silently doing nothing.

## Translation

```jennifer
bioseq.translate("ATGACGTTACATGGATAA");   # "MTLHG*"
```

RNA is accepted directly - `U` is read as `T` for the codon lookup - so the two
spellings of one gene give the same protein.

Three behaviours are worth knowing:

**A trailing partial codon is dropped.** One or two leftover bases cannot form a
codon, so `translate("ATGAC")` is `"M"`.

**A stop is reported, not obeyed.** Translation reads the whole sequence and
emits `*` at each stop rather than halting, because an ORF search needs to see
where all of them fall:

```jennifer
bioseq.translate("ATGTAAATG");   # "M*M"
```

**Ambiguity resolves when it can.** A codon holding an ambiguity code is
expanded, and if every codon it admits gives the same residue, that residue is
the answer - the ambiguity is real but the translation is not ambiguous:

```jennifer
bioseq.translate("CGN");   # "R" - every CGN codon is arginine
bioseq.translate("GGN");   # "G"
bioseq.translate("TTY");   # "F"
bioseq.translate("ATN");   # "X" - genuinely could be I or M
bioseq.translate("NNN");   # "X"
```

### Other genetic codes

`translateWith` takes a codon table. This matters more often than it looks: a
mitochondrial gene read under the standard code stops early at every `TGA`, and
a ciliate gene never stops at all.

```jennifer
bioseq.translate("TGA");                                        # "*"
bioseq.translateWith("TGA", bioseq.VERTEBRATE_MITOCHONDRIAL);   # "W"
bioseq.translateWith("TAATAG", bioseq.CILIATE_NUCLEAR);         # "QQ"
```

The nine tables the deck carries are listed in
[Genetic codes and residues](tables.md).

Note the naming pattern, which recurs across the deck: `translate` is the common
call, and `translateWith` exists for the case where an alternative table is the
actual point. It is not a performance device - see
[Design notes](design.md#tables-are-read-not-passed).

## Reading frames

Frames `1`, `2` and `3` start at the first, second and third base of the forward
strand. Frames `-1`, `-2` and `-3` do the same on the reverse complement, so
frame `-1` starts at the last base and reads back.

```jennifer
bioseq.translateFrame($gene, 2);
bioseq.translateFrame($gene, -1) == bioseq.translate(bioseq.reverseComplement($gene));  # true
```

`sixFrame` returns all six as `Frame{frame, protein}` records, in the order
`1, 2, 3, -1, -2, -3`:

```jennifer
for (def f in bioseq.sixFrame($gene)) {
    io.printf("%d|sign=always  %s\n", $f.frame, $f.protein);
}
```

This is what you read when a sequence's frame is unknown - a raw read, a contig,
an unannotated region - looking for the frame without stops. For the same job
with start and stop codons taken into account, use
[`orfs`](search.md#open-reading-frames).

`sixFrameWith` and `translateFrameWith` take a codon table, the same way
`translateWith` does.

## Residue spellings

```jennifer
bioseq.oneToThree("MKV");         # "MetLysVal"
bioseq.threeToOne("Met-Lys-Val"); # "MKV"
bioseq.threeToOne("MetLysVal");   # "MKV"
```

Codes may be run together or separated by spaces or hyphens, which is how they
are written in prose. The ambiguity codes and the stop symbol have three-letter
spellings too, so anything `translate` produces can be rewritten and read back:

```jennifer
bioseq.oneToThree("*");     # "Ter"
bioseq.oneToThree("X");     # "Xaa"
bioseq.oneToThree("BZJ");   # "AsxGlxXle"
```

## The surface

| | |
| - | - |
| `complement(s)` | pair each base, keep the order |
| `reverseComplement(s)` | the other strand, read 5′ to 3′ |
| `reverse(s)` | reverse the residues, any sequence |
| `transcribe(dna)` / `backTranscribe(rna)` | T ↔ U |
| `translate(nt)` | under the standard code |
| `translateWith(nt, table)` | under a named code |
| `translateFrame(nt, frame)` | one of the six frames |
| `translateFrameWith(nt, frame, table)` | |
| `sixFrame(nt)` → `list of Frame` | all six |
| `sixFrameWith(nt, table)` | |
| `splitCodons(nt)` → `list of string` | three-base codons in reading order |
| `oneToThree(p)` / `threeToOne(p)` | residue spellings |
