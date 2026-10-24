# Transforms

Runnable: `jennifer run examples/transforms.j`.

## Strands

```jennifer
bioseq.complement("ATGC");          # "TACG" - reads 3' to 5'
bioseq.reverseComplement("ATGC");   # "GCAT" - the other strand, 5' to 3'
bioseq.reverse("ATGC");             # "CGTA"
```

`reverseComplement` applied twice returns the original, and equals
`reverse(complement(s))`. `reverse` works on protein too.

The alphabet is preserved: a sequence containing `U` complements with `U`; one
with neither `T` nor `U` is DNA. Ambiguity codes complement as sets (`R` = A/G
becomes `Y` = C/T); `S`, `W` and `N` are their own complements. A sequence
holding both `T` and `U` raises.

```jennifer
bioseq.complement("AUGC");     # "UACG"
bioseq.complement("ACGACG");   # "TGCTGC"
bioseq.complement("RYSWKM");   # "YRSWMK"
bioseq.complement("BDHV");     # "VHDB"
bioseq.complement("SWN");      # "SWN"
```

## Transcription

```jennifer
bioseq.transcribe("ATGC");       # "AUGC"
bioseq.backTranscribe("AUGC");   # "ATGC"
```

Coding-strand convention: `T` and `U` are swapped and nothing is complemented.
`transcribe` requires DNA and `backTranscribe` RNA; the wrong one raises.

## Translation

```jennifer
bioseq.translate("ATGACGTTACATGGATAA");   # "MTLHG*"
```

- RNA is accepted; `U` is read as `T`.
- A trailing one or two bases are dropped: `translate("ATGAC")` is `"M"`.
- A stop becomes `*` and translation continues: `translate("ATGTAAATG")` is
  `"M*M"`.
- A codon with an ambiguity code gives the residue all its expansions agree on,
  otherwise `X`:

```jennifer
bioseq.translate("CGN");   # "R"
bioseq.translate("GGN");   # "G"
bioseq.translate("TTY");   # "F"
bioseq.translate("ATN");   # "X" - I or M
bioseq.translate("NNN");   # "X"
```

### Other genetic codes

`translateWith` takes a codon table. The nine tables are listed in
[Genetic codes](tables.md#genetic-codes).

```jennifer
bioseq.translate("TGA");                                        # "*"
bioseq.translateWith("TGA", bioseq.VERTEBRATE_MITOCHONDRIAL);   # "W"
bioseq.translateWith("TAATAG", bioseq.CILIATE_NUCLEAR);         # "QQ"
```

A `With` suffix always means an explicit alternative to a default (Jennifer has
no default arguments), never a performance option.

## Reading frames

Frames `1`, `2`, `3` start at the first, second and third base. Frames `-1`,
`-2`, `-3` do the same on the reverse complement, so frame `-1` starts at the
last base. A frame starting past the end translates to `""`; any other frame
number raises.

```jennifer
bioseq.translateFrame($gene, 2);
bioseq.translateFrame($gene, -1) == bioseq.translate(bioseq.reverseComplement($gene));  # true
bioseq.translateFrameWith($mito, 1, bioseq.VERTEBRATE_MITOCHONDRIAL);
```

`sixFrame` (and `sixFrameWith`, which takes a table) returns six `Frame`
records in the order `1, 2, 3, -1, -2, -3`:

```jennifer
for (def f in bioseq.sixFrame($gene)) {
    io.printf("%d|sign=always  %s\n", $f.frame, $f.protein);
}
```

For start and stop codons, use [`orfs`](search.md#open-reading-frames).
`splitCodons` returns the codons in reading order, dropping a trailing partial
codon.

## Residue spellings

```jennifer
bioseq.oneToThree("MKV");         # "MetLysVal"
bioseq.threeToOne("Met-Lys-Val"); # "MKV"
bioseq.threeToOne("MetLysVal");   # "MKV"
bioseq.oneToThree("*");           # "Ter"
bioseq.oneToThree("X");           # "Xaa"
bioseq.oneToThree("BZJ");         # "AsxGlxXle"
```

`threeToOne` accepts codes run together or separated by spaces or hyphens.
Ambiguity codes and the stop have three-letter spellings (`THREE_LETTER`), so
any `translate` output round-trips.

Signatures: [cheatsheet](cheatsheet.md#transforms).
