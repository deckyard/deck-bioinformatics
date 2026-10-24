# tools

Generators for the two reference-data files. They are not part of the deck - a
consumer never runs them, and `vendor/` never sees them - but the files they
emit carry a `GENERATED` header, so the generators have to live somewhere a
maintainer can find them.

**Generated data and hand-written API are separate files on purpose.**
`tabledata.j` and `enzymedata.j` hold nothing but `def const` declarations and
are overwritten wholesale; `tables.j` and `enzymes.j` are hand-written
accessors over them and are never generated. Re-running a generator cannot
clobber an accessor, which is exactly what happened once before the split.

```
jennifer run tools/gentables.j     # -> src/topics/tabledata.j
jennifer run tools/genenzymes.j    # -> src/topics/enzymedata.j
jennifer fmt -w src/topics/tabledata.j src/topics/enzymedata.j
jennifer test src/bioinformatics_test.j
```

The generators are Jennifer too. There is nothing in them a second language
was needed for - string handling, maps and a file write - and keeping them in
Jennifer means the repository needs exactly one toolchain, CI included.

They emit through raw `'...'` strings. A generator writes Jennifer source, so
nearly every emitted line contains double quotes and braces; in a raw string
both are literal, with no escape processing and no `{expr}` interpolation slot
to escape around. Only lines carrying a real newline escape, or an apostrophe
of their own, stay cooked - a raw string ends at the first `'`.

Run from the repository root; both write relative to the working directory.

The reasoning behind this split is written up in
[Design notes](https://deckyard.github.io/deck-bioinformatics/design.html#generated-data-hand-written-api).

## Why generate at all

Both files are wide, repetitive data where a transcription slip is invisible on
review and wrong forever afterwards. The generators make the *source* of each
table the thing under review:

- `gentables.j` holds each alternative genetic code as the standard code plus
  the codons NCBI documents it reassigning, which is a four-line diff rather
  than a 64-line blob. It asserts the standard table's own anchors before
  emitting anything.
- `genenzymes.j` asserts that no enzyme name is duplicated, that every
  recognition site is spelled in the IUPAC alphabet, that both cut offsets fall
  inside the site, and that every site is its own reverse complement - with no
  exception list, since the degenerate sites are palindromic too (`N`
  complements to `N`, `W` to `W`, and the specificity lives in the outer bases,
  which pair up). A typo in a new row fails the generator, not the review.

The emitted `.j` still holds full literal maps: the benchmark in
`examples/bench.j` showed a const read costs nothing while a per-call
construction does.

## The tests are the second gate

`tables_test.j` re-states the genetic-code difference lists independently and
fails if a generated table stops matching one. `enzymes_test.j` re-derives the
overhang label from the cut offsets and re-checks the palindromes. So the data
is verified on both sides of the generator, and editing the `.j` by hand
without updating the generator is caught rather than silently overwritten.
