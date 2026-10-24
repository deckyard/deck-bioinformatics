# bioinformatics examples

One runnable program per topic, each the runnable form of a page in the
[documentation](https://deckyard.github.io/deck-bioinformatics/). None needs a
network or an argument. Run them from the repository root, because those that
read `testdata/` use a path relative to the working directory:

```
jennifer run examples/<topic>.j
```

| Example | What it shows | Reads `testdata/` |
| ------- | ------------- | ----------------- |
| [transforms.j](transforms.j) | complement, transcribe, translate, six frames, alternative genetic codes; start here | no |
| [formats.j](formats.j) | FASTA and FASTQ, and the I/O boundary: `fs` is used by the caller, never inside `src/` | yes |
| [composition.j](composition.j) | GC content, mass, melting temperature, pI, GRAVY, codon usage | yes |
| [search.j](search.j) | IUPAC and fuzzy patterns, ORFs, CpG islands, restriction digest | yes |
| [manipulate.j](manipulate.j) | random sequences, shuffle, mutation, windows, and a shuffle-based null model | yes |
| [compare.j](compare.j) | global and local alignment, substitution matrices, Hamming vs edit distance | yes |
| [tables.j](tables.j) | genetic codes, IUPAC codes, residue properties, enzymes | no |
| [bench.j](bench.j) | the measurements behind the implementation conventions | no |

## testdata/

- `sequences.fasta`: three records. The first 72 nt of the human insulin CDS; a
  synthetic cloning construct (labelled synthetic in the file) carrying six
  known restriction sites, which `search.j` digests; and the protein the first
  record translates to.
- `reads.fastq`: four reads off the insulin fragment with realistic quality
  profiles, one decaying toward the 3' end as Illumina data does and one poor
  throughout.
