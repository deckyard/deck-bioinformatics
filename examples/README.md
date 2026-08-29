# bioinformatics examples

One runnable tour per topic. Nothing here needs a network or an argument.

```
jennifer run examples/<topic>.j
```

Run them **from the repository root** - every example except `transforms.j`,
`tables.j` and `bench.j` reads `testdata/` by a path relative to the working
directory.

| Example | What it shows | Reads a file? |
| ------- | ------------- | ------------- |
| [transforms.j](transforms.j) | complement, transcribe, translate, six frames, and why the genetic code matters | no |
| [formats.j](formats.j) | FASTA and FASTQ, and the deck's I/O boundary | yes - `testdata/` |
| [composition.j](composition.j) | GC content, mass, melting temperature, pI, GRAVY, codon usage | yes - `testdata/` |
| [search.j](search.j) | IUPAC and fuzzy patterns, open reading frames, CpG islands, restriction digest | yes - `testdata/` |
| [manipulate.j](manipulate.j) | random sequences, composition-preserving shuffle, mutation, windows - and the null model they exist for | yes - `testdata/` |
| [compare.j](compare.j) | global and local alignment, why protein needs a substitution matrix, Hamming vs edit distance | yes - `testdata/` |
| [tables.j](tables.j) | the reference data: genetic codes, IUPAC codes, residue properties, enzymes | no |
| [bench.j](bench.j) | the measurements behind the deck's implementation conventions | no |

Each example is the runnable form of one documentation page - see the
[docs site](https://jennifer-language.github.io/deck-bioinformatics/).

## Where to start

`transforms.j` - it is the deck's front door, and every verb in it takes and
returns a plain string.

## The I/O boundary

`formats.j` is the only example that touches the filesystem, and that is the
point. Nothing under `src/` opens a file: FASTA and FASTQ are parsed from and
written to **strings**, and `fs` is used here in the example, by the caller.
That is what keeps the deck free of a capability pragma and runnable on
`jennifer-tiny`.

## testdata/

- `sequences.fasta` - three records: the first 72 nt of the human insulin CDS,
  a synthetic cloning construct carrying six known restriction sites (labelled
  synthetic in the file, because it is not a real sequence), and the protein
  the first record translates to. `search.j` digests the construct.
- `reads.fastq` - four reads off the insulin fragment with realistic quality
  profiles, including one that decays toward the 3' end the way Illumina data
  does and one that is poor throughout.
