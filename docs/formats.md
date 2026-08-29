# Formats

FASTA and FASTQ, over strings.

Runnable version: `jennifer run examples/formats.j`.

## The I/O boundary

**Nothing in this deck opens a file.** These verbs take and return text, and
reading the file is the caller's job:

```jennifer
use fs;

def records as list of bioseq.Record init
    bioseq.parseFasta(fs.readString("sequences.fasta"));
```

That is a deliberate line. It keeps the deck free of a filesystem capability -
so it runs on `jennifer-tiny` - and it makes every parser testable without a
fixture on disk. `examples/formats.j` is the only example that touches `fs`, and
it does so as the caller.

## FASTA

```jennifer
def records as list of bioseq.Record init bioseq.parseFasta($text);
$records[0].id;            # "sp|P01308|INS_HUMAN"
$records[0].description;   # "Insulin OS=Homo sapiens OX=9606"
$records[0].sequence;      # the residues, wrapped lines joined
```

The parser handles the whole conventional shape: wrapped sequence lines, blank
lines between records, Windows line endings, and headers with a description
after the id.

The **id is the first whitespace-delimited token**, which is why a UniProt-style
id keeps its pipes and colons intact and only the prose after the first space
becomes the description.

### Only structure is checked

A FASTA file may hold DNA, RNA or protein, and the format does not say which, so
the residues are not validated against any alphabet. What *is* checked is
structural: data before the first `>` header raises, and so does a header with
no id.

A record with a header but no sequence parses to an empty `sequence` rather than
raising. A parser's job is to report what the file says; `len($r.sequence)` is
your guard.

Both parsers stop at `MAX_RECORDS` - a million, which no real file approaches.
It exists because a header line with no sequence is three bytes of input and a
whole `Record` in memory, so a hostile file of nothing but `>` lines expands
about a hundredfold. Bounding the count bounds that.

### Case is preserved

This is the deck's one departure from
[uppercase-out](sequences.md#case-in-uppercase-out). Lowercase in a FASTA body is
**soft-masking** - it marks repeats and low-complexity regions - and discarding
it would throw away real annotation:

```jennifer
bioseq.parseFasta(">m\nATGCatgcATGC\n")[0].sequence;   # "ATGCatgcATGC"
bioseq.filterDna($r.sequence);                        # "ATGCATGCATGC" - opt in
```

### Writing

```jennifer
bioseq.writeFasta($records);                  # wrapped at 60 columns
bioseq.writeFastaWrapped($records, 70);       # another width
bioseq.wrap($sequence, 60);                   # just the wrapping
```

`FASTA_WIDTH` is 60, which is what NCBI and EMBL emit. Parsing what you wrote
returns what you had, at any width.

### Pasted text

```jennifer
bioseq.toPlain(">x\n1 atg acg\n");   # "atgacg"
```

`toPlain` reduces any pasted sequence text to bare residues - dropping headers,
position numbers and all whitespace. Unlike `filterDna` it does not judge the
residues, and unlike `parseFasta` it does not need a header to be present. It is
for the shape a sequence arrives in when copied out of a paper or a database
page.

## FASTQ

```jennifer
def reads as list of bioseq.Read init bioseq.parseFastq($text);
$reads[0].sequence;
$reads[0].quality;   # exactly as long as the sequence
```

The parser reads the four-line form every sequencer emits, and also the wrapped
form the specification permits.

### The trap it avoids

`@` is a legal Phred+33 quality character - it is Q31. So a quality line can
begin with one, and a parser that finds records by looking for `@` at the start
of a line will desynchronise on real data.

This parser reads sequence lines until the `+` separator, then quality lines
**until the quality is as long as the sequence**. That length agreement is the
one structural invariant FASTQ has, and it is enforced on both parse and write.

### Quality scores

```jennifer
bioseq.qualityScores($read);                          # Phred+33
bioseq.qualityScoresWith($read, bioseq.PHRED64);      # archived Illumina
bioseq.meanQuality($read);
bioseq.errorProbability(30);                          # 0.001
```

Phred is logarithmic: `P = 10^(-Q/10)`, so Q20 is one error in a hundred and Q30
one in a thousand.

**The offset is explicit rather than guessed.** Reading a Phred+33 file as
Phred+64 shifts every score by 31, and rather than report quietly wrong numbers,
the deck raises on the resulting negative scores.

**A `Read` is checked before its quality is read.** The struct is exported, so
one can be built by hand rather than parsed - and the quality verbs now verify
that the quality is as long as the sequence, the same invariant the parser and
the writer enforce. Without it a ragged read reported a confident mean over
whichever bases happened to have a score.

> `meanQuality` is the **arithmetic mean of the scores** - the figure every tool
> reports and every quality plot shows. It is *not* the score corresponding to
> the read's mean error probability. Because Phred is logarithmic, averaging
> understates the effect of the worst bases. Use `errorProbability` on the
> individual scores when that distinction matters.

### Dropping to FASTA

```jennifer
bioseq.fastqToFasta($reads);   # list of Record, qualities discarded
```

## The surface

| | |
| - | - |
| `Record` | `id`, `description`, `sequence` |
| `Read` | `id`, `description`, `sequence`, `quality` |
| `parseFasta(text)` | → `list of Record` |
| `writeFasta(records)` / `writeFastaWrapped(records, width)` | |
| `wrap(s, width)` | fixed-width lines |
| `toPlain(text)` | strip everything but residues |
| `parseFastq(text)` | → `list of Read` |
| `writeFastq(reads)` | four-line form |
| `fastqToFasta(reads)` | → `list of Record` |
| `qualityScores(r)` / `qualityScoresWith(r, offset)` | → `list of int` |
| `meanQuality(r)` / `meanQualityWith(r, offset)` | |
| `errorProbability(score)` | |
| `FASTA_WIDTH`, `PHRED33`, `PHRED64` | constants |
