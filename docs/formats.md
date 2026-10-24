# Formats

Runnable: `jennifer run examples/formats.j`.

## The I/O boundary

These verbs take and return strings. Reading the file is the caller's job:

```jennifer
use fs;

def records as list of bioseq.Record init
    bioseq.parseFasta(fs.readString("sequences.fasta"));
```

This keeps the deck free of a filesystem capability (so it runs on
`jennifer-tiny`) and every parser testable without files.

## FASTA

```jennifer
def records as list of bioseq.Record init bioseq.parseFasta($text);
$records[0].id;            # "sp|P01308|INS_HUMAN"
$records[0].description;   # "Insulin OS=Homo sapiens OX=9606"
$records[0].sequence;      # wrapped lines joined
```

- Wrapped sequence lines, blank lines, CRLF endings and descriptions are
  accepted. Spaces and tabs inside sequence lines are dropped.
- The id is the first whitespace-delimited token after `>`, so pipes and colons
  stay in the id; the trimmed rest of the line is the description.
- Only structure is checked: data before the first `>` and a header with no id
  raise, as does text with no record. Residues are not validated against any
  alphabet.
- A header with no sequence lines gives an empty `sequence`.
- Parsing stops with an error past `MAX_RECORDS` (1 000 000) records. A bare `>`
  line is three bytes of input and a whole `Record` in memory, so a hostile file
  of them expands about a hundredfold; the cap bounds that.

### Case is preserved

Lowercase in a FASTA body is soft-masking, so the parsers keep case, unlike the
rest of the deck. Normalise explicitly:

```jennifer
bioseq.parseFasta(">m\nATGCatgcATGC\n")[0].sequence;   # "ATGCatgcATGC"
bioseq.filterDna($r.sequence);                        # "ATGCATGCATGC"
```

### Writing

```jennifer
bioseq.writeFasta($records);                  # wrapped at FASTA_WIDTH = 60
bioseq.writeFastaWrapped($records, 70);
bioseq.wrap($sequence, 60);                   # no trailing newline
```

Output ends in a newline. A record with an empty sequence is written as its
header alone. Parsing the output returns the same records at any width. An
empty list, a record with no id, or a width below 1 raises.

### Pasted text

```jennifer
bioseq.toPlain(">x\n1 atg acg\n");   # "atgacg"
```

`toPlain` drops lines starting with `>` or `;`, digits, and every byte at or
below space (whitespace and control characters). It keeps case and does not
check the alphabet. Use it for sequences copied from a paper or web page.

## FASTQ

```jennifer
def reads as list of bioseq.Read init bioseq.parseFastq($text);
$reads[0].sequence;
$reads[0].quality;   # exactly as long as the sequence
```

Both the four-line form and the wrapped form are read. Sequence lines run to the
`+` line; quality lines are then read until the quality is as long as the
sequence. Quality is not delimited by `@` because `@` is a valid Phred+33
character (Q31). A length mismatch, a missing `+` line or a record not starting
with `@` raises, as does exceeding `MAX_RECORDS`.

`writeFastq` writes the four-line form with a bare `+` line. `fastqToFasta`
drops the qualities.

### Quality scores

```jennifer
bioseq.qualityScores($read);                          # Phred+33
bioseq.qualityScoresWith($read, bioseq.PHRED64);      # Illumina 1.3-1.7
bioseq.meanQuality($read);
bioseq.meanQualityWith($read, bioseq.PHRED64);
bioseq.errorProbability(30);                          # 0.001
```

`P = 10^(-Q/10)`: Q20 is one error in 100, Q30 one in 1000.

The offset is never guessed. Reading Phred+33 data as Phred+64 shifts every
score by 31; the resulting negative scores raise. A `Read` built by hand is
checked too: quality and sequence lengths must match, as the parser and writer
require.

`meanQuality` is the arithmetic mean of the scores, the figure quality plots
show. Phred is logarithmic, so it understates the worst bases; use
`errorProbability` per score when that matters. A read with no bases raises.
