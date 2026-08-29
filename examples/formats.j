#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - reading and writing FASTA and FASTQ.
#
#   jennifer run examples/formats.j
#
# Run it from the repository root - the testdata/ paths are relative to the
# working directory, not to this file.
#
# This is the example that shows the deck's I/O boundary. bioinformatics parses and
# writes *text*; opening the file is the caller's job, and `fs` is used here,
# not inside src/. That is what keeps the deck free of a filesystem capability
# and runnable on jennifer-tiny.

use io;
use fs;
use math;
use strings;
use convert;

import "../src/bioinformatics.j" as bioseq;

# --------------------------------------------------------------------- FASTA

def text as string init fs.readString("testdata/sequences.fasta");
def records as list of bioseq.Record init bioseq.parseFasta($text);

io.printf("%s: %d records\n\n", "testdata/sequences.fasta", len($records));

for (def r in $records) {
    io.printf(
        "  %s|pad=20|align=left %d|pad=4 residues  %s\n",
        $r.id,
        len($r.sequence),
        $r.description);
}

# The id is the first whitespace-delimited token, so a UniProt-style id keeps
# its pipes and only the prose after the first space is the description.
io.printf("\n  first record\n");
io.printf("    id           %s\n", $records[0].id);
io.printf("    description  %s\n", $records[0].description);
io.printf("    sequence     %s\n", $records[0].sequence);

# ------------------------------------------------------- records feed the deck

# A parsed record is just a string in a struct, so every other verb takes it
# straight away.
def gene as bioseq.Record init $records[0];
io.printf("\n  translated   %s\n", bioseq.translate($gene.sequence));
io.printf(
    "  matches the protein record: %t\n",
    bioseq.translate($gene.sequence) == $records[2].sequence);

# ------------------------------------------------------------- writing it back

# Writing wraps at the conventional sixty columns; writeFastaWrapped picks
# another width. Parsing what you wrote returns what you had.
def oneRecord as list of bioseq.Record init [$gene];
io.printf("\n  rewritten at 30 columns\n\n");
for (def line in strings.split(bioseq.writeFastaWrapped($oneRecord, 30), "\n")) {
    if (len($line) > 0) {
        io.printf("    %s\n", $line);
    }
}

io.printf(
    "\n  round-trips: %t\n",
    bioseq.parseFasta(bioseq.writeFasta($records))[1].sequence == $records[1].sequence);

# --------------------------------------------------------------------- FASTQ

def reads as list of bioseq.Read init bioseq.parseFastq(fs.readString("testdata/reads.fastq"));

io.printf("\n\n%s: %d reads\n\n", "testdata/reads.fastq", len($reads));
io.printf(
    "  %s|pad=8|align=left %s|pad=5 %s|pad=8 %s|pad=10 %s\n",
    "read",
    "len",
    "mean Q",
    "err rate",
    "note");

for (def r in $reads) {
    def mean as float init bioseq.meanQuality($r);
    def note as string init "";
    if ($mean < 20.0) {
        $note = "below Q20 - would be trimmed";
    }
    io.printf(
        "  %s|pad=8|align=left %d|pad=5 %f|prec=1|pad=8 %f|prec=4|pad=10 %s\n",
        $r.id,
        len($r.sequence),
        $mean,
        bioseq.errorProbability(convert.toInt(math.round($mean))),
        $note);
}

# Per-base scores are what a quality plot draws. Read 2 is the classic Illumina
# shape: clean at the 5' end, decaying toward the 3'.
io.printf("\n  read2, base by base\n\n    ");
for (def q in bioseq.qualityScores($reads[1])) {
    if ($q >= 30) {
        io.printf("#");
    } elseif ($q >= 20) {
        io.printf("+");
    } elseif ($q >= 10) {
        io.printf("-");
    } else {
        io.printf(".");
    }
}
io.printf("\n    # Q30+   + Q20+   - Q10+   . below Q10\n");

# The offset is explicit rather than guessed: reading a Phred+33 file as
# Phred+64 shifts every score by 31, and the deck raises rather than reporting
# quietly wrong numbers.
io.printf(
    "\n  first base of read1 at Phred+33  Q%d\n",
    bioseq.qualityScoresWith($reads[0], bioseq.PHRED33)[0]);
io.printf(
    "  the same character at Phred+64   Q%d\n",
    bioseq.qualityScoresWith($reads[0], bioseq.PHRED64)[0]);

# --------------------------------------------------------------- FASTQ to FASTA

def asFasta as list of bioseq.Record init bioseq.fastqToFasta($reads);
io.printf("\n  dropped to FASTA: %d records, first is %s\n", len($asFasta), $asFasta[0].id);
