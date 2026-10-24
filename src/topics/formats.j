# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - formats: FASTA and FASTQ, over strings.
# Spliced into bioinformatics.j via include - not a standalone module.
#
# Nothing here opens a file. These verbs take and return text, and reading the
# file is the caller's job - which is what keeps the deck free of a filesystem
# capability and every parser testable without a fixture on disk.
#
# Two deliberate departures from the rest of the deck:
#
#   1. Case is preserved. Everywhere else residues come back uppercase, but
#      lowercase in a FASTA body is soft-masking - repeats and low-complexity
#      regions - and discarding it would throw away real annotation. Pass the
#      sequence through `filterDna` when you want the residues normalised.
#   2. The alphabet is not checked. A FASTA file holds DNA, RNA or protein, and
#      the format does not say which; only structure is validated here.

/**
 * The FASTA line width `writeFasta` uses, as NCBI and EMBL emit.
 */
export def const FASTA_WIDTH as int init 60;

/**
 * The most records one `parseFasta` or `parseFastq` call builds.
 *
 * An amplification guard: a header line with no sequence is three bytes of input
 * and a whole `Record` in memory, so a file of bare `>` lines expands about a
 * hundredfold. A proteome has tens of thousands of records, a transcriptome a few
 * hundred thousand.
 */
export def const MAX_RECORDS as int init 1000000;

/**
 * Phred+33 quality offset: Sanger, Illumina 1.8 and later. The deck's default.
 */
export def const PHRED33 as int init 33;

/**
 * Phred+64 quality offset: Illumina 1.3 to 1.7, archived data only.
 */
export def const PHRED64 as int init 64;

/**
 * One FASTA record.
 *
 * @field id {string} the first whitespace-delimited token after `>`
 * @field description {string} whatever followed it on the header line, "" when
 *        the header carried only an id
 * @field sequence {string} the residues, joined across the wrapped lines with
 *        their case preserved
 */
export def struct Record {
    id as string,
    description as string,
    sequence as string
};

/**
 * One FASTQ read: a sequence and its per-base quality string.
 *
 * @field id {string} the first whitespace-delimited token after `@`
 * @field description {string} the rest of the header line, "" when absent
 * @field sequence {string} the read
 * @field quality {string} the encoded quality, exactly as long as `sequence`
 */
export def struct Read {
    id as string,
    description as string,
    sequence as string,
    quality as string
};

# ----------------------------------------------------------------- line helpers

/**
 * Split text into lines, tolerating both Unix and Windows endings.
 *
 * @param text {string} the text to split
 * @return {list of string} the lines, without their terminators
 * @internal
 */
func splitLines(text as string) {
    def out as list of string;
    for (def line in strings.split($text, "\n")) {
        if (strings.endsWith($line, "\r")) {
            $out[] = $line[0..len($line) - 1];
        } else {
            $out[] = $line;
        }
    }
    return $out;
}

/**
 * Drop every space and tab from a line, keeping everything else.
 *
 * @param line {string} the line to compact
 * @return {string} the line without internal whitespace
 * @internal
 */
func stripSpace(line as string) {
    def buf as bytes init convert.bytesFromString($line, "utf-8");
    def n as int init len($buf);
    def out as bytes;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def code as int init $buf[$i];
        if ($code == 32 or $code == 9) {
            continue;
        }
        $out[] = $code;
    }
    return convert.stringFromBytes($out, "utf-8");
}

/**
 * Split a header line into its id and description.
 *
 * The id is everything up to the first space or tab; the description is the
 * remainder, trimmed. This is the convention every FASTA consumer follows, and
 * it is why an id may itself contain pipes and colons - `sp|P01308|INS_HUMAN`
 * is one token.
 *
 * @param header {string} the header line, without its `>` or `@` marker
 * @param marker {string} the marker character, for the error message
 * @return {list of string} exactly two entries: id, then description
 * @throws {Error} kind "bioinformatics" when the header carries no id
 * @internal
 */
func splitHeader(header as string, marker as string) {
    def trimmed as string init strings.trim($header);
    if (len($trimmed) == 0) {
        raiseError("a '" + $marker + "' header line carries no identifier");
    }
    def cut as int init -1;
    def buf as bytes init convert.bytesFromString($trimmed, "utf-8");
    def n as int init len($buf);
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        if ($buf[$i] == 32 or $buf[$i] == 9) {
            $cut = $i;
            break;
        }
    }
    def out as list of string;
    if ($cut < 0) {
        $out[] = $trimmed;
        $out[] = "";
        return $out;
    }
    $out[] = $trimmed[0..$cut];
    $out[] = strings.trim($trimmed[$cut..]);
    return $out;
}

# ------------------------------------------------------------------------ FASTA

/**
 * Parse FASTA text into records.
 *
 * Accepts wrapped sequence lines, blank lines between records, CRLF line endings
 * and a description after the id. Spaces and tabs inside sequence lines are
 * dropped. Only structure is checked: residues are not validated against an
 * alphabet, and case is preserved (lowercase is soft-masking).
 *
 * A header with no sequence lines gives an empty `sequence`, not an error. At
 * most `MAX_RECORDS` records are built.
 *
 * @param text {string} the contents of a FASTA file
 * @return {list of Record} the records, in file order
 * @throws {Error} kind "bioinformatics" when there is no record, when data appears
 *         before the first header, when a header carries no id, or past `MAX_RECORDS`
 * @example bioseq.parseFasta(">seq1 first\nATGC\n")[0].id == "seq1"
 */
export func parseFasta(text as string) {
    def out as list of Record;
    def id as string init "";
    def description as string init "";
    def sequence as string init "";
    def inRecord as bool init false;

    for (def line in splitLines($text)) {
        if (len(strings.trim($line)) == 0) {
            continue;
        }
        if (strings.startsWith($line, ">")) {
            if ($inRecord) {
                requireRoom(len($out) + 1, "FASTA records");
                $out[] = Record{id: $id, description: $description, sequence: $sequence};
            }
            def parts as list of string init splitHeader($line[1..], ">");
            $id = $parts[0];
            $description = $parts[1];
            $sequence = "";
            $inRecord = true;
        } else {
            if (not $inRecord) {
                raiseError("FASTA data before the first '>' header line" +
                    " (is this a plain sequence? use toPlain)");
            }
            $sequence = $sequence + stripSpace($line);
        }
    }
    if ($inRecord) {
        requireRoom(len($out) + 1, "FASTA records");
        $out[] = Record{id: $id, description: $description, sequence: $sequence};
    }
    if (len($out) == 0) {
        raiseError("no FASTA records found (a record starts with '>')");
    }
    return $out;
}

/**
 * Render records as FASTA text, wrapped at `FASTA_WIDTH` (60) columns.
 *
 * A record with an empty sequence is written as its header line alone.
 *
 * @param records {list of Record} the records to write
 * @return {string} FASTA text, ending in a newline
 * @throws {Error} kind "bioinformatics" when the list is empty or a record has no id
 * @example bioseq.writeFasta([$r]) == ">seq1 first\nATGC\n"
 */
export func writeFasta(records as list of Record) {
    return writeFastaWrapped($records, FASTA_WIDTH);
}

/**
 * Render records as FASTA text at a chosen line width.
 *
 * @param records {list of Record} the records to write
 * @param width {int} the wrap column; must be positive
 * @return {string} FASTA text, ending in a newline
 * @throws {Error} kind "bioinformatics" when the list is empty, the width is not
 *         positive, or a record has no id
 */
export func writeFastaWrapped(records as list of Record, width as int) {
    if (len($records) == 0) {
        raiseError("no records to write");
    }
    if ($width <= 0) {
        raiseError("FASTA line width must be positive, got " + convert.toString($width));
    }
    def out as string init "";
    for (def r in $records) {
        if (len(strings.trim($r.id)) == 0) {
            raiseError("a record has no id, so it cannot be written as FASTA");
        }
        $out = $out + ">" + $r.id;
        if (len($r.description) > 0) {
            $out = $out + " " + $r.description;
        }
        $out = $out + "\n";
        if (len($r.sequence) > 0) {
            $out = $out + wrap($r.sequence, $width) + "\n";
        }
    }
    return $out;
}

/**
 * Break a sequence into fixed-width lines.
 *
 * @param s {string} the sequence to wrap
 * @param width {int} the line width; must be positive
 * @return {string} the wrapped sequence, with no trailing newline
 * @throws {Error} kind "bioinformatics" when the width is not positive
 * @example bioseq.wrap("ATGCATGC", 4) == "ATGC\nATGC"
 */
export func wrap(s as string, width as int) {
    if ($width <= 0) {
        raiseError("wrap width must be positive, got " + convert.toString($width));
    }
    def buf as bytes init convert.bytesFromString($s, "utf-8");
    def n as int init len($buf);
    def out as bytes;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        if ($i > 0 and $i % $width == 0) {
            $out[] = 10;
        }
        $out[] = $buf[$i];
    }
    return convert.stringFromBytes($out, "utf-8");
}

/**
 * Reduce any pasted sequence text to bare residues.
 *
 * Drops lines starting with `>` or `;`, digits, and every byte at or below
 * space (whitespace and control characters alike). Everything else is kept as
 * is: unlike `filterDna` it does not judge the residues, and unlike
 * `parseFasta` it needs no header. Use `filterDna` or a validating verb when the
 * input should be judged rather than cleaned.
 *
 * @param text {string} pasted sequence text
 * @return {string} the residues, in their original case
 * @example bioseq.toPlain(">x\n1 atg acg\n") == "atgacg"
 */
export func toPlain(text as string) {
    def out as string init "";
    for (def line in splitLines($text)) {
        if (strings.startsWith(strings.trim($line), ">")) {
            continue;
        }
        if (strings.startsWith(strings.trim($line), ";")) {
            continue;
        }
        def buf as bytes init convert.bytesFromString($line, "utf-8");
        def kept as bytes;
        def n as int init len($buf);
        for (def i as int init 0; $i < $n; $i = $i + 1) {
            def code as int init $buf[$i];
            if ($code <= 32) {
                continue;
            }
            if ($code >= 48 and $code <= 57) {
                continue;
            }
            $kept[] = $code;
        }
        $out = $out + convert.stringFromBytes($kept, "utf-8");
    }
    return $out;
}

# ------------------------------------------------------------------------ FASTQ

/**
 * Parse FASTQ text into reads.
 *
 * Reads the four-line form and the wrapped form: sequence lines are gathered
 * until a `+` line, then quality lines until the quality is as long as the
 * sequence. Quality is read by length, not by looking for `@`, because `@` is a
 * valid Phred+33 character (Q31). At most `MAX_RECORDS` reads are built.
 *
 * @param text {string} the contents of a FASTQ file
 * @return {list of Read} the reads, in file order
 * @throws {Error} kind "bioinformatics" when there is no record, a record does not
 *         start with `@`, lacks its `+` separator, or has a quality string of the
 *         wrong length, or past `MAX_RECORDS`
 * @example bioseq.parseFastq("@r1\nATGC\n+\nIIII\n")[0].quality == "IIII"
 */
export func parseFastq(text as string) {
    def lines as list of string init splitLines($text);
    def n as int init len($lines);
    def out as list of Read;
    def i as int init 0;

    while ($i < $n) {
        if (len(strings.trim($lines[$i])) == 0) {
            $i = $i + 1;
            continue;
        }
        if (not strings.startsWith($lines[$i], "@")) {
            raiseError("expected a FASTQ record to start with '@' at line " +
                convert.toString($i + 1) + ", found: " + summarise($lines[$i]));
        }
        def parts as list of string init splitHeader($lines[$i][1..], "@");
        $i = $i + 1;

        def sequence as string init "";
        while ($i < $n and not strings.startsWith($lines[$i], "+")) {
            $sequence = $sequence + stripSpace($lines[$i]);
            $i = $i + 1;
        }
        if ($i >= $n) {
            raiseError("FASTQ record '" + $parts[0] + "' has no '+' separator line");
        }
        $i = $i + 1;

        def quality as string init "";
        while ($i < $n and len($quality) < len($sequence)) {
            $quality = $quality + stripSpace($lines[$i]);
            $i = $i + 1;
        }
        if (len($quality) != len($sequence)) {
            raiseError("FASTQ record '" + $parts[0] + "' has " + convert.toString(len($sequence)) +
                " bases but " + convert.toString(len($quality)) + " quality characters");
        }
        requireRoom(len($out) + 1, "FASTQ reads");
        $out[] = Read{
            id: $parts[0],
            description: $parts[1],
            sequence: $sequence,
            quality: $quality
        };
    }
    if (len($out) == 0) {
        raiseError("no FASTQ records found (a record starts with '@')");
    }
    return $out;
}

/**
 * Stop a parse that has produced implausibly many records.
 *
 * @param count {int} how many there would be once the next one is added
 * @param what {string} what is being counted, for the message
 * @throws {Error} kind "bioinformatics" when the count exceeds `MAX_RECORDS`
 * @internal
 */
func requireRoom(count as int, what as string) {
    if ($count > MAX_RECORDS) {
        raiseError("more than " + convert.toString(MAX_RECORDS) + " " + $what +
            " in one parse - refusing to continue, since a file of bare" +
            " header lines expands about a hundredfold in memory");
    }
}

/**
 * Shorten a line for an error message.
 *
 * @param line {string} the offending line
 * @return {string} at most forty characters of it, quoted
 * @internal
 */
func summarise(line as string) {
    if (len($line) <= 40) {
        return "'" + $line + "'";
    }
    return "'" + $line[0..40] + "...'";
}

/**
 * Render reads as four-line FASTQ text.
 *
 * The `+` separator is written bare rather than repeating the id, which is what
 * every current tool emits and what halves the size of a large file.
 *
 * @param reads {list of Read} the reads to write
 * @return {string} FASTQ text, ending in a newline
 * @throws {Error} kind "bioinformatics" when the list is empty, a read has no id, or a
 *         read's quality length disagrees with its sequence length
 */
export func writeFastq(reads as list of Read) {
    if (len($reads) == 0) {
        raiseError("no reads to write");
    }
    def out as string init "";
    for (def r in $reads) {
        if (len(strings.trim($r.id)) == 0) {
            raiseError("a read has no id, so it cannot be written as FASTQ");
        }
        if (len($r.sequence) != len($r.quality)) {
            raiseError("read '" + $r.id + "' has " + convert.toString(len($r.sequence)) +
                " bases but " + convert.toString(len($r.quality)) + " quality characters");
        }
        $out = $out + "@" + $r.id;
        if (len($r.description) > 0) {
            $out = $out + " " + $r.description;
        }
        $out = $out + "\n" + $r.sequence + "\n+\n" + $r.quality + "\n";
    }
    return $out;
}

/**
 * Drop the quality strings, leaving FASTA records.
 *
 * @param reads {list of Read} the reads to convert
 * @return {list of Record} the same sequences, without their qualities
 */
export func fastqToFasta(reads as list of Read) {
    def out as list of Record;
    for (def r in $reads) {
        $out[] = Record{id: $r.id, description: $r.description, sequence: $r.sequence};
    }
    return $out;
}

# --------------------------------------------------------------------- quality

/**
 * Decode a read's quality string into Phred scores, assuming Phred+33.
 *
 * @param r {Read} the read
 * @return {list of int} one score per base
 * @throws {Error} kind "bioinformatics" when the quality and sequence lengths differ,
 *         or a score decodes negative
 * @example bioseq.qualityScores($r)[0] == 40   # 'I' under Phred+33
 */
export func qualityScores(r as Read) {
    return qualityScoresWith($r, PHRED33);
}

/**
 * Decode a read's quality string at a chosen offset.
 *
 * Use `PHRED64` for archived Illumina 1.3-1.7 data; everything current is
 * `PHRED33`. Reading one encoding as the other shifts every score by 31, so the
 * offset is never guessed.
 *
 * @param r {Read} the read
 * @param offset {int} `PHRED33` or `PHRED64`
 * @return {list of int} one score per base
 * @throws {Error} kind "bioinformatics" when the quality and sequence lengths differ,
 *         or a quality character decodes to a negative score (wrong offset)
 */
export func qualityScoresWith(r as Read, offset as int) {
    # A `Read` can arrive from somewhere other than `parseFastq` - built by
    # hand, or carried in from another tool - and the one invariant FASTQ has is
    # that the quality is as long as the sequence. Both the parser and the
    # writer enforce it; without this the quality verbs would quietly report a
    # mean over whichever bases happened to have a score.
    if (len($r.sequence) != len($r.quality)) {
        raiseError("read '" + $r.id + "' has " + convert.toString(len($r.sequence)) +
            " bases but " + convert.toString(len($r.quality)) + " quality characters");
    }
    def buf as bytes init convert.bytesFromString($r.quality, "utf-8");
    def n as int init len($buf);
    def out as list of int;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def score as int init $buf[$i] -$offset;
        if ($score < 0) {
            raiseError("quality character " + describeCode($buf[$i]) + " in read '" + $r.id +
                "' decodes to " + convert.toString($score) + " at offset " +
                convert.toString($offset) + " (wrong Phred offset for this file?)");
        }
        $out[] = $score;
    }
    return $out;
}

/**
 * The mean Phred score across a read, assuming Phred+33.
 *
 * The arithmetic mean of the scores, as quality plots show it. It is not the
 * score of the read's mean error probability: Phred is logarithmic, so the mean
 * score understates the worst bases. Use `errorProbability` on the individual
 * scores when that matters.
 *
 * @param r {Read} the read
 * @return {float} the mean score
 * @throws {Error} kind "bioinformatics" when the read has no bases, its quality and
 *         sequence lengths differ, or a score decodes negative
 */
export func meanQuality(r as Read) {
    return meanQualityWith($r, PHRED33);
}

/**
 * The mean Phred score across a read at a chosen offset.
 *
 * @param r {Read} the read
 * @param offset {int} `PHRED33` or `PHRED64`
 * @return {float} the mean score
 * @throws {Error} kind "bioinformatics" when the read has no bases, its quality and
 *         sequence lengths differ, or a score decodes negative
 */
export func meanQualityWith(r as Read, offset as int) {
    def scores as list of int init qualityScoresWith($r, $offset);
    if (len($scores) == 0) {
        raiseError("read '" + $r.id + "' has no bases to average");
    }
    def total as int init 0;
    for (def s in $scores) {
        $total = $total + $s;
    }
    return $total / len($scores);
}

/**
 * The probability that a base call is wrong, from its Phred score.
 *
 * `P = 10^(-Q/10)`, so Q20 is one error in a hundred and Q30 one in a thousand.
 *
 * @param score {int} a Phred quality score
 * @return {float} the error probability, between 0 and 1
 * @throws {Error} kind "bioinformatics" when the score is negative
 * @example bioseq.errorProbability(30) == 0.001
 */
export func errorProbability(score as int) {
    if ($score < 0) {
        raiseError("a Phred score cannot be negative, got " + convert.toString($score));
    }
    return math.pow(10.0, -$score / 10.0);
}
