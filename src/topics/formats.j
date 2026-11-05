# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

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
    if (not regex.matches('[ \t]', $line)) {
        return $line;
    }
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
    # `indexOf` and string slices both count runes, so a non-ASCII id cuts at
    # the right place.
    def cut as int init strings.indexOf($trimmed, " ");
    def tab as int init strings.indexOf($trimmed, "\t");
    if ($tab >= 0 and ($cut < 0 or $tab < $cut)) {
        $cut = $tab;
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
 * dropped, and lines starting with `;` (comments) are skipped. Only structure is checked: residues are not validated against an
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
    return parseFastaUpTo($text, MAX_RECORDS);
}

/**
 * `parseFasta` with the record cap as a parameter, so a test can reach it.
 *
 * Sequence lines are collected per record and joined once: appending each line
 * to a growing string recopies it and made an 800 kb record take five seconds.
 *
 * @param text {string} the contents of a FASTA file
 * @param limit {int} the most records to build
 * @return {list of Record} the records, in file order
 * @throws {Error} kind "bioinformatics" as for `parseFasta`, past `limit` records
 * @internal
 */
func parseFastaUpTo(text as string, limit as int) {
    def out as list of Record;
    def id as string init "";
    def description as string init "";
    def pieces as list of string;
    def inRecord as bool init false;

    for (def line in splitLines($text)) {
        def trimmed as string init strings.trim($line);
        if (len($trimmed) == 0 or strings.startsWith($trimmed, ";")) {
            continue;
        }
        if (strings.startsWith($line, ">")) {
            if ($inRecord) {
                requireRoom(len($out) + 1, $limit, "FASTA records");
                $out[] = Record{
                    id: $id,
                    description: $description,
                    sequence: strings.join($pieces, "")
                };
            }
            def parts as list of string init splitHeader($line[1..], ">");
            $id = $parts[0];
            $description = $parts[1];
            $pieces = [];
            $inRecord = true;
        } else {
            if (not $inRecord) {
                raiseError("FASTA data before the first '>' header line" +
                    " (is this a plain sequence? use toPlain)");
            }
            $pieces[] = stripSpace($line);
        }
    }
    if ($inRecord) {
        requireRoom(len($out) + 1, $limit, "FASTA records");
        $out[] = Record{id: $id, description: $description, sequence: strings.join($pieces, "")};
    }
    if (len($out) == 0) {
        raiseError("no FASTA records found (a record starts with '>')");
    }
    return $out;
}

/**
 * Render records as FASTA text, wrapped at `FASTA_WIDTH` (60) columns.
 *
 * A record with an empty sequence is written as its header line alone. Fields
 * that would not parse back to the same record raise: see `writeFastaWrapped`.
 *
 * @param records {list of Record} the records to write
 * @return {string} FASTA text, ending in a newline
 * @throws {Error} kind "bioinformatics" when the list is empty or a record has no id,
 *         a field that cannot be written (see `writeFastaWrapped`)
 * @example bioseq.writeFasta([$r]) == ">seq1 first\nATGC\n"
 */
export func writeFasta(records as list of Record) {
    return writeFastaWrapped($records, FASTA_WIDTH);
}

/**
 * Render records as FASTA text at a chosen line width.
 *
 * The id is written trimmed. Output parses back to the same records, so a
 * field that would break that raises: whitespace in the id, a CR or LF in the
 * description, whitespace, `>` or `;` in the sequence (a wrapped line starting
 * with `>` or `;` reads back as a header or a comment).
 *
 * @param records {list of Record} the records to write
 * @param width {int} the wrap column; must be positive
 * @return {string} FASTA text, ending in a newline
 * @throws {Error} kind "bioinformatics" when the list is empty, the width is not
 *         positive, a record has no id, or a field holds a character listed above
 */
export func writeFastaWrapped(records as list of Record, width as int) {
    if (len($records) == 0) {
        raiseError("no records to write");
    }
    if ($width <= 0) {
        raiseError("FASTA line width must be positive, got " + convert.toString($width));
    }
    # Pieces are collected and joined once; `+` onto the growing output
    # recopied it per record (8000 records took eleven seconds).
    def out as list of string;
    for (def r in $records) {
        def id as string init requireWritableHeader($r.id, $r.description, "record", "FASTA");
        if (regex.matches('[\s>;]', $r.sequence)) {
            raiseError("record '" + $id + "' has whitespace, '>' or ';' in its sequence," +
                " which would not parse back as the same record");
        }
        $out[] = ">";
        $out[] = $id;
        if (len($r.description) > 0) {
            $out[] = " ";
            $out[] = $r.description;
        }
        $out[] = "\n";
        if (len($r.sequence) > 0) {
            $out[] = wrap($r.sequence, $width);
            $out[] = "\n";
        }
    }
    return strings.join($out, "");
}

/**
 * Check a record's id and description before writing its header line.
 *
 * @param id {string} the id as given
 * @param description {string} the description
 * @param what {string} "record" or "read", for the message
 * @param format {string} "FASTA" or "FASTQ", for the message
 * @return {string} the trimmed id, the one that is written
 * @throws {Error} kind "bioinformatics" when the id is empty or holds whitespace,
 *         or the description holds a CR or LF
 * @internal
 */
func requireWritableHeader(id as string, description as string, what as string, format as string) {
    def trimmed as string init strings.trim($id);
    if (len($trimmed) == 0) {
        raiseError("a " + $what + " has no id, so it cannot be written as " + $format);
    }
    if (regex.matches('\s', $trimmed)) {
        raiseError("the " + $what + " id " + summarise($trimmed) + " holds whitespace," +
            " which would split it into id and description on parsing");
    }
    if (regex.matches('[\r\n]', $description)) {
        raiseError($what + " '" + $trimmed + "' has a line break in its description");
    }
    return $trimmed;
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
    def n as int init len($s);
    if ($n <= $width) {
        return $s;
    }
    # One native slice per line. ASCII text (every real sequence) is cut on
    # bytes; anything else on runes, so a multi-byte character is never split.
    def lines as list of string;
    def buf as bytes init convert.bytesFromString($s, "utf-8");
    if (len($buf) == $n) {
        for (def i as int init 0; $i < $n; $i = $i + $width) {
            $lines[] = convert.stringFromBytes($buf[$i..math.min($i + $width, $n)], "utf-8");
        }
    } else {
        def chars as list of string init strings.chars($s);
        for (def i as int init 0; $i < $n; $i = $i + $width) {
            $lines[] = strings.join($chars[$i..math.min($i + $width, $n)], "");
        }
    }
    return strings.join($lines, "\n");
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
    def out as list of string;
    for (def line in splitLines($text)) {
        if (strings.startsWith(strings.trim($line), ">")) {
            continue;
        }
        if (strings.startsWith(strings.trim($line), ";")) {
            continue;
        }
        # Most lines hold nothing to drop; keep those whole after one native test.
        if (not regex.matches('[\x00- 0-9]', $line)) {
            $out[] = $line;
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
        $out[] = convert.stringFromBytes($kept, "utf-8");
    }
    return strings.join($out, "");
}

# ------------------------------------------------------------------------ FASTQ

/**
 * Parse FASTQ text into reads.
 *
 * Reads the four-line form and the wrapped form: sequence lines are gathered
 * until a `+` line, then quality lines until the quality is as long as the
 * sequence. Quality is read by length, not by looking for `@`, because `@` is a
 * valid Phred+33 character (Q31). Quality characters must be printable ASCII
 * (`!` to `~`, bytes 33 to 126). At most `MAX_RECORDS` reads are built.
 *
 * @param text {string} the contents of a FASTQ file
 * @return {list of Read} the reads, in file order
 * @throws {Error} kind "bioinformatics" when there is no record, a record does not
 *         start with `@`, lacks its `+` separator, has a quality string of the
 *         wrong length or a quality byte outside 33..126, or past `MAX_RECORDS`
 * @example bioseq.parseFastq("@r1\nATGC\n+\nIIII\n")[0].quality == "IIII"
 */
export func parseFastq(text as string) {
    return parseFastqUpTo($text, MAX_RECORDS);
}

/**
 * `parseFastq` with the record cap as a parameter, so a test can reach it.
 *
 * @param text {string} the contents of a FASTQ file
 * @param limit {int} the most reads to build
 * @return {list of Read} the reads, in file order
 * @throws {Error} kind "bioinformatics" as for `parseFastq`, past `limit` reads
 * @internal
 */
func parseFastqUpTo(text as string, limit as int) {
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

        def seqPieces as list of string;
        while ($i < $n and not strings.startsWith($lines[$i], "+")) {
            $seqPieces[] = stripSpace($lines[$i]);
            $i = $i + 1;
        }
        def sequence as string init strings.join($seqPieces, "");
        if ($i >= $n) {
            raiseError("FASTQ record '" + $parts[0] + "' has no '+' separator line");
        }
        $i = $i + 1;

        # Quality is ASCII once checked, so rune and byte lengths agree.
        def want as int init len($sequence);
        def got as int init 0;
        def qualPieces as list of string;
        while ($i < $n and $got < $want) {
            def piece as string init stripSpace($lines[$i]);
            requireQualityText($piece, $parts[0]);
            $qualPieces[] = $piece;
            $got = $got + len($piece);
            $i = $i + 1;
        }
        def quality as string init strings.join($qualPieces, "");
        if (len($quality) != len($sequence)) {
            raiseError("FASTQ record '" + $parts[0] + "' has " + convert.toString(len($sequence)) +
                " bases but " + convert.toString(len($quality)) + " quality characters");
        }
        requireRoom(len($out) + 1, $limit, "FASTQ reads");
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
 * @param limit {int} the cap, `MAX_RECORDS` outside tests
 * @param what {string} what is being counted, for the message
 * @throws {Error} kind "bioinformatics" when the count exceeds `limit`
 * @internal
 */
func requireRoom(count as int, limit as int, what as string) {
    if ($count > $limit) {
        raiseError("more than " + convert.toString($limit) + " " + $what +
            " in one parse (MAX_RECORDS) - refusing to continue, since a file of bare" +
            " header lines expands about a hundredfold in memory");
    }
}

/**
 * Reject quality text holding a byte outside printable ASCII.
 *
 * Phred+33 and Phred+64 both encode into `!` .. `~`. A byte outside that range
 * is not a score, and a multi-byte character would make the quality's rune
 * length disagree with the number of scores decoded from its bytes.
 *
 * @param quality {string} quality characters
 * @param id {string} the read id, for the message
 * @throws {Error} kind "bioinformatics" on a byte outside 33..126
 * @internal
 */
func requireQualityText(quality as string, id as string) {
    if (regex.matches('^[!-~]*$', $quality)) {
        return;
    }
    def buf as bytes init convert.bytesFromString($quality, "utf-8");
    for (def i as int init 0; $i < len($buf); $i = $i + 1) {
        if ($buf[$i] < 33 or $buf[$i] > 126) {
            raiseError("quality " + describeCode($buf[$i]) + " in read '" + $id +
                "' is not a quality character (printable ASCII '!' to '~', bytes 33 to 126)");
        }
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
 * every current tool emits and what halves the size of a large file. The id is
 * written trimmed. Output parses back to the same reads, so a field that would
 * break that raises: whitespace in the id, a CR or LF in the description,
 * whitespace in the sequence or a leading `+` (read as the separator), a
 * quality byte outside 33..126.
 *
 * @param reads {list of Read} the reads to write
 * @return {string} FASTQ text, ending in a newline
 * @throws {Error} kind "bioinformatics" when the list is empty, a read has no id, a
 *         read's quality length disagrees with its sequence length, or a field holds
 *         a character listed above
 */
export func writeFastq(reads as list of Read) {
    if (len($reads) == 0) {
        raiseError("no reads to write");
    }
    # Pieces are collected and joined once; `+` onto the growing output
    # recopied it per read (8000 reads took 28 seconds).
    def out as list of string;
    for (def r in $reads) {
        def id as string init requireWritableHeader($r.id, $r.description, "read", "FASTQ");
        if (regex.matches('\s', $r.sequence) or strings.startsWith($r.sequence, "+")) {
            raiseError("read '" + $id + "' has whitespace or a leading '+' in its sequence," +
                " which would not parse back as the same read");
        }
        requireQualityText($r.quality, $id);
        if (len($r.sequence) != len($r.quality)) {
            raiseError("read '" + $id + "' has " + convert.toString(len($r.sequence)) +
                " bases but " + convert.toString(len($r.quality)) + " quality characters");
        }
        $out[] = "@";
        $out[] = $id;
        if (len($r.description) > 0) {
            $out[] = " ";
            $out[] = $r.description;
        }
        $out[] = "\n";
        $out[] = $r.sequence;
        $out[] = "\n+\n";
        $out[] = $r.quality;
        $out[] = "\n";
    }
    return strings.join($out, "");
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
 *         a quality byte is outside 33..126, or a score decodes negative
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
 *         a quality byte is outside 33..126, or a quality character decodes to a
 *         negative score (wrong offset)
 */
export func qualityScoresWith(r as Read, offset as int) {
    # A `Read` can arrive from somewhere other than `parseFastq` - built by
    # hand, or carried in from another tool - and the one invariant FASTQ has is
    # that the quality is as long as the sequence. Both the parser and the
    # writer enforce it; without this the quality verbs would quietly report a
    # mean over whichever bases happened to have a score.
    requireQualityText($r.quality, $r.id);
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
 *         sequence lengths differ, a quality byte is outside 33..126, or a score
 *         decodes negative
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
 *         sequence lengths differ, a quality byte is outside 33..126, or a score
 *         decodes negative
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
