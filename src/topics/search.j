# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - search: finding things inside a sequence.
# Spliced into bioinformatics.j via include - not a standalone module.
#
# Offsets are 0-based and half-open, in the coordinates of the sequence the
# caller passed in, including for reverse-strand features (mapped back, and
# tested directly).
#
# Matches overlap: `find("AAA", "AA")` reports two. That is why the IUPAC
# matcher walks the sequence itself instead of delegating to `regex`, whose
# scan is non-overlapping.
#
# Arguments are passed by value, so a helper cannot fill in a list handed to it.
# Every helper here returns what it found and the caller concatenates.
#
# Cost: an exact pattern (67 of the 76 recognition sequences) goes to
# `binary.indexOf`, 200 kb in well under a millisecond. A degenerate pattern is
# expanded into its concrete sequences, each searched natively (`CCWGG` is two
# scans); only one that expands past `MAX_SITE_EXPANSION` falls back to the
# interpreted scan, about 3 microseconds a position. Validation is the floor:
# normalising 200 kb costs about half a second, once per public call.

/**
 * One place a pattern was found.
 *
 * @field start {int} 0-based offset of the first base, in the coordinates of
 *        the sequence searched
 * @field end {int} 0-based offset one past the last base
 * @field text {string} the matched subsequence, as it appears on `strand`
 * @field strand {int} +1 when found on the sequence given, -1 when found on its
 *        reverse complement
 * @field mismatches {int} how many positions differed, always 0 except for
 *        `findFuzzy`
 */
export def struct Match {
    start as int,
    end as int,
    text as string,
    strand as int,
    mismatches as int
};

# ------------------------------------------------------------- exact and fuzzy

/**
 * Find every occurrence of an exact subsequence.
 *
 * Overlapping occurrences are all reported. Only the strand given is searched;
 * reverse-complement the pattern to search the other. Works on protein too.
 *
 * @param s {string} the sequence to search
 * @param pattern {string} the subsequence to look for
 * @return {list of Match} every occurrence, in ascending order of `start`
 * @throws {Error} kind "bioinformatics" when either argument is empty or holds a
 *         non-residue, or the pattern is longer than the sequence
 * @example bioseq.find("AAA", "AA")   # two matches, at 0 and 1
 */
export func find(s as string, pattern as string) {
    def subject as string init requireSequence($s);
    def needle as string init requireSequence($pattern);
    requireFits($subject, $needle);
    return scanExact($subject, $needle, 1);
}

/**
 * Find every occurrence of an IUPAC pattern.
 *
 * The pattern may use ambiguity codes: `GGWCC` matches `GGACC` and `GGTCC`. The
 * relation is one-way: a pattern's `N` admits a subject's `A`, a pattern's `A`
 * does not admit a subject's `N`. `Match.text` is the concrete subsequence
 * found. Overlapping occurrences are all reported.
 *
 * @param s {string} the nucleotide sequence to search
 * @param pattern {string} the pattern, in IUPAC codes
 * @return {list of Match} every occurrence, in ascending order of `start`
 * @throws {Error} kind "bioinformatics" when either argument is empty or not a nucleic
 *         acid, or the pattern is longer than the sequence
 * @example bioseq.findIupac($dna, "GGWCC")
 */
export func findIupac(s as string, pattern as string) {
    def subject as string init requireNucleotide($s);
    def needle as string init requireNucleotide($pattern);
    requireFits($subject, $needle);
    return matchIupac($subject, $needle, 1);
}

/**
 * The IUPAC scan, over already-validated sequences.
 *
 * Each pattern position is expanded once into the bases it admits, then the
 * subject is walked. The expansion is at most four bases, so the inner test is
 * a handful of comparisons and the whole scan is O(n*m) with a tiny m.
 *
 * @param subject {string} a validated nucleotide sequence
 * @param needle {string} a validated IUPAC pattern
 * @param strand {int} the strand label to stamp on each match
 * @return {list of Match} every occurrence
 * @internal
 */
func matchIupac(subject as string, needle as string, strand as int) {
    # An exact pattern is a plain substring, so hand the scan to `binary.indexOf`
    # exactly as `find` does. Sixty-seven of the seventy-six recognition
    # sequences are exact, so this is the path a restriction map normally takes.
    if (isExactSite($needle)) {
        return scanExact($subject, $needle, $strand);
    }
    def expansions as list of bytes;
    for (def code in strings.chars($needle)) {
        $expansions[] = convert.bytesFromString(expandIupac($code), "utf-8");
    }
    def hay as bytes init convert.bytesFromString($subject, "utf-8");
    def m as int init len($needle);
    def anchor as int init anchorByte($expansions[0]);
    def last as int init len($hay) - $m;
    def out as list of Match;
    for (def i as int init 0; $i <= $last; $i = $i + 1) {
        # Reject on the first base with a single integer compare before paying
        # for a call. Nearly every pattern starts with an unambiguous base, so
        # this discards about three quarters of the offsets for the cost of one
        # comparison - and a call per offset was measurably the larger cost.
        if ($anchor >= 0 and $hay[$i] != $anchor) {
            continue;
        }
        if (iupacHitAt($hay, $i, $expansions, $m)) {
            $out[] = Match{
                start: $i,
                end: $i + $m,
                text: convert.stringFromBytes(binary.slice($hay, $i, $i + $m), "utf-8"),
                strand: $strand,
                mismatches: 0
            };
        }
    }
    return $out;
}

/**
 * Find every occurrence of an exact pattern, at native scan speed.
 *
 * The shared engine behind `find` and the exact half of `matchIupac`. Each
 * round searches the remaining suffix and the offset advances one past the hit,
 * so overlapping occurrences are all reported.
 *
 * @param subject {string} a validated sequence
 * @param needle {string} an exact pattern
 * @param strand {int} the strand label to stamp on each match
 * @return {list of Match} every occurrence
 * @internal
 */
func scanExact(subject as string, needle as string, strand as int) {
    def hay as bytes init convert.bytesFromString($subject, "utf-8");
    def pin as bytes init convert.bytesFromString($needle, "utf-8");
    def m as int init len($pin);
    def n as int init len($hay);
    def out as list of Match;
    def offset as int init 0;
    while ($offset + $m <= $n) {
        def rest as bytes init $hay;
        if ($offset > 0) {
            $rest = binary.slice($hay, $offset, $n);
        }
        def at as int init binary.indexOf($rest, $pin);
        if ($at < 0) {
            break;
        }
        def start as int init $offset + $at;
        $out[] = Match{
            start: $start,
            end: $start + $m,
            text: convert.stringFromBytes(binary.slice($hay, $start, $start + $m), "utf-8"),
            strand: $strand,
            mismatches: 0
        };
        $offset = $start + 1;
    }
    return $out;
}

/**
 * Does the pattern match at this offset?
 *
 * @param hay {bytes} the subject
 * @param at {int} the offset to test
 * @param expansions {list of bytes} admitted bases, one entry per pattern position
 * @param m {int} the pattern length
 * @return {bool} true when every position is admitted
 * @internal
 */
func iupacHitAt(hay as bytes, at as int, expansions as list of bytes, m as int) {
    for (def j as int init 0; $j < $m; $j = $j + 1) {
        def allowed as bytes init $expansions[$j];
        def observed as int init $hay[$at + $j];
        def ok as bool init false;
        def choices as int init len($allowed);
        for (def k as int init 0; $k < $choices; $k = $k + 1) {
            if ($allowed[$k] == $observed) {
                $ok = true;
                break;
            }
        }
        if (not $ok) {
            return false;
        }
    }
    return true;
}

/**
 * Find every place a pattern occurs with at most a given number of mismatches.
 *
 * Substitutions only, no insertions or deletions. Each match carries its
 * `mismatches` count; matches come in ascending order of `start`, exact and
 * inexact together. Overlapping matches are all reported. Residues are compared
 * literally (no IUPAC). Use `editDistance` when indels matter.
 *
 * @param s {string} the sequence to search
 * @param pattern {string} the subsequence to look for
 * @param maxMismatches {int} how many differing positions to tolerate
 * @return {list of Match} every occurrence within the tolerance, by `start`
 * @throws {Error} kind "bioinformatics" when an argument is empty or holds a
 *         non-residue, the pattern is longer than the sequence, or the
 *         tolerance is negative
 * @example bioseq.findFuzzy($dna, "GAATTC", 1)
 */
export func findFuzzy(s as string, pattern as string, maxMismatches as int) {
    def subject as string init requireSequence($s);
    def needle as string init requireSequence($pattern);
    requireFits($subject, $needle);
    if ($maxMismatches < 0) {
        raiseError("mismatch tolerance cannot be negative, got " + convert.toString($maxMismatches));
    }
    def hay as bytes init convert.bytesFromString($subject, "utf-8");
    def pin as bytes init convert.bytesFromString($needle, "utf-8");
    def m as int init len($pin);
    def last as int init len($hay) - $m;
    def out as list of Match;
    for (def i as int init 0; $i <= $last; $i = $i + 1) {
        def misses as int init countMismatches($hay, $i, $pin, $m, $maxMismatches);
        if ($misses <= $maxMismatches) {
            $out[] = Match{
                start: $i,
                end: $i + $m,
                text: convert.stringFromBytes(binary.slice($hay, $i, $i + $m), "utf-8"),
                strand: 1,
                mismatches: $misses
            };
        }
    }
    return $out;
}

/**
 * Count differing positions, giving up once the tolerance is exceeded.
 *
 * @param hay {bytes} the subject
 * @param at {int} the offset to compare from
 * @param pin {bytes} the pattern
 * @param m {int} the pattern length
 * @param limit {int} stop counting once this many mismatches are seen
 * @return {int} the mismatch count, or one past `limit` when it was exceeded
 * @internal
 */
func countMismatches(hay as bytes, at as int, pin as bytes, m as int, limit as int) {
    def misses as int init 0;
    for (def j as int init 0; $j < $m; $j = $j + 1) {
        if ($hay[$at + $j] != $pin[$j]) {
            $misses = $misses + 1;
            if ($misses > $limit) {
                return $misses;
            }
        }
    }
    return $misses;
}

/**
 * Reject a pattern longer than the sequence it is searched in.
 *
 * @param subject {string} the sequence
 * @param needle {string} the pattern
 * @throws {Error} kind "bioinformatics" when the pattern cannot fit
 * @internal
 */
func requireFits(subject as string, needle as string) {
    if (len($needle) > len($subject)) {
        raiseError("pattern is " + convert.toString(len($needle)) + " residues, longer than the " +
            convert.toString(len($subject)) + "-residue sequence");
    }
}

# ------------------------------------------------------------ open reading frames

/**
 * One open reading frame.
 *
 * @field frame {int} 1, 2 or 3 on the forward strand; -1, -2 or -3 on the reverse
 * @field strand {int} +1 or -1, the sign of `frame`
 * @field start {int} 0-based offset of the first base of the start codon, in
 *        the coordinates of the sequence searched
 * @field end {int} 0-based offset one past the last base, including the stop
 *        codon when there is one
 * @field length {int} bases spanned, `end` less `start`
 * @field nucleotides {string} the frame's sequence 5' to 3' as read, DNA-spelled,
 *        reverse-complemented for a reverse-strand frame
 * @field protein {string} the translation, without the trailing stop
 * @field stopped {bool} true when a stop codon ended the frame, false when it ran
 *        off the end of the sequence
 */
export def struct Orf {
    frame as int,
    strand as int,
    start as int,
    end as int,
    length as int,
    nucleotides as string,
    protein as string,
    stopped as bool
};

/**
 * Find open reading frames in all six frames, under the standard genetic code.
 *
 * @param nt {string} a DNA or RNA sequence
 * @param minAa {int} the shortest protein to report, in residues
 * @return {list of Orf} frames 1, 2, 3, then -1, -2, -3
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic
 *         acid, or `minAa` is negative
 */
export func orfs(nt as string, minAa as int) {
    return orfsWith($nt, $minAa, 1);
}

/**
 * Find open reading frames under a named genetic code.
 *
 * The code is an NCBI table id, not a codon map, because the search needs the
 * table's start codons: `BACTERIAL` (11) assigns codons as `STANDARD` (1) does
 * but has seven starts instead of three.
 *
 * Within a frame, an ORF begins at the first start codon and runs to the first
 * in-frame stop; the next begins at the next start after that stop. ORFs never
 * overlap within a frame, and nested start codons are not reported separately.
 *
 * @param nt {string} a DNA or RNA sequence
 * @param minAa {int} the shortest protein to report, in residues
 * @param tableId {int} an NCBI genetic code id
 * @return {list of Orf} frames 1, 2, 3, then -1, -2, -3
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic
 *         acid, `minAa` is negative, or the table id is unknown
 * @example bioseq.orfsWith($genome, 50, 11)   # bacterial starts
 */
export func orfsWith(nt as string, minAa as int, tableId as int) {
    def clean as string init requireNucleotide($nt);
    if ($minAa < 0) {
        raiseError("minimum ORF length cannot be negative, got " + convert.toString($minAa));
    }
    def table as map of string to string init codonTable($tableId);
    def starts as list of string init startCodons($tableId);
    def startSet as map of string to bool;
    for (def codon in $starts) {
        $startSet[$codon] = true;
    }
    def forward as string init strings.replace($clean, "U", "T");
    def reverse as string init reverseComplement($forward);
    def out as list of Orf;
    for (def f as int init 1; $f <= 3; $f = $f + 1) {
        for (def o in collectOrfs($forward, $f, $f, $startSet, $table, $minAa)) {
            $out[] = $o;
        }
    }
    for (def f as int init 1; $f <= 3; $f = $f + 1) {
        for (def o in collectOrfs($reverse, $f, -$f, $startSet, $table, $minAa)) {
            $out[] = $o;
        }
    }
    return $out;
}

/**
 * Scan one reading frame of one strand, appending what it finds.
 *
 * `strandSequence` is already oriented 5' to 3' for the strand in question, so
 * the scan is identical either way; only the final coordinate mapping differs.
 *
 * @param strandSequence {string} the strand, oriented for reading
 * @param offset {int} 1, 2 or 3 - how many bases to skip
 * @param frame {int} the frame number to record, signed
 * @param startSet {map of string to bool} the table's start codons
 * @param table {map of string to string} the genetic code
 * @param minAa {int} the shortest protein to report
 * @return {list of Orf} what this frame yielded
 * @internal
 */
func collectOrfs(
    strandSequence as string,
    offset as int,
    frame as int,
    startSet as map of string to bool,
    table as map of string to string,
    minAa as int) {
    def out as list of Orf;
    def strandBytes as bytes init convert.bytesFromString($strandSequence, "utf-8");
    def n as int init len($strandBytes);
    def openAt as int init -1;
    def i as int init $offset - 1;
    while ($i + 3 <= $n) {
        # cut the codon from the byte buffer: a string slice is rune-indexed and
        # would cost O(i), turning the frame scan quadratic
        def kb as bytes;
        $kb[] = $strandBytes[$i];
        $kb[] = $strandBytes[$i + 1];
        $kb[] = $strandBytes[$i + 2];
        def codon as string init convert.stringFromBytes($kb, "utf-8");
        if ($openAt < 0) {
            if (maps.has($startSet, $codon)) {
                $openAt = $i;
            }
        } elseif (maps.has($table, $codon) and $table[$codon] == "*") {
            for (def o in recordOrf($strandSequence, $openAt, $i + 3, $frame, $table, $minAa, true)) {
                $out[] = $o;
            }
            $openAt = -1;
        }
        $i = $i + 3;
    }
    # a frame that runs off the end without a stop is still worth reporting: on
    # a fragment it usually means the gene continues past the sequence you have
    if ($openAt >= 0) {
        def last as int init $openAt + ($n - $openAt) // 3 * 3;
        for (def o in recordOrf($strandSequence, $openAt, $last, $frame, $table, $minAa, false)) {
            $out[] = $o;
        }
    }
    return $out;
}

/**
 * Build one `Orf`, mapping reverse-strand coordinates back to the original.
 *
 * @param strandSequence {string} the strand as read
 * @param from {int} start offset within `strandSequence`
 * @param upto {int} end offset within `strandSequence`, half-open
 * @param frame {int} the signed frame number
 * @param table {map of string to string} the genetic code
 * @param minAa {int} the shortest protein to report
 * @param stopped {bool} whether a stop codon closed the frame
 * @return {list of Orf} one record, or none when the ORF was filtered out
 * @internal
 */
func recordOrf(
    strandSequence as string,
    from as int,
    upto as int,
    frame as int,
    table as map of string to string,
    minAa as int,
    stopped as bool) {
    def out as list of Orf;
    def nucleotides as string init $strandSequence[$from..$upto];
    if (len($nucleotides) < 3) {
        return $out;
    }
    def protein as string init translateClean($nucleotides, $table);
    if ($stopped and len($protein) > 0) {
        $protein = $protein[0..len($protein) - 1];
    }
    if (len($protein) < $minAa) {
        return $out;
    }
    def strand as int init 1;
    def start as int init $from;
    def finish as int init $upto;
    if ($frame < 0) {
        # a position p on the reverse complement of an n-base sequence is
        # position n-1-p on the original, so the half-open span [from, to)
        # becomes [n-to, n-from)
        $strand = -1;
        def n as int init len($strandSequence);
        $start = $n - $upto;
        $finish = $n - $from;
    }
    $out[] = Orf{
        frame: $frame,
        strand: $strand,
        start: $start,
        end: $finish,
        length: $upto - $from,
        nucleotides: $nucleotides,
        protein: $protein,
        stopped: $stopped
    };
    return $out;
}

# ------------------------------------------------------------------ CpG islands

/**
 * A stretch of sequence singled out by a scan.
 *
 * @field start {int} 0-based offset of the first base
 * @field end {int} 0-based offset one past the last base
 * @field length {int} bases spanned
 * @field gc {float} GC percentage across the region
 * @field obsExp {float} observed CpG dinucleotides over the number expected
 *        from the region's C and G counts
 */
export def struct Region {
    start as int,
    end as int,
    length as int,
    gc as float,
    obsExp as float
};

/**
 * Find CpG islands by the Gardiner-Garden and Frommer criteria.
 *
 * A window qualifies when its GC percentage is at least `minGc` and its observed
 * over expected CpG is at least `minObsExp`, expected being
 * `(C count * G count) / length`. The classic values are 200, 50.0 and 0.6.
 * Qualifying windows that touch or overlap are merged, and `gc` and `obsExp`
 * are recomputed over each merged region. A sequence shorter than one window
 * gives an empty list.
 *
 * @param s {string} a DNA sequence
 * @param window {int} the scanning window in bases; must be positive
 * @param minGc {float} the minimum GC percentage, 0 to 100
 * @param minObsExp {float} the minimum observed-over-expected CpG ratio
 * @return {list of Region} the islands, in ascending order of `start`
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not DNA, or a
 *         parameter is out of range
 * @see Gardiner-Garden M, Frommer M (1987) J Mol Biol 196:261-282
 * @example bioseq.cpgIslands($promoter, 200, 50.0, 0.6)
 */
export func cpgIslands(s as string, window as int, minGc as float, minObsExp as float) {
    def clean as string init requireDna($s);
    if ($window <= 0) {
        raiseError("window must be positive, got " + convert.toString($window));
    }
    if ($minGc < 0.0 or $minGc > 100.0) {
        raiseError("minimum GC must be a percentage between 0 and 100");
    }
    if ($minObsExp < 0.0) {
        raiseError("minimum observed/expected ratio cannot be negative");
    }
    def n as int init len($clean);
    def out as list of Region;
    if ($n < $window) {
        return $out;
    }
    def buf as bytes init convert.bytesFromString($clean, "utf-8");

    # prefix sums, so each window's counts are read in constant time rather than
    # recounted - the difference between a linear scan and a quadratic one
    def cAt as list of int;
    def gAt as list of int;
    def cgAt as list of int;
    $cAt[] = 0;
    $gAt[] = 0;
    $cgAt[] = 0;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def isC as int init 0;
        def isG as int init 0;
        if ($buf[$i] == 67) {
            $isC = 1;
        }
        if ($buf[$i] == 71) {
            $isG = 1;
        }
        def isCg as int init 0;
        if ($i > 0 and $buf[$i - 1] == 67 and $buf[$i] == 71) {
            $isCg = 1;
        }
        $cAt[] = $cAt[$i] + $isC;
        $gAt[] = $gAt[$i] + $isG;
        $cgAt[] = $cgAt[$i] + $isCg;
    }

    # collect qualifying windows as a merged span, then describe each span once
    def spanStart as int init -1;
    def spanEnd as int init -1;
    for (def start as int init 0; $start + $window <= $n; $start = $start + 1) {
        def finish as int init $start + $window;
        if (windowQualifies($cAt, $gAt, $cgAt, $start, $finish, $minGc, $minObsExp)) {
            if ($spanStart < 0) {
                $spanStart = $start;
                $spanEnd = $finish;
            } elseif ($start <= $spanEnd) {
                $spanEnd = $finish;
            } else {
                $out[] = describeRegion($cAt, $gAt, $cgAt, $spanStart, $spanEnd);
                $spanStart = $start;
                $spanEnd = $finish;
            }
        }
    }
    if ($spanStart >= 0) {
        $out[] = describeRegion($cAt, $gAt, $cgAt, $spanStart, $spanEnd);
    }
    return $out;
}

/**
 * Does the half-open window meet both thresholds?
 *
 * @param cAt {list of int} prefix sums of C
 * @param gAt {list of int} prefix sums of G
 * @param cgAt {list of int} prefix sums of the CG dinucleotide
 * @param from {int} window start
 * @param upto {int} window end, half-open
 * @param minGc {float} minimum GC percentage
 * @param minObsExp {float} minimum observed/expected ratio
 * @return {bool} true when the window qualifies
 * @internal
 */
func windowQualifies(
    cAt as list of int,
    gAt as list of int,
    cgAt as list of int,
    from as int,
    upto as int,
    minGc as float,
    minObsExp as float) {
    def width as int init $upto - $from;
    def c as int init $cAt[$upto] -$cAt[$from];
    def g as int init $gAt[$upto] -$gAt[$from];
    if (100.0 * ($c + $g) / $width < $minGc) {
        return false;
    }
    if ($c == 0 or $g == 0) {
        return false;
    }
    # a CG straddling the window's left edge belongs to the previous window
    def observed as int init $cgAt[$upto] -$cgAt[$from + 1];
    def expected as float init 1.0 * $c * $g / $width;
    return $observed / $expected >= $minObsExp;
}

/**
 * Describe a merged span as a `Region`.
 *
 * @param cAt {list of int} prefix sums of C
 * @param gAt {list of int} prefix sums of G
 * @param cgAt {list of int} prefix sums of the CG dinucleotide
 * @param from {int} region start
 * @param upto {int} region end, half-open
 * @return {Region} the described region
 * @internal
 */
func describeRegion(
    cAt as list of int,
    gAt as list of int,
    cgAt as list of int,
    from as int,
    upto as int) {
    def width as int init $upto - $from;
    def c as int init $cAt[$upto] -$cAt[$from];
    def g as int init $gAt[$upto] -$gAt[$from];
    def observed as int init $cgAt[$upto] -$cgAt[$from + 1];
    def ratio as float init 0.0;
    if ($c > 0 and $g > 0) {
        $ratio = $observed / (1.0 * $c * $g / $width);
    }
    return Region{
        start: $from,
        end: $upto,
        length: $width,
        gc: 100.0 * ($c + $g) / $width,
        obsExp: $ratio
    };
}

# ------------------------------------------------------------ restriction sites

/**
 * One place an enzyme recognises.
 *
 * @field enzyme {string} the enzyme's name
 * @field start {int} 0-based offset of the first base of the recognition site
 * @field end {int} 0-based offset one past the site
 * @field cut {int} 0-based offset at which the top strand is cut
 * @field site {string} the sequence found, which for a degenerate
 *        recognition sequence is the concrete match rather than the pattern
 * @field strand {int} +1, or -1 for a site found only on the reverse strand
 */
export def struct Site {
    enzyme as string,
    start as int,
    end as int,
    cut as int,
    site as string,
    strand as int
};

/**
 * Find where the named enzymes cut a sequence.
 *
 * Every recognition sequence in the table is palindromic, so one forward scan
 * finds every site. A non-palindromic site would also be searched on the reverse
 * strand and reported with `strand: -1`. The sequence is treated as linear: a
 * site spanning the origin of a circular sequence is not found.
 *
 * @param s {string} a DNA sequence
 * @param names {list of string} enzyme names, any case
 * @return {list of Site} every site, ordered by `cut`, then enzyme name
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not DNA, or an
 *         enzyme is not in the table
 * @example bioseq.restrictionSites($plasmid, ["EcoRI", "BamHI"])
 */
export func restrictionSites(s as string, names as list of string) {
    def clean as string init requireDna($s);
    def found as list of Site;
    for (def name in $names) {
        def e as Enzyme init enzyme($name);
        if (len($e.site) > len($clean)) {
            continue;
        }
        for (def m in matchIupac($clean, $e.site, 1)) {
            $found[] = siteFrom($e, $m);
        }
        def flipped as string init reverseComplement($e.site);
        if ($flipped != $e.site) {
            for (def m in matchIupac($clean, $flipped, -1)) {
                $found[] = siteFrom($e, $m);
            }
        }
    }
    return sortSites($found);
}

/**
 * Build a `Site` from an enzyme and the match that located it.
 *
 * @param e {Enzyme} the enzyme
 * @param m {Match} where its recognition sequence was found
 * @return {Site} the site record
 * @internal
 */
func siteFrom(e as Enzyme, m as Match) {
    return Site{
        enzyme: $e.name,
        start: $m.start,
        end: $m.end,
        cut: $m.start + $e.cutTop,
        site: $m.text,
        strand: $m.strand
    };
}

/**
 * Order sites by cut position, then by enzyme name, via `mergeSortSites`.
 *
 * An earlier insertion sort rebuilt the list on every insertion; an 8 kb digest
 * with the full enzyme table spent about four seconds in it.
 *
 * @param sites {list of Site} the sites to order
 * @return {list of Site} the same sites, ordered
 * @internal
 */
func sortSites(sites as list of Site) {
    if (len($sites) < 2) {
        return $sites;
    }
    return mergeSortSites($sites);
}

/**
 * Merge sort over sites, so ordering stays O(n log n) in the number of sites.
 *
 * A digest of a long sequence with the whole enzyme table produces thousands of
 * sites; the quadratic insertion sort this replaced dominated the digest.
 *
 * @param sites {list of Site} the sites to order
 * @return {list of Site} the same sites, ordered
 * @internal
 */
func mergeSortSites(sites as list of Site) {
    def total as int init len($sites);
    if ($total < 2) {
        return $sites;
    }
    def middle as int init $total // 2;
    def left as list of Site;
    def right as list of Site;
    for (def i as int init 0; $i < $middle; $i = $i + 1) {
        $left[] = $sites[$i];
    }
    for (def i as int init $middle; $i < $total; $i = $i + 1) {
        $right[] = $sites[$i];
    }
    return mergeSiteRuns(mergeSortSites($left), mergeSortSites($right));
}

/**
 * Merge two ordered runs of sites into one.
 *
 * @param left {list of Site} an ordered run
 * @param right {list of Site} another ordered run
 * @return {list of Site} the merged run
 * @internal
 */
func mergeSiteRuns(left as list of Site, right as list of Site) {
    def out as list of Site;
    def i as int init 0;
    def j as int init 0;
    while ($i < len($left) and $j < len($right)) {
        if (sitePrecedes($right[$j], $left[$i])) {
            $out[] = $right[$j];
            $j = $j + 1;
        } else {
            $out[] = $left[$i];
            $i = $i + 1;
        }
    }
    while ($i < len($left)) {
        $out[] = $left[$i];
        $i = $i + 1;
    }
    while ($j < len($right)) {
        $out[] = $right[$j];
        $j = $j + 1;
    }
    return $out;
}

/**
 * Does one site sort before another?
 *
 * By cut position, then by enzyme name so the order is total and the output is
 * reproducible whatever order the enzymes were asked for in.
 *
 * @param a {Site} the site being placed
 * @param b {Site} the site already in the list
 * @return {bool} true when `a` belongs before `b`
 * @internal
 */
func sitePrecedes(a as Site, b as Site) {
    if ($a.cut != $b.cut) {
        return $a.cut < $b.cut;
    }
    return $a.enzyme < $b.enzyme;
}

/**
 * The enzymes in the table that have at least one site in the sequence.
 *
 * Faster than `restrictionSites` over the whole table: it stops at the first
 * site of each enzyme.
 *
 * @param s {string} a DNA sequence
 * @return {list of string} enzyme names, alphabetical
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not DNA
 */
export func enzymesCutting(s as string) {
    def clean as string init requireDna($s);
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def out as list of string;
    for (def e in enzymes()) {
        if (len($e.site) > len($clean)) {
            continue;
        }
        # The exact test is written out here rather than called, because passing
        # the sequence into a helper copies it - Jennifer arguments are values -
        # and doing that once per enzyme dominated everything else. Sixty-seven
        # of the seventy-six sites take this path; only the degenerate ones pay
        # for the call below.
        if (isExactSite($e.site)) {
            if (binary.indexOf($buf, convert.bytesFromString($e.site, "utf-8")) >= 0) {
                $out[] = $e.name;
            }
        } elseif (hasSite($buf, $e.site)) {
            $out[] = $e.name;
        }
    }
    return $out;
}

/**
 * Does this recognition sequence occur at all?
 *
 * The presence-only counterpart of `matchIupac`, and the reason
 * `enzymesCutting` is usable on a real sequence. It stops at the first hit and
 * never builds a `Match`, which for a table of 76 enzymes over a few kilobases
 * is the difference between seconds and milliseconds.
 *
 * @param hay {bytes} the sequence to search, already validated
 * @param needle {string} an IUPAC recognition sequence
 * @return {bool} true when the site occurs at least once
 * @internal
 */
func hasSite(hay as bytes, needle as string) {
    # Sixty-seven of the seventy-six enzymes recognise an exact sequence, and
    # for those the whole question is one native substring search. Only the nine
    # degenerate sites need the interpreted scan below.
    if (isExactSite($needle)) {
        return binary.indexOf($hay, convert.bytesFromString($needle, "utf-8")) >= 0;
    }
    # A degenerate site stands for a small set of concrete ones - `CCWGG` is two
    # sequences, `GANTC` four - and each of those can go to the native search
    # too. Nine of the deck's enzymes take this path, and before it existed they
    # accounted for the entire cost of scanning the table.
    def concrete as list of string init expandSite($needle, MAX_SITE_EXPANSION);
    if (len($concrete) > 0) {
        for (def pattern in $concrete) {
            if (binary.indexOf($hay, convert.bytesFromString($pattern, "utf-8")) >= 0) {
                return true;
            }
        }
        return false;
    }
    def expansions as list of bytes;
    for (def code in strings.chars($needle)) {
        $expansions[] = convert.bytesFromString(expandIupac($code), "utf-8");
    }
    def m as int init len($needle);
    def anchor as int init anchorByte($expansions[0]);
    def last as int init len($hay) - $m;
    for (def i as int init 0; $i <= $last; $i = $i + 1) {
        if ($anchor >= 0 and $hay[$i] != $anchor) {
            continue;
        }
        if (iupacHitAt($hay, $i, $expansions, $m)) {
            return true;
        }
    }
    return false;
}

/**
 * How many concrete sequences a degenerate site may expand into before the
 * expansion is abandoned as not worth it.
 *
 * `GAANNNNTTC` is the widest site in the table at 256, so the cap only bites
 * for something far more ambiguous than any real recognition sequence - which
 * then falls back to the interpreted scan rather than allocating a huge list.
 */
def const MAX_SITE_EXPANSION as int init 1024;

/**
 * Every concrete sequence a degenerate pattern admits.
 *
 * @param needle {string} the pattern, in IUPAC codes
 * @param maxPatterns {int} give up beyond this many
 * @return {list of string} the concrete sequences, or empty when the expansion
 *         would exceed `maxPatterns`
 * @internal
 */
func expandSite(needle as string, maxPatterns as int) {
    def out as list of string;
    $out[] = "";
    for (def code in strings.chars($needle)) {
        def bases as string init expandIupac($code);
        if (len($out) * len($bases) > $maxPatterns) {
            def abandoned as list of string;
            return $abandoned;
        }
        def next as list of string;
        for (def prefix in $out) {
            for (def base in strings.chars($bases)) {
                $next[] = $prefix + $base;
            }
        }
        $out = $next;
    }
    return $out;
}

/**
 * Is every position of this pattern an unambiguous base?
 *
 * When it is, the pattern is a plain substring and the search can be handed to
 * `binary.indexOf` rather than walked here.
 *
 * @param needle {string} the pattern
 * @return {bool} true when the pattern holds only A, C, G, T or U
 * @internal
 */
func isExactSite(needle as string) {
    def buf as bytes init convert.bytesFromString($needle, "utf-8");
    def n as int init len($buf);
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def c as int init $buf[$i];
        if ($c != 65 and $c != 67 and $c != 71 and $c != 84 and $c != 85) {
            return false;
        }
    }
    return true;
}

/**
 * The one byte a pattern's first position demands, or -1 when it admits several.
 *
 * An unambiguous first base turns the scan's rejection test into a single
 * integer compare. An ambiguous one (`RCATGY`, `YACGTR`) gives up the fast path
 * and falls through to the full check, which is correct if slower - and only
 * two enzymes in the table are in that position.
 *
 * @param expansion {bytes} the bases the first pattern position admits
 * @return {int} the sole admitted byte, or -1 when there is more than one
 * @internal
 */
func anchorByte(expansion as bytes) {
    if (len($expansion) == 1) {
        return $expansion[0];
    }
    return -1;
}

/**
 * The enzymes in the table with no site in the sequence: the complement of
 * `enzymesCutting`.
 *
 * @param s {string} a DNA sequence
 * @return {list of string} enzyme names, alphabetical
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not DNA
 */
export func nonCutters(s as string) {
    def cutting as map of string to bool;
    for (def name in enzymesCutting($s)) {
        $cutting[$name] = true;
    }
    def out as list of string;
    for (def e in enzymes()) {
        if (not maps.has($cutting, $e.name)) {
            $out[] = $e.name;
        }
    }
    return $out;
}

# ----------------------------------------------------------------- the digest

/**
 * One piece of a digested sequence.
 *
 * @field start {int} 0-based offset of the fragment's first base in the
 *        original sequence
 * @field end {int} 0-based offset one past its last base; for the fragment that
 *        wraps the origin of a circular sequence this is at most `start`
 * @field length {int} the fragment's length in bases
 * @field sequence {string} the fragment itself
 */
export def struct Fragment {
    start as int,
    end as int,
    length as int,
    sequence as string
};

/**
 * Cut a sequence with the named enzymes and return the pieces.
 *
 * Cuts are made on the top strand at each enzyme's `cutTop` offset; the bottom
 * strand's stagger is described by `Enzyme.overhang`, not reflected in the
 * lengths. Coinciding cuts count once, and a cut at offset 0 or at the end is
 * ignored. Sites are found as by `restrictionSites`, so on a circular sequence a
 * site spanning the origin is missed.
 *
 * Linear: no site gives the whole sequence, n cuts give n + 1 fragments.
 * Circular: no site gives one fragment `0..len`; n cuts give n fragments, the
 * last running from the last cut through the origin, so its `end` is not greater
 * than its `start` (equal when there is one cut, the linearised plasmid).
 *
 * @param s {string} a DNA sequence
 * @param names {list of string} enzyme names, any case
 * @param circular {bool} true to treat the sequence as a closed circle
 * @return {list of Fragment} the pieces, ordered by their start position
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not DNA, or an
 *         enzyme is not in the table
 * @example bioseq.restrictionDigest($plasmid, ["EcoRI"], true)
 */
export func restrictionDigest(s as string, names as list of string, circular as bool) {
    def clean as string init requireDna($s);
    def n as int init len($clean);
    def sites as list of Site init restrictionSites($clean, $names);

    def cuts as list of int;
    def seen as map of int to bool;
    for (def site in $sites) {
        if ($site.cut > 0 and $site.cut < $n and not maps.has($seen, $site.cut)) {
            $seen[$site.cut] = true;
            $cuts[] = $site.cut;
        }
    }

    def out as list of Fragment;
    if (len($cuts) == 0) {
        $out[] = Fragment{start: 0, end: $n, length: $n, sequence: $clean};
        return $out;
    }
    if (not $circular) {
        def from as int init 0;
        for (def cut in $cuts) {
            $out[] = Fragment{
                start: $from,
                end: $cut,
                length: $cut - $from,
                sequence: $clean[$from..$cut]
            };
            $from = $cut;
        }
        $out[] = Fragment{start: $from, end: $n, length: $n - $from, sequence: $clean[$from..$n]};
        return $out;
    }
    # circular: each cut runs to the next, and the last wraps to the first
    for (def i as int init 0; $i < len($cuts); $i = $i + 1) {
        def from as int init $cuts[$i];
        if ($i + 1 < len($cuts)) {
            def upto as int init $cuts[$i + 1];
            $out[] = Fragment{
                start: $from,
                end: $upto,
                length: $upto - $from,
                sequence: $clean[$from..$upto]
            };
        } else {
            def upto as int init $cuts[0];
            $out[] = Fragment{
                start: $from,
                end: $upto,
                length: $n - $from + $upto,
                sequence: $clean[$from..$n] + $clean[0..$upto]
            };
        }
    }
    return $out;
}
