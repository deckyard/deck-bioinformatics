# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - compare: pairwise alignment and distance.
# Spliced into bioinformatics.j via include - not a standalone module.
#
# The deck's one O(n*m) surface. The DP matrix is a flat `list of int` indexed
# `i * width + j`, and the substitution scores for the two sequences are
# precomputed into a map keyed by the two bytes, so the inner loop builds no
# string and calls no helper.
#
# Measured (examples/bench.j): about 380 000 cells a second; 200 x 200 about an
# eighth of a second, 400 x 400 under half a second. Memory binds first: about
# one kilobyte per cell, because every entry is a boxed interpreter value.
# 700 x 700 measured 570 MB resident, so 1000 x 1000 is about a gigabyte and
# 2000 x 2000 an out-of-memory kill. Hence `AlignOptions.maxCells`, on by default. Whole-genome
# throughput is out of scope.
#
# Gap scoring is linear; affine gaps are not implemented.

/**
 * The default `AlignOptions.maxCells`: the largest DP matrix an alignment builds.
 *
 * The matrix has `(len(a) + 1) * (len(b) + 1)` cells, so 999 x 999 residues is
 * the largest square alignment under this limit. At about a kilobyte per cell
 * that is about a gigabyte resident.
 */
export def const DEFAULT_MAX_ALIGNMENT_CELLS as int init 1000000;

/**
 * How two residues are scored against each other, and what a gap costs.
 *
 * @field matchScore {int} score for identical residues, used when `matrix` is
 *        empty
 * @field mismatchScore {int} score for differing residues, likewise
 * @field gapScore {int} score for one gap position, normally negative
 * @field matrix {map of string to int} a substitution matrix such as
 *        `BLOSUM62`; when non-empty it replaces `matchScore` and
 *        `mismatchScore` entirely
 * @field maxCells {int} largest DP matrix to build, `(len(a) + 1) * (len(b) + 1)`
 *        cells; 0 or less removes the limit (about 1 KB per cell)
 */
export def struct AlignOptions {
    matchScore as int,
    mismatchScore as int,
    gapScore as int,
    matrix as map of string to int,
    maxCells as int
};

/**
 * Scoring for DNA or RNA: +2 match, -1 mismatch, -2 per gap position, no matrix,
 * `maxCells` = `DEFAULT_MAX_ALIGNMENT_CELLS`.
 *
 * @return {AlignOptions} the defaults
 */
export func dnaOptions() {
    def none as map of string to int;
    return AlignOptions{
        matchScore: 2,
        mismatchScore: -1,
        gapScore: -2,
        matrix: $none,
        maxCells: DEFAULT_MAX_ALIGNMENT_CELLS
    };
}

/**
 * Scoring for protein: `BLOSUM62`, -8 per gap position,
 * `maxCells` = `DEFAULT_MAX_ALIGNMENT_CELLS`. `matchScore` and `mismatchScore`
 * are set to 1 and -1 but unused while the matrix is present.
 *
 * @return {AlignOptions} the defaults
 */
export func proteinOptions() {
    return AlignOptions{
        matchScore: 1,
        mismatchScore: -1,
        gapScore: -8,
        matrix: BLOSUM62,
        maxCells: DEFAULT_MAX_ALIGNMENT_CELLS
    };
}

/**
 * Look one residue pair up in a substitution matrix.
 *
 * Tries `a + b`, then `b + a`, because `BLOSUM62` stores one order of each pair.
 * A pair the matrix does not cover scores the matrix's `XX` entry (-1 in
 * `BLOSUM62`), or 0 if it has none; it does not raise.
 *
 * @param matrix {map of string to int} the matrix, e.g. `BLOSUM62`
 * @param a {string} one residue
 * @param b {string} the other
 * @return {int} the substitution score
 * @example bioseq.substitutionScore(bioseq.BLOSUM62, "W", "W") == 11
 */
export func substitutionScore(matrix as map of string to int, a as string, b as string) {
    def forward as string init $a + $b;
    if (maps.has($matrix, $forward)) {
        return $matrix[$forward];
    }
    def backward as string init $b + $a;
    if (maps.has($matrix, $backward)) {
        return $matrix[$backward];
    }
    if (maps.has($matrix, "XX")) {
        return $matrix["XX"];
    }
    return 0;
}

/**
 * One pairwise alignment.
 *
 * @field a {string} the first sequence with gaps inserted, as aligned
 * @field b {string} the second, the same length as `a`
 * @field markup {string} one character per column: `|` identical, `:` different
 *        but scoring positively, space otherwise (gaps included)
 * @field score {int} the alignment's total score under the options used
 * @field length {int} columns in the alignment, gaps included
 * @field identities {int} columns where the residues are identical
 * @field similarities {int} columns scoring positively, identities included
 * @field gaps {int} columns holding a gap in one sequence or the other
 * @field startA {int} 0-based offset in the first sequence where the alignment
 *        begins, 0 for a global alignment
 * @field startB {int} the same for the second sequence
 */
export def struct Alignment {
    a as string,
    b as string,
    markup as string,
    score as int,
    length as int,
    identities as int,
    similarities as int,
    gaps as int,
    startA as int,
    startB as int
};

/**
 * Align two sequences end to end - the Needleman-Wunsch algorithm.
 *
 * Every residue of both sequences appears in the result. Terminal gaps are
 * penalised like any other, so for a short sequence against a long one use
 * `alignLocal`.
 *
 * @param a {string} the first sequence
 * @param b {string} the second sequence
 * @param opts {AlignOptions} scoring, from `dnaOptions` or `proteinOptions`
 * @return {Alignment} the alignment
 * @throws {Error} kind "bioinformatics" when either sequence is empty or holds a
 *         non-residue, or the matrix would exceed `opts.maxCells`
 * @see Needleman SB, Wunsch CD (1970) J Mol Biol 48:443-453
 * @example bioseq.alignGlobal($one, $other, bioseq.dnaOptions())
 */
export func alignGlobal(a as string, b as string, opts as AlignOptions) {
    return align($a, $b, $opts, false);
}

/**
 * Find the best-matching stretch within two sequences - the Smith-Waterman
 * algorithm.
 *
 * No cell scores below zero, the traceback starts at the highest-scoring cell
 * and stops at the first zero. The result is the single best local match;
 * `startA` and `startB` locate it. Boundaries are soft: with +2/-1 scoring,
 * random flanking sequence still scores above zero on average, so the span can
 * reach past the region of interest.
 *
 * @param a {string} the first sequence
 * @param b {string} the second sequence
 * @param opts {AlignOptions} scoring, from `dnaOptions` or `proteinOptions`
 * @return {Alignment} the best local alignment
 * @throws {Error} kind "bioinformatics" when either sequence is empty or holds a
 *         non-residue, or the matrix would exceed `opts.maxCells`
 * @see Smith TF, Waterman MS (1981) J Mol Biol 147:195-197
 */
export func alignLocal(a as string, b as string, opts as AlignOptions) {
    return align($a, $b, $opts, true);
}

/**
 * The shared alignment engine.
 *
 * @param a {string} the first sequence
 * @param b {string} the second sequence
 * @param opts {AlignOptions} the scoring scheme
 * @param local {bool} true for Smith-Waterman, false for Needleman-Wunsch
 * @return {Alignment} the alignment
 * @internal
 */
func align(a as string, b as string, opts as AlignOptions, local as bool) {
    def left as string init requireSequence($a);
    def right as string init requireSequence($b);
    def ax as bytes init convert.bytesFromString($left, "utf-8");
    def bx as bytes init convert.bytesFromString($right, "utf-8");
    def n as int init len($ax);
    def m as int init len($bx);
    requireAlignable($n, $m, $opts.maxCells);
    def width as int init $m + 1;
    def gap as int init $opts.gapScore;
    def scores as map of int to int init pairScores($left, $right, $opts);

    # Row 0 is laid down here; each later row is appended cell by cell below,
    # so the matrix is never filled with zeroes only to be overwritten.
    def dp as list of int;
    for (def j as int init 0; $j <= $m; $j = $j + 1) {
        if ($local) {
            $dp[] = 0;
        } else {
            $dp[] = $j * $gap;
        }
    }

    def bestScore as int init 0;
    def bestI as int init $n;
    def bestJ as int init $m;
    for (def i as int init 1; $i <= $n; $i = $i + 1) {
        def row as int init $i * $width;
        def previous as int init $row - $width;
        # The first sequence's residue is fixed for the whole row, so its half
        # of the packed score key is computed once here rather than per cell.
        #
        # A fuller version of this - a per-row table of scores indexed by a
        # compact alphabet - was written and measured at 2% on a 400 x 400
        # alignment, which is inside this machine's noise. It was removed again:
        # seventy lines of machinery for nothing measurable is a bad trade, and
        # the DP loop is otherwise already at the interpreter's floor.
        def fixed as int init $ax[$i - 1] * 256;
        if ($local) {
            $dp[] = 0;
        } else {
            $dp[] = $i * $gap;
        }
        for (def j as int init 1; $j <= $m; $j = $j + 1) {
            def best as int init $dp[$previous + $j - 1] + $scores[$fixed + $bx[$j - 1]];
            def up as int init $dp[$previous + $j] + $gap;
            def sideways as int init $dp[$row + $j - 1] + $gap;
            if ($up > $best) {
                $best = $up;
            }
            if ($sideways > $best) {
                $best = $sideways;
            }
            if ($local and $best < 0) {
                $best = 0;
            }
            $dp[] = $best;
            if ($local and $best > $bestScore) {
                $bestScore = $best;
                $bestI = $i;
                $bestJ = $j;
            }
        }
    }
    if (not $local) {
        $bestScore = $dp[$n * $width + $m];
    }
    return traceback($left, $right, $dp, $width, $bestI, $bestJ, $bestScore, $opts, $local);
}

/**
 * Refuse an alignment whose matrix would not fit in memory.
 *
 * The matrix is `(n + 1) * (m + 1)` boxed values at roughly a kilobyte each, so
 * the product is what matters rather than either length. Raising rather than
 * trying and being killed lets a caller fall back - to a local alignment on a
 * windowed region, or to a distance instead.
 *
 * Passing 0 removes the ceiling entirely. That is deliberate - a caller who has
 * the memory and means it should not be overruled by the library - but it is
 * the one way to make this deck exhaust the host's memory, and two sequences of
 * five thousand residues would need about twenty-five gigabytes. Remove the
 * limit only for inputs whose size you already know.
 *
 * @param n {int} length of the first sequence
 * @param m {int} length of the second
 * @param maxCells {int} the ceiling, or 0 for none
 * @throws {Error} kind "bioinformatics" when the matrix would be too large
 * @internal
 */
func requireAlignable(n as int, m as int, maxCells as int) {
    if ($maxCells <= 0) {
        return;
    }
    def cells as int init ($n + 1) * ($m + 1);
    if ($cells > $maxCells) {
        raiseError("aligning " + convert.toString($n) + " against " + convert.toString($m) +
            " residues needs " + convert.toString($cells) + " matrix cells, over the limit of " +
            convert.toString($maxCells) +
            " (roughly a kilobyte each, so this would exhaust memory)" +
            " - raise AlignOptions.maxCells, or align a smaller region");
    }
}

/**
 * Precompute every substitution score the two sequences can need.
 *
 * Keyed by `byteA * 256 + byteB`, so the inner loop does one integer-keyed
 * lookup instead of building a two-character string per cell. At most a few
 * hundred entries, built once.
 *
 * @param a {string} the first sequence
 * @param b {string} the second sequence
 * @param opts {AlignOptions} the scoring scheme
 * @return {map of int to int} packed residue pair -> score
 * @internal
 */
func pairScores(a as string, b as string, opts as AlignOptions) {
    def useMatrix as bool init len($opts.matrix) > 0;
    def out as map of int to int;
    for (def x in distinctChars($a)) {
        for (def y in distinctChars($b)) {
            def key as int init convert.bytesFromString($x, "utf-8")[0] * 256 +
                convert.bytesFromString($y, "utf-8")[0];
            if (not maps.has($out, $key)) {
                $out[$key] = pairScore($x, $y, $opts, $useMatrix);
            }
        }
    }
    return $out;
}

/**
 * The distinct characters of a sequence, in first-seen order.
 *
 * @param s {string} the sequence
 * @return {list of string} each character once
 */
func distinctChars(s as string) {
    def seen as map of string to bool init {};
    def out as list of string init [];
    for (def c in strings.chars($s)) {
        if (not maps.has($seen, $c)) {
            $seen[$c] = true;
            $out[] = $c;
        }
    }
    return $out;
}

/**
 * Score one residue pair under the given options.
 *
 * @param x {string} one residue
 * @param y {string} the other
 * @param opts {AlignOptions} the scoring scheme
 * @param useMatrix {bool} whether a substitution matrix was supplied
 * @return {int} the score
 * @internal
 */
func pairScore(x as string, y as string, opts as AlignOptions, useMatrix as bool) {
    if ($useMatrix) {
        return substitutionScore($opts.matrix, $x, $y);
    }
    if ($x == $y) {
        return $opts.matchScore;
    }
    return $opts.mismatchScore;
}

/**
 * Walk the matrix back from a cell, building the aligned strings.
 *
 * The predecessor is recovered by recomputing which of the three moves produced
 * the cell's value, rather than by keeping a second direction matrix - the same
 * answer for half the memory.
 *
 * @param a {string} the first sequence
 * @param b {string} the second sequence
 * @param dp {list of int} the filled score matrix
 * @param width {int} the matrix row width, `len(b) + 1`
 * @param fromI {int} the cell to start from, row
 * @param fromJ {int} the cell to start from, column
 * @param score {int} the alignment score to report
 * @param opts {AlignOptions} the scoring scheme
 * @param local {bool} true to stop at the first zero-scoring cell
 * @return {Alignment} the assembled alignment
 * @internal
 */
func traceback(
    a as string,
    b as string,
    dp as list of int,
    width as int,
    fromI as int,
    fromJ as int,
    score as int,
    opts as AlignOptions,
    local as bool) {
    def useMatrix as bool init len($opts.matrix) > 0;
    def gap as int init $opts.gapScore;
    def topRow as bytes;
    def bottomRow as bytes;
    def markupRow as bytes;
    def identities as int init 0;
    def similarities as int init 0;
    def gaps as int init 0;
    def i as int init $fromI;
    def j as int init $fromJ;

    while ($i > 0 or $j > 0) {
        if ($local and $dp[$i * $width + $j] == 0) {
            break;
        }
        def x as string init "";
        def y as string init "";
        if ($i > 0 and $j > 0) {
            $x = $a[$i - 1..$i];
            $y = $b[$j - 1..$j];
        }
        def sub as int init 0;
        if (len($x) > 0) {
            $sub = pairScore($x, $y, $opts, $useMatrix);
        }
        if ($i > 0 and $j > 0 and $dp[$i * $width + $j] == $dp[($i - 1) * $width + $j - 1] + $sub) {
            $topRow[] = convert.bytesFromString($x, "utf-8")[0];
            $bottomRow[] = convert.bytesFromString($y, "utf-8")[0];
            if ($x == $y) {
                $markupRow[] = 124;
                $identities = $identities + 1;
                $similarities = $similarities + 1;
            } elseif ($sub > 0) {
                $markupRow[] = 58;
                $similarities = $similarities + 1;
            } else {
                $markupRow[] = 32;
            }
            $i = $i - 1;
            $j = $j - 1;
        } elseif ($i > 0 and $dp[$i * $width + $j] == $dp[($i - 1) * $width + $j] + $gap) {
            $topRow[] = convert.bytesFromString($a[$i - 1..$i], "utf-8")[0];
            $bottomRow[] = 45;
            $markupRow[] = 32;
            $gaps = $gaps + 1;
            $i = $i - 1;
        } else {
            $topRow[] = 45;
            $bottomRow[] = convert.bytesFromString($b[$j - 1..$j], "utf-8")[0];
            $markupRow[] = 32;
            $gaps = $gaps + 1;
            $j = $j - 1;
        }
    }

    def alignedA as string init reverseBytes($topRow);
    def alignedB as string init reverseBytes($bottomRow);
    return Alignment{
        a: $alignedA,
        b: $alignedB,
        markup: reverseBytes($markupRow),
        score: $score,
        length: len($alignedA),
        identities: $identities,
        similarities: $similarities,
        gaps: $gaps,
        startA: $i,
        startB: $j
    };
}

/**
 * Reverse a byte buffer into a string.
 *
 * The traceback walks backwards, so every row is built in reverse.
 *
 * @param buf {bytes} the buffer to reverse
 * @return {string} the reversed contents
 * @internal
 */
func reverseBytes(buf as bytes) {
    def out as bytes;
    for (def i as int init len($buf) - 1; $i >= 0; $i = $i - 1) {
        $out[] = $buf[$i];
    }
    return convert.stringFromBytes($out, "utf-8");
}

# ---------------------------------------------------------------- reading it

/**
 * The percentage of aligned columns holding identical residues.
 *
 * Measured over the alignment's length, gap columns included, so a global
 * alignment of sequences of different lengths reports a low identity even
 * when the shorter matches perfectly.
 *
 * @param al {Alignment} the alignment
 * @return {float} percent identity, 0 to 100
 * @throws {Error} kind "bioinformatics" when the alignment is empty
 */
export func percentIdentity(al as Alignment) {
    requireAligned($al);
    return 100.0 * $al.identities / $al.length;
}

/**
 * The percentage of aligned columns scoring positively.
 *
 * Identities plus conservative substitutions. Without a substitution matrix
 * only identities score positively, so this equals `percentIdentity`.
 *
 * @param al {Alignment} the alignment
 * @return {float} percent similarity, 0 to 100
 * @throws {Error} kind "bioinformatics" when the alignment is empty
 */
export func percentSimilarity(al as Alignment) {
    requireAligned($al);
    return 100.0 * $al.similarities / $al.length;
}

/**
 * The percentage of aligned columns holding a gap.
 *
 * @param al {Alignment} the alignment
 * @return {float} percent gaps, 0 to 100
 * @throws {Error} kind "bioinformatics" when the alignment is empty
 */
export func percentGaps(al as Alignment) {
    requireAligned($al);
    return 100.0 * $al.gaps / $al.length;
}

/**
 * Reject an alignment that is empty or internally inconsistent.
 *
 * `Alignment` is an exported struct, so one can arrive from somewhere other
 * than `alignGlobal` - built by hand, or carried in from another aligner. A
 * record whose `length` disagrees with its rows used to reach the slicing in
 * `formatAlignment` and fail with a raw interpreter error rather than this
 * deck's, which broke the one-error-shape contract for anyone catching on
 * `kind`. Checking the invariant here keeps that promise.
 *
 * @param al {Alignment} the alignment
 * @throws {Error} kind "bioinformatics" when it is empty or its fields disagree
 * @internal
 */
func requireAligned(al as Alignment) {
    if ($al.length == 0) {
        raiseError("the alignment has no columns to measure");
    }
    if (len($al.a) != $al.length or len($al.b) != $al.length or len($al.markup) != $al.length) {
        raiseError("the alignment is inconsistent: length is " + convert.toString($al.length) +
            " but the rows are " + convert.toString(len($al.a)) + ", " +
            convert.toString(len($al.b)) + " and " + convert.toString(len($al.markup)) + " long");
    }
}

/**
 * Render an alignment as numbered blocks, the way EMBOSS prints one.
 *
 * Each block shows the first sequence, the markup line and the second, with the
 * position each row starts at. Positions are 1-based residue positions in the
 * original sequences (gaps not counted).
 *
 * @param al {Alignment} the alignment
 * @param width {int} residues per block; must be positive
 * @return {string} the formatted alignment, ending in a newline
 * @throws {Error} kind "bioinformatics" when the alignment is empty or the width is not
 *         positive
 */
export func formatAlignment(al as Alignment, width as int) {
    requireAligned($al);
    if ($width <= 0) {
        raiseError("block width must be positive, got " + convert.toString($width));
    }
    def out as string init "";
    def posA as int init $al.startA + 1;
    def posB as int init $al.startB + 1;
    def at as int init 0;
    while ($at < $al.length) {
        def stop as int init $at + $width;
        if ($stop > $al.length) {
            $stop = $al.length;
        }
        def chunkA as string init $al.a[$at..$stop];
        def chunkB as string init $al.b[$at..$stop];
        def usedA as int init countResidues2($chunkA);
        def usedB as int init countResidues2($chunkB);
        $out = $out + io.sprintf("%d|pad=6 %s\n", $posA, $chunkA);
        $out = $out + io.sprintf("%s|pad=6 %s\n", "", $al.markup[$at..$stop]);
        $out = $out + io.sprintf("%d|pad=6 %s\n\n", $posB, $chunkB);
        $posA = $posA + $usedA;
        $posB = $posB + $usedB;
        $at = $stop;
    }
    return $out;
}

/**
 * Count the non-gap characters in an aligned chunk.
 *
 * @param chunk {string} a slice of an aligned sequence
 * @return {int} how many real residues it holds
 * @internal
 */
func countResidues2(chunk as string) {
    def n as int init 0;
    def buf as bytes init convert.bytesFromString($chunk, "utf-8");
    def width as int init len($buf);
    for (def i as int init 0; $i < $width; $i = $i + 1) {
        if ($buf[$i] != 45) {
            $n = $n + 1;
        }
    }
    return $n;
}

# ------------------------------------------------------------------ distances

/**
 * The number of positions at which two equal-length sequences differ.
 *
 * Substitutions only, so the lengths must match; use `editDistance` when they
 * may differ or an indel is plausible. Residues are compared literally.
 *
 * @param a {string} the first sequence
 * @param b {string} the second sequence
 * @return {int} the number of differing positions
 * @throws {Error} kind "bioinformatics" when either sequence is empty, holds a
 *         non-residue, or the two differ in length
 * @example bioseq.hammingDistance("ATGC", "ATCC") == 1
 */
export func hammingDistance(a as string, b as string) {
    def left as string init requireSequence($a);
    def right as string init requireSequence($b);
    if (len($left) != len($right)) {
        raiseError("Hamming distance needs equal lengths, got " + convert.toString(len($left)) +
            " and " + convert.toString(len($right)) + " (use editDistance when the lengths differ)");
    }
    def ax as bytes init convert.bytesFromString($left, "utf-8");
    def bx as bytes init convert.bytesFromString($right, "utf-8");
    def differences as int init 0;
    def n as int init len($ax);
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        if ($ax[$i] != $bx[$i]) {
            $differences = $differences + 1;
        }
    }
    return $differences;
}

/**
 * The fewest single-residue edits that turn one sequence into the other.
 *
 * Levenshtein distance: insertions, deletions and substitutions each cost one.
 * Two rows of `len(b) + 1` entries instead of a full matrix, so memory is O(m)
 * in the second sequence's length (pass the shorter one second); time is O(n*m).
 * `maxCells` does not apply.
 *
 * @param a {string} the first sequence
 * @param b {string} the second sequence
 * @return {int} the edit distance
 * @throws {Error} kind "bioinformatics" when either sequence is empty or holds a
 *         non-residue
 * @example bioseq.editDistance("ATGC", "AGC") == 1
 */
export func editDistance(a as string, b as string) {
    def left as string init requireSequence($a);
    def right as string init requireSequence($b);
    def ax as bytes init convert.bytesFromString($left, "utf-8");
    def bx as bytes init convert.bytesFromString($right, "utf-8");
    def n as int init len($ax);
    def m as int init len($bx);
    def previous as list of int;
    for (def j as int init 0; $j <= $m; $j = $j + 1) {
        $previous[] = $j;
    }
    for (def i as int init 1; $i <= $n; $i = $i + 1) {
        def current as list of int;
        $current[] = $i;
        for (def j as int init 1; $j <= $m; $j = $j + 1) {
            def cost as int init 1;
            if ($ax[$i - 1] == $bx[$j - 1]) {
                $cost = 0;
            }
            def best as int init $previous[$j - 1] + $cost;
            if ($previous[$j] + 1 < $best) {
                $best = $previous[$j] + 1;
            }
            if ($current[$j - 1] + 1 < $best) {
                $best = $current[$j - 1] + 1;
            }
            $current[] = $best;
        }
        $previous = $current;
    }
    return $previous[$m];
}
