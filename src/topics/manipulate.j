# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - manipulate: generating, shuffling, mutating and slicing sequences.
# Spliced into bioinformatics.j via include - not a standalone module.
#
# The random verbs draw from `math.rand`, so `math.randSeed` makes a whole run
# reproducible. The deck keeps no generator of its own: one seed per program,
# no hidden second stream.

/**
 * The four unambiguous DNA bases, the alphabet the random and mutation verbs
 * draw from. No ambiguity codes: an `N` would be missing data, not sequence.
 */
export def const DNA_BASES as string init "ACGT";

/**
 * The four unambiguous RNA bases.
 */
export def const RNA_BASES as string init "ACGU";

/**
 * The twenty standard amino acids, in the conventional one-letter order.
 */
export def const STANDARD_RESIDUES as string init "ACDEFGHIKLMNPQRSTVWY";

# ------------------------------------------------------------------- generation

/**
 * A random DNA sequence of the given length, with the four bases equally likely.
 *
 * No upper bound: `randomDna(2000000000)` tries to allocate two gigabytes.
 *
 * @param n {int} the length in bases; must be positive
 * @return {string} the sequence
 * @throws {Error} kind "bioinformatics" when the length is not positive
 * @example bioseq.randomDna(20)
 */
export func randomDna(n as int) {
    return randomFrom($n, DNA_BASES);
}

/**
 * A random RNA sequence of the given length.
 *
 * @param n {int} the length in bases; must be positive
 * @return {string} the sequence
 * @throws {Error} kind "bioinformatics" when the length is not positive
 */
export func randomRna(n as int) {
    return randomFrom($n, RNA_BASES);
}

/**
 * A random protein of the given length, with the twenty standard residues
 * equally likely.
 *
 * Real proteins are far from uniform; use `randomProteinWith` with a measured
 * composition for a realistic null model.
 *
 * @param n {int} the length in residues; must be positive
 * @return {string} the sequence
 * @throws {Error} kind "bioinformatics" when the length is not positive
 */
export func randomProtein(n as int) {
    return randomFrom($n, STANDARD_RESIDUES);
}

/**
 * A random DNA sequence with a chosen base composition.
 *
 * Weights are relative and normalised, so `{"A": 3.0, "T": 1.0}` and
 * `{"A": 75.0, "T": 25.0}` are the same distribution. A base left out of the map
 * (or weighted 0) never appears.
 *
 * @param n {int} the length in bases; must be positive
 * @param weights {map of string to float} base -> relative frequency
 * @return {string} the sequence
 * @throws {Error} kind "bioinformatics" when the length is not positive, a key is not
 *         one of `ACGT`, a weight is negative, or every weight is zero
 * @example bioseq.randomDnaWith(100, {"A": 30.0, "T": 30.0, "C": 20.0, "G": 20.0})
 */
export func randomDnaWith(n as int, weights as map of string to float) {
    return randomWeighted($n, $weights, DNA_BASES);
}

/**
 * A random RNA sequence with a chosen base composition. Weights as for
 * `randomDnaWith`, keys from `ACGU`.
 *
 * @param n {int} the length in bases; must be positive
 * @param weights {map of string to float} base -> relative frequency
 * @return {string} the sequence
 * @throws {Error} kind "bioinformatics" when the length is not positive, a key is not
 *         one of `ACGU`, a weight is negative, or every weight is zero
 */
export func randomRnaWith(n as int, weights as map of string to float) {
    return randomWeighted($n, $weights, RNA_BASES);
}

/**
 * A random protein with a chosen residue composition.
 *
 * Weights as for `randomDnaWith`, keys from `STANDARD_RESIDUES`.
 * `residueFrequency` of a real protein is a valid weight map.
 *
 * @param n {int} the length in residues; must be positive
 * @param weights {map of string to float} residue -> relative frequency
 * @return {string} the sequence
 * @throws {Error} kind "bioinformatics" when the length is not positive, a key is not
 *         one of the twenty standard residues, a weight is negative, or every
 *         weight is zero
 */
export func randomProteinWith(n as int, weights as map of string to float) {
    return randomWeighted($n, $weights, STANDARD_RESIDUES);
}

/**
 * Draw `n` residues uniformly from an alphabet.
 *
 * @param n {int} how many residues to draw
 * @param alphabet {string} the residues to draw from
 * @return {string} the sequence
 * @throws {Error} kind "bioinformatics" when the length is not positive
 * @internal
 */
func randomFrom(n as int, alphabet as string) {
    requirePositiveLength($n);
    def pool as bytes init convert.bytesFromString($alphabet, "utf-8");
    def top as int init len($pool) - 1;
    def out as bytes;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        $out[] = $pool[math.randInt(0, $top)];
    }
    return convert.stringFromBytes($out, "utf-8");
}

/**
 * Draw `n` residues from an alphabet according to relative weights.
 *
 * Builds a cumulative distribution once, then does one `math.rand` and a linear
 * walk per residue. The alphabet is at most twenty entries, so the walk is
 * cheaper than the arithmetic a binary search would need.
 *
 * @param n {int} how many residues to draw
 * @param weights {map of string to float} residue -> relative frequency
 * @param alphabet {string} the residues the weights may name
 * @return {string} the sequence
 * @throws {Error} kind "bioinformatics" when a parameter is out of range
 * @internal
 */
func randomWeighted(n as int, weights as map of string to float, alphabet as string) {
    requirePositiveLength($n);
    def total as float init 0.0;
    for (def residue in $weights) {
        if (not strings.contains($alphabet, $residue)) {
            raiseError("'" + $residue + "' is not one of the residues this verb" + " generates (" +
                $alphabet + ")");
        }
        if ($weights[$residue] < 0.0) {
            raiseError("weight for '" + $residue + "' is negative");
        }
        $total = $total + $weights[$residue];
    }
    if ($total <= 0.0) {
        raiseError("at least one weight must be greater than zero");
    }

    # cumulative distribution, parallel lists so the draw is two indexed reads
    def codes as list of int;
    def cumulative as list of float;
    def running as float init 0.0;
    for (def residue in $weights) {
        if ($weights[$residue] <= 0.0) {
            continue;
        }
        $running = $running + $weights[$residue] / $total;
        $codes[] = convert.bytesFromString($residue, "utf-8")[0];
        $cumulative[] = $running;
    }

    def out as bytes;
    def last as int init len($codes) - 1;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def draw as float init math.rand();
        def picked as int init $codes[$last];
        for (def j as int init 0; $j <= $last; $j = $j + 1) {
            if ($draw < $cumulative[$j]) {
                $picked = $codes[$j];
                break;
            }
        }
        $out[] = $picked;
    }
    return convert.stringFromBytes($out, "utf-8");
}

/**
 * Reject a non-positive requested length.
 *
 * @param n {int} the requested length
 * @throws {Error} kind "bioinformatics" when it is not positive
 * @internal
 */
func requirePositiveLength(n as int) {
    if ($n <= 0) {
        raiseError("length must be positive, got " + convert.toString($n));
    }
}

# --------------------------------------------------------------------- shuffle

/**
 * Reorder a sequence at random, keeping its composition exactly.
 *
 * Fisher-Yates: every permutation equally likely, residue counts unchanged.
 * The null model for "is this motif surprising?" that keeps composition fixed.
 *
 * @param s {string} any sequence
 * @return {string} the same residues in a random order
 * @throws {Error} kind "bioinformatics" when the sequence is empty or holds a
 *         non-residue
 * @example bioseq.shuffle($promoter)
 */
export func shuffle(s as string) {
    def clean as string init requireSequence($s);
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    for (def i as int init len($buf) - 1; $i > 0; $i = $i - 1) {
        def j as int init math.randInt(0, $i);
        def swap as int init $buf[$i];
        $buf[$i] = $buf[$j];
        $buf[$j] = $swap;
    }
    return convert.stringFromBytes($buf, "utf-8");
}

# ---------------------------------------------------------------------- mutate

/**
 * Substitute residues at random, at a given per-position rate.
 *
 * Substitutions only. Each position mutates with probability `rate`, so the
 * number of changes varies between calls. A mutated position always changes:
 * the replacement is drawn from the alphabet minus the current residue
 * (`DNA_BASES`, `RNA_BASES` or `STANDARD_RESIDUES`). An ambiguity code chosen
 * for mutation is replaced by a draw from the whole alphabet.
 *
 * @param s {string} the sequence to mutate
 * @param rate {float} the per-position substitution probability, 0 to 1
 * @param kind {int} `bioseq.DNA`, `bioseq.RNA` or `bioseq.PROTEIN`, which
 *        decides the alphabet substitutions are drawn from
 * @return {string} the mutated sequence, the same length as the input
 * @throws {Error} kind "bioinformatics" when the rate is outside 0 to 1, the kind is
 *         unknown, or the sequence is empty or of the wrong alphabet
 * @example bioseq.mutate($gene, 0.01, bioseq.DNA)
 */
export func mutate(s as string, rate as float, kind as int) {
    requireKind($kind);
    if ($rate < 0.0 or $rate > 1.0) {
        raiseError("mutation rate must be between 0 and 1, got " + convert.toString($rate));
    }
    def clean as string init "";
    def alphabet as string init DNA_BASES;
    if ($kind == DNA) {
        $clean = requireDna($s);
    } elseif ($kind == RNA) {
        $clean = requireRna($s);
        $alphabet = RNA_BASES;
    } else {
        $clean = requireProtein($s);
        $alphabet = STANDARD_RESIDUES;
    }
    def pool as bytes init convert.bytesFromString($alphabet, "utf-8");
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def out as bytes;
    def n as int init len($buf);
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        if (math.rand() < $rate) {
            $out[] = differentResidue($pool, $buf[$i]);
        } else {
            $out[] = $buf[$i];
        }
    }
    return convert.stringFromBytes($out, "utf-8");
}

/**
 * Draw a residue from the pool that differs from the one given.
 *
 * Drawing from the whole alphabet and rejecting the original would loop; this
 * draws from the alphabet less one position instead, so it always terminates in
 * one step. A residue outside the pool - an ambiguity code that slipped past a
 * lenient alphabet - is replaced by an unbiased draw rather than left alone.
 *
 * @param pool {bytes} the alphabet to draw from
 * @param current {int} the residue being replaced
 * @return {int} a different residue
 * @internal
 */
func differentResidue(pool as bytes, current as int) {
    def at as int init -1;
    def choices as int init len($pool);
    for (def i as int init 0; $i < $choices; $i = $i + 1) {
        if ($pool[$i] == $current) {
            $at = $i;
            break;
        }
    }
    if ($at < 0) {
        return $pool[math.randInt(0, len($pool) - 1)];
    }
    def drawn as int init math.randInt(0, len($pool) - 2);
    if ($drawn >= $at) {
        $drawn = $drawn + 1;
    }
    return $pool[$drawn];
}

# ---------------------------------------------------------- slicing and windows

/**
 * Extract a subsequence by 1-based inclusive position.
 *
 * The deck's one 1-based verb: "residues 25 to 40" is `subrange($s, 25, 40)`,
 * the same as `$s[24..40]`. Positions count residues after whitespace and digits
 * are removed.
 *
 * @param s {string} the sequence
 * @param from {int} the first position, counting from 1
 * @param through {int} the last position, included
 * @return {string} the subsequence
 * @throws {Error} kind "bioinformatics" when the sequence is empty or holds a
 *         non-residue, the range falls outside it, or runs backwards
 * @example bioseq.subrange("ATGCATGC", 3, 5) == "GCA"
 */
export func subrange(s as string, from as int, through as int) {
    def clean as string init requireSequence($s);
    if ($from < 1) {
        raiseError("positions count from 1, got " + convert.toString($from));
    }
    if ($through < $from) {
        raiseError("range runs backwards: " + convert.toString($from) + " to " +
            convert.toString($through));
    }
    if ($through > len($clean)) {
        raiseError("position " + convert.toString($through) + " is past the end of a " +
            convert.toString(len($clean)) + "-residue sequence");
    }
    return $clean[$from - 1..$through];
}

/**
 * Cut a sequence into windows of `size` residues starting every `step`.
 *
 * Windows overlap when `step < size`. The last partial window is dropped, as in
 * `gcWindows`. Works on any sequence.
 *
 * @param s {string} the sequence
 * @param size {int} the window width; must be positive
 * @param step {int} how far to advance between windows; must be positive
 * @return {list of string} the windows, in order
 * @throws {Error} kind "bioinformatics" when the sequence is shorter than one window,
 *         or a parameter is not positive
 * @example bioseq.windows("ATGCATGC", 4, 4) == ["ATGC", "ATGC"]
 */
export func windows(s as string, size as int, step as int) {
    def clean as string init requireSequence($s);
    if ($size <= 0) {
        raiseError("window size must be positive, got " + convert.toString($size));
    }
    if ($step <= 0) {
        raiseError("window step must be positive, got " + convert.toString($step));
    }
    if (len($clean) < $size) {
        raiseError("sequence is " + convert.toString(len($clean)) + " residues, shorter than the " +
            convert.toString($size) + "-residue window");
    }
    def out as list of string;
    def start as int init 0;
    while ($start + $size <= len($clean)) {
        $out[] = $clean[$start..$start + $size];
        $start = $start + $step;
    }
    return $out;
}
