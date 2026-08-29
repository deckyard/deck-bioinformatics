# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - core: alphabets, validation, and the byte-loop conventions every
# other topic is built on. Spliced into bioinformatics.j via include - not a standalone
# module.
#
# Two conventions are set here and followed everywhere else in the deck, both
# measured rather than assumed (examples/bench.j reproduces the numbers):
#
#   1. Per-residue loops walk a `bytes` buffer, never the string. Jennifer
#      strings are rune-indexed, so `$s[$i..$i + 1]` rescans from the start of
#      the string on every access - a 30 kb sequence takes ~17x longer that way.
#      `convert.bytesFromString` once, index the buffer, convert back at the end.
#
#   2. Results accrete into a `bytes` buffer with the `$out[] = code;` append
#      sugar, which is amortised O(N). `lists.push` in a loop is O(N^2) and, at
#      20 000 residues, roughly 3000x slower.

# ---------------------------------------------------------------- molecule kind

/**
 * Molecule kind: deoxyribonucleic acid. Passed where a verb's answer depends on
 * which molecule it is looking at, such as `molecularWeight`.
 */
export def const DNA as int init 1;

/**
 * Molecule kind: ribonucleic acid.
 */
export def const RNA as int init 2;

/**
 * Molecule kind: protein.
 */
export def const PROTEIN as int init 3;

# -------------------------------------------------------------------- alphabets

/**
 * The DNA alphabet: the four bases `ACGT` followed by the eleven IUPAC
 * ambiguity codes `RYSWKMBDHVN`.
 */
export def const DNA_ALPHABET as string init "ACGTRYSWKMBDHVN";

/**
 * The RNA alphabet: `ACGU` plus the same eleven IUPAC ambiguity codes.
 */
export def const RNA_ALPHABET as string init "ACGURYSWKMBDHVN";

/**
 * The protein alphabet: the twenty standard residues, the ambiguity codes
 * `BZXJ`, the rare residues `U` (selenocysteine) and `O` (pyrrolysine), and
 * `*` for a stop, which is what `translate` emits.
 */
export def const PROTEIN_ALPHABET as string init "ACDEFGHIKLMNPQRSTVWYBZXJUO*";

# Byte-code membership sets for the alphabets above, keyed by ASCII code so a
# per-residue loop can test a byte without materialising a one-character string.
# `testAlphabetsAgree` in core_test.j derives these from the alphabet strings and
# fails if the two ever drift apart, so the codes below cannot go stale silently.
#
#   A 65   B 66   C 67   D 68   E 69   F 70   G 71   H 72   I 73   J 74
#   K 75   L 76   M 77   N 78   O 79   P 80   Q 81   R 82   S 83   T 84
#   U 85   V 86   W 87   X 88   Y 89   Z 90   * 42

def const DNA_CODES as map of int to bool init {
    65: true,
    67: true,
    71: true,
    84: true,
    82: true,
    89: true,
    83: true,
    87: true,
    75: true,
    77: true,
    66: true,
    68: true,
    72: true,
    86: true,
    78: true
};

def const RNA_CODES as map of int to bool init {
    65: true,
    67: true,
    71: true,
    85: true,
    82: true,
    89: true,
    83: true,
    87: true,
    75: true,
    77: true,
    66: true,
    68: true,
    72: true,
    86: true,
    78: true
};

def const PROTEIN_CODES as map of int to bool init {
    65: true,
    67: true,
    68: true,
    69: true,
    70: true,
    71: true,
    72: true,
    73: true,
    75: true,
    76: true,
    77: true,
    78: true,
    80: true,
    81: true,
    82: true,
    83: true,
    84: true,
    86: true,
    87: true,
    89: true,
    66: true,
    90: true,
    88: true,
    74: true,
    85: true,
    79: true,
    42: true
};

# The union of DNA_CODES and RNA_CODES, for verbs that accept either nucleic
# acid (gcContent, find, the window helpers).
def const NUCLEOTIDE_CODES as map of int to bool init {
    65: true,
    67: true,
    71: true,
    84: true,
    85: true,
    82: true,
    89: true,
    83: true,
    87: true,
    75: true,
    77: true,
    66: true,
    68: true,
    72: true,
    86: true,
    78: true
};

# ----------------------------------------------------------------- error shape

/**
 * Raise the deck's one error shape.
 *
 * @param message {string} human-readable failure description
 * @throws {Error} always, kind "bioinformatics"
 * @internal
 */
func raiseError(message as string) {
    throw Error{kind: "bioinformatics", message: $message, file: "", line: 0, col: 0};
}

/**
 * Render one byte as a printable character for an error message.
 *
 * @param code {int} an ASCII byte value
 * @return {string} the character itself, or a `byte NN` note when unprintable
 * @internal
 */
func describeCode(code as int) {
    if ($code < 32 or $code > 126) {
        return "byte " + convert.toString($code);
    }
    def one as bytes;
    $one[] = $code;
    return "'" + convert.stringFromBytes($one, "utf-8") + "'";
}

# ------------------------------------------------------------ filters (lenient)

/**
 * Keep only the DNA alphabet, dropping every other character.
 *
 * The lenient counterpart to the validating verbs: where `reverseComplement`
 * raises on a stray character, this silently strips it. Lowercase input is
 * uppercased. `U` is *not* kept - a sequence spelled with `U` is RNA, so use
 * `filterRna`, or `backTranscribe` to move it to DNA.
 *
 * @param s {string} any text - pasted sequence, numbered listing, FASTA body
 * @return {string} the retained residues, uppercased
 * @example bioseq.filterDna("1 atg acg > x") == "ATGACG"
 */
export func filterDna(s as string) {
    return filterBy($s, DNA_ADMITS);
}

/**
 * Keep only the RNA alphabet, dropping every other character.
 *
 * @param s {string} any text
 * @return {string} the retained residues, uppercased
 */
export func filterRna(s as string) {
    return filterBy($s, RNA_ADMITS);
}

/**
 * Keep only the protein alphabet, dropping every other character.
 *
 * @param s {string} any text
 * @return {string} the retained residues, uppercased
 */
export func filterProtein(s as string) {
    return filterBy($s, PROTEIN_ADMITS);
}

/**
 * Shared filter loop: uppercase each byte, keep the ones the set admits.
 *
 * @param s {string} the text to filter
 * @param admits {list of bool} membership by ASCII code
 * @return {string} the retained residues
 * @internal
 */
func filterBy(s as string, admits as list of bool) {
    def buf as bytes init convert.bytesFromString($s, "utf-8");
    def n as int init len($buf);
    def out as bytes;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def code as int init $buf[$i];
        if ($code >= 97 and $code <= 122) {
            $code = $code - 32;
        }
        if ($code < 128 and $admits[$code]) {
            $out[] = $code;
        }
    }
    return convert.stringFromBytes($out, "utf-8");
}

# ------------------------------------------------------------------ predicates

/**
 * Is every character of `s` in the DNA alphabet?
 *
 * Case-insensitive. An empty string is not a sequence, so it answers false.
 *
 * @param s {string} the candidate sequence
 * @return {bool} true when non-empty and wholly within `DNA_ALPHABET`
 */
export func isDna(s as string) {
    return isAll($s, DNA_ADMITS);
}

/**
 * Is every character of `s` in the RNA alphabet?
 *
 * @param s {string} the candidate sequence
 * @return {bool} true when non-empty and wholly within `RNA_ALPHABET`
 */
export func isRna(s as string) {
    return isAll($s, RNA_ADMITS);
}

/**
 * Is every character of `s` in the protein alphabet?
 *
 * Note that `ACGT` satisfies both this and `isDna` - a four-residue peptide and
 * a four-base oligo are spelled the same way. Ask the user, or ask `isDna`
 * first, when the distinction matters.
 *
 * @param s {string} the candidate sequence
 * @return {bool} true when non-empty and wholly within `PROTEIN_ALPHABET`
 */
export func isProtein(s as string) {
    return isAll($s, PROTEIN_ADMITS);
}

/**
 * Shared predicate loop.
 *
 * @param s {string} the candidate sequence
 * @param admits {list of bool} membership by ASCII code
 * @return {bool} true when non-empty and wholly admitted
 * @internal
 */
func isAll(s as string, admits as list of bool) {
    def buf as bytes init convert.bytesFromString($s, "utf-8");
    def n as int init len($buf);
    if ($n == 0) {
        return false;
    }
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def code as int init $buf[$i];
        if ($code >= 97 and $code <= 122) {
            $code = $code - 32;
        }
        if ($code >= 128 or not $admits[$code]) {
            return false;
        }
    }
    return true;
}

# ---------------------------------------------------------- validation (strict)

/**
 * Normalise a sequence for a verb that must not guess.
 *
 * Uppercases, drops whitespace and digits - so a pasted, numbered, line-wrapped
 * listing works - and raises on anything else, naming the character and its
 * position. Gap characters (`-`, `.`) raise too: dropping them would silently
 * shift every downstream coordinate, so an aligned sequence has to be filtered
 * deliberately.
 *
 * @param s {string} the sequence as the caller supplied it
 * @param admits {list of bool} membership by ASCII code
 * @param what {string} alphabet name, used in the error message
 * @return {string} the cleaned, uppercased sequence
 * @throws {Error} kind "bioinformatics" when empty or carrying a foreign character
 * @internal
 */
func normalizeBy(s as string, admits as list of bool, what as string) {
    def buf as bytes init convert.bytesFromString($s, "utf-8");
    def n as int init len($buf);
    def out as bytes;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def code as int init $buf[$i];
        # The overwhelmingly common case is an uppercase residue, so test for it
        # first and let everything else - whitespace, digits, lowercase - fall
        # through to the slower path below. Membership is a list index rather
        # than a map lookup, which measured cheaper and runs once per residue of
        # every input the deck is given.
        if ($code >= 65 and $code <= 90) {
            if ($admits[$code]) {
                $out[] = $code;
                continue;
            }
            raiseError(foreignCharacter($code, $i, $what));
        }
        if ($code == 32 or $code == 9 or $code == 10 or $code == 13) {
            continue;
        }
        if ($code >= 48 and $code <= 57) {
            continue;
        }
        if ($code >= 97 and $code <= 122) {
            $code = $code - 32;
        }
        # Anything still standing is either an admitted non-letter - `*` for a
        # stop is the one that matters - or genuinely foreign.
        if ($code < 128 and $admits[$code]) {
            $out[] = $code;
            continue;
        }
        raiseError(foreignCharacter($code, $i, $what));
    }
    if (len($out) == 0) {
        raiseError("empty " + $what + " sequence");
    }
    return convert.stringFromBytes($out, "utf-8");
}

/**
 * The message for a character that is not a residue of the expected alphabet.
 *
 * @param code {int} the offending byte
 * @param at {int} where it was found
 * @param what {string} the alphabet name
 * @return {string} the message
 * @internal
 */
func foreignCharacter(code as int, at as int, what as string) {
    return "not a " + $what + " sequence: " + describeCode($code) + " at position " +
        convert.toString($at) + " (use filterDna / filterRna / filterProtein to strip it)";
}

/**
 * Validate and normalise a DNA sequence.
 *
 * @param s {string} the sequence as the caller supplied it
 * @return {string} cleaned, uppercased DNA
 * @throws {Error} kind "bioinformatics" when empty or not DNA
 * @internal
 */
func requireDna(s as string) {
    return normalizeBy($s, DNA_ADMITS, "DNA");
}

/**
 * Validate and normalise an RNA sequence.
 *
 * @param s {string} the sequence as the caller supplied it
 * @return {string} cleaned, uppercased RNA
 * @throws {Error} kind "bioinformatics" when empty or not RNA
 * @internal
 */
func requireRna(s as string) {
    return normalizeBy($s, RNA_ADMITS, "RNA");
}

/**
 * Validate and normalise a nucleic-acid sequence of either kind.
 *
 * Used by the verbs that do not care which one they were handed - GC content,
 * pattern search, the window helpers.
 *
 * @param s {string} the sequence as the caller supplied it
 * @return {string} cleaned, uppercased DNA or RNA
 * @throws {Error} kind "bioinformatics" when empty or not a nucleic acid
 * @internal
 */
func requireNucleotide(s as string) {
    return normalizeBy($s, NUCLEOTIDE_ADMITS, "nucleotide");
}

/**
 * Validate and normalise a protein sequence.
 *
 * @param s {string} the sequence as the caller supplied it
 * @return {string} cleaned, uppercased protein
 * @throws {Error} kind "bioinformatics" when empty or not protein
 * @internal
 */
func requireProtein(s as string) {
    return normalizeBy($s, PROTEIN_ADMITS, "protein");
}

/**
 * Validate and normalise a sequence of any kind.
 *
 * The protein alphabet is a superset of both nucleic-acid alphabets, so this
 * admits any residue the deck recognises. Used by the verbs that genuinely do
 * not care - `reverse` is a string operation on residues, not a strand
 * operation.
 *
 * @param s {string} the sequence as the caller supplied it
 * @return {string} cleaned, uppercased residues
 * @throws {Error} kind "bioinformatics" when empty or holding a non-residue
 * @internal
 */
func requireSequence(s as string) {
    return normalizeBy($s, PROTEIN_ADMITS, "sequence");
}

/**
 * Reject a molecule kind that is not one of `DNA` / `RNA` / `PROTEIN`.
 *
 * @param kind {int} the kind to check
 * @throws {Error} kind "bioinformatics" when the value is not a molecule kind
 * @internal
 */
func requireKind(kind as int) {
    if ($kind != DNA and $kind != RNA and $kind != PROTEIN) {
        raiseError("unknown molecule kind " + convert.toString($kind) +
            " (expected bioseq.DNA, bioseq.RNA or bioseq.PROTEIN)");
    }
}

/**
 * Name a molecule kind, for messages and reports.
 *
 * @param kind {int} `DNA`, `RNA` or `PROTEIN`
 * @return {string} "DNA", "RNA" or "protein"
 * @throws {Error} kind "bioinformatics" when the value is not a molecule kind
 */
export func kindName(kind as int) {
    requireKind($kind);
    if ($kind == DNA) {
        return "DNA";
    }
    if ($kind == RNA) {
        return "RNA";
    }
    return "protein";
}
