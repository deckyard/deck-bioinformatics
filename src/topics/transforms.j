# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - transforms: turning one sequence into another.
# Spliced into bioinformatics.j via include - not a standalone module.
#
# Every loop here follows the conventions core.j sets out: walk a `bytes`
# buffer, accrete into a `bytes` buffer, keep the loop body flat. The codon
# lookup rebuilds a three-character key per codon and reads the table map, which
# the benchmark showed beats both string slicing and a packed integer index.

# ---------------------------------------------------------------- complement

/**
 * Complement each base, leaving the sequence in its original 5'-to-3' order.
 *
 * This is the strand-pairing operation on its own: it is almost never what you
 * want by itself, because a complementary strand read in the same direction is
 * written 3'-to-5'. `reverseComplement` is the one to reach for.
 *
 * The alphabet is preserved: an RNA sequence complements with `U`, a DNA
 * sequence with `T`. Ambiguity codes complement as sets, so `R` (A or G)
 * becomes `Y` (C or T), and the self-complementary `S`, `W` and `N` are
 * unchanged.
 *
 * @param s {string} a DNA or RNA sequence
 * @return {string} the complement, in the same alphabet and the same order
 * @throws {Error} kind "bioinformatics" when the sequence is empty, not a nucleic acid,
 *         or mixes `T` and `U` so the alphabet cannot be told
 * @example bioseq.complement("ATGC") == "TACG"
 */
export func complement(s as string) {
    def clean as string init requireNucleotide($s);
    return complementBytes($clean, isRnaSpelling($clean), false);
}

/**
 * Complement each base and reverse the sequence - the other strand, read the
 * way sequences are always read.
 *
 * The workhorse of the topic: a restriction site, a primer, or a feature found
 * on the minus strand is this operation away from the strand you have.
 * Applying it twice returns the original.
 *
 * @param s {string} a DNA or RNA sequence
 * @return {string} the reverse complement, in the same alphabet
 * @throws {Error} kind "bioinformatics" when the sequence is empty, not a nucleic acid,
 *         or mixes `T` and `U`
 * @example bioseq.reverseComplement("ATGACGTTACATGGATAA") == "TTATCCATGTAACGTCAT"
 */
export func reverseComplement(s as string) {
    def clean as string init requireNucleotide($s);
    return complementBytes($clean, isRnaSpelling($clean), true);
}

/**
 * Reverse the residues, without complementing them.
 *
 * Works on any sequence, protein included - it is a string operation on
 * residues, not a strand operation. Rarely biologically meaningful on its own;
 * it is here because SMS has it and because it makes `reverseComplement`
 * decomposable.
 *
 * @param s {string} any sequence
 * @return {string} the residues in reverse order
 * @throws {Error} kind "bioinformatics" when the sequence is empty or holds a character
 *         that is not a residue
 * @example bioseq.reverse("ATGC") == "CGTA"
 */
export func reverse(s as string) {
    def clean as string init requireSequence($s);
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def out as bytes;
    for (def i as int init len($buf) - 1; $i >= 0; $i = $i - 1) {
        $out[] = $buf[$i];
    }
    return convert.stringFromBytes($out, "utf-8");
}

/**
 * Does this sequence spell itself with `U` rather than `T`?
 *
 * A sequence with neither is read as DNA, which is the convention every tool
 * follows. One holding both cannot be assigned an alphabet, and rather than
 * guess, this raises.
 *
 * @param clean {string} an already-validated nucleotide sequence
 * @return {bool} true when the sequence is RNA-spelled
 * @throws {Error} kind "bioinformatics" when the sequence holds both `T` and `U`
 * @internal
 */
func isRnaSpelling(clean as string) {
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    def sawT as bool init false;
    def sawU as bool init false;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        if ($buf[$i] == 84) {
            $sawT = true;
        }
        if ($buf[$i] == 85) {
            $sawU = true;
        }
    }
    if ($sawT and $sawU) {
        raiseError("sequence mixes T and U, so it is neither DNA nor RNA" +
            " (use transcribe or backTranscribe to settle on one)");
    }
    return $sawU;
}

/**
 * The shared complement loop, forward or reversed.
 *
 * @param clean {string} an already-validated nucleotide sequence
 * @param rna {bool} true to complement in the RNA alphabet
 * @param backwards {bool} true to emit the result reversed
 * @return {string} the complemented sequence
 * @internal
 */
func complementBytes(clean as string, rna as bool, backwards as bool) {
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    def out as bytes;
    if ($backwards) {
        for (def i as int init $n - 1; $i >= 0; $i = $i - 1) {
            if ($rna) {
                $out[] = COMPLEMENT_RNA_CODES[$buf[$i]];
            } else {
                $out[] = COMPLEMENT_DNA_CODES[$buf[$i]];
            }
        }
    } else {
        for (def i as int init 0; $i < $n; $i = $i + 1) {
            if ($rna) {
                $out[] = COMPLEMENT_RNA_CODES[$buf[$i]];
            } else {
                $out[] = COMPLEMENT_DNA_CODES[$buf[$i]];
            }
        }
    }
    return convert.stringFromBytes($out, "utf-8");
}

# -------------------------------------------------------------- transcription

/**
 * Transcribe DNA to RNA: every `T` becomes `U`.
 *
 * This is the coding-strand shorthand, not the biological event - it does not
 * complement, because the mRNA sequence matches the coding strand and only the
 * thymines are spelled differently.
 *
 * @param dna {string} a DNA sequence
 * @return {string} the RNA spelling
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not DNA
 * @example bioseq.transcribe("ATGC") == "AUGC"
 */
export func transcribe(dna as string) {
    def clean as string init requireDna($dna);
    return swapCode($clean, 84, 85);
}

/**
 * Reverse-transcribe RNA to DNA: every `U` becomes `T`.
 *
 * @param rna {string} an RNA sequence
 * @return {string} the DNA spelling
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not RNA
 * @example bioseq.backTranscribe("AUGC") == "ATGC"
 */
export func backTranscribe(rna as string) {
    def clean as string init requireRna($rna);
    return swapCode($clean, 85, 84);
}

/**
 * Replace one byte value with another across a sequence.
 *
 * @param clean {string} an already-validated sequence
 * @param oldCode {int} the ASCII code to replace
 * @param newCode {int} the ASCII code to write instead
 * @return {string} the rewritten sequence
 * @internal
 */
func swapCode(clean as string, oldCode as int, newCode as int) {
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    def out as bytes;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        if ($buf[$i] == $oldCode) {
            $out[] = $newCode;
        } else {
            $out[] = $buf[$i];
        }
    }
    return convert.stringFromBytes($out, "utf-8");
}

# --------------------------------------------------------------- translation

/**
 * Translate a nucleotide sequence under the standard genetic code.
 *
 * Reads codons from the first base. A trailing one or two bases cannot form a
 * codon and are dropped. A stop codon becomes `*` and does **not** end the
 * translation - the whole sequence is read, because a caller looking for open
 * reading frames needs to see where the stops fall. A codon holding an
 * ambiguity code translates to `X` unless every base it could stand for gives
 * the same residue, in which case that residue is used.
 *
 * RNA is accepted directly; `U` is read as `T` for the lookup.
 *
 * @param nt {string} a DNA or RNA sequence
 * @return {string} the protein, stops as `*`
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic acid
 * @example bioseq.translate("ATGACGTTACATGGATAA") == "MTLHG*"
 */
export func translate(nt as string) {
    return translateWith($nt, STANDARD);
}

/**
 * Translate under a named genetic code.
 *
 * The alternative codes matter more often than they look: a mitochondrial gene
 * translated under the standard code stops early at every `TGA`, and a ciliate
 * gene never stops at all. Pass the table for the organism, or `codonTable(id)`
 * when the id is carried in data.
 *
 * @param nt {string} a DNA or RNA sequence
 * @param table {map of string to string} a genetic code, e.g. `bioseq.BACTERIAL`
 * @return {string} the protein, stops as `*`
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic acid
 * @example bioseq.translateWith($mito, bioseq.VERTEBRATE_MITOCHONDRIAL)
 */
export func translateWith(nt as string, table as map of string to string) {
    def clean as string init requireNucleotide($nt);
    return translateClean($clean, $table);
}

/**
 * Translate one reading frame.
 *
 * Frames `1`, `2` and `3` start at the first, second and third base of the
 * forward strand. Frames `-1`, `-2` and `-3` do the same on the reverse
 * complement, so frame `-1` starts at the last base and reads back.
 *
 * @param nt {string} a DNA or RNA sequence
 * @param frame {int} one of 1, 2, 3, -1, -2, -3
 * @return {string} the protein for that frame
 * @throws {Error} kind "bioinformatics" when the frame is not one of the six, or the
 *         sequence is empty or not a nucleic acid
 * @example bioseq.translateFrame("GATGACGTTACATGGATAA", 2) == "MTLHG*"
 */
export func translateFrame(nt as string, frame as int) {
    return translateFrameWith($nt, $frame, STANDARD);
}

/**
 * Translate one reading frame under a named genetic code.
 *
 * @param nt {string} a DNA or RNA sequence
 * @param frame {int} one of 1, 2, 3, -1, -2, -3
 * @param table {map of string to string} a genetic code
 * @return {string} the protein for that frame
 * @throws {Error} kind "bioinformatics" when the frame is not one of the six, or the
 *         sequence is empty or not a nucleic acid
 */
export func translateFrameWith(nt as string, frame as int, table as map of string to string) {
    if ($frame == 0 or $frame < -3 or $frame > 3) {
        raiseError("frame must be one of 1, 2, 3, -1, -2, -3, got " + convert.toString($frame));
    }
    def clean as string init requireNucleotide($nt);
    def strand as string init $clean;
    def offset as int init $frame - 1;
    if ($frame < 0) {
        $strand = complementBytes($clean, isRnaSpelling($clean), true);
        $offset = -$frame - 1;
    }
    if ($offset >= len($strand)) {
        return "";
    }
    return translateClean($strand[$offset..], $table);
}

/**
 * Translate all six reading frames under the standard genetic code.
 *
 * Three forward frames then three reverse, in the order `1, 2, 3, -1, -2, -3`.
 * This is what an ORF search reads, and what you look at when a sequence's
 * frame is unknown.
 *
 * @param nt {string} a DNA or RNA sequence
 * @return {list of Frame} one record per frame
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic acid
 */
export func sixFrame(nt as string) {
    return sixFrameWith($nt, STANDARD);
}

/**
 * Translate all six reading frames under a named genetic code.
 *
 * @param nt {string} a DNA or RNA sequence
 * @param table {map of string to string} a genetic code
 * @return {list of Frame} one record per frame
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic acid
 */
export func sixFrameWith(nt as string, table as map of string to string) {
    def clean as string init requireNucleotide($nt);
    def rna as bool init isRnaSpelling($clean);
    def rc as string init complementBytes($clean, $rna, true);
    def out as list of Frame;
    for (def f as int init 1; $f <= 3; $f = $f + 1) {
        $out[] = frameRecord($f, $clean, $f - 1, $table);
    }
    for (def f as int init 1; $f <= 3; $f = $f + 1) {
        $out[] = frameRecord(-$f, $rc, $f - 1, $table);
    }
    return $out;
}

/**
 * Build one `Frame` record from an already-oriented strand.
 *
 * @param frame {int} the frame number to record
 * @param strand {string} the forward or reverse-complemented sequence
 * @param offset {int} how many bases to skip before reading codons
 * @param table {map of string to string} a genetic code
 * @return {Frame} the record
 * @internal
 */
func frameRecord(frame as int, strand as string, offset as int, table as map of string to string) {
    if ($offset >= len($strand)) {
        return Frame{frame: $frame, protein: ""};
    }
    return Frame{frame: $frame, protein: translateClean($strand[$offset..], $table)};
}

/**
 * One reading frame's translation.
 *
 * @field frame {int} 1, 2 or 3 for the forward strand; -1, -2 or -3 for the
 *        reverse complement
 * @field protein {string} the translated residues, stops as `*`
 */
export def struct Frame {
    frame as int,
    protein as string
};

/**
 * The codon loop, over an already-validated sequence.
 *
 * @param clean {string} a validated nucleotide sequence
 * @param table {map of string to string} a genetic code
 * @return {string} the protein
 * @internal
 */
func translateClean(clean as string, table as map of string to string) {
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf) - len($buf) % 3;
    def out as bytes;
    for (def i as int init 0; $i < $n; $i = $i + 3) {
        def kb as bytes;
        # U reads as T so one table serves both alphabets
        def a as int init $buf[$i];
        def b as int init $buf[$i + 1];
        def c as int init $buf[$i + 2];
        if ($a == 85) {
            $a = 84;
        }
        if ($b == 85) {
            $b = 84;
        }
        if ($c == 85) {
            $c = 84;
        }
        $kb[] = $a;
        $kb[] = $b;
        $kb[] = $c;
        def codon as string init convert.stringFromBytes($kb, "utf-8");
        if (maps.has($table, $codon)) {
            def aa as bytes init convert.bytesFromString($table[$codon], "utf-8");
            $out[] = $aa[0];
        } else {
            def resolved as string init resolveAmbiguous($codon, $table);
            def aa as bytes init convert.bytesFromString($resolved, "utf-8");
            $out[] = $aa[0];
        }
    }
    return convert.stringFromBytes($out, "utf-8");
}

/**
 * Translate a codon holding ambiguity codes.
 *
 * Expands each position and reads every codon the ambiguity admits. When they
 * all give the same residue the answer is unambiguous despite the ambiguity -
 * `MGN` is arginine in every expansion, and `GGN` is glycine - so that residue
 * is returned rather than `X`. Otherwise the result is genuinely unknown.
 *
 * The expansion is at most 4x4x4, so this stays cheap even in the worst case,
 * and it is only reached for codons that missed the table outright.
 *
 * @param codon {string} a three-character codon, possibly ambiguous
 * @param table {map of string to string} a genetic code
 * @return {string} the residue, or "X" when the codon admits more than one
 * @internal
 */
func resolveAmbiguous(codon as string, table as map of string to string) {
    def chars as list of string init strings.chars($codon);
    if (not maps.has(IUPAC_DNA, $chars[0]) or not maps.has(IUPAC_DNA, $chars[1]) or
        not maps.has(IUPAC_DNA, $chars[2])) {
        return "X";
    }
    def first as string init "";
    for (def resolved in expandCodon($chars)) {
        if (not maps.has($table, $resolved)) {
            return "X";
        }
        def aa as string init $table[$resolved];
        if ($first == "") {
            $first = $aa;
        } elseif ($aa != $first) {
            return "X";
        }
    }
    return $first;
}

/**
 * Every unambiguous codon an ambiguous one admits.
 *
 * At most 4x4x4 = 64 codons, and only reached for a codon that missed the
 * table outright, so the cost never shows up in a real translation.
 *
 * @param chars {list of string} the codon's three IUPAC characters
 * @return {list of string} the concrete codons
 * @internal
 */
func expandCodon(chars as list of string) {
    def out as list of string;
    for (def a in strings.chars(IUPAC_DNA[$chars[0]])) {
        for (def b in strings.chars(IUPAC_DNA[$chars[1]])) {
            for (def c in strings.chars(IUPAC_DNA[$chars[2]])) {
                $out[] = $a + $b + $c;
            }
        }
    }
    return $out;
}

/**
 * Split a nucleotide sequence into codons.
 *
 * A trailing one or two bases are dropped, the same rule `translate` follows.
 *
 * @param nt {string} a DNA or RNA sequence
 * @return {list of string} three-base codons in reading order
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic acid
 * @example bioseq.splitCodons("ATGACGT") == ["ATG", "ACG"]
 */
export func splitCodons(nt as string) {
    def clean as string init requireNucleotide($nt);
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf) - len($buf) % 3;
    def out as list of string;
    for (def i as int init 0; $i < $n; $i = $i + 3) {
        def kb as bytes;
        $kb[] = $buf[$i];
        $kb[] = $buf[$i + 1];
        $kb[] = $buf[$i + 2];
        $out[] = convert.stringFromBytes($kb, "utf-8");
    }
    return $out;
}

# ---------------------------------------------------------- residue spellings

/**
 * Rewrite a protein in three-letter codes.
 *
 * `MKV` becomes `MetLysVal`. The ambiguity codes and the stop symbol have
 * spellings too - `X` is `Xaa`, `*` is `Ter` - so any sequence `translate`
 * produces can be rewritten.
 *
 * @param protein {string} a protein sequence in one-letter codes
 * @return {string} the same sequence in three-letter codes, unseparated
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a protein
 * @example bioseq.oneToThree("MKV") == "MetLysVal"
 */
export func oneToThree(protein as string) {
    def clean as string init requireProtein($protein);
    def out as bytes;
    for (def one in strings.chars($clean)) {
        def spelled as bytes init convert.bytesFromString(THREE_LETTER[$one], "utf-8");
        $out[] = $spelled[0];
        $out[] = $spelled[1];
        $out[] = $spelled[2];
    }
    return convert.stringFromBytes($out, "utf-8");
}

/**
 * Read a three-letter protein back into one-letter codes.
 *
 * The inverse of `oneToThree`. Codes may be run together (`MetLysVal`) or
 * separated by spaces or hyphens (`Met-Lys-Val`), which is how they are usually
 * written in prose.
 *
 * @param protein {string} a protein sequence in three-letter codes
 * @return {string} the one-letter spelling
 * @throws {Error} kind "bioinformatics" when the length is not a multiple of three or a
 *         group is not a residue code
 * @example bioseq.threeToOne("Met-Lys-Val") == "MKV"
 */
export func threeToOne(protein as string) {
    # Walk the bytes rather than slicing the string every three characters: a
    # string slice is rune-indexed, so taking one at offset i costs O(i) and the
    # whole loop goes quadratic. Separators are dropped in the same pass.
    def buf as bytes init convert.bytesFromString($protein, "utf-8");
    def n as int init len($buf);
    def kept as bytes;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def code as int init $buf[$i];
        if ($code == 45 or $code == 32 or $code == 9 or $code == 10 or $code == 13) {
            continue;
        }
        if ($code >= 97 and $code <= 122) {
            $code = $code - 32;
        }
        $kept[] = $code;
    }
    def total as int init len($kept);
    if ($total == 0) {
        raiseError("empty protein sequence");
    }
    if ($total % 3 != 0) {
        raiseError("a three-letter protein has a length divisible by three, got " +
            convert.toString($total));
    }
    def out as bytes;
    for (def i as int init 0; $i < $total; $i = $i + 3) {
        def group as bytes;
        $group[] = $kept[$i];
        $group[] = $kept[$i + 1];
        $group[] = $kept[$i + 2];
        def code as string init convert.stringFromBytes($group, "utf-8");
        if (not maps.has(ONE_LETTER, $code)) {
            raiseError("not a three-letter residue code: '" + $code + "' at position " +
                convert.toString($i));
        }
        def one as bytes init convert.bytesFromString(ONE_LETTER[$code], "utf-8");
        $out[] = $one[0];
    }
    return convert.stringFromBytes($out, "utf-8");
}
