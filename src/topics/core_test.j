# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - white-box tests for the core topic.
# Spliced into bioinformatics_test.j via include - run with:
#   jennifer test src/bioinformatics_test.j

# The byte-code membership sets in core.j are written as ASCII literals for
# speed. This derives them back from the human-readable alphabet strings and
# fails if the two ever drift apart, so the magic numbers cannot go stale.
func codesOf(alphabet as string) {
    def buf as bytes init convert.bytesFromString($alphabet, "utf-8");
    def out as map of int to bool;
    for (def i as int init 0; $i < len($buf); $i = $i + 1) {
        $out[$buf[$i]] = true;
    }
    return $out;
}

func sameCodes(want as map of int to bool, got as map of int to bool) {
    if (len($want) != len($got)) {
        return false;
    }
    for (def k in $want) {
        if (not maps.has($got, $k)) {
            return false;
        }
    }
    return true;
}

func testAlphabetsAgree() {
    testing.assertTrue(sameCodes(codesOf(DNA_ALPHABET), DNA_CODES));
    testing.assertTrue(sameCodes(codesOf(RNA_ALPHABET), RNA_CODES));
    testing.assertTrue(sameCodes(codesOf(PROTEIN_ALPHABET), PROTEIN_CODES));

    # NUCLEOTIDE_CODES is exactly the union of the two nucleic-acid alphabets.
    def union as map of int to bool init codesOf(DNA_ALPHABET + RNA_ALPHABET);
    testing.assertTrue(sameCodes($union, NUCLEOTIDE_CODES));
}

func testAlphabetShapes() {
    # four bases plus eleven IUPAC ambiguity codes
    testing.assertEqual(len(DNA_ALPHABET), 15);
    testing.assertEqual(len(RNA_ALPHABET), 15);
    # twenty standard residues, BZXJ, UO, and the stop character
    testing.assertEqual(len(PROTEIN_ALPHABET), 27);
    testing.assertTrue(strings.contains(DNA_ALPHABET, "T"));
    testing.assertFalse(strings.contains(DNA_ALPHABET, "U"));
    testing.assertTrue(strings.contains(RNA_ALPHABET, "U"));
    testing.assertFalse(strings.contains(RNA_ALPHABET, "T"));
    testing.assertTrue(strings.contains(PROTEIN_ALPHABET, "*"));
}

# ------------------------------------------------------------------- filtering

func testFilterDnaStripsAndUppercases() {
    testing.assertEqual(filterDna("atgc"), "ATGC");
    testing.assertEqual(filterDna("1 atg acg > x"), "ATGACG");
    testing.assertEqual(filterDna("ATG\nACG\r\n"), "ATGACG");
    # U is RNA, so filterDna drops it rather than guessing
    testing.assertEqual(filterDna("AUGC"), "AGC");
    # nothing recognisable is an empty result, not an error
    testing.assertEqual(filterDna("...???..."), "");
}

func testFilterKeepsAmbiguityCodes() {
    testing.assertEqual(filterDna("ACGTRYSWKMBDHVN"), "ACGTRYSWKMBDHVN");
    testing.assertEqual(filterRna("acgurysWKMBDHVN"), "ACGURYSWKMBDHVN");
}

func testFilterRna() {
    testing.assertEqual(filterRna("augc"), "AUGC");
    # T is DNA, so filterRna drops it
    testing.assertEqual(filterRna("ATGC"), "AGC");
}

func testFilterProtein() {
    testing.assertEqual(filterProtein("mkvl"), "MKVL");
    testing.assertEqual(filterProtein("MKV*"), "MKV*");
    # the standard twenty survive intact
    testing.assertEqual(filterProtein("ACDEFGHIKLMNPQRSTVWY"), "ACDEFGHIKLMNPQRSTVWY");
    # a gap is not a residue
    testing.assertEqual(filterProtein("MK-VL"), "MKVL");
}

func testFilterIsIdempotent() {
    def messy as string init "  1 atg-acg 60 nnn  ";
    def once as string init filterDna($messy);
    testing.assertEqual(filterDna($once), $once);
}

# ------------------------------------------------------------------ predicates

func testIsDna() {
    testing.assertTrue(isDna("ATGC"));
    testing.assertTrue(isDna("atgc"));
    testing.assertTrue(isDna("ACGTRYSWKMBDHVN"));
    testing.assertFalse(isDna("AUGC"));
    testing.assertFalse(isDna("ATG C"));
    # an empty string is not a sequence
    testing.assertFalse(isDna(""));
}

func testIsRna() {
    testing.assertTrue(isRna("AUGC"));
    testing.assertFalse(isRna("ATGC"));
    testing.assertFalse(isRna(""));
}

func testIsProtein() {
    testing.assertTrue(isProtein("MKVLA"));
    testing.assertTrue(isProtein("MKV*"));
    testing.assertFalse(isProtein("MKV-"));
    testing.assertFalse(isProtein(""));
}

func testAcgtSatisfiesBothAlphabets() {
    # documented overlap: a four-residue peptide is spelled like an oligo
    testing.assertTrue(isDna("ACGT"));
    testing.assertTrue(isProtein("ACGT"));
}

func testFilterAndPredicateAgree() {
    # a non-empty filter result always satisfies the matching predicate
    def messy as string init "atg-acg 12 ttacatggataa";
    testing.assertTrue(isDna(filterDna($messy)));
    testing.assertTrue(isProtein(filterProtein("mkv--la*")));
}

# ------------------------------------------------------------------ validation

func testRequireDnaTolerates() {
    # whitespace, digits and case are the pasted-listing case, not an error
    testing.assertEqual(requireDna("  atg acg\n  60 tta  "), "ATGACGTTA");
    testing.assertEqual(requireDna("ATGC"), "ATGC");
}

func testRequireDnaRejects() {
    testing.assertThrows("requireDnaOnRna", "bioinformatics");
    testing.assertThrows("requireDnaOnGap", "bioinformatics");
    testing.assertThrows("requireDnaOnEmpty", "bioinformatics");
    testing.assertThrows("requireDnaOnDigitsOnly", "bioinformatics");
}

func requireDnaOnRna() {
    requireDna("AUGC");
}

func requireDnaOnGap() {
    requireDna("ATG-C");
}

func requireDnaOnEmpty() {
    requireDna("");
}

func requireDnaOnDigitsOnly() {
    # digits are dropped, which leaves nothing - an empty sequence, not a
    # sequence of length zero the caller meant
    requireDna("12345");
}

func testErrorMessageNamesTheCharacter() {
    def caught as string init "";
    try {
        requireDna("ATG?C");
    } catch (e) {
        $caught = $e.message;
    }
    testing.assertTrue(strings.contains($caught, "'?'"));
    testing.assertTrue(strings.contains($caught, "position 3"));
    testing.assertTrue(strings.contains($caught, "DNA"));
}

func testUnprintableCharacterIsDescribedByCode() {
    def caught as string init "";
    try {
        requireProtein("MK\u0007V");
    } catch (e) {
        $caught = $e.message;
    }
    testing.assertTrue(strings.contains($caught, "byte 7"));
}

func testRequireNucleotideAcceptsEither() {
    testing.assertEqual(requireNucleotide("ATGC"), "ATGC");
    testing.assertEqual(requireNucleotide("AUGC"), "AUGC");
    # and even the mixed spelling, which neither strict validator would take
    testing.assertEqual(requireNucleotide("ATGCU"), "ATGCU");
}

func testRequireProtein() {
    testing.assertEqual(requireProtein("mkv la"), "MKVLA");
    testing.assertThrows("requireProteinOnGap", "bioinformatics");
}

func requireProteinOnGap() {
    requireProtein("MK-VL");
}

# ---------------------------------------------------------------- molecule kind

func testKindName() {
    testing.assertEqual(kindName(DNA), "DNA");
    testing.assertEqual(kindName(RNA), "RNA");
    testing.assertEqual(kindName(PROTEIN), "protein");
}

func testKindsAreDistinct() {
    testing.assertTrue(DNA != RNA);
    testing.assertTrue(RNA != PROTEIN);
    testing.assertTrue(DNA != PROTEIN);
}

func testKindNameRejectsUnknown() {
    testing.assertThrows("kindNameOnZero", "bioinformatics");
    testing.assertThrows("kindNameOnHigh", "bioinformatics");
}

func kindNameOnZero() {
    kindName(0);
}

func kindNameOnHigh() {
    kindName(99);
}
