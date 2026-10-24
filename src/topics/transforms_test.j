# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - white-box tests for the transforms topic.
# Spliced into bioinformatics_test.j via include - run with:
#   jennifer test src/bioinformatics_test.j

# INS_CDS and INS_PROTEIN are the shared golden fixtures declared in
# bioinformatics_test.j - a real gene and its published translation.

# ---------------------------------------------------------------- complement

func testComplement() {
    testing.assertEqual(complement("ATGC"), "TACG");
    testing.assertEqual(complement("atgc"), "TACG");
    # order is preserved - this is not the reverse complement
    testing.assertEqual(complement("AAAT"), "TTTA");
}

func testComplementKeepsAlphabet() {
    # a U anywhere makes it RNA, and the complement of A is U
    testing.assertEqual(complement("AUGC"), "UACG");
    testing.assertEqual(complement("ATGC"), "TACG");
    # neither T nor U reads as DNA, the universal convention
    testing.assertEqual(complement("ACGACG"), "TGCTGC");
}

func testComplementOfAmbiguityCodes() {
    # ambiguity codes complement as sets
    testing.assertEqual(complement("RYSWKM"), "YRSWMK");
    testing.assertEqual(complement("BDHV"), "VHDB");
    # the self-complementary ones are fixed points
    testing.assertEqual(complement("SWN"), "SWN");
}

func testReverseComplement() {
    testing.assertEqual(reverseComplement("ATGACGTTACATGGATAA"), "TTATCCATGTAACGTCAT");
    testing.assertEqual(reverseComplement("ATGC"), "GCAT");
    testing.assertEqual(reverseComplement("AUGC"), "GCAU");
}

func testReverseComplementIsAnInvolution() {
    def seqs as list of string init ["ATGACGTTACATGGATAA", "ACGT", "N", "RYSWKMBDHVN", INS_CDS];
    for (def s in $seqs) {
        testing.assertEqual(reverseComplement(reverseComplement($s)), $s);
    }
}

func testReverseComplementIsReverseThenComplement() {
    def s as string init "ATGACGTTACATGGATAA";
    testing.assertEqual(reverseComplement($s), reverse(complement($s)));
    testing.assertEqual(reverseComplement($s), complement(reverse($s)));
}

func testPalindromicSitesAreTheirOwnReverseComplement() {
    # which is exactly why type II enzymes recognise them
    testing.assertEqual(reverseComplement("GAATTC"), "GAATTC");
    testing.assertEqual(reverseComplement("GGATCC"), "GGATCC");
    testing.assertEqual(reverseComplement("GCGGCCGC"), "GCGGCCGC");
}

func testReverse() {
    testing.assertEqual(reverse("ATGC"), "CGTA");
    # works on protein too - it is a residue operation, not a strand operation
    testing.assertEqual(reverse("MKVL"), "LVKM");
    testing.assertEqual(reverse(reverse(INS_CDS)), INS_CDS);
}

func testMixedAlphabetRaises() {
    testing.assertThrows("complementMixed", "bioinformatics");
    testing.assertThrows("reverseComplementMixed", "bioinformatics");
}

func complementMixed() {
    complement("ATGCU");
}

func reverseComplementMixed() {
    reverseComplement("ATGCU");
}

func testComplementRejectsJunk() {
    testing.assertThrows("complementProtein", "bioinformatics");
    testing.assertThrows("complementEmpty", "bioinformatics");
}

func complementProtein() {
    complement("MKVL");
}

func complementEmpty() {
    complement("");
}

# ------------------------------------------------------------- transcription

func testTranscribe() {
    testing.assertEqual(transcribe("ATGC"), "AUGC");
    testing.assertEqual(transcribe("TTTT"), "UUUU");
    # only T changes; nothing is complemented
    testing.assertEqual(transcribe("ACG"), "ACG");
}

func testBackTranscribe() {
    testing.assertEqual(backTranscribe("AUGC"), "ATGC");
    testing.assertEqual(backTranscribe("UUUU"), "TTTT");
}

func testTranscriptionRoundTrips() {
    testing.assertEqual(backTranscribe(transcribe(INS_CDS)), INS_CDS);
    testing.assertEqual(transcribe(backTranscribe("AUGCAUGC")), "AUGCAUGC");
}

func testTranscribeRejectsWrongAlphabet() {
    testing.assertThrows("transcribeRna", "bioinformatics");
    testing.assertThrows("backTranscribeDna", "bioinformatics");
}

func transcribeRna() {
    transcribe("AUGC");
}

func backTranscribeDna() {
    backTranscribe("ATGC");
}

# --------------------------------------------------------------- translation

func testTranslateKnownGene() {
    # the published protein for a real coding sequence
    testing.assertEqual(translate(INS_CDS), INS_PROTEIN);
}

func testTranslateBasics() {
    testing.assertEqual(translate("ATGACGTTACATGGATAA"), "MTLHG*");
    testing.assertEqual(translate("ATG"), "M");
    testing.assertEqual(translate("TGG"), "W");
}

func testTranslateDropsPartialCodon() {
    # one or two trailing bases cannot form a codon
    testing.assertEqual(translate("ATGA"), "M");
    testing.assertEqual(translate("ATGAC"), "M");
    testing.assertEqual(translate("ATGACG"), "MT");
}

func testTranslateDoesNotStopAtAStop() {
    # a stop is reported, not obeyed - an ORF search needs to see them all
    testing.assertEqual(translate("ATGTAAATG"), "M*M");
    testing.assertEqual(len(translate("TAATAATAA")), 3);
}

func testTranslateAcceptsRna() {
    testing.assertEqual(translate("AUGACGUUACAUGGAUAA"), "MTLHG*");
    # and agrees with the DNA spelling of the same sequence
    testing.assertEqual(translate(transcribe(INS_CDS)), translate(INS_CDS));
}

func testTranslateAgreesWithSplitCodons() {
    def codons as list of string init splitCodons(INS_CDS);
    def rebuilt as string init "";
    for (def codon in $codons) {
        $rebuilt = $rebuilt + STANDARD[$codon];
    }
    testing.assertEqual($rebuilt, translate(INS_CDS));
}

func testTranslateAmbiguityResolvesWhenItCan() {
    # every codon MGN admits is arginine, so the answer is not ambiguous
    testing.assertEqual(translate("CGN"), "R");
    testing.assertEqual(translate("GGN"), "G");
    testing.assertEqual(translate("CCN"), "P");
    # but these genuinely could be several residues
    testing.assertEqual(translate("NNN"), "X");
    testing.assertEqual(translate("ATN"), "X");
    # a fully ambiguous third position on a two-fold codon still resolves
    testing.assertEqual(translate("TTY"), "F");
    testing.assertEqual(translate("TTR"), "L");
}

func testTranslateWithAlternativeCode() {
    # TGA is a stop under the standard code and tryptophan in mitochondria
    testing.assertEqual(translate("TGA"), "*");
    testing.assertEqual(translateWith("TGA", VERTEBRATE_MITOCHONDRIAL), "W");
    # AGA codes arginine normally, stops in vertebrate mitochondria
    testing.assertEqual(translate("AGA"), "R");
    testing.assertEqual(translateWith("AGA", VERTEBRATE_MITOCHONDRIAL), "*");
    # the ciliate code reads two of the three stops as glutamine
    testing.assertEqual(translateWith("TAATAG", CILIATE_NUCLEAR), "QQ");
    # bacterial assigns exactly as standard does
    testing.assertEqual(translateWith(INS_CDS, BACTERIAL), translate(INS_CDS));
}

func testSplitCodons() {
    testing.assertEqual(len(splitCodons("ATGACGT")), 2);
    testing.assertEqual(splitCodons("ATGACG")[0], "ATG");
    testing.assertEqual(splitCodons("ATGACG")[1], "ACG");
    # the codon count is the residue count
    testing.assertEqual(len(splitCodons(INS_CDS)), len(translate(INS_CDS)));
}

# -------------------------------------------------------------------- frames

func testTranslateFrameForward() {
    # padding the gene by one base moves its reading frame to 2
    def padded as string init "G" + INS_CDS;
    testing.assertEqual(translateFrame($padded, 2), INS_PROTEIN);
    testing.assertEqual(translateFrame(INS_CDS, 1), INS_PROTEIN);
}

func testTranslateFrameReverse() {
    # frame -1 of a sequence is frame 1 of its reverse complement
    testing.assertEqual(translateFrame(INS_CDS, -1), translate(reverseComplement(INS_CDS)));
}

func testFrameOffsets() {
    # frame N skips N-1 bases, so the ATG's index picks the frame that finds it
    def gene as string init "ATGACGTTACATGGATAA";
    testing.assertEqual(translateFrame("G" + $gene, 2), "MTLHG*");
    testing.assertEqual(translateFrame("GG" + $gene, 3), "MTLHG*");
    # and the wrong frame reads straight past it
    testing.assertTrue(translateFrame("GG" + $gene, 2) != "MTLHG*");
}

func testInvalidFrameRaises() {
    testing.assertThrows("frameZero", "bioinformatics");
    testing.assertThrows("frameFour", "bioinformatics");
    testing.assertThrows("frameMinusFour", "bioinformatics");
}

func frameZero() {
    translateFrame("ATGC", 0);
}

func frameFour() {
    translateFrame("ATGC", 4);
}

func frameMinusFour() {
    translateFrame("ATGC", -4);
}

func testTranslateFrameWith() {
    # TGA is a stop in the standard code and tryptophan in vertebrate mitochondria
    testing.assertEqual(translateFrameWith("GATGTGA", 2, STANDARD), "M*");
    testing.assertEqual(translateFrameWith("GATGTGA", 2, VERTEBRATE_MITOCHONDRIAL), "MW");
    # the reverse complement of TCACAT is ATGTGA
    testing.assertEqual(translateFrameWith("TCACAT", -1, VERTEBRATE_MITOCHONDRIAL), "MW");
    # a frame that starts past the end translates to nothing
    testing.assertEqual(translateFrameWith("AT", 3, STANDARD), "");
}

func testTranslateFrameWithRejectsBadFrame() {
    testing.assertThrows("frameWithZero", "bioinformatics");
}

func frameWithZero() {
    translateFrameWith("ATGC", 0, VERTEBRATE_MITOCHONDRIAL);
}

func testSixFrame() {
    def frames as list of Frame init sixFrame(INS_CDS);
    testing.assertEqual(len($frames), 6);
    # reported in the order 1, 2, 3, -1, -2, -3
    testing.assertEqual($frames[0].frame, 1);
    testing.assertEqual($frames[2].frame, 3);
    testing.assertEqual($frames[3].frame, -1);
    testing.assertEqual($frames[5].frame, -3);
}

func testSixFrameAgreesWithTranslateFrame() {
    for (def f in sixFrame(INS_CDS)) {
        testing.assertEqual($f.protein, translateFrame(INS_CDS, $f.frame));
    }
}

func testSixFrameFrameOneIsPlainTranslation() {
    testing.assertEqual(sixFrame(INS_CDS)[0].protein, translate(INS_CDS));
}

func testSixFrameOnAShortSequence() {
    # two bases cannot fill any frame's first codon
    def frames as list of Frame init sixFrame("AT");
    testing.assertEqual(len($frames), 6);
    for (def f in $frames) {
        testing.assertEqual($f.protein, "");
    }
}

func testSixFrameWithAlternativeCode() {
    def frames as list of Frame init sixFrameWith(INS_CDS, VERTEBRATE_MITOCHONDRIAL);
    testing.assertEqual(len($frames), 6);
    testing.assertEqual($frames[0].protein, translateWith(INS_CDS, VERTEBRATE_MITOCHONDRIAL));
}

# --------------------------------------------------------- residue spellings

func testOneToThree() {
    testing.assertEqual(oneToThree("MKV"), "MetLysVal");
    testing.assertEqual(oneToThree("M"), "Met");
    # the stop symbol and the ambiguity codes have spellings too
    testing.assertEqual(oneToThree("*"), "Ter");
    testing.assertEqual(oneToThree("X"), "Xaa");
    testing.assertEqual(oneToThree("BZJ"), "AsxGlxXle");
}

func testThreeToOne() {
    testing.assertEqual(threeToOne("MetLysVal"), "MKV");
    # separators are how three-letter codes are usually written in prose
    testing.assertEqual(threeToOne("Met-Lys-Val"), "MKV");
    testing.assertEqual(threeToOne("Met Lys Val"), "MKV");
    testing.assertEqual(threeToOne("MET LYS VAL"), "MKV");
    testing.assertEqual(threeToOne("Ter"), "*");
}

func testSpellingsRoundTrip() {
    testing.assertEqual(threeToOne(oneToThree(INS_PROTEIN)), INS_PROTEIN);
    testing.assertEqual(threeToOne(oneToThree("MKV*X")), "MKV*X");
    # and every code in the alphabet survives the trip
    for (def one in strings.chars(PROTEIN_ALPHABET)) {
        testing.assertEqual(threeToOne(oneToThree($one)), $one);
    }
}

func testThreeToOneRejectsBadInput() {
    testing.assertThrows("threeToOneRagged", "bioinformatics");
    testing.assertThrows("threeToOneUnknown", "bioinformatics");
    testing.assertThrows("threeToOneEmpty", "bioinformatics");
}

func threeToOneRagged() {
    threeToOne("MetLy");
}

func threeToOneUnknown() {
    threeToOne("MetQqqVal");
}

func threeToOneEmpty() {
    threeToOne("   ");
}

func testTranslateThenSpellIsChainable() {
    # the composition a user actually writes
    testing.assertEqual(oneToThree(translate("ATGACGTTACATGGATAA")), "MetThrLeuHisGlyTer");
}
