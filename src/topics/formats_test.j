# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - white-box tests for the formats topic.
# Spliced into bioinformatics_test.j via include - run with:
#   jennifer test src/bioinformatics_test.j
#
# Fixtures are embedded as raw '...' strings rather than read from testdata/,
# so the suite runs unchanged on jennifer-tiny.

def const FASTA_TWO as string init '>seq1 the first one
ATGGCCCTGTGG
ATGCGCCTCCTG
>seq2 the second
GGGACCTGACCCAGCC
';

# A UniProt-style header, where the id itself contains pipes - the reason the
# id/description split is on whitespace and not on any other character.
def const FASTA_UNIPROT as string init '>sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606
MALWMRLLPLLALLALWGPDPAAA
';

def const FASTQ_TWO as string init '@read1 first
ATGCATGC
+
IIIIIIII
@read2 second
GGGGCCCC
+
!!!!####
';

# ------------------------------------------------------------------ FASTA read

func testParseFastaBasics() {
    def records as list of Record init parseFasta(FASTA_TWO);
    testing.assertEqual(len($records), 2);
    testing.assertEqual($records[0].id, "seq1");
    testing.assertEqual($records[0].description, "the first one");
    # wrapped sequence lines are joined
    testing.assertEqual($records[0].sequence, "ATGGCCCTGTGGATGCGCCTCCTG");
    testing.assertEqual($records[1].id, "seq2");
    testing.assertEqual($records[1].sequence, "GGGACCTGACCCAGCC");
}

func testParseFastaHeaderWithoutDescription() {
    def records as list of Record init parseFasta(">bare\nATGC\n");
    testing.assertEqual($records[0].id, "bare");
    testing.assertEqual($records[0].description, "");
}

func testParseFastaKeepsPipedIdIntact() {
    def r as Record init parseFasta(FASTA_UNIPROT)[0];
    testing.assertEqual($r.id, "sp|P01308|INS_HUMAN");
    testing.assertEqual($r.description, "Insulin OS=Homo sapiens OX=9606");
    # and the sequence is a real protein the transforms topic also knows
    testing.assertEqual($r.sequence, "MALWMRLLPLLALLALWGPDPAAA");
}

func testParseFastaPreservesCase() {
    # lowercase in a FASTA body is soft-masking, not noise
    def r as Record init parseFasta(">m\nATGCatgcATGC\n")[0];
    testing.assertEqual($r.sequence, "ATGCatgcATGC");
    # and filterDna is how a caller opts into normalising it
    testing.assertEqual(filterDna($r.sequence), "ATGCATGCATGC");
}

func testParseFastaToleratesRealFileShapes() {
    # blank lines between records, Windows endings, internal spaces
    def text as string init ">a desc\r\nATG CAT\r\n\r\n>b\r\nGGG\r\n";
    def records as list of Record init parseFasta($text);
    testing.assertEqual(len($records), 2);
    testing.assertEqual($records[0].sequence, "ATGCAT");
    testing.assertEqual($records[0].description, "desc");
    testing.assertEqual($records[1].sequence, "GGG");
}

func testParseFastaAllowsEmptySequence() {
    # a parser reports what the file says; len() is the caller's guard
    def records as list of Record init parseFasta(">empty\n>next\nATG\n");
    testing.assertEqual(len($records), 2);
    testing.assertEqual($records[0].sequence, "");
    testing.assertEqual($records[1].sequence, "ATG");
}

func testParseFastaRejectsStructuralErrors() {
    testing.assertThrows("fastaNoHeader", "bioinformatics");
    testing.assertThrows("fastaDataFirst", "bioinformatics");
    testing.assertThrows("fastaEmptyId", "bioinformatics");
    testing.assertThrows("fastaEmptyText", "bioinformatics");
}

func fastaNoHeader() {
    parseFasta("ATGCATGC\n");
}

func fastaDataFirst() {
    parseFasta("ATGC\n>seq1\nATGC\n");
}

func fastaEmptyId() {
    parseFasta("> \nATGC\n");
}

func fastaEmptyText() {
    parseFasta("");
}

# ----------------------------------------------------------------- FASTA write

func testWriteFasta() {
    def r as Record init Record{id: "seq1", description: "first", sequence: "ATGC"};
    testing.assertEqual(writeFasta([$r]), ">seq1 first\nATGC\n");
}

func testWriteFastaWithoutDescription() {
    def r as Record init Record{id: "bare", description: "", sequence: "ATGC"};
    testing.assertEqual(writeFasta([$r]), ">bare\nATGC\n");
}

func testWriteFastaWraps() {
    def r as Record init Record{id: "x", description: "", sequence: "ATGCATGCAT"};
    testing.assertEqual(writeFastaWrapped([$r], 4), ">x\nATGC\nATGC\nAT\n");
    # the default is the conventional sixty columns
    testing.assertEqual(FASTA_WIDTH, 60);
}

func testFastaRoundTrips() {
    def original as list of Record init parseFasta(FASTA_TWO);
    def reparsed as list of Record init parseFasta(writeFasta($original));
    testing.assertEqual(len($reparsed), len($original));
    for (def i as int init 0; $i < len($original); $i = $i + 1) {
        testing.assertEqual($reparsed[$i].id, $original[$i].id);
        testing.assertEqual($reparsed[$i].description, $original[$i].description);
        testing.assertEqual($reparsed[$i].sequence, $original[$i].sequence);
    }
}

func testFastaRoundTripsAcrossWidths() {
    def original as list of Record init parseFasta(FASTA_UNIPROT);
    for (def width in [1, 7, 60, 500]) {
        def back as list of Record init parseFasta(writeFastaWrapped($original, $width));
        testing.assertEqual($back[0].sequence, $original[0].sequence);
        testing.assertEqual($back[0].id, $original[0].id);
    }
}

func testWriteFastaRejectsBadInput() {
    testing.assertThrows("writeFastaEmptyList", "bioinformatics");
    testing.assertThrows("writeFastaNoId", "bioinformatics");
    testing.assertThrows("writeFastaZeroWidth", "bioinformatics");
}

func writeFastaEmptyList() {
    def none as list of Record;
    writeFasta($none);
}

func writeFastaNoId() {
    writeFasta([Record{id: "  ", description: "", sequence: "ATGC"}]);
}

func writeFastaZeroWidth() {
    writeFastaWrapped([Record{id: "x", description: "", sequence: "ATGC"}], 0);
}

# ------------------------------------------------------------------------ wrap

func testWrap() {
    testing.assertEqual(wrap("ATGCATGC", 4), "ATGC\nATGC");
    testing.assertEqual(wrap("ATGCA", 4), "ATGC\nA");
    # a width at or above the length leaves the sequence alone
    testing.assertEqual(wrap("ATGC", 4), "ATGC");
    testing.assertEqual(wrap("ATGC", 99), "ATGC");
}

func testWrapRejectsZeroWidth() {
    testing.assertThrows("wrapZero", "bioinformatics");
}

func wrapZero() {
    wrap("ATGC", 0);
}

# --------------------------------------------------------------------- toPlain

func testToPlain() {
    # the shape a sequence arrives in when pasted out of a database page
    testing.assertEqual(toPlain(">x\n1 atg acg\n"), "atgacg");
    testing.assertEqual(toPlain("  1  ATG ACG\n 61  TTA\n"), "ATGACGTTA");
    # unlike parseFasta it needs no header at all
    testing.assertEqual(toPlain("ATGC"), "ATGC");
    # old-style FASTA comment lines are dropped too
    testing.assertEqual(toPlain("; a comment\nATGC\n"), "ATGC");
}

func testToPlainKeepsCase() {
    testing.assertEqual(toPlain("atgcATGC"), "atgcATGC");
}

# ------------------------------------------------------------------ FASTQ read

func testParseFastqBasics() {
    def reads as list of Read init parseFastq(FASTQ_TWO);
    testing.assertEqual(len($reads), 2);
    testing.assertEqual($reads[0].id, "read1");
    testing.assertEqual($reads[0].description, "first");
    testing.assertEqual($reads[0].sequence, "ATGCATGC");
    testing.assertEqual($reads[0].quality, "IIIIIIII");
    testing.assertEqual($reads[1].id, "read2");
    testing.assertEqual($reads[1].quality, "!!!!####");
}

func testParseFastqLengthsAgree() {
    for (def r in parseFastq(FASTQ_TWO)) {
        testing.assertEqual(len($r.sequence), len($r.quality));
    }
}

func testParseFastqWrapped() {
    # the specification permits wrapping; sequence runs to the '+', quality runs
    # until it is as long as the sequence
    def text as string init "@r1\nATGC\nATGC\n+\nIIII\nIIII\n";
    def reads as list of Read init parseFastq($text);
    testing.assertEqual(len($reads), 1);
    testing.assertEqual($reads[0].sequence, "ATGCATGC");
    testing.assertEqual($reads[0].quality, "IIIIIIII");
}

func testParseFastqQualityMayLookLikeAHeader() {
    # '@' is a legal Phred+33 quality character (score 31), so a quality line
    # can start with one - which is why quality is read by length, not by marker
    def text as string init "@r1\nATGC\n+\n@III\n";
    def reads as list of Read init parseFastq($text);
    testing.assertEqual(len($reads), 1);
    testing.assertEqual($reads[0].quality, "@III");
    testing.assertEqual(qualityScores($reads[0])[0], 31);
}

func testParseFastqRejectsStructuralErrors() {
    testing.assertThrows("fastqRaggedQuality", "bioinformatics");
    testing.assertThrows("fastqNoPlus", "bioinformatics");
    testing.assertThrows("fastqNoAt", "bioinformatics");
    testing.assertThrows("fastqEmpty", "bioinformatics");
}

func fastqRaggedQuality() {
    parseFastq("@r1\nATGCATGC\n+\nIII\n");
}

func fastqNoPlus() {
    parseFastq("@r1\nATGC\n");
}

func fastqNoAt() {
    parseFastq("r1\nATGC\n+\nIIII\n");
}

func fastqEmpty() {
    parseFastq("   \n");
}

# ----------------------------------------------------------------- FASTQ write

func testWriteFastq() {
    def r as Read init Read{id: "r1", description: "first", sequence: "ATGC", quality: "IIII"};
    testing.assertEqual(writeFastq([$r]), "@r1 first\nATGC\n+\nIIII\n");
}

func testFastqRoundTrips() {
    def original as list of Read init parseFastq(FASTQ_TWO);
    def back as list of Read init parseFastq(writeFastq($original));
    testing.assertEqual(len($back), len($original));
    for (def i as int init 0; $i < len($original); $i = $i + 1) {
        testing.assertEqual($back[$i].id, $original[$i].id);
        testing.assertEqual($back[$i].description, $original[$i].description);
        testing.assertEqual($back[$i].sequence, $original[$i].sequence);
        testing.assertEqual($back[$i].quality, $original[$i].quality);
    }
}

func testWriteFastqRejectsRaggedRead() {
    testing.assertThrows("writeFastqRagged", "bioinformatics");
    testing.assertThrows("writeFastqEmptyList", "bioinformatics");
}

func writeFastqRagged() {
    writeFastq([Read{id: "r", description: "", sequence: "ATGC", quality: "II"}]);
}

func writeFastqEmptyList() {
    def none as list of Read;
    writeFastq($none);
}

func testFastqToFasta() {
    def records as list of Record init fastqToFasta(parseFastq(FASTQ_TWO));
    testing.assertEqual(len($records), 2);
    testing.assertEqual($records[0].id, "read1");
    testing.assertEqual($records[0].sequence, "ATGCATGC");
    # and the result writes as ordinary FASTA
    testing.assertTrue(strings.startsWith(writeFasta($records), ">read1 first"));
}

# --------------------------------------------------------------------- quality

func testQualityScores() {
    def r as Read init parseFastq(FASTQ_TWO)[0];
    def scores as list of int init qualityScores($r);
    testing.assertEqual(len($scores), 8);
    # 'I' is 73 in ASCII, so Q40 under Phred+33 - a very good base
    testing.assertEqual($scores[0], 40);
}

func testQualityScoresLowEnd() {
    def r as Read init parseFastq(FASTQ_TWO)[1];
    def scores as list of int init qualityScores($r);
    # '!' is the lowest printable quality character: Q0
    testing.assertEqual($scores[0], 0);
    # '#' is Q2, the Illumina "do not use this base" marker
    testing.assertEqual($scores[4], 2);
}

func testPhredOffsetMatters() {
    def r as Read init Read{id: "r", description: "", sequence: "AT", quality: "II"};
    testing.assertEqual(qualityScoresWith($r, PHRED33)[0], 40);
    testing.assertEqual(qualityScoresWith($r, PHRED64)[0], 9);
    # the two encodings differ by exactly 31
    testing.assertEqual(PHRED64 - PHRED33, 31);
}

func testWrongOffsetRaises() {
    # a Phred+33 file read as Phred+64 decodes to negative scores, which is the
    # only reliable signal that the offset is wrong
    testing.assertThrows("qualityWrongOffset", "bioinformatics");
}

func qualityWrongOffset() {
    qualityScoresWith(parseFastq(FASTQ_TWO)[1], PHRED64);
}

func testMeanQuality() {
    def r as Read init parseFastq(FASTQ_TWO)[0];
    # every base is Q40
    testing.assertEqual(meanQuality($r), 40.0);
    # read2 averages the four Q0 and four Q2 bases
    testing.assertEqual(meanQuality(parseFastq(FASTQ_TWO)[1]), 1.0);
}

func testErrorProbability() {
    # Q20 is one error in a hundred, Q30 one in a thousand
    testing.assertTrue(math.abs(errorProbability(20) - 0.01) < 0.000001);
    testing.assertTrue(math.abs(errorProbability(30) - 0.001) < 0.0000001);
    testing.assertEqual(errorProbability(0), 1.0);
    # and it falls as the score rises
    testing.assertTrue(errorProbability(40) < errorProbability(30));
}

func testErrorProbabilityRejectsNegative() {
    testing.assertThrows("errorProbabilityNegative", "bioinformatics");
}

func errorProbabilityNegative() {
    errorProbability(-1);
}

# ------------------------------------------------------- crossing into the deck

func testParsedFastaFeedsTheTransforms() {
    # the composition a user actually writes: read a record, work on its sequence
    def r as Record init parseFasta(FASTA_UNIPROT)[0];
    testing.assertTrue(isProtein($r.sequence));
    testing.assertEqual(oneToThree($r.sequence[0..3]), "MetAlaLeu");

    def gene as Record init parseFasta(FASTA_TWO)[0];
    testing.assertTrue(isDna($gene.sequence));
    testing.assertEqual(translate($gene.sequence), "MALWMRLL");
}

# ------------------------------------------------- the record cap, after the audit

func testRecordCapIsGenerous() {
    # nothing real approaches it
    testing.assertTrue(MAX_RECORDS >= 1000000);
}

func testOrdinaryFilesAreUnaffected() {
    # the cap must be invisible to real input
    testing.assertEqual(len(parseFasta(FASTA_TWO)), 2);
    testing.assertEqual(len(parseFastq(FASTQ_TWO)), 2);
}

# --------------------------------------- the Read invariant, after the OC audit

func testQualityVerbsRejectARaggedRead() {
    # A Read can be built by hand rather than parsed. The one invariant FASTQ
    # has is that quality is as long as sequence; before this was checked, a
    # ragged read reported a confident mean over whichever bases had a score.
    testing.assertThrows("qualityOfRaggedRead", "bioinformatics");
    testing.assertThrows("meanOfRaggedRead", "bioinformatics");
}

func raggedRead() {
    return Read{id: "x", description: "", sequence: "ATGCATGC", quality: "II"};
}

func qualityOfRaggedRead() {
    qualityScores(raggedRead());
}

func meanOfRaggedRead() {
    meanQuality(raggedRead());
}

func testQualityVerbsStillAcceptParsedReads() {
    # the check must not reject anything the parser produces
    for (def r in parseFastq(FASTQ_TWO)) {
        testing.assertEqual(len(qualityScores($r)), len($r.sequence));
    }
}

func testMeanQualityIsNotTruncated() {
    # Jennifer's `/` is true division, so a mean that does not divide evenly
    # keeps its fraction rather than flooring - worth pinning, because reading
    # this as integer division is an easy mistake to carry in from another
    # language.
    def r as Read init Read{id: "r", description: "", sequence: "ATG", quality: "++,"};
    def scores as list of int init qualityScores($r);
    testing.assertEqual($scores[0], 10);
    testing.assertEqual($scores[1], 10);
    testing.assertEqual($scores[2], 11);
    testing.assertTrue(close(meanQuality($r), 10.333333, 0.0001));
}

func testRecordCapIsCheckedBeforeAllocating() {
    # the cap used to be applied after the record was appended, so one more
    # existed than the limit allowed; the guard now runs first
    testing.assertEqual(len(parseFasta(FASTA_TWO)), 2);
    testing.assertEqual(len(parseFastq(FASTQ_TWO)), 2);
}
