# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - white-box tests for the compare topic.
# Spliced into bioinformatics_test.j via include - run with:
#   jennifer test src/bioinformatics_test.j
#
# Alignment is easy to get subtly wrong and hard to eyeball, so the checks here
# lean on invariants that must hold for any correct implementation - the aligned
# rows are the originals with gaps removed, the reported score is the score of
# the alignment returned, identity is symmetric - rather than on remembered
# output. Where an exact result is asserted it is one small enough to verify by
# hand, and the comment says what it should be and why.

func ungapped(s as string) {
    return strings.replace($s, "-", "");
}

# Recompute an alignment's score from its own rows, under the same options.
# If the engine's bookkeeping and its matrix ever disagree, this catches it.
func rescore(al as Alignment, opts as AlignOptions) {
    def total as int init 0;
    def useMatrix as bool init len($opts.matrix) > 0;
    for (def i as int init 0; $i < $al.length; $i = $i + 1) {
        def x as string init $al.a[$i..$i + 1];
        def y as string init $al.b[$i..$i + 1];
        if ($x == "-" or $y == "-") {
            $total = $total + $opts.gapScore;
        } elseif ($useMatrix) {
            $total = $total + substitutionScore($opts.matrix, $x, $y);
        } elseif ($x == $y) {
            $total = $total + $opts.matchScore;
        } else {
            $total = $total + $opts.mismatchScore;
        }
    }
    return $total;
}

# --------------------------------------------------------------- the matrix

func testBlosumAnchors() {
    # the values everyone who has read a matrix knows
    testing.assertEqual(substitutionScore(BLOSUM62, "A", "A"), 4);
    testing.assertEqual(substitutionScore(BLOSUM62, "W", "W"), 11);
    testing.assertEqual(substitutionScore(BLOSUM62, "C", "C"), 9);
    testing.assertEqual(substitutionScore(BLOSUM62, "H", "H"), 8);
    testing.assertEqual(substitutionScore(BLOSUM62, "P", "P"), 7);
}

func testBlosumIsSymmetric() {
    for (def x in strings.chars("ARNDCQEGHILKMFPSTWYV")) {
        for (def y in strings.chars("ARNDCQEGHILKMFPSTWYV")) {
            testing.assertEqual(
                substitutionScore(BLOSUM62, $x, $y),
                substitutionScore(BLOSUM62, $y, $x));
        }
    }
}

func testBlosumRewardsConservativeSubstitution() {
    # the whole reason a matrix beats match/mismatch for protein: some
    # substitutions barely change a protein and some destroy it
    testing.assertTrue(substitutionScore(BLOSUM62, "L", "I") > 0);
    testing.assertTrue(substitutionScore(BLOSUM62, "K", "R") > 0);
    testing.assertTrue(substitutionScore(BLOSUM62, "F", "Y") > 0);
    testing.assertTrue(substitutionScore(BLOSUM62, "L", "D") < 0);
    testing.assertTrue(substitutionScore(BLOSUM62, "W", "P") < 0);
}

func testIdentityAlwaysBeatsSubstitution() {
    for (def x in strings.chars("ARNDCQEGHILKMFPSTWYV")) {
        for (def y in strings.chars("ARNDCQEGHILKMFPSTWYV")) {
            if ($x != $y) {
                testing.assertTrue(substitutionScore(BLOSUM62, $x, $x) >
                    substitutionScore(BLOSUM62, $x, $y));
            }
        }
    }
}

# ------------------------------------------------------------ global alignment

func testIdenticalSequencesAlignPerfectly() {
    def al as Alignment init alignGlobal(INS_CDS, INS_CDS, dnaOptions());
    testing.assertEqual($al.a, INS_CDS);
    testing.assertEqual($al.b, INS_CDS);
    testing.assertEqual($al.gaps, 0);
    testing.assertEqual($al.identities, len(INS_CDS));
    testing.assertEqual($al.length, len(INS_CDS));
    testing.assertTrue(close(percentIdentity($al), 100.0, 0.000001));
    # every column is a match, so the score is length times the match score
    testing.assertEqual($al.score, 2 * len(INS_CDS));
}

func testGlobalAlignmentKeepsBothSequencesWhole() {
    # THE INVARIANT: strip the gaps and you have the originals back
    def al as Alignment init alignGlobal("ATGCATGCTT", "ATGATGCTT", dnaOptions());
    testing.assertEqual(ungapped($al.a), "ATGCATGCTT");
    testing.assertEqual(ungapped($al.b), "ATGATGCTT");
    testing.assertEqual(len($al.a), len($al.b));
    testing.assertEqual(len($al.markup), len($al.a));
    testing.assertEqual($al.length, len($al.a));
}

func testGlobalScoreMatchesItsOwnRows() {
    # the engine's running score must equal the score of what it returned
    def opts as AlignOptions init dnaOptions();
    for (def pair in [
        ["ATGCATGC", "ATGCATGC"],
        ["ATGCATGC", "ATGAATGC"],
        ["ATGCATGCTT", "ATGATGCTT"],
        ["AAAA", "TTTT"],
        ["ATGCATGCATGC", "ATGC"]
    ]) {
        def al as Alignment init alignGlobal($pair[0], $pair[1], $opts);
        testing.assertEqual($al.score, rescore($al, $opts));
    }
}

func testSingleGapIsPlacedNotSpread() {
    # one deletion should cost one gap, not several
    def al as Alignment init alignGlobal("ATGCATGC", "ATGATGC", dnaOptions());
    testing.assertEqual($al.gaps, 1);
    testing.assertEqual($al.length, 8);
    testing.assertEqual($al.identities, 7);
}

func testGlobalIsSymmetricInItsScore() {
    # swapping the arguments swaps the rows but not the score
    def opts as AlignOptions init dnaOptions();
    def forward as Alignment init alignGlobal("ATGCATGCTT", "ATGATGCTT", $opts);
    def backward as Alignment init alignGlobal("ATGATGCTT", "ATGCATGCTT", $opts);
    testing.assertEqual($forward.score, $backward.score);
    testing.assertEqual($forward.identities, $backward.identities);
    testing.assertEqual($forward.a, $backward.b);
    testing.assertEqual($forward.b, $backward.a);
}

func testGlobalAlignsEveryResidue() {
    # a short sequence against a long one is padded to the long one's length,
    # which is exactly why this is the wrong verb for that question
    def al as Alignment init alignGlobal(INS_CDS, INS_CDS[0..12], dnaOptions());
    testing.assertEqual(ungapped($al.b), INS_CDS[0..12]);
    testing.assertEqual($al.gaps, len(INS_CDS) - 12);
    testing.assertTrue(percentIdentity($al) < 20.0);
}

func testMarkupAgreesWithTheRows() {
    def al as Alignment init alignGlobal(INS_PROTEIN, "MALWMRLLPLLGLLGLWGPDPAAA", proteinOptions());
    def identities as int init 0;
    def similarities as int init 0;
    for (def i as int init 0; $i < $al.length; $i = $i + 1) {
        def mark as string init $al.markup[$i..$i + 1];
        def x as string init $al.a[$i..$i + 1];
        def y as string init $al.b[$i..$i + 1];
        if ($mark == "|") {
            testing.assertEqual($x, $y);
            $identities = $identities + 1;
            $similarities = $similarities + 1;
        } elseif ($mark == ":") {
            testing.assertTrue($x != $y);
            testing.assertTrue(substitutionScore(BLOSUM62, $x, $y) > 0);
            $similarities = $similarities + 1;
        }
    }
    testing.assertEqual($identities, $al.identities);
    testing.assertEqual($similarities, $al.similarities);
}

# ------------------------------------------------------------- local alignment

func testLocalFindsTheSharedStretch() {
    # a common core between two otherwise unrelated flanks
    def core as string init "ATGCATGCATGCATGC";
    def one as string init "TTTTTTTT" + $core + "TTTTTTTT";
    def other as string init "AAAAAA" + $core + "AAAAAA";
    def al as Alignment init alignLocal($one, $other, dnaOptions());
    # the local alignment is the core, with no gaps and nothing outside it
    testing.assertEqual($al.a, $core);
    testing.assertEqual($al.b, $core);
    testing.assertEqual($al.gaps, 0);
    testing.assertEqual($al.identities, len($core));
    testing.assertEqual($al.startA, 8);
    testing.assertEqual($al.startB, 6);
}

func testLocalStartOffsetsIndexTheOriginals() {
    # THE COORDINATE CONTRACT: startA and startB locate the alignment in the
    # sequences that were passed in
    def one as string init "TTTTTTTTATGCATGCATGCATGCTTTTTTTT";
    def other as string init "AAAAAAATGCATGCATGCATGCAAAAAA";
    def al as Alignment init alignLocal($one, $other, dnaOptions());
    testing.assertEqual($one[$al.startA..$al.startA + len(ungapped($al.a))], ungapped($al.a));
    testing.assertEqual($other[$al.startB..$al.startB + len(ungapped($al.b))], ungapped($al.b));
}

func testLocalNeverScoresBelowZero() {
    # PROPERTY: Smith-Waterman cannot return a negative score, because it can
    # always choose to align nothing
    def al as Alignment init alignLocal("AAAAAAAA", "TTTTTTTT", dnaOptions());
    testing.assertTrue($al.score >= 0);
}

func testLocalBeatsGlobalOnEmbeddedMatch() {
    # the case that motivates the algorithm: a short sequence inside a long one
    def opts as AlignOptions init dnaOptions();
    def haystack as string init "TTTTTTTTTTTTATGCATGCTTTTTTTTTTTT";
    def needle as string init "ATGCATGC";
    def local as Alignment init alignLocal($haystack, $needle, $opts);
    def global as Alignment init alignGlobal($haystack, $needle, $opts);
    testing.assertTrue($local.score > $global.score);
    testing.assertTrue(percentIdentity($local) > percentIdentity($global));
    testing.assertEqual(percentIdentity($local), 100.0);
}

func testLocalScoreMatchesItsOwnRows() {
    def opts as AlignOptions init dnaOptions();
    def al as Alignment init alignLocal("TTTTATGCATGCTTTT", "AAATGCATGCAAA", $opts);
    testing.assertEqual($al.score, rescore($al, $opts));
}

# ---------------------------------------------------------------- the readouts

func testPercentagesAreConsistent() {
    def al as Alignment init alignGlobal(INS_PROTEIN, "MALWMRLLPLLGLLGLWGPDPAAA", proteinOptions());
    # identities are a subset of similarities, both within the alignment
    testing.assertTrue($al.identities <= $al.similarities);
    testing.assertTrue($al.similarities <= $al.length);
    testing.assertTrue(percentIdentity($al) <= percentSimilarity($al));
    testing.assertTrue(percentIdentity($al) >= 0.0 and percentIdentity($al) <= 100.0);
    testing.assertTrue(close(percentGaps($al), 100.0 * $al.gaps / $al.length, 0.000001));
}

func testSimilarityEqualsIdentityWithoutAMatrix() {
    # with match/mismatch scoring nothing but an identity can score positively
    def al as Alignment init alignGlobal("ATGCATGC", "ATGAATGC", dnaOptions());
    testing.assertEqual($al.identities, $al.similarities);
    testing.assertTrue(close(percentIdentity($al), percentSimilarity($al), 0.000001));
}

func testSimilarityExceedsIdentityWithAMatrix() {
    # L -> I and K -> R are conservative, so they count as similar but not identical
    def al as Alignment init alignGlobal("MLKAW", "MIRAW", proteinOptions());
    testing.assertEqual($al.identities, 3);
    testing.assertEqual($al.similarities, 5);
    testing.assertTrue(percentSimilarity($al) > percentIdentity($al));
}

func testFormatAlignment() {
    def al as Alignment init alignGlobal("ATGCATGCAT", "ATGCATGCAT", dnaOptions());
    def rendered as string init formatAlignment($al, 5);
    # two blocks of five, each three lines plus a blank
    testing.assertTrue(strings.contains($rendered, "ATGCA"));
    testing.assertTrue(strings.contains($rendered, "|||||"));
    # positions count from 1 and advance by the residues consumed
    testing.assertTrue(strings.contains($rendered, "1 "));
    testing.assertTrue(strings.contains($rendered, "6 "));
}

func testReadoutsRejectAnEmptyAlignment() {
    testing.assertThrows("identityOfNothing", "bioinformatics");
    testing.assertThrows("formatZeroWidth", "bioinformatics");
}

func identityOfNothing() {
    percentIdentity(Alignment{
        a: "",
        b: "",
        markup: "",
        score: 0,
        length: 0,
        identities: 0,
        similarities: 0,
        gaps: 0,
        startA: 0,
        startB: 0
    });
}

func formatZeroWidth() {
    formatAlignment(alignGlobal("ATGC", "ATGC", dnaOptions()), 0);
}

func testAlignRejectsEmptySequences() {
    testing.assertThrows("alignEmpty", "bioinformatics");
}

func alignEmpty() {
    alignGlobal("", "ATGC", dnaOptions());
}

# ------------------------------------------------------------------ distances

func testHammingDistance() {
    testing.assertEqual(hammingDistance("ATGC", "ATGC"), 0);
    testing.assertEqual(hammingDistance("ATGC", "ATCC"), 1);
    testing.assertEqual(hammingDistance("AAAA", "TTTT"), 4);
}

func testHammingNeedsEqualLengths() {
    # it has no way to express an insertion, which is why editDistance exists
    testing.assertThrows("hammingRagged", "bioinformatics");
}

func hammingRagged() {
    hammingDistance("ATGC", "ATG");
}

func testEditDistance() {
    testing.assertEqual(editDistance("ATGC", "ATGC"), 0);
    # one substitution
    testing.assertEqual(editDistance("ATGC", "ATCC"), 1);
    # one deletion
    testing.assertEqual(editDistance("ATGC", "AGC"), 1);
    # one insertion
    testing.assertEqual(editDistance("AGC", "ATGC"), 1);
    # the textbook example
    testing.assertEqual(editDistance("kitten", "sitting"), 3);
}

func testEditDistanceIsSymmetric() {
    for (def pair in [
        ["ATGCATGC", "ATGAATGC"],
        ["ATGC", "AGC"],
        ["AAAA", "TTTTTT"],
        [INS_CDS, INS_CDS[0..40]]
    ]) {
        testing.assertEqual(editDistance($pair[0], $pair[1]), editDistance($pair[1], $pair[0]));
    }
}

func testEditDistanceObeysItsBounds() {
    # PROPERTY: at least the length difference, at most the longer length
    for (def pair in [
        ["ATGCATGC", "ATGAATGC"],
        ["ATGC", "AGC"],
        ["AAAA", "TTTTTT"],
        [INS_CDS, INS_CDS[0..40]]
    ]) {
        def d as int init editDistance($pair[0], $pair[1]);
        def shorter as int init len($pair[0]);
        def longer as int init len($pair[1]);
        if ($shorter > $longer) {
            $shorter = len($pair[1]);
            $longer = len($pair[0]);
        }
        testing.assertTrue($d >= $longer - $shorter);
        testing.assertTrue($d <= $longer);
    }
}

func testEditDistanceNeverExceedsHamming() {
    # PROPERTY: on equal-length sequences an edit distance may use indels, so it
    # can only do better than counting substitutions
    for (def pair in [["ATGCATGC", "ATGAATGC"], ["AAAA", "TTTT"], ["ATGCATGCAT", "TGCATGCATG"]]) {
        testing.assertTrue(editDistance($pair[0], $pair[1]) <= hammingDistance($pair[0], $pair[1]));
    }
}

func testDistancesAgreeWithAlignment() {
    # a gapless global alignment of equal-length sequences has exactly
    # hammingDistance mismatched columns
    def one as string init "ATGCATGCATGC";
    def other as string init "ATGAATGCATCC";
    def al as Alignment init alignGlobal($one, $other, dnaOptions());
    testing.assertEqual($al.gaps, 0);
    testing.assertEqual($al.length - $al.identities, hammingDistance($one, $other));
}

# ------------------------------------------------- guards added after the audit

func testAlignmentRefusesAnOversizedMatrix() {
    # The matrix costs about a kilobyte per cell, so two large sequences would
    # exhaust memory rather than merely take a while. The guard raises instead.
    testing.assertThrows("alignTooBig", "bioinformatics");
}

func alignTooBig() {
    def opts as AlignOptions init dnaOptions();
    # 2000 x 2000 is four million cells, well past the default million
    alignGlobal(randomDna(2000), randomDna(2000), $opts);
}

func testAlignmentLimitIsRaisable() {
    # a caller who has the memory and means it can lift the ceiling
    def opts as AlignOptions init dnaOptions();
    testing.assertEqual($opts.maxCells, DEFAULT_MAX_ALIGNMENT_CELLS);
    $opts.maxCells = 0;
    def al as Alignment init alignGlobal(randomDna(60), randomDna(60), $opts);
    testing.assertTrue($al.length >= 60);
}

func testAlignmentLimitCountsCellsNotLength() {
    # a short sequence against a long one is fine; it is the product that binds
    def opts as AlignOptions init dnaOptions();
    def al as Alignment init alignGlobal(randomDna(4), randomDna(2000), $opts);
    testing.assertTrue($al.length >= 2000);
}

func testInconsistentAlignmentRaisesBioseqNotRuntime() {
    # Alignment is exported, so one can arrive built by hand. A record whose
    # length disagrees with its rows used to reach the slicing in
    # formatAlignment and fail with a raw interpreter error, breaking the
    # one-error-shape contract for anyone catching on kind.
    testing.assertThrows("formatInconsistent", "bioinformatics");
    testing.assertThrows("identityInconsistent", "bioinformatics");
}

func brokenAlignment() {
    return Alignment{
        a: "",
        b: "",
        markup: "",
        score: 0,
        length: 5,
        identities: 0,
        similarities: 0,
        gaps: 0,
        startA: 0,
        startB: 0
    };
}

func formatInconsistent() {
    formatAlignment(brokenAlignment(), 3);
}

func identityInconsistent() {
    percentIdentity(brokenAlignment());
}

func testWellFormedAlignmentsStillPass() {
    # the stricter precondition must not reject anything the engine produces
    for (def pair in [
        ["ATGCATGC", "ATGAATGC"],
        ["ATGCATGCTT", "ATGATGCTT"],
        ["AAAA", "TTTT"],
        [INS_PROTEIN, INS_PROTEIN]
    ]) {
        def al as Alignment init alignGlobal($pair[0], $pair[1], dnaOptions());
        testing.assertTrue(percentIdentity($al) >= 0.0);
        testing.assertTrue(len(formatAlignment($al, 10)) > 0);
    }
}
