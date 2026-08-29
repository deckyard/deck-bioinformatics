# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - white-box tests for the manipulate topic.
# Spliced into bioinformatics_test.j via include - run with:
#   jennifer test src/bioinformatics_test.j
#
# Random output cannot be compared against a remembered string without pinning
# the interpreter's own generator, which is not this deck's to pin. So these
# test three things instead: the invariants that hold for every draw (length,
# alphabet, composition), the reproducibility that `math.randSeed` is supposed
# to give, and the statistical shape over enough draws that a broken
# distribution cannot pass. Tolerances are loose enough not to flake and tight
# enough to catch a real fault.

func countOf(s as string, residue as string) {
    def n as int init 0;
    for (def c in strings.chars($s)) {
        if ($c == $residue) {
            $n = $n + 1;
        }
    }
    return $n;
}

func sameComposition(a as string, b as string) {
    def left as map of string to int init residueCount($a);
    def right as map of string to int init residueCount($b);
    if (len($left) != len($right)) {
        return false;
    }
    for (def residue in $left) {
        if (not maps.has($right, $residue)) {
            return false;
        }
        if ($left[$residue] != $right[$residue]) {
            return false;
        }
    }
    return true;
}

# ------------------------------------------------------------------- generation

func testRandomDnaShape() {
    def s as string init randomDna(200);
    testing.assertEqual(len($s), 200);
    testing.assertTrue(isDna($s));
    # only the four unambiguous bases - a generator emitting N would be
    # producing missing data, not sequence
    for (def base in strings.chars($s)) {
        testing.assertTrue(strings.contains(DNA_BASES, $base));
    }
}

func testRandomRnaAndProteinShape() {
    def rna as string init randomRna(100);
    testing.assertEqual(len($rna), 100);
    testing.assertTrue(isRna($rna));
    testing.assertFalse(strings.contains($rna, "T"));

    def protein as string init randomProtein(100);
    testing.assertEqual(len($protein), 100);
    testing.assertTrue(isProtein($protein));
    for (def residue in strings.chars($protein)) {
        testing.assertTrue(hasAminoAcid($residue));
    }
}

func testSeedingIsReproducible() {
    # the property a published analysis depends on
    math.randSeed(20260828);
    def first as string init randomDna(500);
    math.randSeed(20260828);
    def second as string init randomDna(500);
    testing.assertEqual($first, $second);

    # and a different seed gives different output, so seeding is doing something
    math.randSeed(1);
    testing.assertTrue(randomDna(500) != $first);
}

func testRandomDnaIsRoughlyUniform() {
    # over 4000 bases each base should land near 25%. The bound is wide enough
    # that a fair generator will not flake and narrow enough that one stuck on a
    # subset of the alphabet cannot pass.
    math.randSeed(7);
    def s as string init randomDna(4000);
    for (def base in strings.chars(DNA_BASES)) {
        def share as float init 100.0 * countOf($s, $base) / 4000.0;
        testing.assertTrue($share > 20.0 and $share < 30.0);
    }
}

func testWeightedComposition() {
    # 70% A and 30% T, and nothing else at all
    math.randSeed(11);
    def s as string init randomDnaWith(2000, {"A": 70.0, "T": 30.0});
    testing.assertEqual(len($s), 2000);
    testing.assertEqual(countOf($s, "G"), 0);
    testing.assertEqual(countOf($s, "C"), 0);
    def shareA as float init 100.0 * countOf($s, "A") / 2000.0;
    testing.assertTrue($shareA > 65.0 and $shareA < 75.0);
}

func testWeightsAreNormalised() {
    # weights are relative, so these two describe the same distribution
    math.randSeed(3);
    def asCounts as string init randomDnaWith(2000, {"A": 3.0, "T": 1.0});
    math.randSeed(3);
    def asPercent as string init randomDnaWith(2000, {"A": 75.0, "T": 25.0});
    testing.assertEqual($asCounts, $asPercent);
}

func testWeightedCanExcludeResidues() {
    # a base left out never appears - the way to get a purine-only sequence
    math.randSeed(5);
    def purines as string init randomDnaWith(300, {"A": 1.0, "G": 1.0});
    testing.assertEqual(countOf($purines, "T"), 0);
    testing.assertEqual(countOf($purines, "C"), 0);
    testing.assertEqual(countOf($purines, "A") + countOf($purines, "G"), 300);
}

func testGenerationRejectsBadInput() {
    testing.assertThrows("randomZeroLength", "bioinformatics");
    testing.assertThrows("randomNegativeLength", "bioinformatics");
    testing.assertThrows("weightsWrongAlphabet", "bioinformatics");
    testing.assertThrows("weightsNegative", "bioinformatics");
    testing.assertThrows("weightsAllZero", "bioinformatics");
}

func randomZeroLength() {
    randomDna(0);
}

func randomNegativeLength() {
    randomProtein(-5);
}

func weightsWrongAlphabet() {
    randomDnaWith(10, {"A": 1.0, "X": 1.0});
}

func weightsNegative() {
    randomDnaWith(10, {"A": 1.0, "T": -1.0});
}

func weightsAllZero() {
    randomDnaWith(10, {"A": 0.0, "T": 0.0});
}

# --------------------------------------------------------------------- shuffle

func testShufflePreservesComposition() {
    # the defining property, and the reason a shuffle is the honest null model
    math.randSeed(99);
    for (def i as int init 0; $i < 20; $i = $i + 1) {
        def shuffled as string init shuffle(INS_CDS);
        testing.assertEqual(len($shuffled), len(INS_CDS));
        testing.assertTrue(sameComposition($shuffled, INS_CDS));
        # GC content follows from composition, so it survives too
        testing.assertTrue(close(gcContent($shuffled), gcContent(INS_CDS), 0.000001));
    }
}

func testShuffleActuallyReorders() {
    # over twenty shuffles of a 72-base sequence, at least one must differ from
    # the original - the probability of the identity permutation twenty times
    # over is far beyond vanishing
    math.randSeed(4);
    def moved as bool init false;
    for (def i as int init 0; $i < 20; $i = $i + 1) {
        if (shuffle(INS_CDS) != INS_CDS) {
            $moved = true;
        }
    }
    testing.assertTrue($moved);
}

func testShuffleWorksOnProtein() {
    math.randSeed(8);
    def shuffled as string init shuffle(INS_PROTEIN);
    testing.assertTrue(sameComposition($shuffled, INS_PROTEIN));
    # composition-derived measures are unchanged; order-derived ones may not be
    testing.assertTrue(close(
        molecularWeight($shuffled, PROTEIN),
        molecularWeight(INS_PROTEIN, PROTEIN),
        0.000001));
    testing.assertTrue(close(gravy($shuffled), gravy(INS_PROTEIN), 0.000001));
}

func testShuffleOfOneResidue() {
    testing.assertEqual(shuffle("A"), "A");
}

func testShuffleRejectsEmpty() {
    testing.assertThrows("shuffleEmpty", "bioinformatics");
}

func shuffleEmpty() {
    shuffle("");
}

# ---------------------------------------------------------------------- mutate

func testMutateKeepsLengthAndAlphabet() {
    math.randSeed(13);
    def mutated as string init mutate(INS_CDS, 0.1, DNA);
    testing.assertEqual(len($mutated), len(INS_CDS));
    testing.assertTrue(isDna($mutated));
}

func testZeroRateChangesNothing() {
    math.randSeed(2);
    testing.assertEqual(mutate(INS_CDS, 0.0, DNA), INS_CDS);
}

func testFullRateChangesEverything() {
    # at rate 1.0 every position is substituted, and a substitution always
    # changes - so no position may still hold its original residue
    math.randSeed(6);
    def mutated as string init mutate(INS_CDS, 1.0, DNA);
    testing.assertEqual(len($mutated), len(INS_CDS));
    for (def i as int init 0; $i < len(INS_CDS); $i = $i + 1) {
        testing.assertTrue($mutated[$i..$i + 1] != INS_CDS[$i..$i + 1]);
    }
}

func testMutationRateIsApproximatelyHonoured() {
    # a thousand bases at 10% should take roughly a hundred substitutions
    math.randSeed(17);
    def original as string init randomDna(1000);
    def mutated as string init mutate($original, 0.1, DNA);
    def changed as int init 0;
    for (def i as int init 0; $i < len($original); $i = $i + 1) {
        if ($mutated[$i..$i + 1] != $original[$i..$i + 1]) {
            $changed = $changed + 1;
        }
    }
    testing.assertTrue($changed > 60 and $changed < 145);
}

func testMutateProteinUsesTheProteinAlphabet() {
    math.randSeed(21);
    def mutated as string init mutate(INS_PROTEIN, 1.0, PROTEIN);
    for (def residue in strings.chars($mutated)) {
        testing.assertTrue(strings.contains(STANDARD_RESIDUES, $residue));
    }
}

func testMutateRejectsBadInput() {
    testing.assertThrows("mutateRateTooHigh", "bioinformatics");
    testing.assertThrows("mutateRateNegative", "bioinformatics");
    testing.assertThrows("mutateBadKind", "bioinformatics");
    testing.assertThrows("mutateWrongAlphabet", "bioinformatics");
}

func mutateRateTooHigh() {
    mutate("ATGC", 1.5, DNA);
}

func mutateRateNegative() {
    mutate("ATGC", -0.1, DNA);
}

func mutateBadKind() {
    mutate("ATGC", 0.1, 42);
}

func mutateWrongAlphabet() {
    mutate("AUGC", 0.1, DNA);
}

# --------------------------------------------------------------------- subrange

func testSubrangeIsOneBasedInclusive() {
    # the deck's one 1-based verb, so a paper's "residues 3 to 5" transcribes
    # directly
    testing.assertEqual(subrange("ATGCATGC", 3, 5), "GCA");
    testing.assertEqual(subrange("ATGCATGC", 1, 1), "A");
    testing.assertEqual(subrange("ATGCATGC", 1, 8), "ATGCATGC");
    testing.assertEqual(subrange("ATGCATGC", 8, 8), "C");
}

func testSubrangeAgreesWithNativeSlicing() {
    # the translation between the two conventions, stated as a test so the
    # off-by-one cannot drift
    def s as string init INS_CDS;
    for (def from as int init 1; $from <= 10; $from = $from + 1) {
        testing.assertEqual(subrange($s, $from, $from + 5), $s[$from - 1..$from + 5]);
    }
}

func testSubrangeRejectsBadRanges() {
    testing.assertThrows("subrangeZero", "bioinformatics");
    testing.assertThrows("subrangeBackwards", "bioinformatics");
    testing.assertThrows("subrangePastEnd", "bioinformatics");
}

func subrangeZero() {
    subrange("ATGC", 0, 2);
}

func subrangeBackwards() {
    subrange("ATGC", 3, 2);
}

func subrangePastEnd() {
    subrange("ATGC", 2, 99);
}

# ---------------------------------------------------------------------- windows

func testWindows() {
    testing.assertEqual(len(windows("ATGCATGC", 4, 4)), 2);
    testing.assertEqual(windows("ATGCATGC", 4, 4)[0], "ATGC");
    testing.assertEqual(windows("ATGCATGC", 4, 4)[1], "ATGC");
    # overlapping windows, stepping one at a time
    testing.assertEqual(len(windows("ATGCATGC", 4, 1)), 5);
}

func testWindowsDropThePartialTail() {
    # every window is exactly `size` long, so the series stays comparable
    for (def w in windows(INS_CDS, 10, 7)) {
        testing.assertEqual(len($w), 10);
    }
    # ten residues, windows of four stepping four: two full windows, tail dropped
    testing.assertEqual(len(windows("ATGCATGCAT", 4, 4)), 2);
}

func testWindowsAgreeWithGcWindows() {
    # the two window helpers must cut the sequence the same way
    def sliced as list of string init windows(INS_CDS, 18, 9);
    def values as list of float init gcWindows(INS_CDS, 18, 9);
    testing.assertEqual(len($sliced), len($values));
    for (def i as int init 0; $i < len($sliced); $i = $i + 1) {
        testing.assertTrue(close(gcContent($sliced[$i]), $values[$i], 0.000001));
    }
}

func testWindowsRejectBadParameters() {
    testing.assertThrows("windowsTooBig", "bioinformatics");
    testing.assertThrows("windowsBadSize", "bioinformatics");
    testing.assertThrows("windowsBadStep", "bioinformatics");
}

func windowsTooBig() {
    windows("ATGC", 99, 1);
}

func windowsBadSize() {
    windows("ATGC", 0, 1);
}

func windowsBadStep() {
    windows("ATGC", 2, 0);
}

# ------------------------------------------------- the null model, end to end

func testShuffleDestroysMotifsButNotComposition() {
    # The point of the topic. A restriction site in a real sequence is either
    # meaningful or a coincidence; comparing against shuffles of the same
    # sequence is how that question gets answered, because a shuffle differs
    # only in order.
    def construct as string init POLYLINKER;
    testing.assertEqual(len(restrictionSites($construct, ["EcoRI"])), 1);

    math.randSeed(23);
    def kept as int init 0;
    for (def i as int init 0; $i < 30; $i = $i + 1) {
        def shuffled as string init shuffle($construct);
        testing.assertTrue(sameComposition($shuffled, $construct));
        $kept = $kept + len(restrictionSites($shuffled, ["EcoRI"]));
    }
    # a six-base site survives shuffling only by chance, so across thirty
    # shuffles of a 124-base sequence it should be rare
    testing.assertTrue($kept < 10);
}
