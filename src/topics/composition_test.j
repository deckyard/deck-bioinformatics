# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - white-box tests for the composition topic.
# Spliced into bioinformatics_test.j via include - run with:
#   jennifer test src/bioinformatics_test.j
#
# This topic produces numbers, so the tests come in three kinds, and it is worth
# being clear which is which:
#
#   1. External anchors - a value published independently of this deck.
#      Bradykinin's molecular weight is the strongest one here: it is a mass
#      spectrometry standard, and reproducing it validates the residue masses
#      and the residues-plus-water convention together.
#   2. Cross-implementation checks - a value computed separately from the same
#      published parameter set. These catch coding errors, not parameter errors,
#      and are labelled as such rather than dressed up as external validation.
#   3. Property checks - monotonicity and invariants that must hold whatever the
#      parameters are. These are the ones that survive a parameter-set change.

# INS_CDS, INS_PROTEIN and close() are declared in bioinformatics_test.j.

# ----------------------------------------------------------------- base counts

func testBaseCount() {
    def counts as map of string to int init baseCount("AATTGGCC");
    testing.assertEqual($counts["A"], 2);
    testing.assertEqual($counts["T"], 2);
    testing.assertEqual($counts["G"], 2);
    testing.assertEqual($counts["C"], 2);
    # a base that does not occur has no entry, rather than a zero
    testing.assertFalse(maps.has($counts, "N"));
}

func testBaseCountTotalsTheLength() {
    def counts as map of string to int init baseCount(INS_CDS);
    def total as int init 0;
    for (def base in $counts) {
        $total = $total + $counts[$base];
    }
    testing.assertEqual($total, len(INS_CDS));
}

func testBaseFrequencySumsToOneHundred() {
    def freq as map of string to float init baseFrequency(INS_CDS);
    def total as float init 0.0;
    for (def base in $freq) {
        $total = $total + $freq[$base];
    }
    testing.assertTrue(close($total, 100.0, 0.000001));
}

func testResidueCount() {
    def counts as map of string to int init residueCount(INS_PROTEIN);
    testing.assertEqual($counts["L"], 8);
    testing.assertEqual($counts["A"], 6);
    testing.assertEqual($counts["W"], 2);
}

# ------------------------------------------------------------------ GC content

func testGcContent() {
    testing.assertEqual(gcContent("GGCC"), 100.0);
    testing.assertEqual(gcContent("AATT"), 0.0);
    testing.assertEqual(gcContent("ATGC"), 50.0);
    testing.assertTrue(close(gcContent("ATGCGC"), 66.6666666, 0.0001));
}

func testGcAndAtAgreeWithBaseCount() {
    def counts as map of string to int init baseCount(INS_CDS);
    def gc as float init 100.0 * ($counts["G"] + $counts["C"]) / len(INS_CDS);
    testing.assertTrue(close(gcContent(INS_CDS), $gc, 0.000001));
    # on an unambiguous sequence the two partition the whole length
    testing.assertTrue(close(gcContent(INS_CDS) + atContent(INS_CDS), 100.0, 0.000001));
}

func testAmbiguityCountsTowardNeither() {
    # N is a base but is neither GC nor AT, so the two no longer sum to 100 -
    # which is the behaviour that keeps them from over-counting
    testing.assertEqual(gcContent("GGNN"), 50.0);
    testing.assertEqual(atContent("GGNN"), 0.0);
}

func testGcContentAcceptsRna() {
    testing.assertEqual(gcContent("AUGC"), 50.0);
    testing.assertEqual(atContent("AUGC"), 50.0);
}

func testGcWindows() {
    # four windows of four, stepping two, over a sixteen-base sequence
    def windows as list of float init gcWindows("GGGGAAAAGGGGAAAA", 4, 4);
    testing.assertEqual(len($windows), 4);
    testing.assertEqual($windows[0], 100.0);
    testing.assertEqual($windows[1], 0.0);
    testing.assertEqual($windows[2], 100.0);
    testing.assertEqual($windows[3], 0.0);
}

func testGcWindowsDropsThePartialTail() {
    # ten bases, windows of four stepping four: two full windows, and the two
    # leftover bases are dropped so every value covers the same span
    testing.assertEqual(len(gcWindows("GGGGAAAAGG", 4, 4)), 2);
    # a step of one slides across every position
    testing.assertEqual(len(gcWindows("GGGGAAAAGG", 4, 1)), 7);
}

func testGcWindowsRejectsBadParameters() {
    testing.assertThrows("windowsTooLong", "bioinformatics");
    testing.assertThrows("windowsZeroSize", "bioinformatics");
    testing.assertThrows("windowsZeroStep", "bioinformatics");
}

func windowsTooLong() {
    gcWindows("ATGC", 100, 1);
}

func windowsZeroSize() {
    gcWindows("ATGC", 0, 1);
}

func windowsZeroStep() {
    gcWindows("ATGC", 2, 0);
}

# ------------------------------------------------------------ molecular weight

func testBradykininMass() {
    # EXTERNAL ANCHOR. Bradykinin is a mass spectrometry calibration standard
    # with a published average mass of about 1060.2 Da. Reproducing it validates
    # the residue masses and the residues-plus-water convention at once.
    testing.assertTrue(close(molecularWeight("RPPGFSPFR", PROTEIN), 1060.2, 0.1));
}

func testProteinMassConvention() {
    # a peptide is the sum of its residues plus one water
    def gg as float init molecularWeight("GG", PROTEIN);
    testing.assertTrue(close($gg, 2.0 * aminoAcid("G").avgMass + WATER_AVG_MASS, 0.000001));
    # so a single residue weighs one residue plus that same water
    testing.assertTrue(close(
        molecularWeight("G", PROTEIN),
        aminoAcid("G").avgMass + WATER_AVG_MASS,
        0.000001));
}

func testProteinMassGrowsWithLength() {
    testing.assertTrue(molecularWeight(INS_PROTEIN, PROTEIN) > molecularWeight(
        INS_PROTEIN[0..12],
        PROTEIN));
}

func testDnaMass() {
    # cross-implementation check against the standard single-stranded formula
    testing.assertTrue(close(molecularWeight("ATGC", DNA), 1173.83, 0.01));
    # and the convention: residues, plus a water, less a metaphosphate
    def sum as float init DNA_RESIDUE_MASS["A"] + DNA_RESIDUE_MASS["T"] + DNA_RESIDUE_MASS["G"] +
        DNA_RESIDUE_MASS["C"];
    testing.assertTrue(close(
        molecularWeight("ATGC", DNA),
        $sum + WATER_AVG_MASS - PHOSPHATE_MASS,
        0.000001));
}

func testRnaIsHeavierThanDna() {
    # ribose carries one more oxygen per residue than deoxyribose
    testing.assertTrue(molecularWeight("AGC", RNA) > molecularWeight("AGC", DNA));
}

func testMassRejectsAmbiguity() {
    # a base that might be any of four has no mass, and averaging the
    # possibilities would look precise and mean nothing
    testing.assertThrows("massOfN", "bioinformatics");
    testing.assertThrows("massOfX", "bioinformatics");
    testing.assertThrows("massBadKind", "bioinformatics");
}

func massOfN() {
    molecularWeight("ATGN", DNA);
}

func massOfX() {
    molecularWeight("MKVX", PROTEIN);
}

func massBadKind() {
    molecularWeight("ATGC", 42);
}

# --------------------------------------------------------------- melting point

func testWallaceRule() {
    # Tm = 2(A+T) + 4(G+C), exactly
    testing.assertEqual(meltingTemp("AAAA"), 8.0);
    testing.assertEqual(meltingTemp("GGGG"), 16.0);
    testing.assertEqual(meltingTemp("ATGCATGCATGCATGC"), 48.0);
}

func testWallaceIsBlindToOrder() {
    # the rule counts bases and nothing else, which is precisely its weakness
    testing.assertEqual(meltingTemp("GCGCGC"), meltingTemp("GGGCCC"));
}

func testNearestNeighbourSeesOrder() {
    # CROSS-IMPLEMENTATION CHECK against the SantaLucia (1998) formula computed
    # independently. Same composition, twelve degrees apart - the entire reason
    # the nearest-neighbour model exists.
    testing.assertTrue(close(meltingTempNn("GCGCGC"), 22.767, 0.05));
    testing.assertTrue(close(meltingTempNn("GGGCCC"), 10.387, 0.05));
    testing.assertTrue(meltingTempNn("GCGCGC") > meltingTempNn("GGGCCC") + 10.0);
}

func testNearestNeighbourOnATypicalPrimer() {
    # a 50% GC 20-mer at standard PCR conditions lands in the usual primer
    # window, which is the sanity check a bench scientist would apply
    def tm as float init meltingTempNn("ATGCATGCATGCATGCATGC");
    testing.assertTrue(close($tm, 57.683, 0.05));
    testing.assertTrue($tm > 50.0 and $tm < 65.0);
}

func testMeltingTempRisesWithGc() {
    # PROPERTY: three hydrogen bonds beat two, whatever the parameter set
    testing.assertTrue(meltingTempNn("GCGCGCGCGC") > meltingTempNn("ATATATATAT"));
    testing.assertTrue(meltingTemp("GCGCGCGCGC") > meltingTemp("ATATATATAT"));
}

func testMeltingTempRisesWithLength() {
    testing.assertTrue(meltingTempNn("ATGCATGCATGCATGCATGC") > meltingTempNn("ATGCATGC"));
}

func testMeltingTempRisesWithSaltAndConcentration() {
    def base as float init meltingTempNnWith("ATGCATGCATGCATGCATGC", 250.0, 50.0);
    # salt screens the backbone's charge, stabilising the duplex
    testing.assertTrue(meltingTempNnWith("ATGCATGCATGCATGCATGC", 250.0, 100.0) > $base);
    testing.assertTrue(close(meltingTempNnWith("ATGCATGCATGCATGCATGC", 250.0, 100.0), 61.051, 0.05));
    # more strands means more annealing
    testing.assertTrue(meltingTempNnWith("ATGCATGCATGCATGCATGC", 1000.0, 50.0) > $base);
}

func testSelfComplementaryIsHandled() {
    # a palindrome anneals to a copy of itself, so it takes the symmetry
    # correction and the undivided strand concentration
    testing.assertEqual(reverseComplement("GAATTC"), "GAATTC");
    testing.assertTrue(close(meltingTempNn("GAATTC"), -15.116, 0.05));
}

func testMeltingTempRejectsBadInput() {
    testing.assertThrows("tmAmbiguous", "bioinformatics");
    testing.assertThrows("tmTooShort", "bioinformatics");
    testing.assertThrows("tmZeroSalt", "bioinformatics");
}

func tmAmbiguous() {
    meltingTempNn("ATGN");
}

func tmTooShort() {
    meltingTempNn("A");
}

func tmZeroSalt() {
    meltingTempNnWith("ATGCATGC", 250.0, 0.0);
}

# ------------------------------------------------------------ protein chemistry

func testNetChargeDirection() {
    # PROPERTY: charge falls monotonically as pH rises - the fact the pI
    # bisection depends on
    def previous as float init netCharge(INS_PROTEIN, 0.0);
    for (def i as int init 1; $i <= 28; $i = $i + 1) {
        def current as float init netCharge(INS_PROTEIN, 0.5 * $i);
        testing.assertTrue($current < $previous);
        $previous = $current;
    }
}

func testNetChargeOfPolyLysine() {
    # four lysines and an amino terminus are all protonated at neutral pH
    testing.assertTrue(close(netCharge("KKKK", 7.0), 3.9753, 0.001));
    # four glutamates are all deprotonated
    testing.assertTrue(close(netCharge("EEEE", 7.0), -4.0191, 0.001));
}

func testIsoelectricPoint() {
    # CROSS-IMPLEMENTATION CHECK against the same EMBOSS iep pKa set computed
    # independently. Compare against `iep`, not ProtParam: the sets differ.
    testing.assertTrue(close(isoelectricPoint("KKKK"), 11.2783, 0.001));
    testing.assertTrue(close(isoelectricPoint("EEEE"), 3.3709, 0.001));
    testing.assertTrue(close(isoelectricPoint(INS_PROTEIN), 6.3369, 0.001));
}

func testIsoelectricPointIsWhereChargeIsZero() {
    # PROPERTY: the definition itself, which holds for any pKa set
    for (def p in ["KKKK", "EEEE", INS_PROTEIN, "MKWVTFISLLLLFSSAYS"]) {
        testing.assertTrue(math.abs(netCharge($p, isoelectricPoint($p))) < 0.001);
    }
}

func testBasicProteinsHaveHighPi() {
    # PROPERTY: more basic residues, higher pI
    testing.assertTrue(isoelectricPoint("KKKK") > isoelectricPoint("KKKE"));
    testing.assertTrue(isoelectricPoint("KKKE") > isoelectricPoint("EEEE"));
}

func testGravy() {
    # CROSS-IMPLEMENTATION CHECK on the Kyte-Doolittle sum
    testing.assertTrue(close(gravy(INS_PROTEIN), 1.25, 0.000001));
    # a run of one residue is that residue's index exactly
    testing.assertTrue(close(gravy("IIII"), 4.5, 0.000001));
    testing.assertTrue(close(gravy("RRRR"), -4.5, 0.000001));
}

func testGravySignSeparatesHydrophobic() {
    # PROPERTY: the insulin signal peptide is strikingly hydrophobic, which is
    # what a signal peptide is for
    testing.assertTrue(gravy(INS_PROTEIN) > 0.0);
    testing.assertTrue(gravy("DEKRNQ") < 0.0);
}

func testAromaticity() {
    testing.assertTrue(close(aromaticity("FWYA"), 0.75, 0.000001));
    testing.assertTrue(close(aromaticity("AAAA"), 0.0, 0.000001));
}

func testProteinChemistryRejectsAmbiguity() {
    testing.assertThrows("gravyOfX", "bioinformatics");
    testing.assertThrows("chargeOfX", "bioinformatics");
}

func gravyOfX() {
    gravy("MKVX");
}

func chargeOfX() {
    netCharge("MKVX", 7.0);
}

# -------------------------------------------------------------------- summaries

func testDnaStats() {
    def stats as DnaStats init dnaStats(INS_CDS);
    testing.assertEqual($stats.length, 72);
    testing.assertFalse($stats.ambiguous);
    # the summary agrees with the individual verbs
    testing.assertTrue(close($stats.gc, gcContent(INS_CDS), 0.000001));
    testing.assertTrue(close($stats.at, atContent(INS_CDS), 0.000001));
    testing.assertTrue(close($stats.molecularWeight, molecularWeight(INS_CDS, DNA), 0.000001));
    testing.assertTrue(close($stats.meltingTemp, meltingTemp(INS_CDS), 0.000001));
    testing.assertEqual($stats.counts["G"], baseCount(INS_CDS)["G"]);
}

func testDnaStatsOnAnAmbiguousSequence() {
    # composition still works; the numbers that have no defined value are
    # reported as 0.0 with the flag set, rather than raising or being invented
    def stats as DnaStats init dnaStats("ATGCNNNN");
    testing.assertTrue($stats.ambiguous);
    testing.assertEqual($stats.length, 8);
    testing.assertEqual($stats.molecularWeight, 0.0);
    testing.assertEqual($stats.meltingTemp, 0.0);
    testing.assertEqual($stats.gc, 25.0);
}

func testProteinStats() {
    def stats as ProteinStats init proteinStats(INS_PROTEIN);
    testing.assertEqual($stats.length, 24);
    testing.assertFalse($stats.ambiguous);
    testing.assertTrue(close(
        $stats.molecularWeight,
        molecularWeight(INS_PROTEIN, PROTEIN),
        0.000001));
    testing.assertTrue(close($stats.isoelectricPoint, isoelectricPoint(INS_PROTEIN), 0.000001));
    testing.assertTrue(close($stats.gravy, gravy(INS_PROTEIN), 0.000001));
}

func testProteinStatsOnASequenceWithAStop() {
    # a stop has no mass or pKa, so the same flag-and-zero rule applies
    def stats as ProteinStats init proteinStats("MKV*");
    testing.assertTrue($stats.ambiguous);
    testing.assertEqual($stats.length, 4);
    testing.assertEqual($stats.molecularWeight, 0.0);
}

# ----------------------------------------------------------------- codon usage

func testCodonUsage() {
    # six codons: ATG once, ACG once, TTA once, and GGA three times
    def usage as list of CodonUse init codonUsage("ATGACGTTAGGAGGAGGA");
    testing.assertEqual(len($usage), 4);
    def total as int init 0;
    for (def u in $usage) {
        $total = $total + $u.count;
    }
    testing.assertEqual($total, 6);
}

func testCodonUsageFractionIsShareOfSynonyms() {
    # GGA and GGG both encode glycine, so each takes half of the glycine share
    def usage as list of CodonUse init codonUsage("GGAGGG");
    testing.assertEqual(len($usage), 2);
    for (def u in $usage) {
        testing.assertEqual($u.aminoAcid, "G");
        testing.assertTrue(close($u.fraction, 0.5, 0.000001));
        testing.assertTrue(close($u.perThousand, 500.0, 0.000001));
    }
}

func testCodonUsageFractionsSumToOnePerResidue() {
    # PROPERTY: each residue's synonymous fractions partition 1.0
    def totals as map of string to float;
    for (def u in codonUsage(INS_CDS)) {
        if (maps.has($totals, $u.aminoAcid)) {
            $totals[$u.aminoAcid] = $totals[$u.aminoAcid] + $u.fraction;
        } else {
            $totals[$u.aminoAcid] = $u.fraction;
        }
    }
    for (def residue in $totals) {
        testing.assertTrue(close($totals[$residue], 1.0, 0.000001));
    }
}

func testCodonUsagePerThousandSumsToOneThousand() {
    def total as float init 0.0;
    for (def u in codonUsage(INS_CDS)) {
        $total = $total + $u.perThousand;
    }
    testing.assertTrue(close($total, 1000.0, 0.000001));
}

func testCodonUsageAgreesWithTranslation() {
    # every codon's residue is the one translate would emit
    for (def u in codonUsage(INS_CDS)) {
        testing.assertEqual($u.aminoAcid, translate($u.codon));
    }
}

func testCodonUsageAcceptsRna() {
    # U is folded to T so one table serves both spellings
    def fromRna as list of CodonUse init codonUsage("AUGACG");
    testing.assertEqual($fromRna[0].codon, "ATG");
    testing.assertEqual($fromRna[0].aminoAcid, "M");
}

func testCodonUsageUnderAnAlternativeCode() {
    def usage as list of CodonUse init codonUsageWith("TGATGA", VERTEBRATE_MITOCHONDRIAL);
    testing.assertEqual($usage[0].aminoAcid, "W");
    testing.assertEqual(codonUsage("TGATGA")[0].aminoAcid, "*");
}

func testCodonUsageRejectsShortSequence() {
    testing.assertThrows("usageTooShort", "bioinformatics");
}

func usageTooShort() {
    codonUsage("AT");
}
