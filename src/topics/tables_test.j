# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - white-box tests for the tables topic.
# Spliced into bioinformatics_test.j via include - run with:
#   jennifer test src/bioinformatics_test.j
#
# tabledata.j is generated, so what is worth testing is not that 64 lines were
# transcribed correctly but that each alternative genetic code reassigns exactly
# the codons NCBI says it does. The difference lists below are the reviewable
# form of the data; if a generated table stops matching one, this fails.

# The four unambiguous bases, in the NCBI layout order.
def const NCBI_BASES as string init "TCAG";

func allCodons() {
    def out as list of string;
    def b as list of string init strings.chars(NCBI_BASES);
    for (def i in $b) {
        for (def j in $b) {
            for (def k in $b) {
                $out[] = $i + $j + $k;
            }
        }
    }
    return $out;
}

# Is `alt` exactly STANDARD with `overrides` applied - no more, no less?
func standardExcept(alt as map of string to string, overrides as map of string to string) {
    if (len($alt) != len(STANDARD)) {
        return false;
    }
    for (def codon in STANDARD) {
        def want as string init STANDARD[$codon];
        if (maps.has($overrides, $codon)) {
            $want = $overrides[$codon];
        }
        if (not maps.has($alt, $codon)) {
            return false;
        }
        if ($alt[$codon] != $want) {
            return false;
        }
    }
    return true;
}

# Every entry in a difference list must actually differ from STANDARD, so a
# stale or copy-pasted override cannot sit there unnoticed.
func allOverridesAreRealChanges(overrides as map of string to string) {
    for (def codon in $overrides) {
        if (not maps.has(STANDARD, $codon)) {
            return false;
        }
        if (STANDARD[$codon] == $overrides[$codon]) {
            return false;
        }
    }
    return true;
}

# ------------------------------------------------------------ the standard code

func testStandardTableIsComplete() {
    testing.assertEqual(len(STANDARD), 64);
    for (def codon in allCodons()) {
        testing.assertTrue(maps.has(STANDARD, $codon));
        testing.assertEqual(len(STANDARD[$codon]), 1);
    }
}

func testStandardAssignments() {
    # the anchors every biologist checks first
    testing.assertEqual(STANDARD["ATG"], "M");
    testing.assertEqual(STANDARD["TGG"], "W");
    testing.assertEqual(STANDARD["TAA"], "*");
    testing.assertEqual(STANDARD["TAG"], "*");
    testing.assertEqual(STANDARD["TGA"], "*");
    # the sole tryptophan and methionine codons
    testing.assertEqual(countResidue(STANDARD, "W"), 1);
    testing.assertEqual(countResidue(STANDARD, "M"), 1);
    # three stops, six each for leucine, serine and arginine
    testing.assertEqual(countResidue(STANDARD, "*"), 3);
    testing.assertEqual(countResidue(STANDARD, "L"), 6);
    testing.assertEqual(countResidue(STANDARD, "S"), 6);
    testing.assertEqual(countResidue(STANDARD, "R"), 6);
}

func countResidue(table as map of string to string, residue as string) {
    def n as int init 0;
    for (def codon in $table) {
        if ($table[$codon] == $residue) {
            $n = $n + 1;
        }
    }
    return $n;
}

func testEveryTableCoversEveryCodon() {
    for (def id in codonTableIds()) {
        def table as map of string to string init codonTable($id);
        testing.assertEqual(len($table), 64);
        for (def codon in allCodons()) {
            testing.assertTrue(maps.has($table, $codon));
        }
    }
}

# -------------------------------------------- alternative codes, as differences

func testVertebrateMitochondrial() {
    # AGA/AGG stop rather than code for arginine; ATA is methionine; TGA is
    # tryptophan rather than a stop
    def diff as map of string to string init {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"};
    testing.assertTrue(allOverridesAreRealChanges($diff));
    testing.assertTrue(standardExcept(VERTEBRATE_MITOCHONDRIAL, $diff));
    testing.assertTrue(standardExcept(codonTable(2), $diff));
    # four stops become two
    testing.assertEqual(countResidue(VERTEBRATE_MITOCHONDRIAL, "*"), 4);
}

func testYeastMitochondrial() {
    # the whole CTN block reads threonine, not leucine
    def diff as map of string to string init {
        "ATA": "M",
        "CTT": "T",
        "CTC": "T",
        "CTA": "T",
        "CTG": "T",
        "TGA": "W"
    };
    testing.assertTrue(allOverridesAreRealChanges($diff));
    testing.assertTrue(standardExcept(YEAST_MITOCHONDRIAL, $diff));
    testing.assertEqual(countResidue(YEAST_MITOCHONDRIAL, "L"), 2);
}

func testMoldMitochondrial() {
    def diff as map of string to string init {"TGA": "W"};
    testing.assertTrue(allOverridesAreRealChanges($diff));
    testing.assertTrue(standardExcept(MOLD_MITOCHONDRIAL, $diff));
}

func testInvertebrateMitochondrial() {
    def diff as map of string to string init {"AGA": "S", "AGG": "S", "ATA": "M", "TGA": "W"};
    testing.assertTrue(allOverridesAreRealChanges($diff));
    testing.assertTrue(standardExcept(INVERTEBRATE_MITOCHONDRIAL, $diff));
    # serine gains the two arginine codons
    testing.assertEqual(countResidue(INVERTEBRATE_MITOCHONDRIAL, "S"), 8);
    testing.assertEqual(countResidue(INVERTEBRATE_MITOCHONDRIAL, "R"), 4);
}

func testCiliateNuclear() {
    # the famous one: two of the three stops read glutamine
    def diff as map of string to string init {"TAA": "Q", "TAG": "Q"};
    testing.assertTrue(allOverridesAreRealChanges($diff));
    testing.assertTrue(standardExcept(CILIATE_NUCLEAR, $diff));
    testing.assertEqual(countResidue(CILIATE_NUCLEAR, "*"), 1);
}

func testEchinodermMitochondrial() {
    def diff as map of string to string init {"AAA": "N", "AGA": "S", "AGG": "S", "TGA": "W"};
    testing.assertTrue(allOverridesAreRealChanges($diff));
    testing.assertTrue(standardExcept(ECHINODERM_MITOCHONDRIAL, $diff));
}

func testEuplotidNuclear() {
    def diff as map of string to string init {"TGA": "C"};
    testing.assertTrue(allOverridesAreRealChanges($diff));
    testing.assertTrue(standardExcept(EUPLOTID_NUCLEAR, $diff));
}

func testBacterialMatchesStandard() {
    # identical assignments - the table exists for its start codons alone
    def none as map of string to string;
    testing.assertTrue(standardExcept(BACTERIAL, $none));
    testing.assertTrue(standardExcept(codonTable(11), $none));
    # which is exactly what distinguishes it
    testing.assertTrue(len(startCodons(11)) > len(startCodons(1)));
}

# -------------------------------------------------------------- the dispatchers

func testCodonTableIds() {
    def ids as list of int init codonTableIds();
    testing.assertEqual(len($ids), 9);
    testing.assertTrue(lists.contains($ids, 1));
    testing.assertTrue(lists.contains($ids, 11));
    # NCBI numbering has gaps - 7 and 8 are retired
    testing.assertFalse(lists.contains($ids, 7));
    testing.assertFalse(lists.contains($ids, 8));
}

func testCodonTableNames() {
    testing.assertEqual(codonTableName(1), "Standard");
    testing.assertEqual(codonTableName(2), "Vertebrate Mitochondrial");
    # every advertised id names a table
    for (def id in codonTableIds()) {
        testing.assertTrue(len(codonTableName($id)) > 0);
    }
}

func testUnknownTableRaises() {
    testing.assertThrows("codonTableSeven", "bioinformatics");
    testing.assertThrows("codonTableNameSeven", "bioinformatics");
    testing.assertThrows("startCodonsSeven", "bioinformatics");
}

func codonTableSeven() {
    codonTable(7);
}

func codonTableNameSeven() {
    codonTableName(7);
}

func startCodonsSeven() {
    startCodons(7);
}

# --------------------------------------------------------------- start codons

func testStartCodonsAreValidCodons() {
    for (def id in codonTableIds()) {
        def starts as list of string init startCodons($id);
        testing.assertTrue(len($starts) > 0);
        for (def codon in $starts) {
            testing.assertEqual(len($codon), 3);
            # a start codon has to be a codon the table actually assigns
            testing.assertTrue(maps.has(codonTable($id), $codon));
        }
        # ATG initiates under every code the deck carries
        testing.assertTrue(lists.contains($starts, "ATG"));
    }
}

func testStandardStartCodons() {
    def starts as list of string init startCodons(1);
    testing.assertEqual(len($starts), 3);
    testing.assertTrue(lists.contains($starts, "ATG"));
    testing.assertTrue(lists.contains($starts, "TTG"));
    testing.assertTrue(lists.contains($starts, "CTG"));
}

# ---------------------------------------------------------------------- IUPAC

func testIupacExpansions() {
    testing.assertEqual(expandIupac("A"), "A");
    testing.assertEqual(expandIupac("R"), "AG");
    testing.assertEqual(expandIupac("Y"), "CT");
    testing.assertEqual(expandIupac("N"), "ACGT");
    # case-insensitive
    testing.assertEqual(expandIupac("r"), "AG");
    # the RNA-only code resolves too
    testing.assertEqual(expandIupac("U"), "U");
}

func testIupacCodeCount() {
    # four bases plus eleven ambiguity codes
    testing.assertEqual(len(IUPAC_DNA), 15);
    testing.assertEqual(len(IUPAC_RNA), 15);
}

func testIupacExpansionsAreOrderedAndValid() {
    for (def code in IUPAC_DNA) {
        def expansion as string init IUPAC_DNA[$code];
        testing.assertTrue(len($expansion) >= 1);
        testing.assertTrue(len($expansion) <= 4);
        # an expansion only ever names unambiguous bases
        for (def base in strings.chars($expansion)) {
            testing.assertTrue(strings.contains("ACGT", $base));
        }
    }
}

func testIupacClass() {
    testing.assertEqual(iupacClass("R"), "[AG]");
    testing.assertEqual(iupacClass("N"), "[ACGT]");
    testing.assertEqual(iupacClass("A"), "[A]");
}

func testMatchesIupac() {
    testing.assertTrue(matchesIupac("R", "A"));
    testing.assertTrue(matchesIupac("R", "G"));
    testing.assertFalse(matchesIupac("R", "C"));
    testing.assertTrue(matchesIupac("N", "T"));
    # the relation is one-way: N admits A, A does not admit N
    testing.assertFalse(matchesIupac("A", "N"));
    testing.assertTrue(matchesIupac("a", "g") == false);
}

func testIupacRejectsNonCodes() {
    testing.assertThrows("expandBadCode", "bioinformatics");
    testing.assertThrows("matchesBadBase", "bioinformatics");
}

func expandBadCode() {
    expandIupac("Z");
}

func matchesBadBase() {
    matchesIupac("R", "GG");
}

# ----------------------------------------------------------------- complements

func testComplementIsAnInvolution() {
    for (def code in COMPLEMENT_DNA) {
        def back as string init COMPLEMENT_DNA[COMPLEMENT_DNA[$code]];
        testing.assertEqual($back, $code);
    }
    for (def code in COMPLEMENT_RNA) {
        def back as string init COMPLEMENT_RNA[COMPLEMENT_RNA[$code]];
        testing.assertEqual($back, $code);
    }
}

func testComplementPairs() {
    testing.assertEqual(COMPLEMENT_DNA["A"], "T");
    testing.assertEqual(COMPLEMENT_DNA["G"], "C");
    testing.assertEqual(COMPLEMENT_RNA["A"], "U");
    # ambiguity codes complement as sets: R (A or G) -> Y (C or T)
    testing.assertEqual(COMPLEMENT_DNA["R"], "Y");
    testing.assertEqual(COMPLEMENT_DNA["B"], "V");
    # the self-complementary codes
    testing.assertEqual(COMPLEMENT_DNA["S"], "S");
    testing.assertEqual(COMPLEMENT_DNA["W"], "W");
    testing.assertEqual(COMPLEMENT_DNA["N"], "N");
}

func testComplementCoversTheAlphabet() {
    for (def base in strings.chars(DNA_ALPHABET)) {
        testing.assertTrue(maps.has(COMPLEMENT_DNA, $base));
    }
    for (def base in strings.chars(RNA_ALPHABET)) {
        testing.assertTrue(maps.has(COMPLEMENT_RNA, $base));
    }
}

# ---------------------------------------------------------------- amino acids

func testAminoAcidLookupBothSpellings() {
    def k as AminoAcid init aminoAcid("K");
    testing.assertEqual($k.three, "Lys");
    testing.assertEqual($k.name, "Lysine");
    testing.assertEqual($k.charge, 1);
    # the three-letter spelling is the same lookup, in any case
    testing.assertEqual(aminoAcid("Lys").one, "K");
    testing.assertEqual(aminoAcid("lys").one, "K");
    testing.assertEqual(aminoAcid("k").one, "K");
}

func testAminoAcidCount() {
    # the standard twenty, plus selenocysteine and pyrrolysine
    testing.assertEqual(len(aminoAcids()), 22);
}

func testAminoAcidRecordsAreWellFormed() {
    for (def aa in aminoAcids()) {
        testing.assertEqual(len($aa.one), 1);
        testing.assertEqual(len($aa.three), 3);
        testing.assertTrue(len($aa.name) > 0);
        testing.assertTrue($aa.monoMass > 0.0);
        testing.assertTrue($aa.avgMass > 0.0);
        # Kyte-Doolittle runs from -4.5 (arginine) to 4.5 (isoleucine)
        testing.assertTrue($aa.hydropathy >= -4.5 and $aa.hydropathy <= 4.5);
        testing.assertTrue($aa.charge >= -1 and $aa.charge <= 1);
        # a side chain ionises exactly when it carries a charge
        if ($aa.charge == 0) {
            testing.assertEqual($aa.pKa, 0.0);
        } else {
            testing.assertTrue($aa.pKa > 0.0);
        }
        # the record agrees with the key it was found under
        testing.assertEqual(aminoAcid($aa.one).three, $aa.three);
    }
}

func testMonoisotopicIsLighterExceptForSelenocysteine() {
    for (def aa in aminoAcids()) {
        if ($aa.one == "U") {
            # selenium is the documented exception: its most abundant isotope
            # (80Se) is heavier than the element's average mass, so the
            # monoisotopic residue outweighs the average one
            testing.assertTrue($aa.monoMass > $aa.avgMass);
        } else {
            testing.assertTrue($aa.monoMass < $aa.avgMass);
        }
    }
}

func testIsobaricResidues() {
    # leucine and isoleucine are indistinguishable by mass - the reason mass
    # spectrometry reports "Leu/Ile"
    testing.assertEqual(aminoAcid("L").monoMass, aminoAcid("I").monoMass);
    testing.assertEqual(aminoAcid("L").avgMass, aminoAcid("I").avgMass);
    # and glycine is the lightest residue
    for (def aa in aminoAcids()) {
        testing.assertTrue($aa.monoMass >= aminoAcid("G").monoMass);
    }
}

func testKnownMasses() {
    # ExPASy residue masses, to the digits the deck stores
    testing.assertEqual(aminoAcid("G").monoMass, 57.02146);
    testing.assertEqual(aminoAcid("W").monoMass, 186.07931);
    testing.assertEqual(aminoAcid("A").avgMass, 71.0788);
}

func testChargedResidues() {
    # three basic, four acidic among the standard twenty
    def basic as int init 0;
    def acidic as int init 0;
    for (def aa in aminoAcids()) {
        if ($aa.one == "U" or $aa.one == "O") {
            continue;
        }
        if ($aa.charge > 0) {
            $basic = $basic + 1;
        }
        if ($aa.charge < 0) {
            $acidic = $acidic + 1;
        }
    }
    testing.assertEqual($basic, 3);
    testing.assertEqual($acidic, 4);
}

func testHasAminoAcid() {
    testing.assertTrue(hasAminoAcid("K"));
    testing.assertTrue(hasAminoAcid("Lys"));
    # the ambiguity codes are deliberately absent
    testing.assertFalse(hasAminoAcid("X"));
    testing.assertFalse(hasAminoAcid("B"));
    testing.assertFalse(hasAminoAcid("Z"));
    testing.assertFalse(hasAminoAcid("J"));
}

func testUnknownAminoAcidRaises() {
    testing.assertThrows("aminoAcidX", "bioinformatics");
    testing.assertThrows("aminoAcidNonsense", "bioinformatics");
}

func aminoAcidX() {
    aminoAcid("X");
}

func aminoAcidNonsense() {
    aminoAcid("Zzz");
}

# --------------------------------------------------------- terminus constants

func testTerminusPkaValues() {
    # the EMBOSS iep parameter set, named in the docs because pI is set-dependent
    testing.assertEqual(PKA_N_TERMINUS, 8.6);
    testing.assertEqual(PKA_C_TERMINUS, 3.6);
    # the amino terminus is the more basic of the two
    testing.assertTrue(PKA_N_TERMINUS > PKA_C_TERMINUS);
}

func testWaterMasses() {
    testing.assertTrue(WATER_MONO_MASS < WATER_AVG_MASS);
    testing.assertTrue(WATER_AVG_MASS > 18.0 and WATER_AVG_MASS < 18.1);
}
