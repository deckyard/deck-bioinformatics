# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - white-box tests for the enzymes topic.
# Spliced into bioinformatics_test.j via include - run with:
#   jennifer test src/bioinformatics_test.j
#
# The table is generated, so these test the invariants a hand-added row could
# violate - a site outside the IUPAC alphabet, a cut offset past the end of the
# site, an overhang label that disagrees with the offsets it is derived from -
# plus the published geometry of the enzymes everyone actually recognises.

# Reverse-complement over the IUPAC alphabet, written out here rather than
# borrowed from transforms.j so this test does not depend on a later topic.
func rcIupac(site as string) {
    def out as string init "";
    def chars as list of string init strings.chars($site);
    for (def i as int init len($chars) - 1; $i >= 0; $i = $i - 1) {
        $out = $out + COMPLEMENT_DNA[$chars[$i]];
    }
    return $out;
}

func testTableIsPopulated() {
    testing.assertTrue(len(enzymes()) >= 70);
    testing.assertEqual(len(enzymes()), len(enzymeNames()));
}

func testEveryRowIsWellFormed() {
    for (def e in enzymes()) {
        testing.assertTrue(len($e.name) > 0);
        testing.assertTrue(len($e.site) >= 4);

        # a site is spelled in the DNA IUPAC alphabet and nothing else
        for (def base in strings.chars($e.site)) {
            testing.assertTrue(strings.contains(DNA_ALPHABET, $base));
        }

        # both cuts land within the recognition site
        testing.assertTrue($e.cutTop >= 0 and $e.cutTop <= len($e.site));
        testing.assertTrue($e.cutBottom >= 0 and $e.cutBottom <= len($e.site));

        # the overhang label is derived from the offsets, so it must agree
        if ($e.cutTop == $e.cutBottom) {
            testing.assertEqual($e.overhang, "blunt");
            testing.assertEqual($e.overhangLength, 0);
            testing.assertTrue(isBlunt($e));
        } elseif ($e.cutTop < $e.cutBottom) {
            testing.assertEqual($e.overhang, "5'");
            testing.assertEqual($e.overhangLength, $e.cutBottom - $e.cutTop);
            testing.assertFalse(isBlunt($e));
        } else {
            testing.assertEqual($e.overhang, "3'");
            testing.assertEqual($e.overhangLength, $e.cutTop - $e.cutBottom);
            testing.assertFalse(isBlunt($e));
        }
    }
}

func testLookupRoundTrips() {
    for (def name in enzymeNames()) {
        testing.assertTrue(hasEnzyme($name));
        testing.assertEqual(enzyme($name).name, $name);
    }
}

func testLookupIsCaseInsensitive() {
    testing.assertEqual(enzyme("ecori").name, "EcoRI");
    testing.assertEqual(enzyme("ECORI").name, "EcoRI");
    testing.assertEqual(enzyme("  EcoRI  ").name, "EcoRI");
    testing.assertTrue(hasEnzyme("hindiii"));
}

func testKnownGeometry() {
    # G^AATTC - the canonical 5' overhang
    def eco as Enzyme init enzyme("EcoRI");
    testing.assertEqual($eco.site, "GAATTC");
    testing.assertEqual($eco.cutTop, 1);
    testing.assertEqual($eco.cutBottom, 5);
    testing.assertEqual($eco.overhang, "5'");
    testing.assertEqual($eco.overhangLength, 4);

    # CTGCA^G - a 3' overhang, the mirror case
    def pst as Enzyme init enzyme("PstI");
    testing.assertEqual($pst.cutTop, 5);
    testing.assertEqual($pst.cutBottom, 1);
    testing.assertEqual($pst.overhang, "3'");
    testing.assertEqual($pst.overhangLength, 4);

    # CCC^GGG - blunt
    def sma as Enzyme init enzyme("SmaI");
    testing.assertEqual($sma.cutTop, 3);
    testing.assertEqual($sma.cutBottom, 3);
    testing.assertTrue(isBlunt($sma));

    # GC^GGCCGC - the eight-cutter used to find rare sites
    testing.assertEqual(enzyme("NotI").site, "GCGGCCGC");
    testing.assertEqual(len(enzyme("NotI").site), 8);
}

func testDegenerateSitesAreSpelledInIupac() {
    # EcoRII reads CCWGG - W is A or T, and the site is stored that way rather
    # than expanded into two rows
    testing.assertEqual(enzyme("EcoRII").site, "CCWGG");
    testing.assertTrue(strings.contains(enzyme("HinfI").site, "N"));
    testing.assertTrue(strings.contains(enzyme("AflIII").site, "R"));
    testing.assertTrue(strings.contains(enzyme("XmnI").site, "NNNN"));
}

func testEverySiteIsPalindromic() {
    # A type II enzyme binds as a homodimer, so its recognition site reads the
    # same on both strands. This holds for the degenerate sites too - N
    # complements to N and W to W, so the specificity in the outer bases pairs
    # up - which is why there is no exception list here. A site that fails this
    # is a typo.
    for (def e in enzymes()) {
        testing.assertEqual(rcIupac($e.site), $e.site);
    }
}

func testIsoschizomers() {
    # same site, same cut - true isoschizomers
    testing.assertTrue(lists.contains(isoschizomers("HpaII"), "MspI"));
    testing.assertTrue(lists.contains(isoschizomers("MspI"), "HpaII"));
    testing.assertEqual(enzyme("HpaII").cutTop, enzyme("MspI").cutTop);

    # same site, different cut - neoschizomers, which is why isoschizomers()
    # returns them and the caller reads cutTop to tell them apart
    testing.assertTrue(lists.contains(isoschizomers("SmaI"), "XmaI"));
    testing.assertTrue(enzyme("SmaI").cutTop != enzyme("XmaI").cutTop);
    testing.assertEqual(enzyme("SmaI").overhang, "blunt");
    testing.assertEqual(enzyme("XmaI").overhang, "5'");

    # Sau3AI and MboI are the classic pair for partial digests
    testing.assertTrue(lists.contains(isoschizomers("Sau3AI"), "MboI"));

    # an enzyme is never its own isoschizomer
    for (def name in enzymeNames()) {
        testing.assertFalse(lists.contains(isoschizomers($name), $name));
    }
}

func testIsoschizomersAreSymmetric() {
    for (def e in enzymes()) {
        for (def other in isoschizomers($e.name)) {
            testing.assertTrue(lists.contains(isoschizomers($other), $e.name));
        }
    }
}

func testSiteLengthsAreSensible() {
    # four-, six- and eight-cutters; nothing shorter than four
    for (def e in enzymes()) {
        testing.assertTrue(len($e.site) >= 4 and len($e.site) <= 10);
    }
    testing.assertEqual(len(enzyme("AluI").site), 4);
    testing.assertEqual(len(enzyme("BamHI").site), 6);
    testing.assertEqual(len(enzyme("AscI").site), 8);
}

func testUnknownEnzymeRaises() {
    testing.assertThrows("enzymeNonsense", "bioinformatics");
    testing.assertThrows("isoschizomersNonsense", "bioinformatics");
    testing.assertFalse(hasEnzyme("BsaI"));
}

func enzymeNonsense() {
    enzyme("NotAnEnzyme");
}

func isoschizomersNonsense() {
    isoschizomers("NotAnEnzyme");
}
