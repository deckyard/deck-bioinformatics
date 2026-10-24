# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - white-box tests for the search topic.
# Spliced into bioinformatics_test.j via include - run with:
#   jennifer test src/bioinformatics_test.j
#
# Two things here are worth more than the rest: that overlapping matches are
# reported, and that a feature found on the reverse strand comes back in the
# coordinates of the sequence the caller passed in. Both are silent failures
# otherwise - a plausible answer that is simply wrong - so both are checked
# against the sequence itself rather than against a remembered number.

# The synthetic polylinker from testdata/sequences.fasta, carrying one site each
# for EcoRI, KpnI, BamHI, HindIII, PstI and ApaI, between two fillers.
def const POLYLINKER as string init 'TTAGCTAACGCTAGCTAAGGCATTACGATCGGATTACCGTAACGGAATTCGGTACCGGATCCAAGCTTCTGCAGGGGCCCCCATTAGGCATAGCATTAGCGGATCAATTGCCATAGGCATTACC';

# ------------------------------------------------------------------ exact find

func testFindReportsOverlaps() {
    # three overlapping occurrences, not one
    def hits as list of Match init find("AAAA", "AA");
    testing.assertEqual(len($hits), 3);
    testing.assertEqual($hits[0].start, 0);
    testing.assertEqual($hits[1].start, 1);
    testing.assertEqual($hits[2].start, 2);
}

func testFindCoordinatesSliceBack() {
    # the contract: $s[$m.start..$m.end] is the match
    for (def m in find(POLYLINKER, "GGATCC")) {
        testing.assertEqual(POLYLINKER[$m.start..$m.end], $m.text);
        testing.assertEqual($m.text, "GGATCC");
        testing.assertEqual($m.end - $m.start, 6);
        testing.assertEqual($m.strand, 1);
        testing.assertEqual($m.mismatches, 0);
    }
}

func testFindMisses() {
    testing.assertEqual(len(find("ATGCATGC", "TTTT")), 0);
}

func testFindWorksOnProtein() {
    # a residue search, not a strand operation
    def hits as list of Match init find(INS_PROTEIN, "LLPLL");
    testing.assertEqual(len($hits), 1);
    testing.assertEqual(INS_PROTEIN[$hits[0].start..$hits[0].end], "LLPLL");
}

func testFindRejectsBadInput() {
    testing.assertThrows("findPatternTooLong", "bioinformatics");
    testing.assertThrows("findEmptyPattern", "bioinformatics");
}

func findPatternTooLong() {
    find("ATG", "ATGCATGC");
}

func findEmptyPattern() {
    find("ATGC", "");
}

# ------------------------------------------------------------------ IUPAC find

func testFindIupac() {
    # W is A or T, so GGWCC matches both spellings
    def hits as list of Match init findIupac("GGACCTTGGTCC", "GGWCC");
    testing.assertEqual(len($hits), 2);
    testing.assertEqual($hits[0].text, "GGACC");
    testing.assertEqual($hits[1].text, "GGTCC");
}

func testFindIupacReportsTheConcreteMatch() {
    # the pattern is degenerate; what comes back is what was actually found
    def hits as list of Match init findIupac("ATGCAAGCTT", "AAGCTT");
    testing.assertEqual(len($hits), 1);
    testing.assertEqual($hits[0].text, "AAGCTT");
    testing.assertEqual($hits[0].text, POLYLINKER[62..68]);
}

func testFindIupacNMatchesAnything() {
    testing.assertEqual(len(findIupac("ATGC", "NNNN")), 1);
    testing.assertEqual(len(findIupac("ATGC", "N")), 4);
    testing.assertEqual(len(findIupac("ATGCATGC", "GANTC")), 0);
}

func testIupacIsOneWay() {
    # a pattern's N admits a subject's A; a pattern's A does not admit an N
    testing.assertEqual(len(findIupac("AAAA", "N")), 4);
    testing.assertEqual(len(findIupac("NNNN", "A")), 0);
}

func testFindIupacAgreesWithFindOnConcretePatterns() {
    # with no ambiguity code the two must agree exactly
    def exact as list of Match init find(POLYLINKER, "GGATCC");
    def iupac as list of Match init findIupac(POLYLINKER, "GGATCC");
    testing.assertEqual(len($exact), len($iupac));
    for (def i as int init 0; $i < len($exact); $i = $i + 1) {
        testing.assertEqual($exact[$i].start, $iupac[$i].start);
    }
}

# ------------------------------------------------------------------ fuzzy find

func testFindFuzzyExactIsZeroMismatches() {
    def hits as list of Match init findFuzzy("ATGCATGC", "ATGC", 0);
    testing.assertEqual(len($hits), 2);
    for (def m in $hits) {
        testing.assertEqual($m.mismatches, 0);
    }
}

func testFindFuzzyTolerance() {
    # ATTC differs from ATGC at one position
    def none as list of Match init findFuzzy("TTATTCTT", "ATGC", 0);
    testing.assertEqual(len($none), 0);
    def one as list of Match init findFuzzy("TTATTCTT", "ATGC", 1);
    testing.assertEqual(len($one), 1);
    testing.assertEqual($one[0].mismatches, 1);
    testing.assertEqual($one[0].text, "ATTC");
}

func testFindFuzzyToleranceIsMonotone() {
    # PROPERTY: allowing more mismatches can only find more, never fewer
    def previous as int init 0;
    for (def tolerance as int init 0; $tolerance <= 6; $tolerance = $tolerance + 1) {
        def count as int init len(findFuzzy(POLYLINKER, "GAATTC", $tolerance));
        testing.assertTrue($count >= $previous);
        $previous = $count;
    }
}

func testFindFuzzyAtFullToleranceMatchesEverywhere() {
    # tolerating as many mismatches as there are positions matches at every offset
    def hits as list of Match init findFuzzy("ATGCATGC", "ATGC", 4);
    testing.assertEqual(len($hits), 5);
}

func testFindFuzzyRejectsNegativeTolerance() {
    testing.assertThrows("fuzzyNegative", "bioinformatics");
}

func fuzzyNegative() {
    findFuzzy("ATGC", "AT", -1);
}

# --------------------------------------------------------------------- ORFs

func testOrfFindsAGene() {
    # a clean construct: start, eight codons, stop
    def dna as string init "TTTATGGCCCTGTGGATGCGCCTCCTGTAATTT";
    def found as list of Orf init orfs($dna, 3);
    testing.assertTrue(len($found) >= 1);
    def first as Orf init $found[0];
    testing.assertEqual($first.strand, 1);
    testing.assertTrue($first.stopped);
    testing.assertEqual($first.protein, "MALWMRLL");
    # the reported span slices back to the ORF, stop codon included
    testing.assertEqual($dna[$first.start..$first.end], $first.nucleotides);
    testing.assertTrue(strings.startsWith($first.nucleotides, "ATG"));
    testing.assertTrue(strings.endsWith($first.nucleotides, "TAA"));
}

func testOrfCoordinatesAlwaysSliceBack() {
    # THE COORDINATE CONTRACT, forward and reverse. A reverse-strand ORF is
    # reported in the original sequence's coordinates, so slicing the original
    # and reverse-complementing must reproduce the ORF as read.
    def dna as string init POLYLINKER;
    for (def o in orfs($dna, 4)) {
        if ($o.strand > 0) {
            testing.assertEqual($dna[$o.start..$o.end], $o.nucleotides);
        } else {
            testing.assertEqual(reverseComplement($dna[$o.start..$o.end]), $o.nucleotides);
        }
        testing.assertEqual($o.length, $o.end - $o.start);
        testing.assertTrue($o.start >= 0 and $o.end <= len($dna));
    }
}

func testOrfProteinMatchesItsNucleotides() {
    for (def o in orfs(POLYLINKER, 3)) {
        def full as string init translate($o.nucleotides);
        if ($o.stopped) {
            testing.assertEqual($full[0..len($full) - 1], $o.protein);
        } else {
            testing.assertEqual($full, $o.protein);
        }
    }
}

func testOrfMinimumLengthFilters() {
    # PROPERTY: raising the threshold can only remove ORFs
    def previous as int init 999;
    for (def minAa as int init 1; $minAa <= 20; $minAa = $minAa + 1) {
        def count as int init len(orfs(POLYLINKER, $minAa));
        testing.assertTrue($count <= $previous);
        $previous = $count;
    }
    # and every survivor is at least that long
    for (def o in orfs(POLYLINKER, 6)) {
        testing.assertTrue(len($o.protein) >= 6);
    }
}

func testOrfUnstoppedIsReportedAsSuch() {
    # the insulin fragment has no in-frame stop - it is a fragment, and the gene
    # continues past it
    def found as list of Orf init orfs(INS_CDS, 20);
    testing.assertTrue(len($found) >= 1);
    def first as Orf init $found[0];
    testing.assertFalse($first.stopped);
    testing.assertEqual($first.protein, INS_PROTEIN);
}

func testOrfsCoverAllSixFrames() {
    def frames as map of int to bool;
    for (def o in orfs(POLYLINKER, 1)) {
        $frames[$o.frame] = true;
    }
    # every frame reported is one of the six, and the sign matches the strand
    for (def frame in $frames) {
        testing.assertTrue($frame >= -3 and $frame <= 3 and $frame != 0);
    }
    for (def o in orfs(POLYLINKER, 1)) {
        if ($o.frame > 0) {
            testing.assertEqual($o.strand, 1);
        } else {
            testing.assertEqual($o.strand, -1);
        }
    }
}

func testOrfStartCodonsFollowTheTable() {
    # GTG starts a gene in bacteria but not under the standard code, so the same
    # sequence yields a longer ORF under table 11 - which is why orfsWith takes
    # a table id rather than a codon map
    def dna as string init "GTGGCCCTGTGGATGCGCCTCCTGTAA";
    def standard as list of Orf init orfsWith($dna, 3, 1);
    def bacterial as list of Orf init orfsWith($dna, 3, 11);
    testing.assertTrue(len($bacterial) > 0);
    testing.assertTrue(len($standard) > 0);
    # bacterial opens on the GTG at the very start
    testing.assertTrue(strings.startsWith($bacterial[0].nucleotides, "GTG"));
    testing.assertEqual($bacterial[0].protein, "VALWMRLL");
    # the standard code does not accept GTG, so its frame-1 ORF opens later and
    # is correspondingly shorter. It opens on CTG rather than the internal ATG,
    # because CTG is a standard start codon too - which is exactly the sort of
    # thing a table lookup gets right and an assumption gets wrong.
    testing.assertFalse(strings.startsWith($standard[0].nucleotides, "GTG"));
    testing.assertTrue(len($standard[0].protein) < len($bacterial[0].protein));
    testing.assertTrue(lists.contains(startCodons(1), "CTG"));
    testing.assertFalse(lists.contains(startCodons(1), "GTG"));
    testing.assertTrue(lists.contains(startCodons(11), "GTG"));
}

func testOrfRejectsBadInput() {
    testing.assertThrows("orfNegativeMin", "bioinformatics");
    testing.assertThrows("orfUnknownTable", "bioinformatics");
}

func orfNegativeMin() {
    orfs("ATGCATGCTAA", -1);
}

func orfUnknownTable() {
    orfsWith("ATGCATGCTAA", 1, 7);
}

# --------------------------------------------------------------- CpG islands

func testCpgIslandFound() {
    # a CG-dense stretch flanked by AT-rich sequence
    def dna as string init "AAAATTTTAAAATTTT" + strings.repeat("CG", 20) + "AAAATTTTAAAATTTT";
    def islands as list of Region init cpgIslands($dna, 20, 50.0, 0.6);
    testing.assertEqual(len($islands), 1);
    def island as Region init $islands[0];
    testing.assertTrue($island.gc > 60.0);
    testing.assertTrue($island.obsExp > 2.0);
    testing.assertEqual($island.length, $island.end - $island.start);
    # It covers the CG block, which sits at [16, 56). The island runs a little
    # wider because a window straddling the boundary is still half CG-dense and
    # so still qualifies - island edges are approximate by nature, and merging
    # overlapping windows is what produces the reported span.
    testing.assertTrue($island.start <= 16);
    testing.assertTrue($island.end >= 56);
}

func testNoIslandInAtRichSequence() {
    testing.assertEqual(len(cpgIslands(strings.repeat("ATAT", 40), 20, 50.0, 0.6)), 0);
}

func testGcRichIsNotEnough() {
    # A GC-rich stretch with no CpG dinucleotide at all fails the ratio test,
    # which is the whole point of the observed/expected criterion. Note how hard
    # such a sequence is to construct: any alternating GC block has a C followed
    # by a G at some junction, so this needs every G before every C.
    def dna as string init strings.repeat("G", 60) + strings.repeat("C", 60);
    testing.assertTrue(gcContent($dna) == 100.0);
    testing.assertEqual(len(find($dna, "CG")), 0);
    testing.assertEqual(len(cpgIslands($dna, 20, 50.0, 0.6)), 0);
}

func testCpgRegionStatisticsAreConsistent() {
    def dna as string init "AAAATTTTAAAATTTT" + strings.repeat("CG", 20) + "AAAATTTTAAAATTTT";
    for (def r in cpgIslands($dna, 20, 50.0, 0.6)) {
        # the reported GC is the GC of the reported span
        testing.assertTrue(close($r.gc, gcContent($dna[$r.start..$r.end]), 0.000001));
    }
}

func testCpgShortSequenceYieldsNothing() {
    testing.assertEqual(len(cpgIslands("CGCG", 200, 50.0, 0.6)), 0);
}

func testCpgRejectsBadParameters() {
    testing.assertThrows("cpgZeroWindow", "bioinformatics");
    testing.assertThrows("cpgBadGc", "bioinformatics");
}

func cpgZeroWindow() {
    cpgIslands("ATGCATGC", 0, 50.0, 0.6);
}

func cpgBadGc() {
    cpgIslands("ATGCATGC", 4, 150.0, 0.6);
}

# ---------------------------------------------------------- restriction sites

func testPolylinkerSitesAreFound() {
    # the construct in testdata/ was built with exactly one site each
    def names as list of string init ["EcoRI", "KpnI", "BamHI", "HindIII", "PstI", "ApaI"];
    def sites as list of Site init restrictionSites(POLYLINKER, $names);
    testing.assertEqual(len($sites), 6);
    for (def site in $sites) {
        # the reported span slices back to the recognition sequence
        testing.assertEqual(POLYLINKER[$site.start..$site.end], $site.site);
        testing.assertEqual($site.site, enzyme($site.enzyme).site);
        testing.assertEqual($site.cut, $site.start + enzyme($site.enzyme).cutTop);
    }
}

func testSitesComeBackOrdered() {
    def names as list of string init ["PstI", "EcoRI", "BamHI", "HindIII", "KpnI", "ApaI"];
    def sites as list of Site init restrictionSites(POLYLINKER, $names);
    # given in scrambled order, returned in cut order
    for (def i as int init 1; $i < len($sites); $i = $i + 1) {
        testing.assertTrue($sites[$i].cut >= $sites[$i - 1].cut);
    }
    # the polylinker was laid out in this order
    testing.assertEqual($sites[0].enzyme, "EcoRI");
    testing.assertEqual($sites[len($sites) - 1].enzyme, "ApaI");
}

func testDegenerateSiteMatches() {
    # EcoRII reads CCWGG, so it cuts both CCAGG and CCTGG
    def dna as string init "AAACCAGGTTTCCTGGAAA";
    def sites as list of Site init restrictionSites($dna, ["EcoRII"]);
    testing.assertEqual(len($sites), 2);
    testing.assertEqual($sites[0].site, "CCAGG");
    testing.assertEqual($sites[1].site, "CCTGG");
}

func testEnzymesCuttingAndNonCutters() {
    def cutting as list of string init enzymesCutting(POLYLINKER);
    def quiet as list of string init nonCutters(POLYLINKER);
    # the two partition the whole table, with no overlap
    testing.assertEqual(len($cutting) + len($quiet), len(enzymes()));
    for (def name in $cutting) {
        testing.assertFalse(lists.contains($quiet, $name));
        testing.assertTrue(len(restrictionSites(POLYLINKER, [$name])) > 0);
    }
    # the six the construct was built around are all in the cutting set
    for (def name in ["EcoRI", "KpnI", "BamHI", "HindIII", "PstI", "ApaI"]) {
        testing.assertTrue(lists.contains($cutting, $name));
    }
}

func testNonCutterReallyDoesNotCut() {
    for (def name in nonCutters(POLYLINKER)) {
        testing.assertEqual(len(restrictionSites(POLYLINKER, [$name])), 0);
    }
}

func testUnknownEnzymeRaisesInSearch() {
    testing.assertThrows("sitesUnknownEnzyme", "bioinformatics");
}

func sitesUnknownEnzyme() {
    restrictionSites("ATGCATGC", ["NotAnEnzyme"]);
}

# --------------------------------------------------------------- the digest

func testLinearDigestReassembles() {
    # THE INVARIANT: the fragments, in order, are the original sequence
    def names as list of string init ["EcoRI", "BamHI", "HindIII"];
    def pieces as list of Fragment init restrictionDigest(POLYLINKER, $names, false);
    testing.assertEqual(len($pieces), 4);
    def rebuilt as string init "";
    def total as int init 0;
    for (def f in $pieces) {
        $rebuilt = $rebuilt + $f.sequence;
        $total = $total + $f.length;
        testing.assertEqual($f.length, len($f.sequence));
    }
    testing.assertEqual($rebuilt, POLYLINKER);
    testing.assertEqual($total, len(POLYLINKER));
}

func testLinearDigestCutCountIsFragmentCount() {
    # PROPERTY: n cuts give n+1 linear fragments
    for (def names in [["EcoRI"], ["EcoRI", "BamHI"], ["EcoRI", "BamHI", "PstI"]]) {
        def cuts as int init len(restrictionSites(POLYLINKER, $names));
        testing.assertEqual(len(restrictionDigest(POLYLINKER, $names, false)), $cuts + 1);
    }
}

func testUncutSequenceIsOneFragment() {
    # DraI reads TTTAAA, which the polylinker does not contain
    testing.assertEqual(len(restrictionSites(POLYLINKER, ["DraI"])), 0);
    def pieces as list of Fragment init restrictionDigest(POLYLINKER, ["DraI"], false);
    testing.assertEqual(len($pieces), 1);
    testing.assertEqual($pieces[0].sequence, POLYLINKER);
}

func testCircularDigestReassembles() {
    # n cuts give n fragments on a circle, and their lengths still total the
    # original - the wrapping fragment is the one whose end precedes its start
    def names as list of string init ["EcoRI", "BamHI", "HindIII"];
    def pieces as list of Fragment init restrictionDigest(POLYLINKER, $names, true);
    testing.assertEqual(len($pieces), 3);
    def total as int init 0;
    def wrapped as int init 0;
    for (def f in $pieces) {
        $total = $total + $f.length;
        testing.assertEqual($f.length, len($f.sequence));
        if ($f.end < $f.start) {
            $wrapped = $wrapped + 1;
        }
    }
    testing.assertEqual($total, len(POLYLINKER));
    testing.assertEqual($wrapped, 1);
}

func testCircularSingleCutLinearises() {
    # cutting a plasmid once gives one linear fragment of the full length,
    # which is exactly what linearising is
    def pieces as list of Fragment init restrictionDigest(POLYLINKER, ["EcoRI"], true);
    testing.assertEqual(len($pieces), 1);
    testing.assertEqual($pieces[0].length, len(POLYLINKER));
}

func testDigestFragmentsSliceBack() {
    for (def f in restrictionDigest(POLYLINKER, ["EcoRI", "PstI"], false)) {
        testing.assertEqual(POLYLINKER[$f.start..$f.end], $f.sequence);
    }
}

# ------------------------------------------------- ordering, after the audit

func testSortSitesOrdersTotally() {
    # the insertion sort was rewritten to shift in place rather than rebuild the
    # list; the ordering it produces must not have changed
    def names as list of string init ["PstI", "ApaI", "EcoRI", "KpnI", "BamHI", "HindIII"];
    def sites as list of Site init restrictionSites(POLYLINKER, $names);
    testing.assertEqual(len($sites), 6);
    for (def i as int init 1; $i < len($sites); $i = $i + 1) {
        # ascending by cut, and by enzyme name where cuts coincide
        testing.assertTrue($sites[$i].cut >= $sites[$i - 1].cut);
        if ($sites[$i].cut == $sites[$i - 1].cut) {
            testing.assertTrue($sites[$i].enzyme >= $sites[$i - 1].enzyme);
        }
    }
}

func testSiteOrderIsIndependentOfEnzymeOrder() {
    # asking in a different order must give the same answer
    def one as list of Site init restrictionSites(POLYLINKER, ["EcoRI", "BamHI", "PstI"]);
    def other as list of Site init restrictionSites(POLYLINKER, ["PstI", "EcoRI", "BamHI"]);
    testing.assertEqual(len($one), len($other));
    for (def i as int init 0; $i < len($one); $i = $i + 1) {
        testing.assertEqual($one[$i].enzyme, $other[$i].enzyme);
        testing.assertEqual($one[$i].cut, $other[$i].cut);
    }
}

func testCoincidentCutsAreOrderedByName() {
    # SmaI and XmaI share a site but cut it differently; HpaII and MspI cut it
    # identically, so their cuts coincide and the name decides
    def dna as string init "AAAACCGGAAAA";
    def sites as list of Site init restrictionSites($dna, ["MspI", "HpaII"]);
    testing.assertEqual(len($sites), 2);
    testing.assertEqual($sites[0].cut, $sites[1].cut);
    testing.assertEqual($sites[0].enzyme, "HpaII");
    testing.assertEqual($sites[1].enzyme, "MspI");
}

func testPresenceScanAgreesWithTheFullScan() {
    # enzymesCutting now uses a presence-only scan that stops at the first hit
    # and never builds a Match. It must agree exactly with the full scan.
    def cutting as list of string init enzymesCutting(POLYLINKER);
    for (def e in enzymes()) {
        def found as bool init len(restrictionSites(POLYLINKER, [$e.name])) > 0;
        testing.assertEqual(lists.contains($cutting, $e.name), $found);
    }
}
