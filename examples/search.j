#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - finding things inside a sequence.
#
#   jennifer run examples/search.j
#
# Run it from the repository root - the testdata/ path is relative to the
# working directory, not to this file.
#
# Pattern search with and without ambiguity, approximate matching, open reading
# frames, CpG islands, and a restriction digest.

use io;
use fs;
use strings;

import "../src/bioinformatics.j" as bioseq;

def fasta as string init fs.readString("testdata/sequences.fasta");
def records as list of bioseq.Record init bioseq.parseFasta($fasta);
def gene as string init $records[0].sequence;
def construct as string init $records[1].sequence;

# --------------------------------------------------------------- exact search

io.printf("%s (%d nt)\n\n", $records[1].id, len($construct));

io.printf("  GGATCC occurs at:");
for (def m in bioseq.find($construct, "GGATCC")) {
    io.printf(" %d", $m.start);
}
io.printf("\n");

# Overlapping occurrences are all reported: a site inside another site is a real
# finding, and a non-overlapping scan would silently drop it.
io.printf("  \"AA\" in \"AAAA\" occurs %d times, not one\n", len(bioseq.find("AAAA", "AA")));

# ------------------------------------------------------------ ambiguous search

# An IUPAC pattern searches for a family of sequences at once. EcoRII reads
# CCWGG - W is A or T - so one pattern finds both spellings.
io.printf("\n  IUPAC pattern CCWGG against CCAGG / CCTGG\n");
for (def m in bioseq.findIupac("AAACCAGGTTTCCTGGAAA", "CCWGG")) {
    io.printf("    at %d|pad=3  found %s\n", $m.start, $m.text);
}

# ------------------------------------------------------------ approximate search

# A primer that no longer matches exactly may still bind. findFuzzy counts the
# substitutions, so exact hits are distinguishable from near misses.
io.printf("\n  GAATTC allowing mismatches\n");
for (def tolerance in [0, 1, 2]) {
    def hits as list of bioseq.Match init bioseq.findFuzzy($construct, "GAATTC", $tolerance);
    io.printf("    up to %d mismatch  %d|pad=3 site", $tolerance, len($hits));
    if (len($hits) != 1) {
        io.printf("s");
    }
    io.printf("\n");
}

# ------------------------------------------------------------------- ORFs

# When a sequence's frame is unknown, all six are scanned. Each ORF comes back
# in the coordinates of the sequence passed in - including the reverse-strand
# ones, which are mapped back rather than reported in reverse-complement space.
io.printf("\n\nopen reading frames in %s\n\n", $records[0].id);
io.printf("  %s|pad=6 %s|pad=12 %s|pad=8 %s\n", "frame", "span", "stop?", "protein");

for (def o in bioseq.orfs($gene, 8)) {
    def closed as string init "runs on";
    if ($o.stopped) {
        $closed = "stop";
    }
    def span as string init io.sprintf("%d-%d", $o.start, $o.end);
    io.printf(
        "  %d|sign=always|pad=6 %s|pad=12|align=left %s|pad=8|align=left %s\n",
        $o.frame,
        $span,
        $closed,
        $o.protein);
}

# The insulin fragment's frame-1 ORF runs off the end without a stop, which is
# exactly right: it is a fragment, and the gene continues past what we have.
io.printf("\n  frame 1 has no in-frame stop - the CDS continues past this fragment\n");

# A reverse-strand ORF's coordinates index the original sequence, so slicing it
# out and reverse-complementing gives back the reading frame.
for (def o in bioseq.orfs($gene, 8)) {
    if ($o.strand < 0) {
        io.printf(
            "\n  frame %d spans [%d, %d) of the forward sequence;\n",
            $o.frame,
            $o.start,
            $o.end);
        io.printf(
            "  reverse-complementing that slice reproduces it: %t\n",
            bioseq.reverseComplement($gene[$o.start..$o.end]) == $o.nucleotides);
    }
}

# ---------------------------------------------------------------- CpG islands

# CpG is scarce in vertebrate genomes - methylated cytosine deaminates to
# thymine, so the dinucleotide has been eroded everywhere except where
# methylation is kept off, which is where promoters sit.
def promoterLike as string init "AAAATTTTAAAATTTT" + strings.repeat("CG", 20) + "AAAATTTTAAAATTTT";

io.printf("\n\nCpG islands in a synthetic CG-dense stretch\n\n");
for (def r in bioseq.cpgIslands($promoterLike, 20, 50.0, 0.6)) {
    def span as string init io.sprintf("[%d, %d)", $r.start, $r.end);
    io.printf(
        "  %s|pad=12|align=left %d bases, GC %f|prec=1%%, observed/expected CpG %f|prec=2\n",
        $span,
        $r.length,
        $r.gc,
        $r.obsExp);
}

# GC-rich is not the same as CpG-rich, and the ratio test is what separates
# them. This sequence is 100% GC and contains no CpG at all.
def gcOnly as string init strings.repeat("G", 60) + strings.repeat("C", 60);
io.printf(
    "\n  a 100%% GC sequence with no CpG dinucleotide: %d islands\n",
    len(bioseq.cpgIslands($gcOnly, 20, 50.0, 0.6)));

# ----------------------------------------------------------- restriction sites

io.printf("\n\nrestriction map of %s\n\n", $records[1].id);

def panel as list of string init ["EcoRI", "KpnI", "BamHI", "HindIII", "PstI", "ApaI"];
io.printf(
    "  %s|pad=10|align=left %s|pad=10|align=left %s|pad=8 %s\n",
    "enzyme",
    "site",
    "at",
    "cuts after");

for (def site in bioseq.restrictionSites($construct, $panel)) {
    io.printf(
        "  %s|pad=10|align=left %s|pad=10|align=left %d|pad=8 %d\n",
        $site.enzyme,
        $site.site,
        $site.start,
        $site.cut);
}

# Which enzymes leave the construct alone is the more useful question when the
# job is to keep an insert intact.
def quiet as list of string init bioseq.nonCutters($construct);
io.printf(
    "\n  %d of the %d enzymes in the table cut it; %d leave it alone\n",
    len(bioseq.enzymesCutting($construct)),
    len(bioseq.enzymes()),
    len($quiet));

# ---------------------------------------------------------------- the digest

io.printf("\n\ndigest with EcoRI + BamHI + HindIII\n\n");

def cutters as list of string init ["EcoRI", "BamHI", "HindIII"];

io.printf("  linear\n");
for (def f in bioseq.restrictionDigest($construct, $cutters, false)) {
    def span as string init io.sprintf("[%d, %d)", $f.start, $f.end);
    io.printf("    %s|pad=14|align=left %d|pad=4 bp\n", $span, $f.length);
}

# The fragments, in order, are the original sequence - the invariant a digest
# has to satisfy.
def rebuilt as string init "";
for (def f in bioseq.restrictionDigest($construct, $cutters, false)) {
    $rebuilt = $rebuilt + $f.sequence;
}
io.printf("    fragments reassemble to the original: %t\n", $rebuilt == $construct);

# On a circle the same three cuts give three fragments rather than four, and one
# of them wraps the origin - its end precedes its start.
io.printf("\n  circular\n");
for (def f in bioseq.restrictionDigest($construct, $cutters, true)) {
    def note as string init "";
    if ($f.end < $f.start) {
        $note = "  <- wraps the origin";
    }
    def span as string init io.sprintf("[%d, %d)", $f.start, $f.end);
    io.printf("    %s|pad=14|align=left %d|pad=4 bp%s\n", $span, $f.length, $note);
}

# Cutting a plasmid once linearises it: one fragment, full length.
def once as list of bioseq.Fragment init bioseq.restrictionDigest($construct, ["EcoRI"], true);
io.printf(
    "\n  one cut on a circle gives %d fragment of %d bp - a linearised plasmid\n",
    len($once),
    $once[0].length);
