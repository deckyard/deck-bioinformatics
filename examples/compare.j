#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - comparing two sequences.
#
#   jennifer run examples/compare.j
#
# Run it from the repository root - the testdata/ path is relative to the
# working directory, not to this file.
#
# Global and local alignment, why a substitution matrix matters for protein,
# and the two distances.

use io;
use fs;
use math;
use strings;

import "../src/bioinformatics.j" as bioseq;

math.randSeed(20260828);

def fasta as string init fs.readString("testdata/sequences.fasta");
def records as list of bioseq.Record init bioseq.parseFasta($fasta);
def gene as string init $records[0].sequence;
def protein as string init $records[2].sequence;

# ------------------------------------------------------------ global alignment

# Two versions of the same gene: the original and a copy carrying substitutions
# and a small deletion. Global alignment is the right question here, because
# both sequences are the same thing throughout.
def variant as string init bioseq.mutate($gene, 0.08, bioseq.DNA);
$variant = $variant[0..30] + $variant[33..];

def dna as bioseq.AlignOptions init bioseq.dnaOptions();
def global as bioseq.Alignment init bioseq.alignGlobal($gene, $variant, $dna);

io.printf("global alignment: a gene against a mutated copy\n\n");
io.printf("%s", bioseq.formatAlignment($global, 60));
io.printf(
    "  score %d, %d columns, %f|prec=1%% identical, %d gap columns\n",
    $global.score,
    $global.length,
    bioseq.percentIdentity($global),
    $global.gaps);

# Strip the gaps and the two rows are the originals again - the invariant a
# global alignment has to satisfy.
def intact as bool init strings.replace($global.a, "-", "") == $gene and
    strings.replace($global.b, "-", "") == $variant;
io.printf("\n  both sequences survive intact: %t\n", $intact);

# ------------------------------------------------------------- local alignment

# A shared stretch buried in two otherwise unrelated sequences. Global alignment
# would spend its whole score on the flanks; local alignment ignores them.
def core as string init "ATGCGCCTCCTGCCCCTGCTG";
def one as string init bioseq.randomDna(30) + $core + bioseq.randomDna(30);
def other as string init bioseq.randomDna(20) + $core + bioseq.randomDna(40);

def local as bioseq.Alignment init bioseq.alignLocal($one, $other, $dna);
def globalOnSame as bioseq.Alignment init bioseq.alignGlobal($one, $other, $dna);

io.printf("\n\nlocal alignment: a shared core inside unrelated flanks\n\n");
io.printf(
    "  local   score %d|pad=5  %f|prec=1%% identical over %d columns\n",
    $local.score,
    bioseq.percentIdentity($local),
    $local.length);
io.printf(
    "  global  score %d|pad=5  %f|prec=1%% identical over %d columns\n",
    $globalOnSame.score,
    bioseq.percentIdentity($globalOnSame),
    $globalOnSame.length);

io.printf("\n  the local alignment found:\n");
io.printf("    %s\n", $local.a);
io.printf("    %s\n", $local.markup);
io.printf("    %s\n", $local.b);

# Its start offsets index the sequences that were passed in, so the match can be
# located in the originals - and the shared core is inside what it found.
io.printf(
    "\n  found at offset %d of the first sequence and %d of the second\n",
    $local.startA,
    $local.startB);
io.printf(
    "  and it contains the planted core: %t\n",
    strings.contains(strings.replace($local.a, "-", ""), $core));

# Note that it reaches past the core, which was planted at offset 30. That is
# Smith-Waterman working correctly, not a bug: it keeps extending outward while
# the running score rises, and with +2 for a match against -1 for a mismatch,
# random sequence breaks even at one match in three. So the boundaries of a
# local alignment are soft unless the flanks are genuinely dissimilar.
io.printf("  (the core was planted at offset 30 - the alignment reaches past it,\n");
io.printf("   because a +2/-1 scheme still scores random flanks above zero)\n");

# With flanks that cannot align at all, the boundaries land exactly.
def cleanOne as string init strings.repeat("T", 20) + $core + strings.repeat("T", 20);
def cleanOther as string init strings.repeat("A", 12) + $core + strings.repeat("A", 12);
def exact as bioseq.Alignment init bioseq.alignLocal($cleanOne, $cleanOther, $dna);
io.printf("\n  with poly-T against poly-A flanks, the boundaries are exact:\n");
io.printf(
    "    found %s at offsets %d and %d - planted at 20 and 12\n",
    $exact.a,
    $exact.startA,
    $exact.startB);

# ------------------------------------------------- why protein needs a matrix

# L against I is a conservative substitution that barely changes a protein.
# L against D is not. A match/mismatch score cannot tell them apart.
io.printf("\n\nwhy a substitution matrix matters for protein\n\n");
io.printf("  %s|pad=10|align=left %s\n", "pair", "BLOSUM62");
for (def pair in [["L", "L"], ["L", "I"], ["K", "R"], ["F", "Y"], ["L", "D"], ["W", "P"]]) {
    def note as string init "";
    def score as int init bioseq.substitutionScore(bioseq.BLOSUM62, $pair[0], $pair[1]);
    if ($pair[0] == $pair[1]) {
        $note = "identical";
    } elseif ($score > 0) {
        $note = "conservative - scores positively";
    } else {
        $note = "disruptive";
    }
    io.printf("  %s -> %s     %d|sign=always|pad=4   %s\n", $pair[0], $pair[1], $score, $note);
}

# The consequence: an alignment scored with BLOSUM62 reports similarity above
# identity, because it can see that some differences hardly matter.
# M->L, R->K, L->I and L->V are all conservative: the protein still works.
def conservative as string init "MALWLKILPLVALLALWGPDPAAA";
def prot as bioseq.Alignment init bioseq.alignGlobal(
    $protein,
    $conservative,
    bioseq.proteinOptions());

io.printf("\n  aligning two related proteins\n\n");
io.printf("    %s\n    %s\n    %s\n", $prot.a, $prot.markup, $prot.b);
io.printf(
    "\n    %f|prec=1%% identical, %f|prec=1%% similar\n",
    bioseq.percentIdentity($prot),
    bioseq.percentSimilarity($prot));
io.printf("    (| identical, : conservative substitution)\n");

# ------------------------------------------------------------------ distances

io.printf("\n\ndistances\n\n");

def a as string init "ATGCATGCATGC";
def b as string init "ATGAATGCATCC";
io.printf("  %s\n  %s\n", $a, $b);
io.printf(
    "    Hamming  %d   substitutions only, so lengths must match\n",
    bioseq.hammingDistance($a, $b));
io.printf("    edit     %d   insertions and deletions allowed too\n", bioseq.editDistance($a, $b));

# Where the two part company: an insertion shifts everything after it, so
# Hamming counts a difference at nearly every position while edit distance sees
# a single event.
def shifted as string init "TGCATGCATGCA";
io.printf("\n  %s\n  %s   (the first sequence, shifted by one)\n", $a, $shifted);
io.printf("    Hamming  %d   every position looks wrong\n", bioseq.hammingDistance($a, $shifted));
io.printf("    edit     %d   one deletion and one insertion\n", bioseq.editDistance($a, $shifted));

io.printf("\n  and Hamming refuses unequal lengths rather than guessing:\n");
try {
    bioseq.hammingDistance("ATGC", "ATG");
} catch (e) {
    io.printf("    %s\n", $e.message);
}
