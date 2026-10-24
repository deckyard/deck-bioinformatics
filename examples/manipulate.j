#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - generating, shuffling, mutating and slicing.
#
#   jennifer run examples/manipulate.j
#
# Run it from the repository root - the testdata/ path is relative to the
# working directory, not to this file.
#
# Random sequences, shuffle, mutation and windows, used to build null models: a
# finding in a real sequence means something only against the same question
# asked of sequence with no structure.

use io;
use fs;
use math;

import "../src/bioinformatics.j" as bioseq;

# One seed makes the whole run reproducible. bioinformatics keeps no generator of its
# own, so this single call governs every random verb below.
math.randSeed(20260828);

def fasta as string init fs.readString("testdata/sequences.fasta");
def records as list of bioseq.Record init bioseq.parseFasta($fasta);
def gene as string init $records[0].sequence;
def construct as string init $records[1].sequence;

# ------------------------------------------------------------------- generation

io.printf("random sequences\n\n");
io.printf("  uniform DNA      %s\n", bioseq.randomDna(40));
io.printf("  uniform RNA      %s\n", bioseq.randomRna(40));
io.printf("  uniform protein  %s\n", bioseq.randomProtein(40));

# Weights are relative, not percentages - they get normalised - and a residue
# left out simply never appears.
io.printf(
    "\n  AT-rich (70/30)  %s\n",
    bioseq.randomDnaWith(40, {"A": 35.0, "T": 35.0, "C": 15.0, "G": 15.0}));
io.printf("  purines only     %s\n", bioseq.randomDnaWith(40, {"A": 1.0, "G": 1.0}));

# A uniform random protein is not what a real one looks like. Matching the
# composition of a real sequence gives a far better null.
def composition as map of string to float init bioseq.residueFrequency($records[2].sequence);
io.printf("\n  real protein     %s\n", $records[2].sequence);
io.printf(
    "  same composition %s\n",
    bioseq.randomProteinWith(len($records[2].sequence), $composition));

# ---------------------------------------------------------------- reproducible

math.randSeed(1);
def first as string init bioseq.randomDna(30);
math.randSeed(1);
def second as string init bioseq.randomDna(30);
io.printf("\n  the same seed gives the same sequence: %t\n", $first == $second);

# --------------------------------------------------------------------- shuffle

# A shuffle keeps composition exactly and destroys order. That is what makes it
# the honest comparison: it differs from the real sequence in the one respect
# being tested.
io.printf("\n\nshuffling keeps composition, destroys order\n\n");
io.printf("  original   %s\n", $gene);
io.printf("  shuffled   %s\n", bioseq.shuffle($gene));
io.printf(
    "\n  GC content   original %f|prec=1%%   shuffled %f|prec=1%%\n",
    bioseq.gcContent($gene),
    bioseq.gcContent(bioseq.shuffle($gene)));

# ------------------------------------------------- is that site worth anything?

# The construct carries one EcoRI site. Is that surprising in a sequence of this
# length and composition? Shuffle it many times and count how often a site
# survives - that is the null distribution, and the answer.
io.printf("\n\nis the EcoRI site in %s surprising?\n\n", $records[1].id);

def trials as int init 500;
def withSite as int init 0;
for (def i as int init 0; $i < $trials; $i = $i + 1) {
    if (len(bioseq.restrictionSites(bioseq.shuffle($construct), ["EcoRI"])) > 0) {
        $withSite = $withSite + 1;
    }
}

io.printf(
    "  the real construct has %d EcoRI site\n",
    len(bioseq.restrictionSites($construct, ["EcoRI"])));
io.printf("  of %d shuffles with identical composition, %d carried one\n", $trials, $withSite);
io.printf(
    "  so a site like this arises by chance about %f|prec=1%% of the time\n",
    100.0 * $withSite / $trials);
io.printf("\n  (the site was placed there deliberately - but the method is the\n");
io.printf("   same one you would use when you do not already know that)\n");

# ---------------------------------------------------------------------- mutate

# Substitutions only, at a per-position rate. A mutated position always changes,
# so the realised rate matches the one asked for.
io.printf("\n\nmutation\n\n");
io.printf("  original     %s\n", $gene);

for (def rate in [0.01, 0.05, 0.2]) {
    def mutated as string init bioseq.mutate($gene, $rate, bioseq.DNA);
    def changed as int init 0;
    for (def i as int init 0; $i < len($gene); $i = $i + 1) {
        if ($mutated[$i..$i + 1] != $gene[$i..$i + 1]) {
            $changed = $changed + 1;
        }
    }
    io.printf("  rate %f|prec=2   %s  %d changed\n", $rate, $mutated, $changed);
}

# What a mutation does to the protein is the interesting part: most
# substitutions in a coding sequence are silent or conservative, and a few are
# not.
io.printf("\n  protein consequences at 5%%\n\n");
io.printf("    original  %s\n", bioseq.translate($gene));
for (def i as int init 0; $i < 3; $i = $i + 1) {
    io.printf("    mutated   %s\n", bioseq.translate(bioseq.mutate($gene, 0.05, bioseq.DNA)));
}

# ------------------------------------------------------------ slicing and windows

# subrange is the deck's one 1-based inclusive verb, named for it, so a paper
# saying "bases 13 to 24" transcribes directly.
io.printf("\n\nextracting\n\n");
io.printf("  bases 13-24 (as a paper would quote them)  %s\n", bioseq.subrange($gene, 13, 24));
io.printf("  the same slice in the deck's own coordinates  %s\n", $gene[12..24]);

io.printf("\n  24-base windows, stepping 12\n\n");
for (def w in bioseq.windows($gene, 24, 12)) {
    io.printf("    %s   GC %f|prec=1%%\n", $w, bioseq.gcContent($w));
}
