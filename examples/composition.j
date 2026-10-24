#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - what a sequence is made of, and what that implies.
#
#   jennifer run examples/composition.j
#
# Run it from the repository root - the testdata/ path is relative to the
# working directory, not to this file.
#
# GC content, molecular weight, melting temperature, isoelectric point, GRAVY
# and codon usage, each with the parameter set it uses named.

use io;
use fs;
use math;
use convert;

import "../src/bioinformatics.j" as bioseq;

def fasta as string init fs.readString("testdata/sequences.fasta");
def records as list of bioseq.Record init bioseq.parseFasta($fasta);
def gene as string init $records[0].sequence;
def protein as string init $records[2].sequence;

# ---------------------------------------------------------------- DNA summary

def stats as bioseq.DnaStats init bioseq.dnaStats($gene);

io.printf("%s\n\n", $records[0].id);
io.printf("  length        %d nt\n", $stats.length);
io.printf("  GC            %f|prec=1%%\n", $stats.gc);
io.printf("  AT            %f|prec=1%%\n", $stats.at);
io.printf("  mass          %f|prec=1 Da  (single strand, 5' OH)\n", $stats.molecularWeight);
# The Wallace rule is calibrated for 14-20 base oligos. Applied to a 72-base
# gene fragment it extrapolates to a number with no physical meaning - shown
# here with the caveat rather than hidden, because a summary struct that
# silently reports 250 C is worse than one that explains itself.
io.printf("  Tm (Wallace)  %f|prec=1 C", $stats.meltingTemp);
if ($stats.length > 30) {
    io.printf("   <- extrapolated well past the rule's range");
}
io.printf("\n");

io.printf("\n  base composition\n");
for (def base in bioseq.baseCount($gene)) {
    io.printf(
        "    %s  %d|pad=3  %f|prec=1%%\n",
        $base,
        bioseq.baseCount($gene)[$base],
        bioseq.baseFrequency($gene)[$base]);
}

# ------------------------------------------------------------ GC along the way

# A sequence's average GC says little; its local variation says a great deal.
# This is the series a GC plot draws.
io.printf("\n  GC in 18-base windows\n\n");
for (def value in bioseq.gcWindows($gene, 18, 9)) {
    io.printf("    %f|prec=1|pad=5%%  ", $value);
    def bars as int init convert.toInt(math.round($value / 5.0));
    for (def i as int init 0; $i < $bars; $i = $i + 1) {
        io.printf("#");
    }
    io.printf("\n");
}

# ----------------------------------------------------- why the model matters

# The Wallace rule counts bases. The nearest-neighbour model reads them in
# order, and stacking energy depends on which base follows which - so two
# sequences of identical composition melt twelve degrees apart.
io.printf("\n\nmelting temperature, two ways\n\n");
io.printf("  %s|pad=10|align=left %s|pad=10 %s\n", "oligo", "Wallace", "nearest-neighbour");

for (def oligo in ["GCGCGC", "GGGCCC", "ATGCATGCATGCATGCATGC"]) {
    io.printf(
        "  %s|pad=10|align=left %f|prec=1|pad=10 %f|prec=1\n",
        $oligo,
        bioseq.meltingTemp($oligo),
        bioseq.meltingTempNn($oligo));
}

io.printf("\n  GCGCGC and GGGCCC have identical composition.\n");
io.printf(
    "  Wallace cannot tell them apart; the NN model puts them %f|prec=1 C apart.\n",
    bioseq.meltingTempNn("GCGCGC") - bioseq.meltingTempNn("GGGCCC"));

# Conditions matter too, which is why they are parameters and not constants.
io.printf("\n  a 20-mer at 250 nM primer\n");
for (def salt in [25.0, 50.0, 100.0, 200.0]) {
    io.printf(
        "    %f|prec=0|pad=4 mM Na+   %f|prec=1 C\n",
        $salt,
        bioseq.meltingTempNnWith("ATGCATGCATGCATGCATGC", 250.0, $salt));
}

# ------------------------------------------------------------ protein chemistry

def pstats as bioseq.ProteinStats init bioseq.proteinStats($protein);

io.printf("\n\n%s\n\n", $records[2].id);
io.printf("  sequence      %s\n", $protein);
io.printf("  length        %d residues\n", $pstats.length);
io.printf("  mass          %f|prec=2 Da\n", $pstats.molecularWeight);
io.printf("  pI            %f|prec=2   (EMBOSS iep pKa set)\n", $pstats.isoelectricPoint);
io.printf("  GRAVY         %f|prec=3\n", $pstats.gravy);
io.printf("  aromaticity   %f|prec=3\n", $pstats.aromaticity);

# GRAVY is positive here, and strikingly so: this fragment is the insulin
# signal peptide, whose job is to be hydrophobic enough to cross a membrane.
if ($pstats.gravy > 0.0) {
    io.printf("\n  positive GRAVY - hydrophobic overall, as a signal peptide should be\n");
}

# Charge against pH is the curve the pI is read off. It crosses zero exactly
# once, which is what makes the bisection safe.
io.printf("\n  net charge against pH\n\n");
for (def step as int init 0; $step <= 14; $step = $step + 2) {
    def ph as float init 1.0 * $step;
    def charge as float init bioseq.netCharge($protein, $ph);
    def marker as string init "";
    if (math.abs($charge) < 0.25) {
        $marker = "  <- near the pI";
    }
    io.printf("    pH %f|prec=1|pad=4   %f|prec=2|pad=7|sign=always%s\n", $ph, $charge, $marker);
}

# --------------------------------------------------------------- codon usage

# A gene leans on some synonymous codons more than others, and that bias is
# what codon-optimisation software rewrites.
io.printf("\n\ncodon usage in %s\n\n", $records[0].id);
io.printf(
    "  %s|pad=7|align=left %s|pad=4 %s|pad=7 %s\n",
    "codon",
    "aa",
    "count",
    "share of that residue's codons");

for (def usage in bioseq.codonUsage($gene)) {
    if ($usage.count > 1) {
        io.printf(
            "  %s|pad=7|align=left %s|pad=4 %d|pad=7 %f|prec=2\n",
            $usage.codon,
            $usage.aminoAcid,
            $usage.count,
            $usage.fraction);
    }
}

io.printf(
    "\n  (only codons used more than once are shown; %d distinct codons in all)\n",
    len(bioseq.codonUsage($gene)));
