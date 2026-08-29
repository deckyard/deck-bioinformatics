#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - turning one sequence into another.
#
#   jennifer run examples/transforms.j
#
# Complement, transcribe, translate, and read all six frames. No network, no
# files, no arguments.

use io;
use strings;

import "../src/bioinformatics.j" as bioseq;

# The start of the human preproinsulin coding sequence, as it might be pasted
# out of a database: wrapped, numbered, and in mixed case. filterDna strips all
# of that; every other verb tolerates the whitespace and digits but would raise
# on a genuine non-base, which is the difference between the two tiers.
def const PASTED as string init '
    1  atggccctgt ggatgcgcct cctgcccctg ctggcgctgc tggccctctg
   51  gggacctgac ccagccgcag cc
';

def gene as string init bioseq.filterDna(PASTED);

io.printf("gene            %s\n", $gene);
io.printf("length          %d nt\n\n", len($gene));

# ------------------------------------------------------------- the two strands

# complement() pairs the bases but leaves them in place, which means the result
# reads 3'-to-5' - almost never what you want on its own. reverseComplement()
# is the strand you would actually write down.
io.printf("complement      %s\n", bioseq.complement($gene));
io.printf("revComplement   %s\n", bioseq.reverseComplement($gene));

# Applying it twice returns the original, which is the quickest sanity check
# there is on a strand operation.
io.printf(
    "rc twice = gene %t\n\n",
    bioseq.reverseComplement(bioseq.reverseComplement($gene)) == $gene);

# ------------------------------------------------------------- transcription

# Transcription here is the coding-strand shorthand: T becomes U and nothing is
# complemented, because the mRNA matches the coding strand.
def rna as string init bioseq.transcribe($gene);
io.printf("transcribed     %s\n", $rna);
io.printf("back again      %t\n\n", bioseq.backTranscribe($rna) == $gene);

# ---------------------------------------------------------------- translation

def protein as string init bioseq.translate($gene);
io.printf("protein         %s\n", $protein);
io.printf("three-letter    %s\n", bioseq.oneToThree($protein));

# RNA translates directly - U is read as T for the codon lookup - so the two
# spellings of one gene give the same protein.
io.printf("from RNA        %t\n\n", bioseq.translate($rna) == $protein);

# ------------------------------------------------------- the genetic code matters

# The same codons under a different code. TGA is a stop in the nucleus and
# tryptophan in vertebrate mitochondria; a mitochondrial gene read under the
# standard code stops early at every one of them.
io.printf("codon TGA reads\n");
for (def id in [1, 2, 6]) {
    io.printf(
        "  %s|pad=28|align=left %s\n",
        bioseq.codonTableName($id),
        bioseq.codonTable($id)["TGA"]);
}

# Ambiguity resolves when it can: every codon CGN admits is arginine, so the
# answer is not ambiguous at all. NNN genuinely is.
io.printf("\nambiguous codons\n");
for (def codon in ["CGN", "GGN", "TTY", "ATN", "NNN"]) {
    io.printf("  %s -> %s\n", $codon, bioseq.translate($codon));
}

# ------------------------------------------------------------------- six frames

# When a sequence's frame is unknown - a raw read, a contig, an unannotated
# region - you read all six and look for the one without stops.
io.printf("\nsix frames\n\n");

for (def f in bioseq.sixFrame($gene)) {
    def stops as int init 0;
    for (def aa in strings.chars($f.protein)) {
        if ($aa == "*") {
            $stops = $stops + 1;
        }
    }
    def label as string init "stops";
    if ($stops == 1) {
        $label = "stop";
    }
    def note as string init "";
    if ($stops == 0) {
        $note = "  <- open";
    }
    io.printf(
        "  %d|sign=always|pad=2 %s|pad=26|align=left %d %s|pad=5|align=left%s\n",
        $f.frame,
        $f.protein,
        $stops,
        $label,
        $note);
}

io.printf(
    "\n  frame 1 is the plain translation: %t\n",
    bioseq.sixFrame($gene)[0].protein == $protein);
io.printf(
    "  frame -1 is frame 1 of the reverse complement: %t\n",
    bioseq.translateFrame($gene, -1) == bioseq.translate(bioseq.reverseComplement($gene)));
