#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - the reference data the deck carries.
#
#   jennifer run examples/tables.j
#
# Genetic codes, IUPAC ambiguity codes, residue properties and restriction
# enzymes. Needs no network, file or argument.

use io;
use strings;

import "../src/bioinformatics.j" as bioseq;

# ------------------------------------------------------------- genetic codes

io.printf("genetic codes\n\n");

for (def id in bioseq.codonTableIds()) {
    def starts as list of string init bioseq.startCodons($id);
    def label as string init "starts";
    if (len($starts) == 1) {
        $label = "start";
    }
    io.printf(
        "  %d|pad=3  %s|pad=44|align=left %d %s\n",
        $id,
        bioseq.codonTableName($id),
        len($starts),
        $label);
}

# What actually separates one code from another is a handful of codons. The
# ciliate code is the famous case: two of the three stops read glutamine, which
# is why a naive translation of a Paramecium gene stops far too early.
io.printf("\n  TAA under the standard code    %s\n", bioseq.STANDARD["TAA"]);
io.printf("  TAA under the ciliate code    %s\n", bioseq.CILIATE_NUCLEAR["TAA"]);
io.printf("  TGA under the standard code   %s\n", bioseq.STANDARD["TGA"]);
io.printf("  TGA in vertebrate mitochondria %s\n", bioseq.VERTEBRATE_MITOCHONDRIAL["TGA"]);

# Bacterial reads identically to standard - it differs only in what may start a
# gene, which is why an ORF search takes a table id rather than a codon map.
io.printf("\n  standard starts   %s\n", strings.join(bioseq.startCodons(1), " "));
io.printf("  bacterial starts  %s\n", strings.join(bioseq.startCodons(11), " "));

# ---------------------------------------------------------------- IUPAC codes

io.printf("\nIUPAC ambiguity codes\n\n");

for (def code in bioseq.IUPAC_DNA) {
    def expansion as string init bioseq.expandIupac($code);
    if (len($expansion) > 1) {
        io.printf(
            "  %s  %s|pad=6|align=left as a regex class %s\n",
            $code,
            $expansion,
            bioseq.iupacClass($code));
    }
}

io.printf("\n  R admits G?  %t\n", bioseq.matchesIupac("R", "G"));
io.printf("  R admits C?  %t\n", bioseq.matchesIupac("R", "C"));

# ------------------------------------------------------------- amino acids

io.printf("\nresidue properties\n\n");
io.printf(
    "  %s|pad=4|align=left %s|pad=16|align=left %s|pad=11 %s|pad=8 %s\n",
    "aa",
    "name",
    "mono mass",
    "pKa",
    "hydropathy");

for (def aa in bioseq.aminoAcids()) {
    io.printf(
        "  %s %s|pad=3|align=left %s|pad=16|align=left %f|prec=5|pad=11 ",
        $aa.one,
        $aa.three,
        $aa.name,
        $aa.monoMass);
    if ($aa.pKa > 0.0) {
        io.printf("%f|prec=1|pad=8", $aa.pKa);
    } else {
        io.printf("%s|pad=8|align=right", "-");
    }
    io.printf(" %f|prec=1\n", $aa.hydropathy);
}

# The two most hydrophobic and most hydrophilic residues, straight off the
# Kyte-Doolittle scale that `gravy` will sum.
io.printf("\n  most hydrophobic  Ile %f|prec=1\n", bioseq.aminoAcid("I").hydropathy);
io.printf("  most hydrophilic  Arg %f|prec=1\n", bioseq.aminoAcid("R").hydropathy);

# Leucine and isoleucine weigh exactly the same - the reason mass spectrometry
# reports them together as "Leu/Ile".
io.printf(
    "  Leu and Ile are isobaric: %t\n",
    bioseq.aminoAcid("L").monoMass == bioseq.aminoAcid("I").monoMass);

# --------------------------------------------------------- restriction enzymes

io.printf("\nrestriction enzymes (%d in the table)\n\n", len(bioseq.enzymes()));

def shown as list of string init ["EcoRI", "PstI", "SmaI", "NotI", "EcoRII", "HinfI"];
for (def name in $shown) {
    def e as bioseq.Enzyme init bioseq.enzyme($name);
    io.printf(
        "  %s|pad=9|align=left %s|pad=12|align=left cut %d/%d  %s|pad=6|align=left",
        $e.name,
        $e.site,
        $e.cutTop,
        $e.cutBottom,
        $e.overhang);
    if ($e.overhangLength > 0) {
        io.printf(" %d nt", $e.overhangLength);
    }
    io.printf("\n");
}

# Same site, different cut: SmaI leaves blunt ends, XmaI a four-base overhang.
# Both are reported as isoschizomers, and cutTop is what tells them apart.
io.printf("\n  SmaI isoschizomers  %s\n", strings.join(bioseq.isoschizomers("SmaI"), ", "));
io.printf("  SmaI CCC^GGG -> %s\n", bioseq.enzyme("SmaI").overhang);
io.printf(
    "  XmaI C^CCGGG -> %s, %d nt\n",
    bioseq.enzyme("XmaI").overhang,
    bioseq.enzyme("XmaI").overhangLength);

# The eight-cutters are the ones that rarely cut at all, which is what makes
# them useful for linearising a plasmid exactly once.
io.printf("\n  eight-base recognition sites:");
for (def e in bioseq.enzymes()) {
    if (len($e.site) == 8) {
        io.printf(" %s", $e.name);
    }
}
io.printf("\n");
