#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# Generate src/topics/tabledata.j.
#
#   jennifer run tools/gentables.j
#
# The alternative genetic codes are expressed here as the STANDARD table plus
# the documented differences, which is how NCBI itself describes them and how a
# reviewer can actually check them. The emitted .j holds full literal maps (the
# benchmark says a const read is free, a construction is not); tables_test.j
# re-derives the differences and fails if the two ever disagree.
#
# Numeric values are carried as strings so they reach the output exactly as
# written here, with no formatting round-trip in between. Run from the
# repository root.
#
# Emitted text uses raw '...' strings. A generator writes Jennifer source, so
# almost every line it emits contains double quotes and braces - and in a raw
# string both are literal, with no escape processing and no `{expr}` slot to
# work around. Only the lines carrying a real escape (a newline) or an
# apostrophe of their own stay cooked, since a raw string ends at the first `'`.

use io;
use fs;
use strings;
use lists;
use maps;
use convert;

# --------------------------------------------------------------------- data

def const BASES as string init 'TCAG';

# NCBI transl_table=1, in the canonical TTT..GGG layout.
def const STANDARD_AAS_HEAD as string init 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRR';
def const STANDARD_AAS_TAIL as string init 'IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';

def struct CodonTable {
    id as int,
    constName as string,
    name as string,
    diffs as string,
    starts as string
};

def tables as list of CodonTable init [
    CodonTable{id: 1, constName: 'STANDARD', name: 'Standard', diffs: '', starts: 'TTG,CTG,ATG'},
    CodonTable{
        id: 2,
        constName: 'VERTEBRATE_MITOCHONDRIAL',
        name: 'Vertebrate Mitochondrial',
        diffs: 'AGA=*,AGG=*,ATA=M,TGA=W',
        starts: 'ATT,ATC,ATA,ATG,GTG'
    },
    CodonTable{
        id: 3,
        constName: 'YEAST_MITOCHONDRIAL',
        name: 'Yeast Mitochondrial',
        diffs: 'ATA=M,CTT=T,CTC=T,CTA=T,CTG=T,TGA=W',
        starts: 'ATA,ATG,GTG'
    },
    CodonTable{
        id: 4,
        constName: 'MOLD_MITOCHONDRIAL',
        name: 'Mold, Protozoan and Coelenterate Mitochondrial; Mycoplasma/Spiroplasma',
        diffs: 'TGA=W',
        starts: 'TTA,TTG,CTG,ATT,ATC,ATA,ATG,GTG'
    },
    CodonTable{
        id: 5,
        constName: 'INVERTEBRATE_MITOCHONDRIAL',
        name: 'Invertebrate Mitochondrial',
        diffs: 'AGA=S,AGG=S,ATA=M,TGA=W',
        starts: 'TTG,ATT,ATC,ATA,ATG,GTG'
    },
    CodonTable{
        id: 6,
        constName: 'CILIATE_NUCLEAR',
        name: 'Ciliate, Dasycladacean and Hexamita Nuclear',
        diffs: 'TAA=Q,TAG=Q',
        starts: 'ATG'
    },
    CodonTable{
        id: 9,
        constName: 'ECHINODERM_MITOCHONDRIAL',
        name: 'Echinoderm and Flatworm Mitochondrial',
        diffs: 'AAA=N,AGA=S,AGG=S,TGA=W',
        starts: 'ATG,GTG'
    },
    CodonTable{
        id: 10,
        constName: 'EUPLOTID_NUCLEAR',
        name: 'Euplotid Nuclear',
        diffs: 'TGA=C',
        starts: 'ATG'
    },
    CodonTable{
        id: 11,
        constName: 'BACTERIAL',
        name: 'Bacterial, Archaeal and Plant Plastid',
        diffs: '',
        starts: 'TTG,CTG,ATT,ATC,ATA,ATG,GTG'
    }
];

# one -> three, name, monoisotopic residue mass, average residue mass,
# side-chain pKa or 0.0, charge sign, Kyte-Doolittle hydropathy.
# Masses: ExPASy residue masses (free amino acid less water).
# pKa: the EMBOSS `iep` set (Epk.dat), named because pI differs between sets.
def struct Residue {
    one as string,
    three as string,
    name as string,
    mono as string,
    avg as string,
    pKa as string,
    charge as int,
    hydropathy as string
};

def residues as list of Residue init [
    Residue{
        one: 'A',
        three: 'Ala',
        name: 'Alanine',
        mono: '71.03711',
        avg: '71.0788',
        pKa: '0.0',
        charge: 0,
        hydropathy: '1.8'
    },
    Residue{
        one: 'R',
        three: 'Arg',
        name: 'Arginine',
        mono: '156.10111',
        avg: '156.1875',
        pKa: '12.5',
        charge: 1,
        hydropathy: '-4.5'
    },
    Residue{
        one: 'N',
        three: 'Asn',
        name: 'Asparagine',
        mono: '114.04293',
        avg: '114.1038',
        pKa: '0.0',
        charge: 0,
        hydropathy: '-3.5'
    },
    Residue{
        one: 'D',
        three: 'Asp',
        name: 'Aspartic acid',
        mono: '115.02694',
        avg: '115.0886',
        pKa: '3.9',
        charge: -1,
        hydropathy: '-3.5'
    },
    Residue{
        one: 'C',
        three: 'Cys',
        name: 'Cysteine',
        mono: '103.00919',
        avg: '103.1388',
        pKa: '8.5',
        charge: -1,
        hydropathy: '2.5'
    },
    Residue{
        one: 'E',
        three: 'Glu',
        name: 'Glutamic acid',
        mono: '129.04259',
        avg: '129.1155',
        pKa: '4.1',
        charge: -1,
        hydropathy: '-3.5'
    },
    Residue{
        one: 'Q',
        three: 'Gln',
        name: 'Glutamine',
        mono: '128.05858',
        avg: '128.1307',
        pKa: '0.0',
        charge: 0,
        hydropathy: '-3.5'
    },
    Residue{
        one: 'G',
        three: 'Gly',
        name: 'Glycine',
        mono: '57.02146',
        avg: '57.0519',
        pKa: '0.0',
        charge: 0,
        hydropathy: '-0.4'
    },
    Residue{
        one: 'H',
        three: 'His',
        name: 'Histidine',
        mono: '137.05891',
        avg: '137.1411',
        pKa: '6.5',
        charge: 1,
        hydropathy: '-3.2'
    },
    Residue{
        one: 'I',
        three: 'Ile',
        name: 'Isoleucine',
        mono: '113.08406',
        avg: '113.1594',
        pKa: '0.0',
        charge: 0,
        hydropathy: '4.5'
    },
    Residue{
        one: 'L',
        three: 'Leu',
        name: 'Leucine',
        mono: '113.08406',
        avg: '113.1594',
        pKa: '0.0',
        charge: 0,
        hydropathy: '3.8'
    },
    Residue{
        one: 'K',
        three: 'Lys',
        name: 'Lysine',
        mono: '128.09496',
        avg: '128.1741',
        pKa: '10.8',
        charge: 1,
        hydropathy: '-3.9'
    },
    Residue{
        one: 'M',
        three: 'Met',
        name: 'Methionine',
        mono: '131.04049',
        avg: '131.1926',
        pKa: '0.0',
        charge: 0,
        hydropathy: '1.9'
    },
    Residue{
        one: 'F',
        three: 'Phe',
        name: 'Phenylalanine',
        mono: '147.06841',
        avg: '147.1766',
        pKa: '0.0',
        charge: 0,
        hydropathy: '2.8'
    },
    Residue{
        one: 'P',
        three: 'Pro',
        name: 'Proline',
        mono: '97.05276',
        avg: '97.1167',
        pKa: '0.0',
        charge: 0,
        hydropathy: '-1.6'
    },
    Residue{
        one: 'S',
        three: 'Ser',
        name: 'Serine',
        mono: '87.03203',
        avg: '87.0782',
        pKa: '0.0',
        charge: 0,
        hydropathy: '-0.8'
    },
    Residue{
        one: 'T',
        three: 'Thr',
        name: 'Threonine',
        mono: '101.04768',
        avg: '101.1051',
        pKa: '0.0',
        charge: 0,
        hydropathy: '-0.7'
    },
    Residue{
        one: 'W',
        three: 'Trp',
        name: 'Tryptophan',
        mono: '186.07931',
        avg: '186.2132',
        pKa: '0.0',
        charge: 0,
        hydropathy: '-0.9'
    },
    Residue{
        one: 'Y',
        three: 'Tyr',
        name: 'Tyrosine',
        mono: '163.06333',
        avg: '163.176',
        pKa: '10.1',
        charge: -1,
        hydropathy: '-1.3'
    },
    Residue{
        one: 'V',
        three: 'Val',
        name: 'Valine',
        mono: '99.06841',
        avg: '99.1326',
        pKa: '0.0',
        charge: 0,
        hydropathy: '4.2'
    },
    Residue{
        one: 'U',
        three: 'Sec',
        name: 'Selenocysteine',
        mono: '150.95364',
        avg: '150.0379',
        pKa: '5.2',
        charge: -1,
        hydropathy: '2.5'
    },
    Residue{
        one: 'O',
        three: 'Pyl',
        name: 'Pyrrolysine',
        mono: '237.14773',
        avg: '237.3018',
        pKa: '0.0',
        charge: 0,
        hydropathy: '-3.9'
    }
];

def const IUPAC_PAIRS as string init 'A=A,C=C,G=G,T=T,R=AG,Y=CT,S=CG,' +
    'W=AT,K=GT,M=AC,B=CGT,D=AGT,' + 'H=ACT,V=ACG,N=ACGT';

def const COMPLEMENT_PAIRS as string init 'A=T,C=G,G=C,T=A,R=Y,Y=R,S=S,' +
    'W=W,K=M,M=K,B=V,V=B,D=H,H=D,N=N';

def const THREE_LETTER_PAIRS as string init 'A=Ala,B=Asx,C=Cys,D=Asp,E=Glu,F=Phe,' +
    'G=Gly,H=His,I=Ile,J=Xle,K=Lys,L=Leu,' + 'M=Met,N=Asn,O=Pyl,P=Pro,Q=Gln,R=Arg,' +
    'S=Ser,T=Thr,U=Sec,V=Val,W=Trp,X=Xaa,' + 'Y=Tyr,Z=Glx,*=Ter';

def const PROTEIN_ALPHABET as string init 'ACDEFGHIKLMNPQRSTVWYBZXJUO*';

# Average residue masses: each is the nucleoside monophosphate less one water,
# so a chain's mass is the sum of its residues plus one water.
def const DNA_MASS_PAIRS as string init 'A=313.21,C=289.18,G=329.21,T=304.19';
def const RNA_MASS_PAIRS as string init 'A=329.21,C=305.18,G=345.21,U=306.17';

# SantaLucia (1998) PNAS 95:1460-1465, unified parameters.
# dH in kcal/mol, dS in cal/(mol K), for the ten unique nearest-neighbour pairs.
def const NN_PAIRS as string init 'AA=-7.9;-22.2,AT=-7.2;-20.4,TA=-7.2;-21.3,' +
    'CA=-8.5;-22.7,GT=-8.4;-22.4,CT=-7.8;-21.0,GA=-8.2;-22.2,' +
    'CG=-10.6;-27.2,GC=-9.8;-24.4,GG=-8.0;-19.9';

# NCBI BLOSUM62, the default substitution matrix for protein alignment.
def const BLOSUM_ORDER as string init 'ARNDCQEGHILKMFPSTWYVBZX*';

def blosumRows as list of string init [
    ' 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4',
    '-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4',
    '-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4',
    '-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4',
    ' 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4',
    '-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4',
    '-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4',
    ' 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4',
    '-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4',
    '-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4',
    '-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4',
    '-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4',
    '-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4',
    '-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4',
    '-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4',
    ' 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4',
    ' 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4',
    '-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4',
    '-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4',
    ' 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4',
    '-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4',
    '-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4',
    ' 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4',
    '-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1'
];

# ------------------------------------------------------------------- helpers

func fail(message as string) {
    io.printf("gentables: %s\n", $message);
    exit 1;
}

# "A=T,C=G" -> a map, preserving the order it was written in.
func parsePairs(spec as string) {
    def out as map of string to string;
    if (len($spec) == 0) {
        return $out;
    }
    for (def item in strings.split($spec, ',')) {
        def at as int init strings.indexOf($item, '=');
        if ($at < 0) {
            fail('malformed pair: ' + $item);
        }
        $out[$item[0..$at]] = $item[$at + 1..];
    }
    return $out;
}

func tokens(line as string) {
    def out as list of string;
    for (def part in strings.split($line, ' ')) {
        if (len($part) > 0) {
            $out[] = $part;
        }
    }
    return $out;
}

# The 64 codons in the canonical NCBI layout order.
func allCodons() {
    def out as list of string;
    def letters as list of string init strings.chars(BASES);
    for (def a in $letters) {
        for (def b in $letters) {
            for (def c in $letters) {
                $out[] = $a + $b + $c;
            }
        }
    }
    return $out;
}

func revcomp(pair as string) {
    def pairMap as map of string to string init parsePairs(COMPLEMENT_PAIRS);
    def flipped as string init '';
    def chars as list of string init strings.chars($pair);
    for (def i as int init len($chars) - 1; $i >= 0; $i = $i - 1) {
        $flipped = $flipped + $pairMap[$chars[$i]];
    }
    return $flipped;
}

# A 128-entry list, `fill` everywhere except the codes named in `values`.
func denseList(values as map of int to string, fill as string) {
    def out as list of string;
    for (def i as int init 0; $i < 128; $i = $i + 1) {
        if (maps.has($values, $i)) {
            $out[] = $values[$i];
        } else {
            $out[] = $fill;
        }
    }
    return $out;
}

func codeOf(one as string) {
    return convert.bytesFromString($one, 'utf-8')[0];
}

def emitted as list of string;

func emit(line as string) {
    $emitted[] = $line;
}

func emitBlock(name as string, kind as string, rows as list of string) {
    emit('def const ' + $name + ' as list of ' + $kind + ' init [');
    emit(strings.join($rows, ",\n"));
    emit('];');
    emit('');
}

# -------------------------------------------------------------- derived data

def standardAas as string init STANDARD_AAS_HEAD + STANDARD_AAS_TAIL;
def codons as list of string init allCodons();
def standard as map of string to string;
for (def i as int init 0; $i < len($codons); $i = $i + 1) {
    $standard[$codons[$i]] = $standardAas[$i..$i + 1];
}

if (len($standard) != 64) {
    fail('the standard table does not have 64 codons');
}
if ($standard['ATG'] != 'M' or $standard['TAA'] != '*' or $standard['TGG'] != 'W') {
    fail("the standard table's anchors are wrong");
}

def iupac as map of string to string init parsePairs(IUPAC_PAIRS);
def complements as map of string to string init parsePairs(COMPLEMENT_PAIRS);
def threeLetter as map of string to string init parsePairs(THREE_LETTER_PAIRS);
def dnaMass as map of string to string init parsePairs(DNA_MASS_PAIRS);
def rnaMass as map of string to string init parsePairs(RNA_MASS_PAIRS);

if (len($threeLetter) != len(PROTEIN_ALPHABET)) {
    fail('three-letter codes must cover the protein alphabet exactly');
}
for (def one in strings.chars(PROTEIN_ALPHABET)) {
    if (not maps.has($threeLetter, $one)) {
        fail("no three-letter code for '" + $one + "'");
    }
}

# A duplex reads the same from either strand, so a dinucleotide and its reverse
# complement share one parameter. Expand the ten published values into all sixteen.
def nnUnique as map of string to string;
for (def item in strings.split(NN_PAIRS, ',')) {
    def at as int init strings.indexOf($item, '=');
    $nnUnique[$item[0..$at]] = $item[$at + 1..];
}
def nnAll as map of string to string;
for (def a in strings.chars('ACGT')) {
    for (def b in strings.chars('ACGT')) {
        def pair as string init $a + $b;
        if (maps.has($nnUnique, $pair)) {
            $nnAll[$pair] = $nnUnique[$pair];
        } else {
            def flipped as string init revcomp($pair);
            if (not maps.has($nnUnique, $flipped)) {
                fail('no nearest-neighbour parameter for ' + $pair);
            }
            $nnAll[$pair] = $nnUnique[$flipped];
        }
    }
}
if (len($nnAll) != 16) {
    fail('the nearest-neighbour table is not complete');
}
for (def pair in ['AT', 'TA', 'CG', 'GC']) {
    if (revcomp($pair) != $pair or $nnAll[$pair] != $nnUnique[$pair]) {
        fail('symmetric pair ' + $pair + ' is inconsistent');
    }
}

# BLOSUM62, keyed by the concatenated residue pair.
def blosum as map of string to int;
def blosumLetters as list of string init strings.chars(BLOSUM_ORDER);
if (len($blosumRows) != len($blosumLetters)) {
    fail('BLOSUM62 is not square');
}
for (def i as int init 0; $i < len($blosumLetters); $i = $i + 1) {
    def cells as list of string init tokens($blosumRows[$i]);
    if (len($cells) != len($blosumLetters)) {
        fail('BLOSUM62 row ' + $blosumLetters[$i] + ' has the wrong width');
    }
    for (def j as int init 0; $j < len($blosumLetters); $j = $j + 1) {
        $blosum[$blosumLetters[$i] + $blosumLetters[$j]] = convert.toInt($cells[$j]);
    }
}
# a substitution matrix must be symmetric
for (def a in $blosumLetters) {
    for (def b in $blosumLetters) {
        if ($blosum[$a + $b] != $blosum[$b + $a]) {
            fail('BLOSUM62 ' + $a + $b + ' is not symmetric');
        }
    }
}
# the anchors everyone knows
def anchors as map of string to int init {
    'AA': 4,
    'WW': 11,
    'CC': 9,
    'HH': 8,
    'PP': 7,
    'RK': 2,
    'FY': 3,
    'IV': 3,
    'LM': 2,
    'WC': -2
};
for (def pair in $anchors) {
    if ($blosum[$pair] != $anchors[$pair]) {
        fail('BLOSUM62[' + $pair + '] should be ' + convert.toString($anchors[$pair]));
    }
}
# identities beat the substitutions of that residue
for (def a in strings.chars('ARNDCQEGHILKMFPSTWYV')) {
    for (def b in strings.chars('ARNDCQEGHILKMFPSTWYV')) {
        if ($a != $b and $blosum[$a + $a] <= $blosum[$a + $b]) {
            fail($a + $a + ' should beat ' + $a + $b);
        }
    }
}

# ------------------------------------------------------------------- emission

emit('# SPDX-License-Identifier: LGPL-3.0-only');
emit('# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>');
emit('# pragma-jennifer-version: >=0.25.0');
emit('');
emit('# bioinformatics - table data: the reference data the rest of the deck reads.');
emit('# Spliced into bioinformatics.j via include - not a standalone module.');
emit('#');
emit('# GENERATED by tools/gentables.j. The alternative genetic codes are emitted');
emit('# as full literal maps because a const read costs nothing while a per-call');
emit('# construction does (examples/bench.j). What a reviewer needs to check is not');
emit('# 64 lines of near-identical data but the handful of codons each table');
emit('# reassigns - tables_test.j holds exactly that difference list and fails if a');
emit('# table here stops matching it.');
emit('');

for (def t in $tables) {
    def diffs as map of string to string init parsePairs($t.diffs);
    def table as map of string to string;
    for (def codon in $codons) {
        if (maps.has($diffs, $codon)) {
            $table[$codon] = $diffs[$codon];
        } else {
            $table[$codon] = $standard[$codon];
        }
    }
    emit('/**');
    emit(' * NCBI genetic code ' + convert.toString($t.id) + ' - ' + $t.name + '.');
    emit(' *');
    if (len($diffs) > 0) {
        def parts as list of string;
        for (def codon in $diffs) {
            $parts[] = '`' + $codon + '` -> `' + $diffs[$codon] + '`';
        }
        $parts = lists.sort($parts);
        def plural as string init 's';
        if (len($diffs) == 1) {
            $plural = '';
        }
        def lead as string init 'Differs from `STANDARD` in ' + convert.toString(len($diffs)) +
            ' codon' + $plural + ': ';
        # wrap at the 100-column lint limit, allowing for the " * " prefix
        def line as string init $lead;
        def wrapped as list of string;
        for (def k as int init 0; $k < len($parts); $k = $k + 1) {
            def piece as string init $parts[$k] + ', ';
            if ($k == len($parts) - 1) {
                $piece = $parts[$k] + '.';
            }
            if (len($line) + len($piece) > 93) {
                $wrapped[] = strings.trimRight($line);
                $line = $piece;
            } else {
                $line = $line + $piece;
            }
        }
        $wrapped[] = strings.trimRight($line);
        for (def wl in $wrapped) {
            emit(' * ' + $wl);
        }
    } elseif ($t.id == 1) {
        emit(" * The eukaryotic nuclear code, and the deck's default everywhere a");
        emit(' * genetic code is not named explicitly.');
    } else {
        emit(' * Identical to `STANDARD` in its codon assignments; it differs only in');
        emit(' * which codons may initiate translation - see `startCodons`.');
    }
    if ($t.id == 3) {
        emit(' *');
        emit(' * NCBI additionally marks `CGA` and `CGC` as absent in this code rather');
        emit(' * than as arginine. A map cannot express "absent", so they are kept as');
        emit(' * `R` here; a sequence relying on that distinction needs a hand check.');
    }
    emit(' *');
    emit(' * @see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG' +
        convert.toString($t.id));
    emit(' */');
    emit('export def const ' + $t.constName + ' as map of string to string init {');
    def rows as list of string;
    for (def codon in $codons) {
        $rows[] = '    "' + $codon + '": "' + $table[$codon] + '"';
    }
    emit(strings.join($rows, ",\n"));
    emit('};');
    emit('');
}

emit('# Codon table id -> the table itself, its NCBI name, and its start codons.');
emit('# Kept as three parallel maps rather than a struct so each stays a plain');
emit('# const the dispatchers below can read without building anything.');
emit('');
emit('def const TABLE_NAMES as map of int to string init {');
def nameRows as list of string;
for (def t in $tables) {
    $nameRows[] = '    ' + convert.toString($t.id) + ': "' + $t.name + '"';
}
emit(strings.join($nameRows, ",\n"));
emit('};');
emit('');
emit('def const TABLE_STARTS as map of int to string init {');
def startRows as list of string;
for (def t in $tables) {
    $startRows[] = '    ' + convert.toString($t.id) + ': "' + $t.starts + '"';
}
emit(strings.join($startRows, ",\n"));
emit('};');
emit('');

emit('/**');
emit(' * IUPAC nucleotide ambiguity codes for DNA: code -> the bases it stands for.');
emit(' *');
emit(' * The four unambiguous bases map to themselves, so a lookup never needs a');
emit(' * special case. `N` is the fully ambiguous code.');
emit(' *');
emit(' * @see https://www.bioinformatics.org/sms2/iupac.html');
emit(' */');
emit('export def const IUPAC_DNA as map of string to string init {');
def iupacRows as list of string;
for (def code in $iupac) {
    $iupacRows[] = '    "' + $code + '": "' + $iupac[$code] + '"';
}
emit(strings.join($iupacRows, ",\n"));
emit('};');
emit('');
emit('/**');
emit(' * IUPAC nucleotide ambiguity codes for RNA - `IUPAC_DNA` with `U` for `T`.');
emit(' */');
emit('export def const IUPAC_RNA as map of string to string init {');
def rnaRows as list of string;
for (def code in $iupac) {
    if ($code == 'T') {
        continue;
    }
    $rnaRows[] = '    "' + $code + '": "' + strings.replace($iupac[$code], 'T', 'U') + '"';
}
$rnaRows[] = '    "U": "U"';
emit(strings.join($rnaRows, ",\n"));
emit('};');
emit('');

emit('/**');
emit(' * Watson-Crick complement of each DNA code, ambiguity codes included.');
emit(' *');
emit(' * The ambiguity codes complement as sets: `R` (A or G) complements to `Y`');
emit(' * (C or T), `B` (not A) to `V` (not T), and the self-complementary `S`, `W`');
emit(' * and `N` map to themselves.');
emit(' */');
emit('export def const COMPLEMENT_DNA as map of string to string init {');
def compRows as list of string;
for (def code in $complements) {
    $compRows[] = '    "' + $code + '": "' + $complements[$code] + '"';
}
emit(strings.join($compRows, ",\n"));
emit('};');
emit('');
emit('/**');
emit(' * Watson-Crick complement of each RNA code - `COMPLEMENT_DNA` with `U`');
emit(' * standing in for `T` on both sides.');
emit(' */');
emit('export def const COMPLEMENT_RNA as map of string to string init {');
def rnaComp as map of string to string;
for (def code in $complements) {
    def from as string init $code;
    if ($from == 'T') {
        $from = 'U';
    }
    $rnaComp[$from] = strings.replace($complements[$code], 'T', 'U');
}
def rnaCompRows as list of string;
for (def code in $rnaComp) {
    $rnaCompRows[] = '    "' + $code + '": "' + $rnaComp[$code] + '"';
}
emit(strings.join($rnaCompRows, ",\n"));
emit('};');
emit('');

emit('/**');
emit(" * One residue's identity, mass and physico-chemical properties.");
emit(' *');
emit(' * @field one {string} the one-letter code, e.g. "K"');
emit(' * @field three {string} the three-letter code, e.g. "Lys"');
emit(' * @field name {string} the full name, e.g. "Lysine"');
emit(' * @field monoMass {float} monoisotopic residue mass in daltons (the free');
emit(' *        amino acid less one water)');
emit(' * @field avgMass {float} average residue mass in daltons');
emit(' * @field pKa {float} side-chain dissociation constant, or 0.0 when the side');
emit(' *        chain does not ionise - the EMBOSS `iep` parameter set');
emit(' * @field charge {int} +1 basic, -1 acidic, 0 neutral side chain');
emit(' * @field hydropathy {float} Kyte-Doolittle hydropathy index, the input to');
emit(' *        `gravy`');
emit(' */');
emit('export def struct AminoAcid {');
emit('    one as string,');
emit('    three as string,');
emit('    name as string,');
emit('    monoMass as float,');
emit('    avgMass as float,');
emit('    pKa as float,');
emit('    charge as int,');
emit('    hydropathy as float');
emit('};');
emit('');
emit('/**');
emit(' * pKa of the free alpha-amino group at the N terminus (EMBOSS `iep`).');
emit(' */');
emit('export def const PKA_N_TERMINUS as float init 8.6;');
emit('');
emit('/**');
emit(' * pKa of the free alpha-carboxyl group at the C terminus (EMBOSS `iep`).');
emit(' */');
emit('export def const PKA_C_TERMINUS as float init 3.6;');
emit('');
emit('/**');
emit(' * Monoisotopic mass of water, added once when residue masses are summed into');
emit(' * a peptide mass.');
emit(' */');
emit('export def const WATER_MONO_MASS as float init 18.010565;');
emit('');
emit('/**');
emit(' * Average mass of water, the counterpart to `WATER_MONO_MASS`.');
emit(' */');
emit('export def const WATER_AVG_MASS as float init 18.01528;');
emit('');
emit('# The residue table, keyed by one-letter code. `U` (selenocysteine) and `O`');
emit('# (pyrrolysine) follow the standard twenty; the ambiguity codes B, Z, X and J');
emit('# have no single mass or pKa and are deliberately absent, so a caller that');
emit('# needs them has to decide what they mean rather than get a silent wrong answer.');
emit('#');
emit('# Masses: ExPASy residue masses. Hydropathy: Kyte & Doolittle (1982),');
emit('# J Mol Biol 157:105-132. pKa: the EMBOSS `iep` set (Epk.dat).');
emit('def const RESIDUES as map of string to AminoAcid init {');
def residueRows as list of string;
for (def r in $residues) {
    $residueRows[] = '    "' + $r.one + '": AminoAcid{one: "' + $r.one + '", three: "' +
        $r.three + '", name: "' + $r.name + '", monoMass: ' + $r.mono + ', avgMass: ' + $r.avg +
        ', pKa: ' + $r.pKa + ', charge: ' + convert.toString($r.charge) + ', hydropathy: ' +
        $r.hydropathy + '}';
}
emit(strings.join($residueRows, ",\n"));
emit('};');
emit('');

emit('# Byte-keyed complements, for the per-residue loops in transforms.j. Same');
emit('# content as COMPLEMENT_DNA / COMPLEMENT_RNA above, keyed by ASCII code so a');
emit('# loop can complement a byte without materialising a one-character string -');
emit('# the convention core.j sets out. tables_test.j checks the two agree.');
emit('def const COMPLEMENT_DNA_CODES as map of int to int init {');
def compCodeRows as list of string;
for (def code in $complements) {
    $compCodeRows[] = '    ' + convert.toString(codeOf($code)) + ': ' +
        convert.toString(codeOf($complements[$code]));
}
emit(strings.join($compCodeRows, ",\n"));
emit('};');
emit('');
emit('def const COMPLEMENT_RNA_CODES as map of int to int init {');
def rnaCodeRows as list of string;
for (def code in $rnaComp) {
    $rnaCodeRows[] = '    ' + convert.toString(codeOf($code)) + ': ' +
        convert.toString(codeOf($rnaComp[$code]));
}
emit(strings.join($rnaCodeRows, ",\n"));
emit('};');
emit('');

emit('/**');
emit(' * One-letter residue code -> its three-letter code.');
emit(' *');
emit(' * Covers the whole protein alphabet, unlike `RESIDUES`: the ambiguity codes');
emit(' * (`Asx`, `Glx`, `Xaa`, `Xle`) and the stop symbol (`Ter`) have a spelling');
emit(' * even though they have no mass or pKa.');
emit(' */');
emit('export def const THREE_LETTER as map of string to string init {');
def tlRows as list of string;
for (def one in $threeLetter) {
    $tlRows[] = '    "' + $one + '": "' + $threeLetter[$one] + '"';
}
emit(strings.join($tlRows, ",\n"));
emit('};');
emit('');
emit('# Three-letter code (uppercased) -> one-letter code, the inverse of');
emit('# THREE_LETTER, for threeToOne.');
emit('def const ONE_LETTER as map of string to string init {');
def olRows as list of string;
for (def one in $threeLetter) {
    $olRows[] = '    "' + strings.upper($threeLetter[$one]) + '": "' + $one + '"';
}
emit(strings.join($olRows, ",\n"));
emit('};');
emit('');
emit('# Three-letter code (uppercased) -> one-letter code, so `aminoAcid` can take');
emit('# either spelling.');
emit('def const THREE_TO_ONE as map of string to string init {');
def ttoRows as list of string;
for (def r in $residues) {
    $ttoRows[] = '    "' + strings.upper($r.three) + '": "' + $r.one + '"';
}
emit(strings.join($ttoRows, ",\n"));
emit('};');
emit('');

emit('/**');
emit(' * Average residue mass of each DNA base, in daltons.');
emit(' *');
emit(' * Each value is the deoxyribonucleoside monophosphate less one water, so a');
emit(" * strand's mass is the sum of its residues plus one water. These are the");
emit(' * coefficients of the standard single-stranded oligo formula.');
emit(' */');
emit('export def const DNA_RESIDUE_MASS as map of string to float init {');
def dnaMassRows as list of string;
for (def base in $dnaMass) {
    $dnaMassRows[] = '    "' + $base + '": ' + $dnaMass[$base];
}
emit(strings.join($dnaMassRows, ",\n"));
emit('};');
emit('');
emit('/**');
emit(' * Average residue mass of each RNA base, in daltons - the ribonucleoside');
emit(' * counterpart of `DNA_RESIDUE_MASS`.');
emit(' */');
emit('export def const RNA_RESIDUE_MASS as map of string to float init {');
def rnaMassRows as list of string;
for (def base in $rnaMass) {
    $rnaMassRows[] = '    "' + $base + '": ' + $rnaMass[$base];
}
emit(strings.join($rnaMassRows, ",\n"));
emit('};');
emit('');
emit('/**');
emit(' * Average mass of a metaphosphate group (HPO3), subtracted to report an');
emit(" * oligo with a free 5' hydroxyl rather than a 5' phosphate.");
emit(' */');
emit('export def const PHOSPHATE_MASS as float init 79.98;');
emit('');

emit('/**');
emit(' * Nearest-neighbour enthalpy dH, in kcal/mol, for every dinucleotide.');
emit(' *');
emit(' * SantaLucia (1998) unified parameters. A duplex reads the same from either');
emit(' * strand, so a dinucleotide and its reverse complement share one parameter;');
emit(' * the ten published values are expanded to all sixteen here.');
emit(' *');
emit(' * @see SantaLucia J (1998) PNAS 95:1460-1465');
emit(' */');
emit('export def const NN_ENTHALPY as map of string to float init {');
def nnKeys as list of string;
for (def pair in $nnAll) {
    $nnKeys[] = $pair;
}
$nnKeys = lists.sort($nnKeys);
def enthalpyRows as list of string;
def entropyRows as list of string;
for (def pair in $nnKeys) {
    def parts as list of string init strings.split($nnAll[$pair], ';');
    $enthalpyRows[] = '    "' + $pair + '": ' + $parts[0];
    $entropyRows[] = '    "' + $pair + '": ' + $parts[1];
}
emit(strings.join($enthalpyRows, ",\n"));
emit('};');
emit('');
emit('/**');
emit(' * Nearest-neighbour entropy dS, in cal/(mol K), the companion to');
emit(' * `NN_ENTHALPY`.');
emit(' */');
emit('export def const NN_ENTROPY as map of string to float init {');
emit(strings.join($entropyRows, ",\n"));
emit('};');
emit('');
emit('# Helix-initiation corrections, applied once per duplex end, and the symmetry');
emit('# correction for a self-complementary strand (SantaLucia 1998).');
emit('def const NN_INIT_GC_ENTHALPY as float init 0.1;');
emit('def const NN_INIT_GC_ENTROPY as float init -2.8;');
emit('def const NN_INIT_AT_ENTHALPY as float init 2.3;');
emit('def const NN_INIT_AT_ENTROPY as float init 4.1;');
emit('def const NN_SYMMETRY_ENTROPY as float init -1.4;');
emit('');
emit('/**');
emit(' * The gas constant in cal/(mol K), as used by the nearest-neighbour model.');
emit(' */');
emit('export def const GAS_CONSTANT as float init 1.987;');
emit('');

emit('/**');
emit(' * BLOSUM62, the default substitution matrix for protein alignment.');
emit(' *');
emit(' * Keyed by the concatenated residue pair, so `BLOSUM62["WW"]` is 11 and');
emit(' * `BLOSUM62["AR"]` is -1. The matrix is symmetric, so only one order of');
emit(' * each pair is stored - read it through `substitutionScore`, which tries');
emit(' * both, rather than indexing it directly.');
emit(' *');
emit(' * Covers the twenty standard residues plus the ambiguity codes `B`, `Z` and');
emit(' * `X` and the stop symbol `*`, because an alignment has to score whatever a');
emit(' * translation produced.');
emit(' *');
emit(' * @see Henikoff S, Henikoff JG (1992) PNAS 89:10915-10919');
emit(' */');
emit('export def const BLOSUM62 as map of string to int init {');
def blosumRowsOut as list of string;
for (def i as int init 0; $i < len($blosumLetters); $i = $i + 1) {
    for (def j as int init $i; $j < len($blosumLetters); $j = $j + 1) {
        def pair as string init $blosumLetters[$i] + $blosumLetters[$j];
        $blosumRowsOut[] = '    "' + $pair + '": ' + convert.toString($blosum[$pair]);
    }
}
emit(strings.join($blosumRowsOut, ",\n"));
emit('};');
emit('');

# ------------------------------------------------- dense by-code lookup tables

emit('# Membership by byte value, for the validation loops in core.j. Same content');
emit('# as the *_CODES maps above; a list index is cheaper than a map lookup and');
emit('# these run over every residue of every input the deck is given.');

def alphabets as map of string to string init {
    'DNA_ADMITS': 'ACGTRYSWKMBDHVN',
    'RNA_ADMITS': 'ACGURYSWKMBDHVN',
    'PROTEIN_ADMITS': PROTEIN_ALPHABET,
    'NUCLEOTIDE_ADMITS': 'ACGTURYSWKMBDHVN'
};
for (def name in $alphabets) {
    def values as map of int to string;
    for (def one in strings.chars($alphabets[$name])) {
        $values[codeOf($one)] = 'true';
    }
    emitBlock($name, 'bool', denseList($values, 'false'));
}

emit('# Residue properties by byte value, for gravy / molecularWeight / netCharge.');
emit('# A residue with no entry reads as the zero value, so callers must still ask');
emit('# hasAminoAcid first - these are the fast path, not the validation.');

def hydro as map of int to string;
def avgMassByCode as map of int to string;
def pkaByCode as map of int to string;
def chargeByCode as map of int to string;
def knownByCode as map of int to string;
for (def r in $residues) {
    def code as int init codeOf($r.one);
    $hydro[$code] = $r.hydropathy;
    $avgMassByCode[$code] = $r.avg;
    $pkaByCode[$code] = $r.pKa;
    $chargeByCode[$code] = convert.toString($r.charge);
    $knownByCode[$code] = 'true';
}
emitBlock('HYDROPATHY_BY_CODE', 'float', denseList($hydro, '0.0'));
emitBlock('AVG_MASS_BY_CODE', 'float', denseList($avgMassByCode, '0.0'));
emitBlock('PKA_BY_CODE', 'float', denseList($pkaByCode, '0.0'));
emitBlock('CHARGE_BY_CODE', 'int', denseList($chargeByCode, '0'));
emitBlock('RESIDUE_KNOWN_BY_CODE', 'bool', denseList($knownByCode, 'false'));

def dnaMassByCode as map of int to string;
for (def base in $dnaMass) {
    $dnaMassByCode[codeOf($base)] = $dnaMass[$base];
}
emitBlock('DNA_MASS_BY_CODE', 'float', denseList($dnaMassByCode, '-1.0'));
def rnaMassByCode as map of int to string;
for (def base in $rnaMass) {
    $rnaMassByCode[codeOf($base)] = $rnaMass[$base];
}
emitBlock('RNA_MASS_BY_CODE', 'float', denseList($rnaMassByCode, '-1.0'));

def dnaCompByCode as map of int to string;
for (def code in $complements) {
    $dnaCompByCode[codeOf($code)] = convert.toString(codeOf($complements[$code]));
}
emitBlock('COMPLEMENT_DNA_BY_CODE', 'int', denseList($dnaCompByCode, '-1'));
def rnaCompByCode as map of int to string;
for (def code in $rnaComp) {
    $rnaCompByCode[codeOf($code)] = convert.toString(codeOf($rnaComp[$code]));
}
emitBlock('COMPLEMENT_RNA_BY_CODE', 'int', denseList($rnaCompByCode, '-1'));

fs.writeString('src/topics/tabledata.j', strings.join($emitted, "\n") + "\n");
io.printf(
    "wrote src/topics/tabledata.j - %d genetic codes, %d residues\n",
    len($tables),
    len($residues));
