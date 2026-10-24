#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# Generate src/topics/enzymedata.j.
#
#   jennifer run tools/genenzymes.j
#
# Recognition sequences and cut positions are hand-curated from published enzyme
# specifications rather than imported in bulk from REBASE, whose data carries its
# own usage terms (see docs/enzymes.md). Cut offsets are measured from the first
# base of the recognition site: `cutTop` is where the top strand is nicked,
# `cutBottom` where the bottom strand is, both counted along the top strand.
#
#     EcoRI  G^AATTC   ->  site "GAATTC", cutTop 1, cutBottom 5, 5' overhang AATT
#     PstI   CTGCA^G   ->  site "CTGCAG", cutTop 5, cutBottom 1, 3' overhang TGCA
#     SmaI   CCC^GGG   ->  site "CCCGGG", cutTop 3, cutBottom 3, blunt
#
# The generator is where the review happens: it asserts every invariant it can
# before writing anything, so a malformed row fails here rather than reaching
# the deck. Run from the repository root.
#
# Emitted text uses raw '...' strings. A generator writes Jennifer source, so
# almost every line it emits contains double quotes and braces - and in a raw
# string both are literal, with no escape processing and no `{expr}` slot to
# work around. Only the lines carrying a real escape (a newline) or an
# apostrophe of their own stay cooked, since a raw string ends at the first `'`.

use io;
use fs;
use strings;
use maps;
use convert;

def struct Enzyme {
    name as string,
    site as string,
    cutTop as int,
    cutBottom as int
};

def enzymes as list of Enzyme init [
    Enzyme{name: 'AatII', site: 'GACGTC', cutTop: 5, cutBottom: 1},
    Enzyme{name: 'AflII', site: 'CTTAAG', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'AflIII', site: 'ACRYGT', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'AgeI', site: 'ACCGGT', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'AluI', site: 'AGCT', cutTop: 2, cutBottom: 2},
    Enzyme{name: 'ApaI', site: 'GGGCCC', cutTop: 5, cutBottom: 1},
    Enzyme{name: 'AscI', site: 'GGCGCGCC', cutTop: 2, cutBottom: 6},
    Enzyme{name: 'AvrII', site: 'CCTAGG', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'BamHI', site: 'GGATCC', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'BclI', site: 'TGATCA', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'BglII', site: 'AGATCT', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'BsaAI', site: 'YACGTR', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'BspEI', site: 'TCCGGA', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'BspHI', site: 'TCATGA', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'BsrGI', site: 'TGTACA', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'BstBI', site: 'TTCGAA', cutTop: 2, cutBottom: 4},
    Enzyme{name: 'BstEII', site: 'GGTNACC', cutTop: 1, cutBottom: 6},
    Enzyme{name: 'ClaI', site: 'ATCGAT', cutTop: 2, cutBottom: 4},
    Enzyme{name: 'DdeI', site: 'CTNAG', cutTop: 1, cutBottom: 4},
    Enzyme{name: 'DraI', site: 'TTTAAA', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'EagI', site: 'CGGCCG', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'Eco47III', site: 'AGCGCT', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'EcoRI', site: 'GAATTC', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'EcoRII', site: 'CCWGG', cutTop: 0, cutBottom: 5},
    Enzyme{name: 'EcoRV', site: 'GATATC', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'FseI', site: 'GGCCGGCC', cutTop: 6, cutBottom: 2},
    Enzyme{name: 'HaeIII', site: 'GGCC', cutTop: 2, cutBottom: 2},
    Enzyme{name: 'HhaI', site: 'GCGC', cutTop: 3, cutBottom: 1},
    Enzyme{name: 'HindIII', site: 'AAGCTT', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'HinfI', site: 'GANTC', cutTop: 1, cutBottom: 4},
    Enzyme{name: 'HpaI', site: 'GTTAAC', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'HpaII', site: 'CCGG', cutTop: 1, cutBottom: 3},
    Enzyme{name: 'KasI', site: 'GGCGCC', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'KpnI', site: 'GGTACC', cutTop: 5, cutBottom: 1},
    Enzyme{name: 'MboI', site: 'GATC', cutTop: 0, cutBottom: 4},
    Enzyme{name: 'MfeI', site: 'CAATTG', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'MluI', site: 'ACGCGT', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'MscI', site: 'TGGCCA', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'MspI', site: 'CCGG', cutTop: 1, cutBottom: 3},
    Enzyme{name: 'NaeI', site: 'GCCGGC', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'NarI', site: 'GGCGCC', cutTop: 2, cutBottom: 4},
    Enzyme{name: 'NcoI', site: 'CCATGG', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'NdeI', site: 'CATATG', cutTop: 2, cutBottom: 4},
    Enzyme{name: 'NheI', site: 'GCTAGC', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'NotI', site: 'GCGGCCGC', cutTop: 2, cutBottom: 6},
    Enzyme{name: 'NruI', site: 'TCGCGA', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'NsiI', site: 'ATGCAT', cutTop: 5, cutBottom: 1},
    Enzyme{name: 'NspI', site: 'RCATGY', cutTop: 5, cutBottom: 1},
    Enzyme{name: 'PacI', site: 'TTAATTAA', cutTop: 5, cutBottom: 3},
    Enzyme{name: 'PciI', site: 'ACATGT', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'PmeI', site: 'GTTTAAAC', cutTop: 4, cutBottom: 4},
    Enzyme{name: 'PmlI', site: 'CACGTG', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'PstI', site: 'CTGCAG', cutTop: 5, cutBottom: 1},
    Enzyme{name: 'PvuI', site: 'CGATCG', cutTop: 4, cutBottom: 2},
    Enzyme{name: 'PvuII', site: 'CAGCTG', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'RsaI', site: 'GTAC', cutTop: 2, cutBottom: 2},
    Enzyme{name: 'SacI', site: 'GAGCTC', cutTop: 5, cutBottom: 1},
    Enzyme{name: 'SacII', site: 'CCGCGG', cutTop: 4, cutBottom: 2},
    Enzyme{name: 'SalI', site: 'GTCGAC', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'Sau3AI', site: 'GATC', cutTop: 0, cutBottom: 4},
    Enzyme{name: 'SbfI', site: 'CCTGCAGG', cutTop: 6, cutBottom: 2},
    Enzyme{name: 'ScaI', site: 'AGTACT', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'SmaI', site: 'CCCGGG', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'SnaBI', site: 'TACGTA', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'SpeI', site: 'ACTAGT', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'SphI', site: 'GCATGC', cutTop: 5, cutBottom: 1},
    Enzyme{name: 'SrfI', site: 'GCCCGGGC', cutTop: 4, cutBottom: 4},
    Enzyme{name: 'SspI', site: 'AATATT', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'StuI', site: 'AGGCCT', cutTop: 3, cutBottom: 3},
    Enzyme{name: 'StyI', site: 'CCWWGG', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'SwaI', site: 'ATTTAAAT', cutTop: 4, cutBottom: 4},
    Enzyme{name: 'TaqI', site: 'TCGA', cutTop: 1, cutBottom: 3},
    Enzyme{name: 'XbaI', site: 'TCTAGA', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'XhoI', site: 'CTCGAG', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'XmaI', site: 'CCCGGG', cutTop: 1, cutBottom: 5},
    Enzyme{name: 'XmnI', site: 'GAANNNNTTC', cutTop: 5, cutBottom: 5}
];

def const IUPAC as string init 'ACGTRYSWKMBDHVN';
def const COMPLEMENTS as string init 'TGCAYRSWMKVDHBN';

func fail(message as string) {
    io.printf("genenzymes: %s\n", $message);
    exit 1;
}

func complementOf(base as string) {
    def at as int init strings.indexOf(IUPAC, $base);
    if ($at < 0) {
        fail("no complement for '" + $base + "'");
    }
    return COMPLEMENTS[$at..$at + 1];
}

func reverseComplement(site as string) {
    def flipped as string init '';
    def chars as list of string init strings.chars($site);
    for (def i as int init len($chars) - 1; $i >= 0; $i = $i - 1) {
        $flipped = $flipped + complementOf($chars[$i]);
    }
    return $flipped;
}

# ------------------------------------------------------------------ validation

def seen as map of string to bool;
for (def e in $enzymes) {
    if (maps.has($seen, $e.name)) {
        fail('duplicate enzyme ' + $e.name);
    }
    $seen[$e.name] = true;

    for (def base in strings.chars($e.site)) {
        if (not strings.contains(IUPAC, $base)) {
            fail($e.name + ": '" + $base + "' is not an IUPAC nucleotide code");
        }
    }
    if ($e.cutTop < 0 or $e.cutTop > len($e.site)) {
        fail($e.name + ': cutTop outside the recognition site');
    }
    if ($e.cutBottom < 0 or $e.cutBottom > len($e.site)) {
        fail($e.name + ': cutBottom outside the recognition site');
    }

    # A type II enzyme binds as a homodimer, so its site reads the same on both
    # strands. This holds for the degenerate sites too, because N complements to
    # N and W to W: the specificity lives in the outer bases, and those pair up.
    # A new row whose site fails this is almost certainly a typo, so it is a hard
    # check with no exception list.
    def flipped as string init reverseComplement($e.site);
    if ($flipped != $e.site) {
        fail($e.name + ': site ' + $e.site + ' is not palindromic (rc ' + $flipped + ')');
    }
}

# ------------------------------------------------------------------- emission

def out as list of string;

func emit(line as string) {
    $out[] = $line;
}

emit('# SPDX-License-Identifier: LGPL-3.0-only');
emit('# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>');
emit('# pragma-jennifer-version: >=0.25.0');
emit('');
emit('# bioinformatics - enzymes: the restriction enzyme table.');
emit('# Spliced into bioinformatics.j via include - not a standalone module.');
emit('#');
emit('# GENERATED by tools/genenzymes.j. Recognition sequences and cut positions');
emit('# are hand-curated from published enzyme specifications rather than imported');
emit('# in bulk from REBASE, whose data carries its own usage terms. Adding an');
emit('# enzyme is a four-field row, which is the point: the table is meant to be');
emit('# extended by whoever needs the enzyme, not to be exhaustive on day one.');
emit('#');
emit('# Cut offsets are counted along the top strand from the first base of the');
emit('# recognition site:');
emit('#');
emit("#     EcoRI  G^AATTC   site \"GAATTC\", cutTop 1, cutBottom 5, 5' overhang");
emit("#     PstI   CTGCA^G   site \"CTGCAG\", cutTop 5, cutBottom 1, 3' overhang");
emit('#     SmaI   CCC^GGG   site "CCCGGG", cutTop 3, cutBottom 3, blunt');
emit('#');
emit('# Type IIS enzymes (BsaI, BsmBI, ...), which cut at a distance outside their');
emit('# recognition site, are deliberately absent: their cut geometry does not fit');
emit('# these two offsets and modelling them properly is a separate job.');
emit('');
emit('/**');
emit(' * One restriction enzyme: what it recognises and where it cuts.');
emit(' *');
emit(' * @field name {string} the conventional name, e.g. "EcoRI"');
emit(' * @field site {string} the recognition sequence in IUPAC codes, so a');
emit(' *        degenerate site like `CCWGG` is written as it is specified');
emit(' * @field cutTop {int} offset into the site where the top strand is nicked');
emit(' * @field cutBottom {int} offset where the bottom strand is nicked, counted');
emit(' *        along the top strand');
emit(" * @field overhang \{string\} \"5'\", \"3'\" or \"blunt\"");
emit(' * @field overhangLength {int} number of unpaired bases left behind, 0 when blunt');
emit(' */');
emit('export def struct Enzyme {');
emit('    name as string,');
emit('    site as string,');
emit('    cutTop as int,');
emit('    cutBottom as int,');
emit('    overhang as string,');
emit('    overhangLength as int');
emit('};');
emit('');
emit('# Keyed by the uppercased name so lookup can be case-insensitive while the');
emit('# `name` field keeps the conventional spelling (EcoRI, not ECORI).');
emit('def const ENZYME_TABLE as map of string to Enzyme init {');

def rows as list of string;
for (def e in $enzymes) {
    def overhang as string init 'blunt';
    def overhangLength as int init 0;
    if ($e.cutTop < $e.cutBottom) {
        $overhang = "5'";
        $overhangLength = $e.cutBottom - $e.cutTop;
    } elseif ($e.cutTop > $e.cutBottom) {
        $overhang = "3'";
        $overhangLength = $e.cutTop - $e.cutBottom;
    }
    $rows[] = '    "' + strings.upper($e.name) + '": Enzyme{name: "' + $e.name +
        '", site: "' + $e.site + '", cutTop: ' + convert.toString($e.cutTop) + ', cutBottom: ' +
        convert.toString($e.cutBottom) + ', overhang: "' + $overhang + '", overhangLength: ' +
        convert.toString($overhangLength) + '}';
}
emit(strings.join($rows, ",\n"));
emit('};');
emit('');

fs.writeString('src/topics/enzymedata.j', strings.join($out, "\n") + "\n");
io.printf("wrote src/topics/enzymedata.j - %d enzymes\n", len($enzymes));
