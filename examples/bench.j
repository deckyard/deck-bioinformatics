#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - the measurements behind the deck's implementation conventions.
#
#   jennifer run examples/bench.j
#
# Each measurement settled a choice about how the deck's loops are written, and
# is kept so a regression shows up as a number. The commentary records the
# numbers from the machine they were first taken on; compare ratios, not
# milliseconds. docs/design.md discusses them.

use io;
use time;
use convert;

import "../src/bioinformatics.j" as bioseq;

def const N as int init 20000;

func makeSeq(n as int) {
    def alphabet as bytes init convert.bytesFromString("ACGT", "utf-8");
    def b as bytes;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        $b[] = $alphabet[($i * 7 + $i // 4) % 4];
    }
    return convert.stringFromBytes($b, "utf-8");
}

func elapsedMs(a as time.Time, b as time.Time) {
    def d as time.Duration init time.sub($b, $a);
    return $d.nanos / 1000000.0;
}

func report(label as string, ms as float) {
    io.printf("  %s|pad=34|align=left %f|prec=2 ms\n", $label, $ms);
}

def seq as string init makeSeq(N);
io.printf("bioinformatics benchmarks - %d nt\n\n", len($seq));

# ---------------------------------------------------------------------------
# 1. Per-residue access. Jennifer strings are rune-indexed, so `$s[$i..$i + 1]`
#    rescans the string on every access and the loop goes quadratic. The gap is
#    around 40x at this N, on every machine it has been run on. Hence: every
#    per-residue loop in this deck walks a `bytes` buffer.
# ---------------------------------------------------------------------------

io.printf("1. per-residue access\n");

def t0 as time.Time init time.now();
def gcA as int init 0;
for (def i as int init 0; $i < N; $i = $i + 1) {
    def c as string init $seq[$i..$i + 1];
    if ($c == "G" or $c == "C") {
        $gcA = $gcA + 1;
    }
}
def t1 as time.Time init time.now();
report("string slice $s[i..i+1]", elapsedMs($t0, $t1));

def t2 as time.Time init time.now();
def buf as bytes init convert.bytesFromString($seq, "utf-8");
def gcB as int init 0;
for (def i as int init 0; $i < N; $i = $i + 1) {
    def code as int init $buf[$i];
    if ($code == 71 or $code == 67) {
        $gcB = $gcB + 1;
    }
}
def t3 as time.Time init time.now();
report("bytes index $b[i]", elapsedMs($t2, $t3));
io.printf("  (counts agree: %t)\n\n", $gcA == $gcB);

# ---------------------------------------------------------------------------
# 2. Accreting a result. `lists.push` returns a fresh list every pass, so the
#    loop copies the whole thing N times. Against the `$out[] =` append sugar
#    that is three to four thousand times slower - the exact multiple moves with
#    the machine, the order of magnitude does not. Hence: results are built into
#    a `bytes` buffer and converted once at the end.
# ---------------------------------------------------------------------------

io.printf("2. building a %d-residue result\n", N);

def t4 as time.Time init time.now();
def outStr as string init "";
for (def i as int init 0; $i < N; $i = $i + 1) {
    $outStr = $outStr + "A";
}
def t5 as time.Time init time.now();
report("string + in a loop", elapsedMs($t4, $t5));

def t6 as time.Time init time.now();
def outBytes as bytes;
for (def i as int init 0; $i < N; $i = $i + 1) {
    $outBytes[] = 65;
}
def built as string init convert.stringFromBytes($outBytes, "utf-8");
def t7 as time.Time init time.now();
report("bytes append then convert", elapsedMs($t6, $t7));
io.printf("  (both %d long: %t)\n\n", N, len($outStr) == N and len($built) == N);

# `lists.push` is left out of the timed section on purpose: at N = 20000 it runs
# for the better part of a minute and would make this file unpleasant to run.
# Uncomment to see it, and note the shape rather than the number.
#
#   def slow as list of string;
#   for (def i as int init 0; $i < N; $i = $i + 1) { $slow = lists.push($slow, "A"); }

# ---------------------------------------------------------------------------
# 3. Table access. The deck was designed on the assumption that reading an
#    exported module `def const` map copies it, and that tables therefore had to
#    be passed in as parameters. Measured, it does not: reading the const inside
#    the loop is indistinguishable from hoisting it into a local first, over the
#    64-entry genetic code as well as over a small alphabet set. Threading it
#    through a call is not measurably worse either - all three sit inside each
#    other's run-to-run spread, which is wide enough that any ranking between
#    them would be noise. Hence: the deck reads its tables directly, and takes a
#    table argument only where an alternative genetic code is genuinely
#    wanted.
# ---------------------------------------------------------------------------

io.printf("3. alphabet membership over %d residues\n", N);

def t8 as time.Time init time.now();
def clean as string init bioseq.filterDna($seq);
def t9 as time.Time init time.now();
report("filterDna (const set read)", elapsedMs($t8, $t9));

def t10 as time.Time init time.now();
def ok as bool init bioseq.isDna($seq);
def t11 as time.Time init time.now();
report("isDna (const set read)", elapsedMs($t10, $t11));
io.printf("  (clean: %t, %d residues)\n", $ok, len($clean));

# ---------------------------------------------------------------------------
# 4. Alignment - the deck's one O(n*m) surface, and the reason the docs state a
#    practical ceiling rather than pretending there is none. The matrix is a
#    flat `list of int` and the substitution scores are precomputed into an
#    integer-keyed map, so the inner loop is indexed reads and one lookup.
#    Watch the time quadruple as the sequences double.
# ---------------------------------------------------------------------------

io.printf("\n4. global alignment, growing square\n");

def opts as bioseq.AlignOptions init bioseq.dnaOptions();
for (def size in [50, 100, 200, 400]) {
    def one as string init $seq[0..$size];
    def other as string init bioseq.mutate($one, 0.1, bioseq.DNA);
    def a0 as time.Time init time.now();
    def al as bioseq.Alignment init bioseq.alignGlobal($one, $other, $opts);
    def a1 as time.Time init time.now();
    def label as string init io.sprintf("%d x %d (%d cells)", $size, $size, $size * $size);
    report($label, elapsedMs($a0, $a1));
    if ($al.length < $size) {
        io.printf("  (unexpectedly short alignment)\n");
    }
}

io.printf("  (a 400 x 400 alignment is 160 000 cells; a chromosome is not\n");
io.printf("   a question this deck answers - see docs/design.md)\n");
