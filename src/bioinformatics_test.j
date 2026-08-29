# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# White-box unit tests for the bioinformatics module.
#
# Runs as the co-located test overlay: `jennifer test src/bioinformatics_test.j`.
# The tests are split into topic files under `topics/`, spliced together here
# with `include` - one test file per implementation topic file.
#
# The whole deck is network-free, subprocess-free and filesystem-free, so there
# is no "needs the live service" tier here: everything the deck does is covered.
# Fixtures are embedded as raw '...' strings rather than read from testdata/, so
# the suite runs unchanged on jennifer-tiny.

use testing;
use strings;
use lists;
use maps;
use convert;
use math;

# Shared golden fixtures. The insulin coding fragment and its published protein
# are used by more than one topic - transforms translates it, composition weighs
# it - so they live here rather than in whichever topic file happened to need
# them first. The test files are spliced into one scope, so a name defined in a
# topic test is visible to every other one; keeping the shared ones here makes
# that deliberate rather than accidental.
#
# Homo sapiens insulin (INS), first 72 nt of the preproinsulin CDS, and the
# protein it translates to. A real gene, so a slip in the codon table shows up
# as a wrong protein rather than as a passing tautology.
def const INS_CDS as string init 'ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTG' +
    'GCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCC';
def const INS_PROTEIN as string init "MALWMRLLPLLALLALWGPDPAAA";

# Float comparison, for the topics that produce measurements rather than strings.
func close(got as float, want as float, tolerance as float) {
    return math.abs($got - $want) <= $tolerance;
}

include "topics/core_test.j";
include "topics/tables_test.j";
include "topics/enzymes_test.j";
include "topics/transforms_test.j";
include "topics/formats_test.j";
include "topics/composition_test.j";
include "topics/search_test.j";
include "topics/manipulate_test.j";
include "topics/compare_test.j";
