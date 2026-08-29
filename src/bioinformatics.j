# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * @jennifer/bioinformatics - sequence manipulation for DNA, RNA and protein.
 *
 * Modelled on the Sequence Manipulation Suite catalogue: transforms
 * (complement, transcribe, translate), composition and physical properties
 * (GC content, molecular weight, melting temperature, isoelectric point),
 * pattern and feature search (IUPAC patterns, ORFs, restriction sites),
 * FASTA / FASTQ, sequence generation, and pairwise alignment.
 *
 * Sequences are plain **strings**. `bioseq.reverseComplement("ATGC")` is
 * `"GCAT"` - no constructor stands between you and the first call. Input is
 * accepted in any case; residues come back uppercase. A molecule kind is an
 * argument only where it changes the answer (`molecularWeight($s, bioseq.DNA)`),
 * never where it does not (`reverseComplement($s)`). Offsets carried in a struct
 * are 0-based and half-open, matching Jennifer's own slicing, so
 * `$s[$m.start..$m.end]` is the match.
 *
 * Nothing here touches the network, a subprocess, or the filesystem: FASTA and
 * FASTQ are parsed from and written to strings, and opening the file is the
 * caller's job. The deck therefore runs unchanged on `jennifer-tiny`, and every
 * function is a pure transformation that tests without fixtures on disk.
 *
 * Everything raises one error shape - `Error{kind: "bioinformatics"}` - with a message
 * written for a human, and validates its arguments before doing any work.
 *
 * The implementation is split into topic files, spliced together here with
 * `include` - the module boundary (and every `export`) is this file. The
 * co-located white-box tests live in `bioinformatics_test.j`.
 *
 * @module bioinformatics
 * @version 0.1.0
 * @license LGPL-3.0-only
 * @see https://www.bioinformatics.org/sms2/
 * @see https://docs.jennifer-lang.dev/milestones.html#m274---bioinformatics
 * @example
 *   import "@jennifer/bioinformatics/" as bioseq;
 *
 *   def gene as string init bioseq.filterDna("atg acg tta cat gga taa");
 *   io.printf("%s\n", bioseq.reverseComplement($gene));   # TTATCCATGTAACGTCAT
 *   io.printf("%s\n", bioseq.translate($gene));           # MTLHG*
 *   io.printf("%f|prec=1%%\n", bioseq.gcContent($gene));  # 33.3
 */

# Each topic adds what it needs as it lands; lint keeps this list honest.
use maps;
use convert;
use binary;
use io;
use strings;
use math;

# Reference data first: the generated *data.j files declare the tables that the
# hand-written topic files immediately above them read.
include "topics/core.j";
include "topics/tabledata.j";
include "topics/tables.j";
include "topics/enzymedata.j";
include "topics/enzymes.j";
include "topics/transforms.j";
include "topics/formats.j";
include "topics/composition.j";
include "topics/search.j";
include "topics/manipulate.j";
include "topics/compare.j";
