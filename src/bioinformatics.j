# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * @deckyard/bioinformatics - sequence manipulation for DNA, RNA and protein.
 *
 * Modelled on the Sequence Manipulation Suite catalogue: transforms, composition
 * and physical properties, pattern/ORF/restriction search, FASTA and FASTQ,
 * random sequences, and pairwise alignment.
 *
 * Sequences are plain strings. Input is accepted in any case; residues come back
 * uppercase. A molecule kind is an argument only where it changes the answer
 * (`molecularWeight($s, bioseq.DNA)`). Offsets in structs are 0-based and
 * half-open, so `$s[$m.start..$m.end]` is the match. Every error is
 * `Error{kind: "bioinformatics"}`, raised before any work is done.
 *
 * Nothing here touches the network, a subprocess or the filesystem: FASTA and
 * FASTQ are parsed from and written to strings. The deck runs on `jennifer-tiny`.
 *
 * The topic files are spliced in with `include`; this file is the module
 * boundary. Tests: `bioinformatics_test.j`.
 *
 * @module bioinformatics
 * @version 0.2.0-dev
 * @license LGPL-3.0-only
 * @see https://www.bioinformatics.org/sms2/
 * @see https://docs.jennifer-lang.dev/milestones.html#m274---bioinformatics
 * @example
 *   import "@deckyard/bioinformatics/" as bioseq;
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
