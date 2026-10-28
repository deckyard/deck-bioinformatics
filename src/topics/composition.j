# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - composition: what a sequence is made of, and what that implies.
# Spliced into bioinformatics.j via include - not a standalone module.
#
# This is the topic where numbers come from published parameter sets rather than
# from structure, and where two sets can disagree about the same sequence. Every
# formula here names its source in its docblock, and the parameters themselves
# are exported constants, so a caller who needs a different set can see exactly
# what this one used.
#
# Ambiguity codes are rejected rather than averaged. There is no defensible mass
# for a base that might be any of four, and quietly returning the mean of the
# possibilities would produce a number that looks authoritative and is not.

# ------------------------------------------------------------------ base counts

/**
 * Count each residue in a nucleotide sequence.
 *
 * Only residues present get an entry: guard a lookup with `maps.has`.
 *
 * @param s {string} a DNA or RNA sequence
 * @return {map of string to int} residue -> count
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic acid
 * @example bioseq.baseCount("ATGC")["A"] == 1
 */
export func baseCount(s as string) {
    return countResidues(requireNucleotide($s));
}

/**
 * Count each residue in a protein sequence.
 *
 * @param p {string} a protein sequence
 * @return {map of string to int} residue -> count
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a protein
 */
export func residueCount(p as string) {
    return countResidues(requireProtein($p));
}

/**
 * Tally the residues of an already-validated sequence.
 *
 * @param clean {string} a validated sequence
 * @return {map of string to int} residue -> count
 * @internal
 */
func countResidues(clean as string) {
    # Validated input is ASCII, so the tally is a list indexed by byte and only
    # the residues seen become map entries.
    def tally as list of int;
    for (def i as int init 0; $i < 128; $i = $i + 1) {
        $tally[] = 0;
    }
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        $tally[$buf[$i]] = $tally[$buf[$i]] + 1;
    }
    def out as map of string to int;
    for (def code as int init 0; $code < 128; $code = $code + 1) {
        if ($tally[$code] > 0) {
            def one as bytes;
            $one[] = $code;
            $out[convert.stringFromBytes($one, "utf-8")] = $tally[$code];
        }
    }
    return $out;
}

/**
 * The proportion of the sequence each residue makes up, as a percentage.
 *
 * @param s {string} a DNA or RNA sequence
 * @return {map of string to float} residue -> percentage
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic acid
 */
export func baseFrequency(s as string) {
    return toPercentages(baseCount($s));
}

/**
 * The proportion of a protein each residue makes up, as a percentage.
 *
 * @param p {string} a protein sequence
 * @return {map of string to float} residue -> percentage
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a protein
 */
export func residueFrequency(p as string) {
    return toPercentages(residueCount($p));
}

/**
 * Convert a count map to percentages of its own total.
 *
 * @param counts {map of string to int} residue -> count
 * @return {map of string to float} residue -> percentage
 * @internal
 */
func toPercentages(counts as map of string to int) {
    def total as int init 0;
    for (def residue in $counts) {
        $total = $total + $counts[$residue];
    }
    def out as map of string to float;
    for (def residue in $counts) {
        $out[$residue] = 100.0 * $counts[$residue] / $total;
    }
    return $out;
}

# ------------------------------------------------------------------ GC content

/**
 * The percentage of the sequence that is G or C.
 *
 * Ambiguity codes count toward the length but not toward G or C, so
 * `gcContent` and `atContent` sum to less than 100 on a sequence carrying `N`
 * (`S` is not counted as GC).
 *
 * @param s {string} a DNA or RNA sequence
 * @return {float} the GC percentage, 0 to 100
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic acid
 * @example bioseq.gcContent("ATGCGC") == 66.66666666666667
 */
export func gcContent(s as string) {
    def clean as string init requireNucleotide($s);
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    def hits as int init 0;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        if ($buf[$i] == 71 or $buf[$i] == 67) {
            $hits = $hits + 1;
        }
    }
    return 100.0 * $hits / $n;
}

/**
 * The percentage of the sequence that is A, T or U. Ambiguity codes count toward
 * the length only.
 *
 * @param s {string} a DNA or RNA sequence
 * @return {float} the AT percentage, 0 to 100
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a nucleic acid
 */
export func atContent(s as string) {
    def clean as string init requireNucleotide($s);
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    def hits as int init 0;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        if ($buf[$i] == 65 or $buf[$i] == 84 or $buf[$i] == 85) {
            $hits = $hits + 1;
        }
    }
    return 100.0 * $hits / $n;
}

/**
 * GC content measured in a sliding window.
 *
 * Windows start at 0, `step`, `2 * step`, ...; the last partial window is
 * dropped, so every value covers exactly `size` bases. `windows` cuts the same
 * way.
 *
 * @param s {string} a DNA or RNA sequence
 * @param size {int} the window width in bases; must be positive
 * @param step {int} how far to advance between windows; must be positive
 * @return {list of float} the GC percentage of each full window
 * @throws {Error} kind "bioinformatics" when the sequence is shorter than one window,
 *         or a parameter is not positive
 * @example bioseq.gcWindows($genome, 100, 50)
 */
export func gcWindows(s as string, size as int, step as int) {
    def clean as string init requireNucleotide($s);
    if ($size <= 0) {
        raiseError("window size must be positive, got " + convert.toString($size));
    }
    if ($step <= 0) {
        raiseError("window step must be positive, got " + convert.toString($step));
    }
    if (len($clean) < $size) {
        raiseError("sequence is " + convert.toString(len($clean)) + " bases, shorter than the " +
            convert.toString($size) + "-base window");
    }
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    def out as list of float;
    def start as int init 0;
    while ($start + $size <= $n) {
        def hits as int init 0;
        for (def i as int init $start; $i < $start + $size; $i = $i + 1) {
            if ($buf[$i] == 71 or $buf[$i] == 67) {
                $hits = $hits + 1;
            }
        }
        $out[] = 100.0 * $hits / $size;
        $start = $start + $step;
    }
    return $out;
}

# ------------------------------------------------------------- molecular weight

/**
 * The average molecular weight of a sequence, in daltons.
 *
 * For `DNA` and `RNA`: one strand with a free 5' hydroxyl (an unphosphorylated
 * synthetic oligo), computed as the residue masses plus `WATER_AVG_MASS` minus
 * `PHOSPHATE_MASS`. Roughly double it for a duplex. For `PROTEIN`: the average
 * residue masses plus one water, the ExPASy convention.
 *
 * Ambiguity codes (and `*`, `B`, `Z`, `X`, `J` in protein) raise rather than
 * being averaged.
 *
 * @param s {string} the sequence
 * @param kind {int} `bioseq.DNA`, `bioseq.RNA` or `bioseq.PROTEIN`
 * @return {float} the mass in daltons
 * @throws {Error} kind "bioinformatics" when the kind is unknown, the sequence is empty
 *         or of the wrong alphabet, or it carries an ambiguity code
 * @example bioseq.molecularWeight("ATGC", bioseq.DNA) == 1173.83
 */
export func molecularWeight(s as string, kind as int) {
    requireKind($kind);
    if ($kind == PROTEIN) {
        def protein as string init requireProtein($s);
        def pbuf as bytes init convert.bytesFromString($protein, "utf-8");
        def pn as int init len($pbuf);
        def ptotal as float init WATER_AVG_MASS;
        for (def i as int init 0; $i < $pn; $i = $i + 1) {
            def code as int init $pbuf[$i];
            if (not RESIDUE_KNOWN_BY_CODE[$code]) {
                raiseError("no molecular weight for " + describeCode($code) +
                    " - resolve it, or drop it, before weighing the sequence");
            }
            $ptotal = $ptotal + AVG_MASS_BY_CODE[$code];
        }
        return $ptotal;
    }
    def clean as string init "";
    def wantDna as bool init $kind == DNA;
    if ($wantDna) {
        $clean = requireDna($s);
    } else {
        $clean = requireRna($s);
    }
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    def total as float init WATER_AVG_MASS - PHOSPHATE_MASS;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def code as int init $buf[$i];
        def mass as float init RNA_MASS_BY_CODE[$code];
        if ($wantDna) {
            $mass = DNA_MASS_BY_CODE[$code];
        }
        if ($mass < 0.0) {
            raiseError("no molecular weight for the ambiguity code " + describeCode($code) +
                " - resolve it, or drop it, before weighing the sequence");
        }
        $total = $total + $mass;
    }
    return $total;
}

# ------------------------------------------------------------- melting point

/**
 * Melting temperature by the Wallace rule, in degrees Celsius.
 *
 * `Tm = 2(A+T) + 4(G+C)`. Calibrated for oligos of about 14 to 20 bases; it
 * ignores order, salt and strand concentration, and on longer sequences returns
 * a number with no physical meaning. Use `meltingTempNn` when the answer
 * matters. RNA is rejected.
 *
 * @param s {string} a DNA sequence
 * @return {float} the melting temperature in Celsius
 * @throws {Error} kind "bioinformatics" when the sequence is empty, not DNA, or carries
 *         an ambiguity code
 * @see Wallace RB et al. (1979) Nucleic Acids Res 6:3543-3557
 * @example bioseq.meltingTemp("ATGCATGCATGCATGC") == 48.0
 */
export func meltingTemp(s as string) {
    def clean as string init requireUnambiguousDna($s);
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def at as int init 0;
    def gc as int init 0;
    for (def i as int init 0; $i < len($buf); $i = $i + 1) {
        if ($buf[$i] == 71 or $buf[$i] == 67) {
            $gc = $gc + 1;
        } else {
            $at = $at + 1;
        }
    }
    return 2.0 * $at + 4.0 * $gc;
}

/**
 * Melting temperature by the nearest-neighbour model at 250 nM primer and 50 mM
 * Na+: `meltingTempNnWith(s, 250.0, 50.0)`.
 *
 * @param s {string} a DNA sequence of at least two bases
 * @return {float} the melting temperature in Celsius
 * @throws {Error} kind "bioinformatics" when the sequence is too short, not DNA, or
 *         carries an ambiguity code
 */
export func meltingTempNn(s as string) {
    return meltingTempNnWith($s, 250.0, 50.0);
}

/**
 * Melting temperature by the nearest-neighbour model, at chosen conditions.
 *
 * Enthalpy and entropy are summed over the dinucleotides, plus helix initiation
 * for each terminal base (G-C or A-T), a salt term of `0.368 * (N - 1) * ln[Na+]`
 * on the entropy, and, for a self-complementary strand, the symmetry term with
 * `Ct` in place of `Ct/4`. Then `Tm = dH / (dS + R ln(Ct/4)) - 273.15`.
 *
 * Parameters: SantaLucia (1998) unified set, exported as `NN_ENTHALPY`,
 * `NN_ENTROPY` and `GAS_CONSTANT`. Other sets differ by a degree or two.
 *
 * @param s {string} a DNA sequence of at least two bases
 * @param primerNm {float} total strand concentration in nanomolar
 * @param saltMm {float} monovalent cation concentration in millimolar
 * @return {float} the melting temperature in Celsius
 * @throws {Error} kind "bioinformatics" when the sequence is too short, not DNA,
 *         carries an ambiguity code, or a concentration is not positive
 * @see SantaLucia J (1998) PNAS 95:1460-1465
 */
export func meltingTempNnWith(s as string, primerNm as float, saltMm as float) {
    def clean as string init requireUnambiguousDna($s);
    if (len($clean) < 2) {
        raiseError("the nearest-neighbour model needs at least two bases, got " +
            convert.toString(len($clean)));
    }
    if ($primerNm <= 0.0) {
        raiseError("primer concentration must be positive");
    }
    if ($saltMm <= 0.0) {
        raiseError("salt concentration must be positive");
    }

    def enthalpy as float init 0.0;
    def entropy as float init 0.0;
    def pairs as int init len($clean) - 1;
    for (def i as int init 0; $i < $pairs; $i = $i + 1) {
        def pair as string init $clean[$i..$i + 2];
        $enthalpy = $enthalpy + NN_ENTHALPY[$pair];
        $entropy = $entropy + NN_ENTROPY[$pair];
    }

    # helix initiation, once for each end, keyed on whether that end is a strong
    # G-C pair or a weak A-T one
    def ends as list of string init [$clean[0..1], $clean[len($clean) - 1..]];
    for (def base in $ends) {
        if ($base == "G" or $base == "C") {
            $enthalpy = $enthalpy + NN_INIT_GC_ENTHALPY;
            $entropy = $entropy + NN_INIT_GC_ENTROPY;
        } else {
            $enthalpy = $enthalpy + NN_INIT_AT_ENTHALPY;
            $entropy = $entropy + NN_INIT_AT_ENTROPY;
        }
    }

    # a self-complementary strand anneals to a copy of itself, which changes
    # both the symmetry term and the effective concentration
    def divisor as float init 4.0;
    if ($clean == reverseComplement($clean)) {
        $entropy = $entropy + NN_SYMMETRY_ENTROPY;
        $divisor = 1.0;
    }

    # salt correction, applied per phosphate
    def sodium as float init $saltMm / 1000.0;
    $entropy = $entropy + 0.368 * (len($clean) - 1) * math.ln($sodium);

    def strands as float init $primerNm / 1000000000.0;
    def denominator as float init $entropy + GAS_CONSTANT * math.ln($strands / $divisor);
    return $enthalpy * 1000.0 / $denominator - 273.15;
}

/**
 * Validate a DNA sequence and reject ambiguity codes.
 *
 * The thermodynamic and mass verbs all need concrete bases; an ambiguity code
 * has no stacking parameter and no mass.
 *
 * @param s {string} the sequence as supplied
 * @return {string} cleaned, uppercased, unambiguous DNA
 * @throws {Error} kind "bioinformatics" when empty, not DNA, or ambiguous
 * @internal
 */
func requireUnambiguousDna(s as string) {
    def clean as string init requireDna($s);
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    for (def i as int init 0; $i < len($buf); $i = $i + 1) {
        def code as int init $buf[$i];
        if ($code != 65 and $code != 67 and $code != 71 and $code != 84) {
            raiseError("ambiguity code " + describeCode($code) + " at position " +
                convert.toString($i) + " has no thermodynamic parameters - resolve it first");
        }
    }
    return $clean;
}

# ------------------------------------------------------------ protein chemistry

/**
 * The net charge of a protein at a given pH.
 *
 * Sums the fractional charge of the two termini and the side chains of D, E, C,
 * Y (acidic) and H, K, R (basic). pKa values are the EMBOSS `iep` set, exported
 * through `aminoAcid` and the `PKA_*` constants. Decreases monotonically with pH.
 *
 * @param p {string} a protein sequence
 * @param ph {float} the pH to evaluate at
 * @return {float} the net charge in elementary charges
 * @throws {Error} kind "bioinformatics" when the sequence is empty, not a protein, or
 *         carries an ambiguity code
 * @example bioseq.netCharge("MKKK", 7.0)
 */
export func netCharge(p as string, ph as float) {
    def clean as string init requireProtein($p);
    def counts as list of int init ionisableCounts($clean);
    return chargeOfCounts($counts, presentCodes($counts), $ph);
}

/**
 * Tally the ionisable residues once, so a pH sweep does not re-walk the protein.
 *
 * `isoelectricPoint` evaluates the charge up to a hundred times while it
 * bisects. Walking every residue each time made the cost proportional to length
 * times iterations; counting first makes the sweep proportional to the seven
 * ionisable residue types instead. On a 20 000-residue protein that was the
 * difference between seven seconds and a few milliseconds.
 *
 * @param clean {string} a validated protein sequence
 * @return {list of int} counts indexed by residue byte value
 * @throws {Error} kind "bioinformatics" when a residue has no pKa
 * @internal
 */
func ionisableCounts(clean as string) {
    def counts as list of int;
    for (def i as int init 0; $i < 128; $i = $i + 1) {
        $counts[] = 0;
    }
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def code as int init $buf[$i];
        if (not RESIDUE_KNOWN_BY_CODE[$code]) {
            raiseError("no pKa for " + describeCode($code) +
                " - resolve the ambiguity code before charging the sequence");
        }
        if (CHARGE_BY_CODE[$code] != 0) {
            $counts[$code] = $counts[$code] + 1;
        }
    }
    return $counts;
}

/**
 * The net charge implied by a tally of ionisable residues.
 *
 * @param counts {list of int} counts indexed by residue byte value
 * @param present {list of int} the codes to charge (`presentCodes`)
 * @param ph {float} the pH to evaluate at
 * @return {float} the net charge
 * @internal
 */
func chargeOfCounts(counts as list of int, present as list of int, ph as float) {
    def nTerminus as float init 1.0 / (1.0 + math.pow(10.0, $ph - PKA_N_TERMINUS));
    def cTerminus as float init 1.0 / (1.0 + math.pow(10.0, PKA_C_TERMINUS - $ph));
    def charge as float init $nTerminus - $cTerminus;
    for (def code in $present) {
        def many as int init $counts[$code];
        def sign as int init CHARGE_BY_CODE[$code];
        def acidity as float init PKA_BY_CODE[$code];
        if ($sign > 0) {
            $charge = $charge + $many / (1.0 + math.pow(10.0, $ph - $acidity));
        } elseif ($sign < 0) {
            $charge = $charge - $many / (1.0 + math.pow(10.0, $acidity - $ph));
        }
    }
    return $charge;
}

/**
 * The residue codes a tally counts at least once.
 *
 * @param counts {list of int} counts indexed by residue byte value
 * @return {list of int} the codes with a non-zero count
 * @internal
 */
func presentCodes(counts as list of int) {
    def out as list of int;
    for (def code as int init 0; $code < 128; $code = $code + 1) {
        if ($counts[$code] > 0) {
            $out[] = $code;
        }
    }
    return $out;
}

/**
 * The pH at which a protein carries no net charge.
 *
 * Bisection on `netCharge` over pH 0 to 14, stopping when the charge is within
 * 1e-4 of zero or after 100 halvings. Uses the EMBOSS `iep` pKa set: compare
 * against `iep`, not ProtParam (pKa sets differ by a few tenths of a unit).
 *
 * @param p {string} a protein sequence
 * @return {float} the isoelectric point
 * @throws {Error} kind "bioinformatics" when the sequence is empty, not a protein, or
 *         carries an ambiguity code
 * @example bioseq.isoelectricPoint("MKWVTFISLLLLFSSAYS")
 */
export func isoelectricPoint(p as string) {
    def clean as string init requireProtein($p);
    def counts as list of int init ionisableCounts($clean);
    def present as list of int init presentCodes($counts);
    def low as float init 0.0;
    def high as float init 14.0;
    def mid as float init 7.0;
    for (def i as int init 0; $i < 100; $i = $i + 1) {
        $mid = ($low + $high) / 2.0;
        def charge as float init chargeOfCounts($counts, $present, $mid);
        if (math.abs($charge) < 0.0001) {
            return $mid;
        }
        if ($charge > 0.0) {
            $low = $mid;
        } else {
            $high = $mid;
        }
    }
    return $mid;
}

/**
 * The grand average of hydropathy - the mean Kyte-Doolittle index.
 *
 * Positive is hydrophobic. Soluble globular proteins typically score -0.5 to
 * -1.5; membrane-spanning ones score positive.
 *
 * @param p {string} a protein sequence
 * @return {float} the mean hydropathy
 * @throws {Error} kind "bioinformatics" when the sequence is empty, not a protein, or
 *         carries an ambiguity code
 * @see Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132
 */
export func gravy(p as string) {
    def clean as string init requireProtein($p);
    def buf as bytes init convert.bytesFromString($clean, "utf-8");
    def n as int init len($buf);
    def total as float init 0.0;
    for (def i as int init 0; $i < $n; $i = $i + 1) {
        def code as int init $buf[$i];
        if (not RESIDUE_KNOWN_BY_CODE[$code]) {
            raiseError("no hydropathy index for " + describeCode($code) +
                " - resolve the ambiguity code first");
        }
        $total = $total + HYDROPATHY_BY_CODE[$code];
    }
    return $total / $n;
}

/**
 * The proportion of aromatic residues - phenylalanine, tryptophan and tyrosine.
 *
 * @param p {string} a protein sequence
 * @return {float} the fraction, 0 to 1
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a protein
 */
export func aromaticity(p as string) {
    def clean as string init requireProtein($p);
    def hits as int init 0;
    for (def residue in strings.chars($clean)) {
        if ($residue == "F" or $residue == "W" or $residue == "Y") {
            $hits = $hits + 1;
        }
    }
    return 1.0 * $hits / len($clean);
}

# ------------------------------------------------------------------- summaries

/**
 * A summary of a DNA sequence.
 *
 * @field length {int} bases
 * @field gc {float} GC percentage
 * @field at {float} AT percentage
 * @field counts {map of string to int} residue -> count
 * @field molecularWeight {float} single-stranded mass in daltons, 0.0 when the
 *        sequence carries ambiguity codes and has no defined mass
 * @field meltingTemp {float} Wallace-rule Tm in Celsius, 0.0 when undefined;
 *        meaningless beyond about 20 bases
 * @field ambiguous {bool} true when the sequence carries any non-ACGT code, in
 *        which case the two numbers above are 0.0 rather than guessed
 */
export def struct DnaStats {
    length as int,
    gc as float,
    at as float,
    counts as map of string to int,
    molecularWeight as float,
    meltingTemp as float,
    ambiguous as bool
};

/**
 * Summarise a DNA sequence in one pass.
 *
 * A sequence carrying ambiguity codes still gets its length, counts, GC and AT;
 * mass and melting temperature are 0.0 and `ambiguous` is set. RNA is rejected.
 *
 * @param s {string} a DNA sequence
 * @return {DnaStats} the summary
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not DNA
 */
export func dnaStats(s as string) {
    def clean as string init requireDna($s);
    def counts as map of string to int init countResidues($clean);
    def ambiguous as bool init false;
    for (def base in $counts) {
        if ($base != "A" and $base != "C" and $base != "G" and $base != "T") {
            $ambiguous = true;
        }
    }
    def mass as float init 0.0;
    def tm as float init 0.0;
    if (not $ambiguous) {
        $mass = molecularWeight($clean, DNA);
        $tm = meltingTemp($clean);
    }
    return DnaStats{
        length: len($clean),
        gc: gcContent($clean),
        at: atContent($clean),
        counts: $counts,
        molecularWeight: $mass,
        meltingTemp: $tm,
        ambiguous: $ambiguous
    };
}

/**
 * A summary of a protein sequence.
 *
 * @field length {int} residues
 * @field molecularWeight {float} average mass in daltons, 0.0 when undefined
 * @field isoelectricPoint {float} pI, 0.0 when undefined
 * @field gravy {float} mean Kyte-Doolittle hydropathy, 0.0 when undefined
 * @field aromaticity {float} fraction of F, W and Y
 * @field counts {map of string to int} residue -> count
 * @field ambiguous {bool} true when the sequence carries a code with no mass or
 *        pKa (`B`, `Z`, `X`, `J`) or a stop, in which case the three numbers
 *        above are 0.0 rather than guessed
 */
export def struct ProteinStats {
    length as int,
    molecularWeight as float,
    isoelectricPoint as float,
    gravy as float,
    aromaticity as float,
    counts as map of string to int,
    ambiguous as bool
};

/**
 * Summarise a protein sequence in one pass.
 *
 * A sequence carrying `B`, `Z`, `X`, `J` or `*` still gets its length, counts and
 * aromaticity; mass, pI and GRAVY are 0.0 and `ambiguous` is set.
 *
 * @param p {string} a protein sequence
 * @return {ProteinStats} the summary
 * @throws {Error} kind "bioinformatics" when the sequence is empty or not a protein
 */
export func proteinStats(p as string) {
    def clean as string init requireProtein($p);
    def counts as map of string to int init countResidues($clean);
    def ambiguous as bool init false;
    for (def residue in $counts) {
        if (not hasAminoAcid($residue)) {
            $ambiguous = true;
        }
    }
    def mass as float init 0.0;
    def pi as float init 0.0;
    def hydropathy as float init 0.0;
    if (not $ambiguous) {
        $mass = molecularWeight($clean, PROTEIN);
        $pi = isoelectricPoint($clean);
        $hydropathy = gravy($clean);
    }
    return ProteinStats{
        length: len($clean),
        molecularWeight: $mass,
        isoelectricPoint: $pi,
        gravy: $hydropathy,
        aromaticity: aromaticity($clean),
        counts: $counts,
        ambiguous: $ambiguous
    };
}

# ----------------------------------------------------------------- codon usage

/**
 * How often one codon is used, and how it compares to its synonyms.
 *
 * @field codon {string} the codon
 * @field aminoAcid {string} the residue it encodes under the table used
 * @field count {int} occurrences in the sequence
 * @field fraction {float} share of the codons encoding this same residue, 0 to 1
 * @field perThousand {float} occurrences per thousand codons in the sequence
 */
export def struct CodonUse {
    codon as string,
    aminoAcid as string,
    count as int,
    fraction as float,
    perThousand as float
};

/**
 * Tabulate codon usage under the standard genetic code.
 *
 * @param nt {string} a coding sequence
 * @return {list of CodonUse} one entry per codon present
 * @throws {Error} kind "bioinformatics" when the sequence is empty, not a nucleic acid,
 *         or shorter than one codon
 */
export func codonUsage(nt as string) {
    return codonUsageWith($nt, STANDARD);
}

/**
 * Tabulate codon usage under a named genetic code.
 *
 * Codons are reported in DNA spelling (`U` read as `T`). A trailing partial codon
 * is dropped. A codon the table does not hold, such as one with an ambiguity
 * code, is counted under residue `X`. Codons absent from the sequence are absent
 * from the result.
 *
 * @param nt {string} a coding sequence
 * @param table {map of string to string} a genetic code
 * @return {list of CodonUse} one entry per codon present, in first-seen
 *         order
 * @throws {Error} kind "bioinformatics" when the sequence is empty, not a nucleic acid,
 *         or shorter than one codon
 */
export func codonUsageWith(nt as string, table as map of string to string) {
    def codons as list of string init splitCodons($nt);
    if (len($codons) == 0) {
        raiseError("sequence is shorter than one codon");
    }
    def counts as map of string to int;
    def order as list of string;
    for (def codon in $codons) {
        def key as string init strings.replace($codon, "U", "T");
        if (maps.has($counts, $key)) {
            $counts[$key] = $counts[$key] + 1;
        } else {
            $counts[$key] = 1;
            $order[] = $key;
        }
    }
    # totals per residue, so each codon's share of its synonyms can be reported
    def perResidue as map of string to int;
    for (def codon in $counts) {
        def residue as string init "X";
        if (maps.has($table, $codon)) {
            $residue = $table[$codon];
        }
        if (maps.has($perResidue, $residue)) {
            $perResidue[$residue] = $perResidue[$residue] + $counts[$codon];
        } else {
            $perResidue[$residue] = $counts[$codon];
        }
    }
    def total as int init len($codons);
    def out as list of CodonUse;
    for (def codon in $order) {
        def residue as string init "X";
        if (maps.has($table, $codon)) {
            $residue = $table[$codon];
        }
        $out[] = CodonUse{
            codon: $codon,
            aminoAcid: $residue,
            count: $counts[$codon],
            fraction: 1.0 * $counts[$codon] / $perResidue[$residue],
            perThousand: 1000.0 * $counts[$codon] / $total
        };
    }
    return $out;
}
