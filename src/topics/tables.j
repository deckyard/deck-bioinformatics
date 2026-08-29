# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - tables: lookup over the reference data.
# Spliced into bioinformatics.j via include - not a standalone module.
#
# The data itself is generated into tabledata.j - genetic codes, IUPAC codes,
# complements, residue properties. This file is hand-written and holds only the
# accessors over it, so regenerating the data never touches the API.

# --------------------------------------------------------------- codon tables

/**
 * Look a genetic code up by its NCBI `transl_table` number.
 *
 * The named constants (`STANDARD`, `VERTEBRATE_MITOCHONDRIAL`, ...) are the
 * direct way to reach a table; this is the indirect one, for code that carries
 * a table number around rather than a choice made at authoring time.
 *
 * @param id {int} an NCBI genetic code id, one of `codonTableIds()`
 * @return {map of string to string} codon -> one-letter residue, stop as `*`
 * @throws {Error} kind "bioinformatics" when the id is not a supported table
 * @example bioseq.codonTable(2)["AGA"] == "*"
 */
export func codonTable(id as int) {
    if ($id == 1) {
        return STANDARD;
    }
    if ($id == 2) {
        return VERTEBRATE_MITOCHONDRIAL;
    }
    if ($id == 3) {
        return YEAST_MITOCHONDRIAL;
    }
    if ($id == 4) {
        return MOLD_MITOCHONDRIAL;
    }
    if ($id == 5) {
        return INVERTEBRATE_MITOCHONDRIAL;
    }
    if ($id == 6) {
        return CILIATE_NUCLEAR;
    }
    if ($id == 9) {
        return ECHINODERM_MITOCHONDRIAL;
    }
    if ($id == 10) {
        return EUPLOTID_NUCLEAR;
    }
    if ($id == 11) {
        return BACTERIAL;
    }
    raiseError("unknown genetic code " + convert.toString($id) +
        " (supported: see bioseq.codonTableIds())");
}

/**
 * The NCBI name of a genetic code.
 *
 * @param id {int} an NCBI genetic code id
 * @return {string} the table's name, e.g. "Vertebrate Mitochondrial"
 * @throws {Error} kind "bioinformatics" when the id is not a supported table
 */
export func codonTableName(id as int) {
    if (not maps.has(TABLE_NAMES, $id)) {
        raiseError("unknown genetic code " + convert.toString($id) +
            " (supported: see bioseq.codonTableIds())");
    }
    return TABLE_NAMES[$id];
}

/**
 * Every genetic code this deck carries, in ascending NCBI order.
 *
 * The numbering has gaps - 7 and 8 were retired - so this is the list to walk
 * rather than a range.
 *
 * @return {list of int} the supported `transl_table` numbers
 */
export func codonTableIds() {
    def ids as list of int;
    for (def id in TABLE_NAMES) {
        $ids[] = $id;
    }
    return $ids;
}

/**
 * The codons that may initiate translation under a given genetic code.
 *
 * This is what separates several tables that agree on every codon assignment:
 * `BACTERIAL` reads identically to `STANDARD` but recognises seven starts
 * rather than three, which is why an ORF search takes a table id and not just
 * a codon map.
 *
 * @param id {int} an NCBI genetic code id
 * @return {list of string} the start codons
 * @throws {Error} kind "bioinformatics" when the id is not a supported table
 * @example bioseq.startCodons(11)   # TTG, CTG, ATT, ATC, ATA, ATG, GTG
 */
export func startCodons(id as int) {
    if (not maps.has(TABLE_STARTS, $id)) {
        raiseError("unknown genetic code " + convert.toString($id) +
            " (supported: see bioseq.codonTableIds())");
    }
    return strings.split(TABLE_STARTS[$id], ",");
}

# ---------------------------------------------------------------------- IUPAC

/**
 * The bases an IUPAC nucleotide code stands for.
 *
 * @param code {string} one IUPAC nucleotide character, any case
 * @return {string} the bases it admits, e.g. "AG" for `R`
 * @throws {Error} kind "bioinformatics" when the character is not an IUPAC code
 * @example bioseq.expandIupac("r") == "AG"
 */
export func expandIupac(code as string) {
    def key as string init strings.upper($code);
    if (maps.has(IUPAC_DNA, $key)) {
        return IUPAC_DNA[$key];
    }
    if (maps.has(IUPAC_RNA, $key)) {
        return IUPAC_RNA[$key];
    }
    raiseError("not an IUPAC nucleotide code: '" + $code + "'");
}

/**
 * Render an IUPAC code as a regular-expression character class.
 *
 * This is how the ambiguity-aware search in `findIupac` is built: an IUPAC
 * pattern becomes a character class per position, which `regex` then matches in
 * linear time.
 *
 * @param code {string} one IUPAC nucleotide character, any case
 * @return {string} a bracketed class, e.g. "[AG]" for `R`
 * @throws {Error} kind "bioinformatics" when the character is not an IUPAC code
 */
export func iupacClass(code as string) {
    return "[" + expandIupac($code) + "]";
}

/**
 * Does an IUPAC code admit a given base?
 *
 * The relation is one-way on purpose: `N` admits `A`, but `A` does not admit
 * `N`. Ask with the pattern code first and the observed base second.
 *
 * @param code {string} the IUPAC code from a pattern
 * @param base {string} the observed base
 * @return {bool} true when `code` admits `base`
 * @throws {Error} kind "bioinformatics" when either character is not an IUPAC code
 * @example bioseq.matchesIupac("R", "G") == true
 */
export func matchesIupac(code as string, base as string) {
    def admitted as string init expandIupac($code);
    def observed as string init strings.upper($base);
    if (len($observed) != 1) {
        raiseError("expected a single base, got '" + $base + "'");
    }
    return strings.contains($admitted, $observed);
}

# ---------------------------------------------------------------- amino acids

/**
 * Look up a residue's identity, mass and properties.
 *
 * Takes either spelling - `aminoAcid("K")` and `aminoAcid("Lys")` are the same
 * lookup. The ambiguity codes `B`, `Z`, `X` and `J` are deliberately absent:
 * they have no single mass or pKa, so a caller that means one of them has to
 * decide what it means rather than receive a plausible wrong number.
 *
 * @param code {string} a one- or three-letter residue code, any case
 * @return {AminoAcid} the residue record
 * @throws {Error} kind "bioinformatics" when the code is not a known residue
 * @example bioseq.aminoAcid("Trp").monoMass == 186.07931
 */
export func aminoAcid(code as string) {
    def key as string init strings.upper(strings.trim($code));
    if (len($key) == 3) {
        if (not maps.has(THREE_TO_ONE, $key)) {
            raiseError("unknown amino acid '" + $code + "'");
        }
        $key = THREE_TO_ONE[$key];
    }
    if (not maps.has(RESIDUES, $key)) {
        raiseError("unknown amino acid '" + $code +
            "' (the ambiguity codes B, Z, X and J have no single mass or pKa)");
    }
    return RESIDUES[$key];
}

/**
 * Every residue this deck carries properties for.
 *
 * The standard twenty in the conventional one-letter order, followed by
 * selenocysteine and pyrrolysine.
 *
 * @return {list of AminoAcid} the residue records
 */
export func aminoAcids() {
    def all as list of AminoAcid;
    for (def one in RESIDUES) {
        $all[] = RESIDUES[$one];
    }
    return $all;
}

/**
 * Is this a residue the deck carries properties for?
 *
 * The cheap guard in front of `aminoAcid` for a caller walking a sequence that
 * may contain ambiguity codes.
 *
 * @param code {string} a one- or three-letter residue code, any case
 * @return {bool} true when `aminoAcid` would succeed
 */
export func hasAminoAcid(code as string) {
    def key as string init strings.upper(strings.trim($code));
    if (len($key) == 3) {
        if (not maps.has(THREE_TO_ONE, $key)) {
            return false;
        }
        $key = THREE_TO_ONE[$key];
    }
    return maps.has(RESIDUES, $key);
}
