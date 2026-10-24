# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# bioinformatics - enzymes: lookup over the restriction enzyme table.
# Spliced into bioinformatics.j via include - not a standalone module.
#
# The table itself is generated into enzymedata.j; this file is hand-written and
# holds only the accessors over it. Regenerating the data never touches this.
#
# Searching a sequence for sites and cutting it lives in search.j, which needs
# the IUPAC matcher - this topic is the catalogue, not the operation.

/**
 * Look up a restriction enzyme by name.
 *
 * Case-insensitive (`enzyme("ecori")` works); the record carries the
 * conventional spelling.
 *
 * @param name {string} the enzyme name, any case
 * @return {Enzyme} the enzyme record
 * @throws {Error} kind "bioinformatics" when the deck does not carry that enzyme
 * @example bioseq.enzyme("EcoRI").site == "GAATTC"
 */
export func enzyme(name as string) {
    def key as string init strings.upper(strings.trim($name));
    if (not maps.has(ENZYME_TABLE, $key)) {
        raiseError("unknown restriction enzyme '" + $name +
            "' (see bioseq.enzymeNames(); the table is meant to be extended)");
    }
    return ENZYME_TABLE[$key];
}

/**
 * Is this an enzyme the deck carries?
 *
 * @param name {string} the enzyme name, any case
 * @return {bool} true when `enzyme` would succeed
 */
export func hasEnzyme(name as string) {
    return maps.has(ENZYME_TABLE, strings.upper(strings.trim($name)));
}

/**
 * Every enzyme in the table, in alphabetical order.
 *
 * @return {list of Enzyme} the enzyme records
 */
export func enzymes() {
    def all as list of Enzyme;
    for (def key in ENZYME_TABLE) {
        $all[] = ENZYME_TABLE[$key];
    }
    return $all;
}

/**
 * The names of every enzyme in the table, in alphabetical order.
 *
 * @return {list of string} conventional enzyme names
 */
export func enzymeNames() {
    def names as list of string;
    for (def key in ENZYME_TABLE) {
        $names[] = ENZYME_TABLE[$key].name;
    }
    return $names;
}

/**
 * Enzymes that recognise the same sequence as the named one.
 *
 * Matches the stored site string exactly, so both true isoschizomers
 * (`HpaII` / `MspI`, same cut) and neoschizomers (`SmaI` / `XmaI`, different cut)
 * are returned; compare `cutTop` to tell them apart.
 *
 * @param name {string} the enzyme name, any case
 * @return {list of string} the other enzymes with the same site, alphabetical
 * @throws {Error} kind "bioinformatics" when the deck does not carry that enzyme
 * @example bioseq.isoschizomers("SmaI")   # ["XmaI"]
 */
export func isoschizomers(name as string) {
    def target as Enzyme init enzyme($name);
    def out as list of string;
    for (def key in ENZYME_TABLE) {
        def other as Enzyme init ENZYME_TABLE[$key];
        if ($other.site == $target.site and $other.name != $target.name) {
            $out[] = $other.name;
        }
    }
    return $out;
}

/**
 * Does this enzyme leave blunt ends?
 *
 * @param e {Enzyme} the enzyme record
 * @return {bool} true when both strands are cut at the same offset
 */
export func isBlunt(e as Enzyme) {
    return $e.cutTop == $e.cutBottom;
}
