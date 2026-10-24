# Cheatsheet

Every export, reached as `bioseq.NAME`. All errors are
`Error{kind: "bioinformatics"}`.

```jennifer
import "@deckyard/bioinformatics/" as bioseq;
```

Conventions ([details](sequences.md)): input in any case, residues out
uppercase (FASTA/FASTQ parsing preserves case); offsets 0-based half-open, in
forward-strand coordinates even for reverse-strand features; `subrange` is the
one 1-based inclusive verb; `filter*` strip silently, everything else raises on
a foreign character; physical quantities raise on ambiguity codes; a `With`
suffix makes a default explicit.

## Alphabets and validation

| Function | Returns | Notes |
| - | - | - |
| `filterDna(s)` | `string` | keep `DNA_ALPHABET` only, uppercased |
| `filterRna(s)` | `string` | keep `RNA_ALPHABET` only |
| `filterProtein(s)` | `string` | keep `PROTEIN_ALPHABET` only |
| `isDna(s)` | `bool` | non-empty and wholly in the alphabet |
| `isRna(s)` | `bool` | |
| `isProtein(s)` | `bool` | |
| `kindName(kind)` | `string` | `"DNA"`, `"RNA"` or `"protein"` |

| Constant | Value | Meaning |
| - | - | - |
| `DNA` | `1` | molecule kind |
| `RNA` | `2` | molecule kind |
| `PROTEIN` | `3` | molecule kind |
| `DNA_ALPHABET` | `"ACGTRYSWKMBDHVN"` | bases plus IUPAC codes |
| `RNA_ALPHABET` | `"ACGURYSWKMBDHVN"` | |
| `PROTEIN_ALPHABET` | `"ACDEFGHIKLMNPQRSTVWYBZXJUO*"` | standard, ambiguous, rare, stop |

## Transforms

| Function | Returns | Notes |
| - | - | - |
| `complement(s)` | `string` | pairs each base, same order |
| `reverseComplement(s)` | `string` | the other strand, 5' to 3' |
| `reverse(s)` | `string` | any sequence |
| `transcribe(dna)` | `string` | T to U |
| `backTranscribe(rna)` | `string` | U to T |
| `translate(nt)` | `string` | standard code, stops as `*` |
| `translateWith(nt, table)` | `string` | codon map, e.g. `VERTEBRATE_MITOCHONDRIAL` |
| `translateFrame(nt, frame)` | `string` | frame 1, 2, 3, -1, -2, -3 |
| `translateFrameWith(nt, frame, table)` | `string` | |
| `sixFrame(nt)` | `list of Frame` | order 1, 2, 3, -1, -2, -3 |
| `sixFrameWith(nt, table)` | `list of Frame` | |
| `splitCodons(nt)` | `list of string` | trailing partial codon dropped |
| `oneToThree(protein)` | `string` | `MKV` to `MetLysVal` |
| `threeToOne(protein)` | `string` | run together or space/hyphen separated |

### `Frame`

| Field | Meaning |
| - | - |
| `frame as int` | 1, 2, 3, or -1, -2, -3 for the reverse complement |
| `protein as string` | the translation, stops as `*` |

## Composition and properties

| Function | Returns | Notes |
| - | - | - |
| `baseCount(s)` | `map of string to int` | present residues only |
| `baseFrequency(s)` | `map of string to float` | percentages |
| `residueCount(p)` | `map of string to int` | |
| `residueFrequency(p)` | `map of string to float` | percentages |
| `gcContent(s)` | `float` | percentage; ambiguity codes count in length only |
| `atContent(s)` | `float` | percentage of A, T, U |
| `gcWindows(s, size, step)` | `list of float` | full windows only |
| `molecularWeight(s, kind)` | `float` | average Da; nucleic acid as 5'-OH single strand |
| `meltingTemp(s)` | `float` | Wallace rule, Celsius, DNA |
| `meltingTempNn(s)` | `float` | nearest-neighbour, 250 nM primer, 50 mM Na+ |
| `meltingTempNnWith(s, primerNm, saltMm)` | `float` | nM, mM |
| `netCharge(p, ph)` | `float` | EMBOSS `iep` pKa set |
| `isoelectricPoint(p)` | `float` | bisection over pH 0 to 14 |
| `gravy(p)` | `float` | mean Kyte-Doolittle hydropathy |
| `aromaticity(p)` | `float` | fraction of F, W, Y |
| `dnaStats(s)` | `DnaStats` | |
| `proteinStats(p)` | `ProteinStats` | |
| `codonUsage(nt)` | `list of CodonUse` | standard code, first-seen order |
| `codonUsageWith(nt, table)` | `list of CodonUse` | |

### `DnaStats`

| Field | Meaning |
| - | - |
| `length as int` | bases |
| `gc as float` | GC percentage |
| `at as float` | AT percentage |
| `counts as map of string to int` | residue counts |
| `molecularWeight as float` | Da, 0.0 when ambiguous |
| `meltingTemp as float` | Wallace Tm, 0.0 when ambiguous |
| `ambiguous as bool` | any code other than A, C, G, T |

### `ProteinStats`

| Field | Meaning |
| - | - |
| `length as int` | residues |
| `molecularWeight as float` | Da, 0.0 when ambiguous |
| `isoelectricPoint as float` | 0.0 when ambiguous |
| `gravy as float` | 0.0 when ambiguous |
| `aromaticity as float` | fraction of F, W, Y |
| `counts as map of string to int` | residue counts |
| `ambiguous as bool` | carries `B`, `Z`, `X`, `J` or `*` |

### `CodonUse`

| Field | Meaning |
| - | - |
| `codon as string` | DNA-spelled codon |
| `aminoAcid as string` | its residue under the table, `X` if not in it |
| `count as int` | occurrences |
| `fraction as float` | share among synonymous codons, 0 to 1 |
| `perThousand as float` | per thousand codons |

## Search

| Function | Returns | Notes |
| - | - | - |
| `find(s, pattern)` | `list of Match` | exact, overlapping |
| `findIupac(s, pattern)` | `list of Match` | IUPAC codes in the pattern |
| `findFuzzy(s, pattern, maxMismatches)` | `list of Match` | substitutions only |
| `orfs(nt, minAa)` | `list of Orf` | six frames, standard code |
| `orfsWith(nt, minAa, tableId)` | `list of Orf` | NCBI table id |
| `cpgIslands(s, window, minGc, minObsExp)` | `list of Region` | classic: 200, 50.0, 0.6 |
| `restrictionSites(s, names)` | `list of Site` | ordered by cut, then enzyme |
| `enzymesCutting(s)` | `list of string` | enzymes with a site |
| `nonCutters(s)` | `list of string` | enzymes without one |
| `restrictionDigest(s, names, circular)` | `list of Fragment` | cuts at `cutTop` |

### `Match`

| Field | Meaning |
| - | - |
| `start as int` | 0-based offset of the first base |
| `end as int` | one past the last base |
| `text as string` | the matched subsequence |
| `strand as int` | +1, or -1 on the reverse complement |
| `mismatches as int` | 0 except for `findFuzzy` |

### `Orf`

| Field | Meaning |
| - | - |
| `frame as int` | 1, 2, 3 forward; -1, -2, -3 reverse |
| `strand as int` | sign of `frame` |
| `start as int` | first base of the start codon, forward coordinates |
| `end as int` | one past the last base, stop included |
| `length as int` | `end - start` |
| `nucleotides as string` | the frame as read, 5' to 3' |
| `protein as string` | translation without the stop |
| `stopped as bool` | false when it ran off the end |

### `Region`

| Field | Meaning |
| - | - |
| `start as int` | 0-based offset |
| `end as int` | one past the last base |
| `length as int` | bases |
| `gc as float` | GC percentage |
| `obsExp as float` | observed over expected CpG |

### `Site`

| Field | Meaning |
| - | - |
| `enzyme as string` | enzyme name |
| `start as int` | first base of the recognition site |
| `end as int` | one past the site |
| `cut as int` | top-strand cut offset |
| `site as string` | concrete sequence found |
| `strand as int` | +1, or -1 if found only on the reverse strand |

### `Fragment`

| Field | Meaning |
| - | - |
| `start as int` | first base in the original sequence |
| `end as int` | one past the last; not greater than `start` for a fragment through a circular origin |
| `length as int` | bases |
| `sequence as string` | the fragment |

## Formats

| Function | Returns | Notes |
| - | - | - |
| `parseFasta(text)` | `list of Record` | case preserved |
| `writeFasta(records)` | `string` | wrapped at `FASTA_WIDTH` |
| `writeFastaWrapped(records, width)` | `string` | |
| `wrap(s, width)` | `string` | no trailing newline |
| `toPlain(text)` | `string` | strip headers, digits, whitespace |
| `parseFastq(text)` | `list of Read` | four-line or wrapped |
| `writeFastq(reads)` | `string` | four-line form |
| `fastqToFasta(reads)` | `list of Record` | |
| `qualityScores(r)` | `list of int` | Phred+33 |
| `qualityScoresWith(r, offset)` | `list of int` | `PHRED33` or `PHRED64` |
| `meanQuality(r)` | `float` | arithmetic mean, Phred+33 |
| `meanQualityWith(r, offset)` | `float` | |
| `errorProbability(score)` | `float` | `10^(-Q/10)` |

| Constant | Value | Meaning |
| - | - | - |
| `FASTA_WIDTH` | `60` | default line width |
| `MAX_RECORDS` | `1000000` | records per parse |
| `PHRED33` | `33` | Sanger, Illumina 1.8+ |
| `PHRED64` | `64` | Illumina 1.3 to 1.7 |

### `Record`

| Field | Meaning |
| - | - |
| `id as string` | first token after `>` |
| `description as string` | rest of the header, or "" |
| `sequence as string` | joined sequence lines |

### `Read`

| Field | Meaning |
| - | - |
| `id as string` | first token after `@` |
| `description as string` | rest of the header, or "" |
| `sequence as string` | the bases |
| `quality as string` | encoded quality, same length as `sequence` |

## Manipulation

| Function | Returns | Notes |
| - | - | - |
| `randomDna(n)` | `string` | uniform over `DNA_BASES` |
| `randomRna(n)` | `string` | uniform over `RNA_BASES` |
| `randomProtein(n)` | `string` | uniform over `STANDARD_RESIDUES` |
| `randomDnaWith(n, weights)` | `string` | relative weights |
| `randomRnaWith(n, weights)` | `string` | |
| `randomProteinWith(n, weights)` | `string` | |
| `shuffle(s)` | `string` | Fisher-Yates, composition kept |
| `mutate(s, rate, kind)` | `string` | substitutions, a mutated position always changes |
| `subrange(s, from, through)` | `string` | **1-based inclusive** |
| `windows(s, size, step)` | `list of string` | full windows only |

| Constant | Value | Meaning |
| - | - | - |
| `DNA_BASES` | `"ACGT"` | generator alphabet |
| `RNA_BASES` | `"ACGU"` | |
| `STANDARD_RESIDUES` | `"ACDEFGHIKLMNPQRSTVWY"` | |

## Comparison

| Function | Returns | Notes |
| - | - | - |
| `alignGlobal(a, b, opts)` | `Alignment` | Needleman-Wunsch |
| `alignLocal(a, b, opts)` | `Alignment` | Smith-Waterman |
| `dnaOptions()` | `AlignOptions` | +2 / -1, gap -2 |
| `proteinOptions()` | `AlignOptions` | `BLOSUM62`, gap -8 |
| `substitutionScore(matrix, a, b)` | `int` | either pair order |
| `percentIdentity(al)` | `float` | over alignment length |
| `percentSimilarity(al)` | `float` | |
| `percentGaps(al)` | `float` | |
| `formatAlignment(al, width)` | `string` | EMBOSS-style blocks |
| `hammingDistance(a, b)` | `int` | equal lengths only |
| `editDistance(a, b)` | `int` | Levenshtein |

| Constant | Value | Meaning |
| - | - | - |
| `DEFAULT_MAX_ALIGNMENT_CELLS` | `1000000` | default `maxCells`, about 1 GB |

### `AlignOptions`

| Field | Meaning |
| - | - |
| `matchScore as int` | identical residues, when `matrix` is empty |
| `mismatchScore as int` | differing residues, when `matrix` is empty |
| `gapScore as int` | per gap position |
| `matrix as map of string to int` | substitution matrix, or empty |
| `maxCells as int` | DP matrix limit; 0 for none |

### `Alignment`

| Field | Meaning |
| - | - |
| `a as string` | first sequence, gapped |
| `b as string` | second sequence, gapped |
| `markup as string` | `\|` identical, `:` similar, space otherwise |
| `score as int` | total score |
| `length as int` | columns |
| `identities as int` | identical columns |
| `similarities as int` | positively scoring columns, identities included |
| `gaps as int` | gap columns |
| `startA as int` | 0-based start in the first sequence |
| `startB as int` | 0-based start in the second |

## Reference data

| Function | Returns | Notes |
| - | - | - |
| `codonTable(id)` | `map of string to string` | NCBI `transl_table` |
| `codonTableName(id)` | `string` | |
| `codonTableIds()` | `list of int` | 1-6, 9-11 |
| `startCodons(id)` | `list of string` | |
| `expandIupac(code)` | `string` | `R` to `AG` |
| `iupacClass(code)` | `string` | `R` to `[AG]` |
| `matchesIupac(code, base)` | `bool` | one-way |
| `aminoAcid(code)` | `AminoAcid` | one- or three-letter code |
| `aminoAcids()` | `list of AminoAcid` | 22 records |
| `hasAminoAcid(code)` | `bool` | |
| `enzyme(name)` | `Enzyme` | case-insensitive |
| `hasEnzyme(name)` | `bool` | |
| `enzymes()` | `list of Enzyme` | alphabetical |
| `enzymeNames()` | `list of string` | alphabetical |
| `isoschizomers(name)` | `list of string` | same site string |
| `isBlunt(e)` | `bool` | `cutTop == cutBottom` |

| Constant | Value | Meaning |
| - | - | - |
| `PKA_N_TERMINUS` | `8.6` | EMBOSS `iep` |
| `PKA_C_TERMINUS` | `3.6` | EMBOSS `iep` |
| `WATER_MONO_MASS` | `18.010565` | Da |
| `WATER_AVG_MASS` | `18.01528` | Da |
| `PHOSPHATE_MASS` | `79.98` | Da, HPO3 |
| `GAS_CONSTANT` | `1.987` | cal/(mol K) |

Map constants:

| Constant | Type | Value |
| - | - | - |
| `STANDARD` | `map of string to string` | NCBI code 1 |
| `VERTEBRATE_MITOCHONDRIAL` | `map of string to string` | code 2 |
| `YEAST_MITOCHONDRIAL` | `map of string to string` | code 3 |
| `MOLD_MITOCHONDRIAL` | `map of string to string` | code 4 |
| `INVERTEBRATE_MITOCHONDRIAL` | `map of string to string` | code 5 |
| `CILIATE_NUCLEAR` | `map of string to string` | code 6 |
| `ECHINODERM_MITOCHONDRIAL` | `map of string to string` | code 9 |
| `EUPLOTID_NUCLEAR` | `map of string to string` | code 10 |
| `BACTERIAL` | `map of string to string` | code 11 |
| `IUPAC_DNA` | `map of string to string` | code to bases |
| `IUPAC_RNA` | `map of string to string` | code to bases |
| `COMPLEMENT_DNA` | `map of string to string` | code to complement |
| `COMPLEMENT_RNA` | `map of string to string` | code to complement |
| `THREE_LETTER` | `map of string to string` | one-letter to three-letter |
| `DNA_RESIDUE_MASS` | `map of string to float` | average residue mass, Da |
| `RNA_RESIDUE_MASS` | `map of string to float` | average residue mass, Da |
| `NN_ENTHALPY` | `map of string to float` | dinucleotide dH, kcal/mol |
| `NN_ENTROPY` | `map of string to float` | dinucleotide dS, cal/(mol K) |
| `BLOSUM62` | `map of string to int` | pair to score, one order per pair |

### `AminoAcid`

| Field | Meaning |
| - | - |
| `one as string` | one-letter code |
| `three as string` | three-letter code |
| `name as string` | full name |
| `monoMass as float` | monoisotopic residue mass, Da |
| `avgMass as float` | average residue mass, Da |
| `pKa as float` | side chain, 0.0 if it does not ionise |
| `charge as int` | +1, -1 or 0 |
| `hydropathy as float` | Kyte-Doolittle |

### `Enzyme`

| Field | Meaning |
| - | - |
| `name as string` | conventional spelling |
| `site as string` | recognition sequence, IUPAC |
| `cutTop as int` | top-strand cut offset into the site |
| `cutBottom as int` | bottom-strand cut, counted along the top strand |
| `overhang as string` | `"5'"`, `"3'"` or `"blunt"` |
| `overhangLength as int` | 0 when blunt |
