# Cheatsheet

The whole surface on one page. Every name is reached as `bioseq.NAME`.

```jennifer
import "@jennifer/bioinformatics/" as bioseq;
```

## Alphabets and validation

| | |
| - | - |
| `filterDna(s)` `filterRna(s)` `filterProtein(s)` | strip everything that is not a residue |
| `isDna(s)` `isRna(s)` `isProtein(s)` | non-empty and wholly in the alphabet |
| `kindName(kind)` | `"DNA"` / `"RNA"` / `"protein"` |
| `DNA` `RNA` `PROTEIN` | molecule kinds |
| `DNA_ALPHABET` `RNA_ALPHABET` `PROTEIN_ALPHABET` | the admitted residues |

## Transforms

| | |
| - | - |
| `complement(s)` | pair each base, keep the order |
| `reverseComplement(s)` | the other strand, 5′ to 3′ |
| `reverse(s)` | reverse the residues |
| `transcribe(dna)` `backTranscribe(rna)` | T ↔ U |
| `translate(nt)` | standard code |
| `translateWith(nt, table)` | named code |
| `translateFrame(nt, frame)` | frame 1..3 or -1..-3 |
| `translateFrameWith(nt, frame, table)` | |
| `sixFrame(nt)` `sixFrameWith(nt, table)` | → `list of Frame` |
| `splitCodons(nt)` | → `list of string` |
| `oneToThree(p)` `threeToOne(p)` | residue spellings |
| `Frame` | `frame`, `protein` |

## Composition and properties

| | |
| - | - |
| `baseCount(s)` `baseFrequency(s)` | nucleotide composition |
| `residueCount(p)` `residueFrequency(p)` | protein composition |
| `gcContent(s)` `atContent(s)` | percentages |
| `gcWindows(s, size, step)` | → `list of float` |
| `molecularWeight(s, kind)` | daltons |
| `meltingTemp(s)` | Wallace rule |
| `meltingTempNn(s)` | nearest-neighbour, 250 nM / 50 mM |
| `meltingTempNnWith(s, primerNm, saltMm)` | |
| `netCharge(p, ph)` `isoelectricPoint(p)` | |
| `gravy(p)` `aromaticity(p)` | |
| `dnaStats(s)` | → `DnaStats` |
| `proteinStats(p)` | → `ProteinStats` |
| `codonUsage(nt)` `codonUsageWith(nt, table)` | → `list of CodonUse` |
| `DnaStats` | `length`, `gc`, `at`, `counts`, `molecularWeight`, `meltingTemp`, `ambiguous` |
| `ProteinStats` | `length`, `molecularWeight`, `isoelectricPoint`, `gravy`, `aromaticity`, `counts`, `ambiguous` |
| `CodonUse` | `codon`, `aminoAcid`, `count`, `fraction`, `perThousand` |

## Search

| | |
| - | - |
| `find(s, pattern)` | exact, overlapping |
| `findIupac(s, pattern)` | ambiguity codes in the pattern |
| `findFuzzy(s, pattern, maxMismatches)` | substitutions only |
| `orfs(nt, minAa)` `orfsWith(nt, minAa, tableId)` | → `list of Orf` |
| `cpgIslands(s, window, minGc, minObsExp)` | → `list of Region` |
| `restrictionSites(s, names)` | → `list of Site` |
| `enzymesCutting(s)` `nonCutters(s)` | → `list of string` |
| `restrictionDigest(s, names, circular)` | → `list of Fragment` |
| `Match` | `start`, `end`, `text`, `strand`, `mismatches` |
| `Orf` | `frame`, `strand`, `start`, `end`, `length`, `nucleotides`, `protein`, `stopped` |
| `Region` | `start`, `end`, `length`, `gc`, `obsExp` |
| `Site` | `enzyme`, `start`, `end`, `cut`, `site`, `strand` |
| `Fragment` | `start`, `end`, `length`, `sequence` |

## Formats

| | |
| - | - |
| `parseFasta(text)` | → `list of Record` |
| `writeFasta(records)` | wrapped at 60 |
| `writeFastaWrapped(records, width)` | |
| `wrap(s, width)` `toPlain(text)` | |
| `parseFastq(text)` | → `list of Read` |
| `writeFastq(reads)` `fastqToFasta(reads)` | |
| `qualityScores(r)` `qualityScoresWith(r, offset)` | → `list of int` |
| `meanQuality(r)` `meanQualityWith(r, offset)` | |
| `errorProbability(score)` | `10^(-Q/10)` |
| `Record` | `id`, `description`, `sequence` |
| `Read` | `id`, `description`, `sequence`, `quality` |
| `FASTA_WIDTH` `PHRED33` `PHRED64` | constants |
| `MAX_RECORDS` | parser record cap |

## Manipulation

| | |
| - | - |
| `randomDna(n)` `randomRna(n)` `randomProtein(n)` | uniform |
| `randomDnaWith(n, weights)` `randomRnaWith(n, weights)` `randomProteinWith(n, weights)` | weighted |
| `shuffle(s)` | composition-preserving |
| `mutate(s, rate, kind)` | point substitutions |
| `subrange(s, from, through)` | **1-based inclusive** |
| `windows(s, size, step)` | → `list of string` |
| `DNA_BASES` `RNA_BASES` `STANDARD_RESIDUES` | generator alphabets |

## Comparison

| | |
| - | - |
| `alignGlobal(a, b, opts)` | Needleman-Wunsch |
| `alignLocal(a, b, opts)` | Smith-Waterman |
| `dnaOptions()` `proteinOptions()` | scoring defaults |
| `substitutionScore(matrix, a, b)` | order-independent |
| `percentIdentity(al)` `percentSimilarity(al)` `percentGaps(al)` | |
| `formatAlignment(al, width)` | numbered blocks |
| `hammingDistance(a, b)` `editDistance(a, b)` | |
| `AlignOptions` | `matchScore`, `mismatchScore`, `gapScore`, `matrix`, `maxCells` |
| `Alignment` | `a`, `b`, `markup`, `score`, `length`, `identities`, `similarities`, `gaps`, `startA`, `startB` |
| `BLOSUM62` | the substitution matrix |
| `DEFAULT_MAX_ALIGNMENT_CELLS` | the matrix-size guard, one million cells |

## Reference data

| | |
| - | - |
| `codonTable(id)` `codonTableName(id)` `codonTableIds()` | genetic codes |
| `startCodons(id)` | per-table start codons |
| `STANDARD` `VERTEBRATE_MITOCHONDRIAL` `YEAST_MITOCHONDRIAL` | codon tables |
| `MOLD_MITOCHONDRIAL` `INVERTEBRATE_MITOCHONDRIAL` `CILIATE_NUCLEAR` | |
| `ECHINODERM_MITOCHONDRIAL` `EUPLOTID_NUCLEAR` `BACTERIAL` | |
| `expandIupac(code)` `iupacClass(code)` `matchesIupac(code, base)` | ambiguity |
| `IUPAC_DNA` `IUPAC_RNA` | code → bases |
| `COMPLEMENT_DNA` `COMPLEMENT_RNA` | code → complement |
| `aminoAcid(code)` `aminoAcids()` `hasAminoAcid(code)` | residue properties |
| `AminoAcid` | `one`, `three`, `name`, `monoMass`, `avgMass`, `pKa`, `charge`, `hydropathy` |
| `THREE_LETTER` | one-letter → three-letter |
| `enzyme(name)` `hasEnzyme(name)` `enzymes()` `enzymeNames()` | enzymes |
| `isoschizomers(name)` `isBlunt(e)` | |
| `Enzyme` | `name`, `site`, `cutTop`, `cutBottom`, `overhang`, `overhangLength` |

### Parameter constants

`PKA_N_TERMINUS` `PKA_C_TERMINUS` `WATER_MONO_MASS` `WATER_AVG_MASS`
`DNA_RESIDUE_MASS` `RNA_RESIDUE_MASS` `PHOSPHATE_MASS` `NN_ENTHALPY`
`NN_ENTROPY` `GAS_CONSTANT`

## Conventions at a glance

- Sequences are **strings**; input any case, residues out uppercase.
- Offsets are **0-based half-open** - `$s[$m.start..$m.end]` is the match.
  `subrange` is the one 1-based inclusive exception.
- Reverse-strand features report **forward-sequence coordinates**.
- `filter*` strips silently; everything else raises on a non-residue, naming it.
- Ambiguity codes raise where a physical measurement is asked for.
- One error shape: `Error{kind: "bioinformatics"}`.
- A `With` suffix means an alternative table or offset, never a performance knob.
