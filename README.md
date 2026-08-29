# @jennifer/bioinformatics

**Sequence manipulation for DNA, RNA and protein**, in the
[Jennifer](https://jennifer-lang.dev/) language.

Transforms, composition and physical properties, pattern and feature search,
FASTA and FASTQ, sequence generation, and pairwise alignment.

```jennifer
import "@jennifer/bioinformatics/" as bioseq;

def gene as string init bioseq.filterDna("1 atg acg tta cat gga taa");

bioseq.reverseComplement($gene);            # TTATCCATGTAACGTCAT
bioseq.translate($gene);                    # MTLHG*
bioseq.gcContent($gene);                    # 33.33
bioseq.molecularWeight($gene, bioseq.DNA);  # 5546.66
```

Sequences are plain strings: no `Seq` type to construct, no wrapper to unwrap.
The package is named for the field and the alias for what the calls act on, so
`bioseq.reverseComplement($s)` reads better in a pipeline than
`bioinformatics.reverseComplement($s)`. Nothing in the deck depends on the name
the caller gives it.

The deck runs on `jennifer-tiny`. Nothing under `src/` opens a socket, spawns a
process or touches the filesystem, so it needs no capability and `[engines]`
lists both builds.

## Documentation

[**jennifer-language.github.io/deck-bioinformatics**](https://jennifer-language.github.io/deck-bioinformatics/)
carries a guide per topic, a cheatsheet of the whole surface, and the design
notes.

The source is [`docs/`](docs/), built with
[Grimoire](https://github.com/jennifer-language/grimoire) and published on every
push to `main`. To read it locally:

```
git clone https://github.com/jennifer-language/grimoire .grimoire
jennifer run .grimoire/bin/grimoire serve
```

## Development

```
jennifer test src/bioinformatics_test.j # the white-box suite
jennifer run examples/compare.j         # global and local alignment, distances
jennifer run examples/search.j          # patterns, ORFs, CpG islands, digest
jennifer run examples/manipulate.j      # random sequences, shuffle, mutate, windows
jennifer run examples/composition.j     # GC, mass, Tm, pI, GRAVY, codon usage
jennifer run examples/formats.j         # FASTA / FASTQ over testdata/
jennifer run examples/transforms.j      # complement, transcribe, translate
jennifer run examples/tables.j          # the reference data the deck carries
jennifer run examples/bench.j           # the measurements behind the conventions

jennifer fmt -l src/topics/*.j          # formatting check
jennifer lint src/bioinformatics.j      # the module, which splices every topic
jennifer lint src/bioinformatics_test.j # the test overlay, likewise
jennifer lint examples/*.j tools/*.j

jennifer run tools/gentables.j          # regenerate src/topics/tabledata.j
jennifer run tools/genenzymes.j         # regenerate src/topics/enzymedata.j
```

Lint reaches every line through four compilable units. A topic file cannot be
linted alone, being `include`d rather than compiled, so the module and the test
overlay are what cover it through the splice.

Two topic files are generated, see [`tools/`](tools/README.md). Generated data
and hand-written accessors live in separate files (`tabledata.j` vs `tables.j`),
so regenerating never clobbers API. The measurements in `examples/bench.j`
settled the implementation conventions and are kept so a regression shows up as
a number; the story is in [Design notes](docs/design.md).

## Continuous integration

| workflow | when | what |
| --- | --- | --- |
| [`test.yml`](.github/workflows/test.yml) | push, PR | fmt and lint over `src/`, `examples/` and `tools/`; the suite; every example; the generators; then every example again on `jennifer-tiny` |
| [`docs.yml`](.github/workflows/docs.yml) | push to `main` | build `docs/` with Grimoire, publish to Pages |
| [`release.yml`](.github/workflows/release.yml) | version tag | re-run the gate, check `deck.toml` matches the tag, cut the release |

## License

LGPL-3.0-only, see [LICENSE.md](LICENSE.md).
