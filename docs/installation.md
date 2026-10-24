# Installation

## With jvc

```sh
jvc add @deckyard/bioinformatics
```

or declare it in your project's `deck.toml` and run `jvc install`:

```toml
[decks]
"@deckyard/bioinformatics" = ">=0.1.0"
```

`jvc` vendors the deck's `src/` to `vendor/deckyard/bioinformatics/`. Import it
by its scoped name; the trailing slash resolves to the entry module:

```jennifer
import "@deckyard/bioinformatics/" as bioseq;
```

## Without jvc

Unpack the `src/` subtree of a release archive
(`deckyard-bioinformatics-<version>.tar.gz`) into
`vendor/deckyard/bioinformatics/`, or import a checkout directly. There is no
build step:

```sh
git clone https://github.com/deckyard/deck-bioinformatics
```

```jennifer
import "./deck-bioinformatics/src/bioinformatics.j" as bioseq;
```

The examples in this repository import `../src/bioinformatics.j` the same way.

## Requirements

| | |
| - | - |
| Jennifer | `>=0.25.0` |
| Deck dependencies | none |
| Capabilities | none (no `net`, `exec` or filesystem) |
| Builds | `jennifer` and `jennifer-tiny`, both listed in `[engines]` |

The floor comes from the `math` functions the thermodynamic and pI formulas need
(`math.ln`, `math.pow`). The deck also uses `strings`, `maps`, `convert`,
`binary` and `io`.

## Checking it works

From a checkout:

```sh
jennifer test src/bioinformatics_test.j   # fixtures are embedded, no files needed
jennifer run examples/transforms.j        # needs no testdata/
```
