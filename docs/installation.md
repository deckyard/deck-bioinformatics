# Installation

`@jennifer/bioinformatics` is a **deck**: a versioned bundle of Jennifer
modules, vendored into the project that uses it rather than installed globally.

## With jvc

Once [jvc](https://jennifer-lang.dev/) is available, declare the dependency in
your project's `deck.toml`:

```toml
[decks]
"@jennifer/bioinformatics" = ">=0.1.0"
```

then

```sh
jvc install
```

which writes the deck into the project-local `vendor/` tree. Import it by its
scoped name:

```jennifer
import "@jennifer/bioinformatics/" as bioseq;
```

The trailing slash expands to the package-named entry point. The alias is not
optional dressing: the package is named for the field, but the calls are named
for the thing they act on, and `bioseq.reverseComplement($s)` reads better in a
pipeline than `bioinformatics.reverseComplement($s)`. Every example in this
documentation binds `bioseq`, so keeping that alias makes the snippets
copy-paste directly. Import it under another name if you prefer - nothing in
the deck depends on what the caller calls it.

## Without jvc

The deck is plain Jennifer source with no build step, so a checkout works
directly:

```sh
git clone https://github.com/jennifer-language/deck-bioinformatics
```

```jennifer
import "./deck-bioinformatics/src/bioinformatics.j" as bioseq;
```

That is how the examples in this repository import it, so they run against a
checkout with nothing installed.

## Requirements

| | |
| - | - |
| Jennifer | `>=0.25.0` |
| Deck dependencies | none |
| Capabilities | none - no `net`, no `exec` |
| Builds | `jennifer` and `jennifer-tiny`, both listed in `[engines]` |

The engine floor comes from the mathematical foundations the thermodynamic and
isoelectric-point formulas need (`exp`, `ln`, `log`). Everything else the deck
uses - `strings`, `lists`, `maps`, `math`, `regex`, `convert` - has been in the
standard library far longer.

## Checking it works

From a checkout:

```sh
jennifer test src/bioinformatics_test.j     # the white-box suite
jennifer run examples/transforms.j  # the front door, no files needed
```

The suite is network-free, subprocess-free and filesystem-free, with every
fixture embedded in the test files, so it runs unchanged under `jennifer-tiny`.
