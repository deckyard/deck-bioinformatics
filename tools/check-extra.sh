#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
#
# tabledata.j and enzymedata.j are generated. Regenerate them here and fail if
# the result differs from what is committed, so the two cannot drift apart
# unnoticed.
#
# Usage: tools/check-extra.sh
set -euo pipefail
cd "$(dirname "$0")/.."

JENNIFER="${JENNIFER:-jennifer}"
generated=(src/topics/tabledata.j src/topics/enzymedata.j)

$JENNIFER run tools/gentables.j
$JENNIFER run tools/genenzymes.j
$JENNIFER fmt -w "${generated[@]}"

if ! git diff --quiet -- "${generated[@]}"; then
    echo "::error::generated data differs from its generators - re-run tools/ and commit" >&2
    git --no-pager diff --stat -- "${generated[@]}" >&2
    exit 1
fi
echo "generated data is in step with its generators"
