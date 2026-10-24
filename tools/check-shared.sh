#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
#
# Guard the shared CI layer against local edits.
#
# The workflows and most of tools/ are written by the template's
# tools/sync-deck.sh and must stay identical across every deck - that is what
# lets one fix reach all of them, and what keeps them from drifting apart again.
# The sync records a sha256 per shared file in tools/shared.lock; this fails
# when a file no longer matches it.
#
# A deck that needs one more gate puts it in tools/check-extra.sh, or in a
# workflow of its own, never in a shared file. A shared file that is wrong is
# fixed in the template and synced to every deck.
#
# The template itself has no lock - it is where the shared files come from.
#
# Usage: tools/check-shared.sh
set -euo pipefail
cd "$(dirname "$0")/.."

lock=tools/shared.lock

if [[ ! -f "$lock" ]]; then
    if [[ -f tools/sync-deck.sh ]]; then
        echo "this is the template - it is the source of the shared files"
        exit 0
    fi
    echo "error: no ${lock} - this deck has never been synced with the template." >&2
    echo "From the template, run: tools/sync-deck.sh <path to this deck>" >&2
    exit 1
fi

source_line="$(grep -m1 '^# source: ' "$lock" | sed 's/^# source: //' || true)"
report="$(grep -v '^#' "$lock" | sha256sum --check --quiet - 2>&1 || true)"

if [[ -n "$report" ]]; then
    echo "error: shared files differ from the template they were synced from:" >&2
    sed -E 's/^/  /; s/: FAILED open or read$/: missing/; s/: FAILED$/: changed/' <<<"$report" \
        | grep -vE 'WARNING|^  sha256sum: ' >&2
    echo "" >&2
    echo "These files are shared by every deck. Undo the local change, and either" >&2
    echo "move deck-specific logic to tools/check-extra.sh or fix the file in the" >&2
    echo "template and sync it (synced from: ${source_line:-unknown})." >&2
    exit 1
fi
echo "all $(grep -vc '^#' "$lock") shared files match the template (${source_line:-unknown})"
