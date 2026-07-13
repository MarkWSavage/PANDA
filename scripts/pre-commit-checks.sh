#!/usr/bin/env bash
set -euo pipefail

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

mapfile -t staged_py < <(git diff --cached --name-only --diff-filter=ACMR -- '*.py')
mapfile -t staged_cc < <(git diff --cached --name-only --diff-filter=ACMR -- '*.cc')

exit_code=0

if [ "${#staged_py[@]}" -gt 0 ]; then
    if [ -x "$repo_root/.lint-venv/bin/ruff" ]; then
        echo "== ruff (blocking) =="
        if ! "$repo_root/.lint-venv/bin/ruff" check "${staged_py[@]}"; then
            exit_code=1
        fi
    else
        echo "warning: .lint-venv/bin/ruff not found -- run: python3 -m venv .lint-venv && .lint-venv/bin/pip install ruff" >&2
    fi
fi

if [ "${#staged_cc[@]}" -gt 0 ]; then
    if [ -f "$repo_root/build/compile_commands.json" ]; then
        echo "== clang-tidy (warnings only, non-blocking) =="
        clang-tidy -p "$repo_root/build" --header-filter='include/' "${staged_cc[@]}" || true
    else
        echo "warning: build/compile_commands.json not found -- skipping clang-tidy (run: cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -S . -B build)" >&2
    fi
fi

exit "$exit_code"
