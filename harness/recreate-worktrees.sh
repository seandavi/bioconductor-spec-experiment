#!/bin/sh
# Recreate the experiment worktrees the eval harness expects (eval/ports.toml),
# checked out at the frozen tags that were scored. Run from the repo root.
#   sh harness/recreate-worktrees.sh            # experiments only (for `just run <port>`)
#   sh harness/recreate-worktrees.sh --replay   # also replay/<exp> (for eval/replay/)
# eval/ports.toml has WORKTREES as an absolute path; edit it if your clone lives elsewhere.
set -eu
root=$(pwd)
wt="$root/../$(basename "$root").worktrees"
mkdir -p "$wt"
for e in e0-control e0b-no-r e1-differential e2-conformance e3-spec-adr e4-primitive-cases; do
  [ -d "$wt/$e" ] || git worktree add --detach "$wt/$e" "frozen/$e"
  if [ "${1:-}" = "--replay" ]; then
    [ -d "$wt/replay-$e" ] || git worktree add "$wt/replay-$e" "replay/$e"
  fi
done
git worktree list
