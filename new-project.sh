#!/usr/bin/env bash
# Scaffold a project folder with the four files.
#
#   ./bin/new-project.sh ~/projects/2026_MyProject
#   ./bin/new-project.sh ~/projects/2026_MyProject --name "Clarkia demography" --code cx-demo
#
# Never overwrites an existing file: it reports and skips. Safe to re-run on a folder that
# already has some of the four, to fill in the rest.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET=""; NAME=""; CODE=""

while [ "$#" -gt 0 ]; do
  case "$1" in
    --name) NAME="${2:-}"; shift 2 ;;
    --code) CODE="${2:-}"; shift 2 ;;
    -h|--help) sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) TARGET="$1"; shift ;;
  esac
done

if [ -z "$TARGET" ]; then
  echo "usage: new-project.sh <folder> [--name \"Full name\"] [--code handle]" >&2
  exit 2
fi

mkdir -p "$TARGET"
TARGET="$(cd "$TARGET" && pwd)"
[ -n "$NAME" ] || NAME="$(basename "$TARGET" | sed 's/^[0-9]\{4\}[_-]//')"
[ -n "$CODE" ] || CODE="$(printf '%s' "$NAME" | tr '[:upper:] ' '[:lower:]-')"

wrote=0; skipped=0
for f in CLAUDE.md PROJECT_INDEX.md TODO.md; do
  if [ -e "$TARGET/$f" ]; then
    echo "  skip   $f (already exists)"; skipped=$((skipped+1)); continue
  fi
  sed -e "s|<Project name>|$NAME|g" -e "s|<short-handle>|$CODE|g" \
      "$REPO/templates/$f" > "$TARGET/$f"
  echo "  wrote  $f"; wrote=$((wrote+1))
done

echo
echo "$TARGET"
echo "  $wrote written, $skipped skipped"
echo
echo "handoff.md is written by the agent at the end of a session, not now."
echo "Next: open CLAUDE.md and fill it in. That file is most of the value."
