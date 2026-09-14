#!/usr/bin/env bash
# Preserve IQ-TREE's exit status in an explicit outcome, not an empty Newick.
set -euo pipefail
if (( $# < 8 )); then
    echo "Usage: run_iqtree_with_fallback.sh OG TIMEOUT_SECONDS KILL_AFTER_SECONDS SOURCE_TREE FINAL_TREE FALLBACK_JSON -- COMMAND..." >&2
    exit 125
fi
og=$1
budget=$2
kill_after=$3
source_tree=$4
final_tree=$5
fallback=$6
shift 6
[[ $1 == -- ]] || exit 125
shift
[[ $og =~ ^[A-Za-z0-9_.-]+$ ]] || exit 125
command -v timeout >/dev/null || exit 125
command -v "$1" >/dev/null || exit 127
started=$SECONDS
status=0
timeout --signal=TERM --kill-after="$kill_after" "$budget" "$@" || status=$?
if (( status == 0 )) && [[ -s "$source_tree" ]]; then
    mv "$source_tree" "$final_tree"
    exit 0
fi
case "$status" in
    125|126|127) exit "$status" ;; # supervisor/executable setup error
    124) reason=timeout ;;
    137) reason=killed ;; # SIGKILL is not proof of OOM or of timeout
    0) reason=missing_tree ;;
    *) reason=tool_exit ;;
esac
printf '{"orthogroup":"%s","reason":"%s","exit_status":%s,"elapsed_seconds":%s}\n' \
    "$og" "$reason" "$status" "$((SECONDS - started))" > "$fallback"
echo "IQ-TREE fallback: $og reason=$reason exit=$status" >&2
