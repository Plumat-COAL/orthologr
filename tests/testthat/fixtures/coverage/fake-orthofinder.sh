#!/bin/sh
# A stand-in for OrthoFinder (POSIX sh), used by test-coverage-tools.R to exercise orthofinder2()'s run path.
#   orthofinder -h                       prints a version line
#   orthofinder -f DIR ... -n NAME ...   creates DIR/OrthoFinder/Results_NAME and records its arguments there
# FAKE_OF_MODE: ok (default) | noresults (exit 0 without a results folder) | fail (exit status 2)
mode="${FAKE_OF_MODE:-ok}"
if [ "$1" = "-h" ]; then
        echo "OrthoFinder version 2.5.5 (fake) Copyright (C) 2014 David Emms"
        exit 0
fi
if [ "$mode" = "fail" ]; then
        echo "fake OrthoFinder: failed on purpose" >&2
        exit 2
fi
args="$*"
dir=""
name=""
while [ $# -gt 0 ]; do
        case "$1" in
                -f) dir="$2"; shift ;;
                -n) name="$2"; shift ;;
        esac
        shift
done
if [ "$mode" = "noresults" ]; then
        exit 0
fi
mkdir -p "$dir/OrthoFinder/Results_$name"
echo "$args" > "$dir/OrthoFinder/Results_$name/args.txt"
exit 0
