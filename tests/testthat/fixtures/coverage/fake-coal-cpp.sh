#!/bin/sh
# A stand-in for COAL's 'coal-cpp' program (POSIX sh), used by test-coverage-coal.R to exercise orthologr's
# "binary" COAL backend without COAL. It understands the three command lines orthologr uses:
#   coal-cpp --version
#   coal-cpp dnds --help
#   coal-cpp dnds AXT -o OUT -m METHOD ...      (writes a one-row result table)
#   coal-cpp ALN NUC -output FORMAT -o OUT      (PAL2NAL's command line; writes NUC to OUT)
# The environment variable FAKE_COAL_MODE chooses what it does:
#   ok (default) | nonative (help without --comeron-flavour and without GY) | nogy (help without GY) |
#   error (row with status 'error') | notconv (row with status 'not_converged') |
#   empty (header only) | fail (exit status 3) | nooutput (exit 0, no output file) | undefined (Comeron note)
mode="${FAKE_COAL_MODE:-ok}"
if [ "$1" = "--version" ]; then
        echo "coal-cpp 0.0.0-fake"
        exit 0
fi
if [ "$1" = "dnds" ] && [ "$2" = "--help" ]; then
        echo "usage: coal-cpp dnds AXT -o OUT [-m METHOD]"
        if [ "$mode" = "nonative" ] || [ "$mode" = "nogy" ]; then
                echo "  -m METHOD      YN (default), NG, LWL, MLWL, LPB, MLPB, MYN; repeat for several"
        else
                echo "  -m METHOD      YN (default), NG, LWL, MLWL, LPB, MLPB, MYN, GY, Comeron (alias C95), Li; repeat for"
                echo "                 several (GY: maximum likelihood, slower)"
        fi
        if [ "$mode" != "nonative" ]; then
                echo "  --comeron-flavour {libsequence,gestimator}"
        fi
        exit 0
fi
if [ "$mode" = "fail" ]; then
        echo "fake coal-cpp: failed on purpose" >&2
        exit 3
fi
out=""
method="YN"
cmd="$1"
in1="$1"
in2="$2"
while [ $# -gt 0 ]; do
        case "$1" in
                -o) out="$2"; shift ;;
                -m) method="$2"; shift ;;
        esac
        shift
done
if [ "$mode" = "nooutput" ]; then
        exit 0
fi
if [ "$cmd" = "dnds" ]; then
        printf 'name\tmethod\tdN\tdS\tdNdS\tstatus\tnote\n' > "$out"
        case "$mode" in
                empty) ;;
                error) printf 'pair\t%s\tNA\tNA\tNA\terror\tsequences of unequal length\n' "$method" >> "$out" ;;
                notconv) printf 'pair\t%s\tNA\tNA\tNA\tnot_converged\tthe iteration did not converge\n' "$method" >> "$out" ;;
                undefined) printf 'pair\t%s\t0.123456\t1.5\t0.082304\tok\tKs: correction undefined\n' "$method" >> "$out" ;;
                *) printf 'pair\t%s\t0.1\t0.5\t0.2\tok\t\n' "$method" >> "$out" ;;
        esac
        exit 0
fi
cat "$in2" > "$out"
exit 0
