#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

ipath="/data/2026-testbeam/actual/physics"
opath="/data/2026-testbeam/process"

DCALIB="${ROOT_DIR}/process/bin/dcalib"
MIN_PAIRS=1000
WRITE_LOGS=0

run=""
period=""

usage()
{
    cat <<EOF_USAGE
usage:
  $0 --run RUN --period PERIOD

required:
  --run, -r RUN          run name/directory
  --period, -p PERIOD    expected pulser period in coarse clock cycles

options:
  --help, -h             show this help message

example:
  $0 --run 12345 --period 10000
EOF_USAGE
}

fail()
{
    echo "ERROR: $*" >&2
    echo >&2
    usage >&2
    exit 1
}

run_job()
{
    local logfile=$1
    shift

    if [ "${WRITE_LOGS}" -eq 1 ]; then
        ( "$@" ) > "${logfile}" 2>&1
    else
        ( "$@" )
    fi
}

while [ $# -gt 0 ]; do
    case "$1" in
        --run|-r)
            [ $# -ge 2 ] || fail "$1 requires RUN"
            run=$2
            shift 2
            ;;
        --period|-p)
            [ $# -ge 2 ] || fail "$1 requires PERIOD"
            period=$2
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            fail "unknown option: $1"
            ;;
    esac
done

[ -n "${run}" ] || fail "missing --run"
[ -n "${period}" ] || fail "missing --period"
[[ "${period}" =~ ^[0-9]+$ ]] || fail "--period must be a positive integer"
[ "${period}" -gt 0 ] || fail "--period must be greater than zero"

if [ ! -x "${DCALIB}" ]; then
    fail "${DCALIB} does not exist or is not executable"
fi

irpath="${ipath}/${run}"
if [ ! -d "${irpath}" ]; then
   fail "${irpath} does not exist"
fi
orpath="${opath}/${run}"
mkdir -p "${orpath}"

pids=()
for idpath in "${irpath}"/kc705* "${irpath}"/rdo*; do

    [ -d "${idpath}" ] || continue

    device=$(basename "${idpath}")
    odpath="${orpath}/${device}"
    mkdir -p "${odpath}"
    echo " --- TDC calibration for device ${device}: ${idpath} "

    decoded_files=("${idpath}"/decoded/alcdaq.fifo_*.root)
    if [ ${#decoded_files[@]} -eq 0 ]; then
        echo " --- no decoded files found for ${device} "
        continue
    fi

    for fpath in "${decoded_files[@]}"; do
        fname=$(basename "${fpath}")
        fifo=${fname#alcdaq.}; fifo=${fifo%.root}

        root_output="${odpath}/dcalib.${fifo}.root"
        txt_output="${odpath}/dcalib.${fifo}.conf"

        run_job "${odpath}/dcalib.${fifo}.log" \
            "${DCALIB}" \
                --input "${fpath}" \
                --output "${root_output}" \
                --calibration-output "${txt_output}" \
                --period "${period}" \
                --min-pairs "${MIN_PAIRS}" &
        pids+=("$!")
    done
done

if [ ${#pids[@]} -eq 0 ]; then
    echo " --- no TDC calibration jobs started "
    exit 1
fi

wait "${pids[@]}"
echo " --- TDC calibration jobs completed "
