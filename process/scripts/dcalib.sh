#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

ipath="/data/2026-testbeam/actual/testpulse"
opath="/data/2026-testbeam/process"

DCALIB="${ROOT_DIR}/process/bin/dcalib"
MIN_PAIRS=1000
WRITE_LOGS=0

run=""
period=""
DEVICE_FILTER=()
FIFO_FILTER=()
DEVICE_ALL=1
FIFO_ALL=1

usage()
{
    cat <<EOF_USAGE
usage:
  $0 --run RUN --period PERIOD [--devices all|DEVICE ...] [--fifos all|FIFO ...]

required:
  --run, -r RUN          run name/directory
  --period, -p PERIOD    expected pulser period in coarse clock cycles

options:
  --devices DEVICE ...   devices to process, default: all
  --fifos FIFO ...       FIFOs to process, default: all
  --help, -h             show this help message

examples:
  $0 --run 12345 --period 10000
  $0 --run 12345 --period 10000 --devices all --fifos all
  $0 --run 12345 --period 10000 --devices 192 196 --fifos 0 4 8
  $0 --run 12345 --period 10000 --devices {192..195} --fifos {0..16}
  $0 --run 12345 --period 10000 --devices "{192..195}" --fifos "{0..16}"
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

add_range_or_value()
{
    local target=$1
    local value=$2

    if [[ "${value}" =~ ^\{([0-9]+)\.\.([0-9]+)\}$ ]]; then
        local first=${BASH_REMATCH[1]}
        local last=${BASH_REMATCH[2]}
        [ "${first}" -le "${last}" ] || fail "invalid range ${value}"
        local v
        for ((v = first; v <= last; ++v)); do
            eval "${target}+=(\"${v}\")"
        done
        return
    fi

    [[ "${value}" =~ ^[0-9]+$ ]] || fail "expected non-negative integer or range, got '${value}'"
    eval "${target}+=(\"${value}\")"
}

parse_filter_values()
{
    local target=$1
    shift
    local values=()

    while [ $# -gt 0 ]; do
        case "$1" in
            --*|-*) break ;;
            *) values+=("$1"); shift ;;
        esac
    done

    [ ${#values[@]} -gt 0 ] || fail "filter option requires at least one value"

    PARSE_FILTER_MODE="list"
    PARSE_FILTER_CONSUMED=$((1 + ${#values[@]}))

    if [ ${#values[@]} -eq 1 ] && [ "${values[0]}" = "all" ]; then
        PARSE_FILTER_MODE="all"
        return
    fi

    for value in "${values[@]}"; do
        [ "${value}" != "all" ] || fail "'all' must be the only value for a filter option"
        add_range_or_value "${target}" "${value}"
    done
}

contains_value()
{
    local value=$1
    shift
    local item
    for item in "$@"; do
        if [ "${item}" = "${value}" ]; then
            return 0
        fi
    done
    return 1
}

extract_last_integer()
{
    local text=$1
    local last=""
    local match
    while read -r match; do
        last=${match}
    done < <(grep -oE '[0-9]+' <<< "${text}" || true)
    [ -n "${last}" ] || return 1
    echo "${last}"
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
        --devices)
            parse_filter_values DEVICE_FILTER "${@:2}"
            if [ "${PARSE_FILTER_MODE}" = "all" ]; then
                DEVICE_ALL=1
                DEVICE_FILTER=()
            else
                DEVICE_ALL=0
            fi
            shift "${PARSE_FILTER_CONSUMED}"
            ;;
        --fifos)
            parse_filter_values FIFO_FILTER "${@:2}"
            if [ "${PARSE_FILTER_MODE}" = "all" ]; then
                FIFO_ALL=1
                FIFO_FILTER=()
            else
                FIFO_ALL=0
            fi
            shift "${PARSE_FILTER_CONSUMED}"
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
    device_id=$(extract_last_integer "${device}") || fail "could not extract numeric device id from ${device}"
    if [ "${DEVICE_ALL}" -ne 1 ] && ! contains_value "${device_id}" "${DEVICE_FILTER[@]}"; then
        continue
    fi

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
        fifo_id=${fname#alcdaq.fifo_}; fifo_id=${fifo_id%.root}
        [[ "${fifo_id}" =~ ^[0-9]+$ ]] || fail "could not extract numeric FIFO id from ${fname}"
        if [ "${FIFO_ALL}" -ne 1 ] && ! contains_value "${fifo_id}" "${FIFO_FILTER[@]}"; then
            continue
        fi

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
