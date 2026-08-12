#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

actual_base="/data/2026-testbeam/actual"
run_type="physics"
opath="/data/2026-testbeam/process"

DECODER="${ROOT_DIR}/process/bin/decoder"

run=""
ALLOWED_SPILL_ERRORS=0
OVERWRITE=0
DEVICE_FILTER=()
FIFO_FILTER=()
DEVICE_ALL=1
FIFO_ALL=1

usage()
{
    cat <<EOF_USAGE
usage:
  $0 --run RUN [--run-type TYPE] [--devices all|DEVICE ...] [--fifos all|FIFO ...] [options]

required:
  --run, -r RUN                  run name/directory

options:
  --run-type TYPE                input run type under /data/2026-testbeam/actual, default: physics
                                  supported: physics, testpulse
  --devices DEVICE ...           device directory names to process, default: all
  --fifos FIFO ...               FIFO numbers to process, default: all
  --allowed-spill-errors N       maximum errors before a spill payload is emptied, default: 0
  --overwrite                    overwrite existing decoded ROOT files
  --help, -h                     show this help message

examples:
  $0 --run 12345
  $0 --run 12345 --run-type testpulse --devices kc705-200 rdo-{192..195} --fifos {0..16}
  $0 --run 12345 --allowed-spill-errors 1
EOF_USAGE
}

fail()
{
    echo "ERROR: $*" >&2
    echo >&2
    usage >&2
    exit 1
}

add_fifo_range_or_value()
{
    local target=$1
    local value=$2

    if [[ "${value}" =~ ^\{([0-9]+)\.\.([0-9]+)\}$ ]]; then
        local first=${BASH_REMATCH[1]}
        local last=${BASH_REMATCH[2]}
        [ "${first}" -le "${last}" ] || fail "invalid FIFO range ${value}"
        local v
        for ((v = first; v <= last; ++v)); do
            eval "${target}+=(\"${v}\")"
        done
        return
    fi

    [[ "${value}" =~ ^[0-9]+$ ]] || fail "expected FIFO number or range, got '${value}'"
    eval "${target}+=(\"${value}\")"
}

add_device_range_or_value()
{
    local target=$1
    local value=$2

    if [[ "${value}" =~ ^([A-Za-z0-9_.-]*)\{([0-9]+)\.\.([0-9]+)\}$ ]]; then
        local prefix=${BASH_REMATCH[1]}
        local first=${BASH_REMATCH[2]}
        local last=${BASH_REMATCH[3]}
        [ "${first}" -le "${last}" ] || fail "invalid device range ${value}"
        local width=${#first}
        local v
        for ((v = first; v <= last; ++v)); do
            printf -v formatted "%0${width}d" "${v}"
            eval "${target}+=(\"${prefix}${formatted}\")"
        done
        return
    fi

    [[ "${value}" =~ ^[A-Za-z0-9_.-]+$ ]] || fail "expected device name or prefixed range, got '${value}'"
    eval "${target}+=(\"${value}\")"
}

parse_filter_values()
{
    local target=$1
    local kind=$2
    shift 2
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
        if [ "${kind}" = "device" ]; then
            add_device_range_or_value "${target}" "${value}"
        else
            add_fifo_range_or_value "${target}" "${value}"
        fi
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

while [ $# -gt 0 ]; do
    case "$1" in
        --run|-r)
            [ $# -ge 2 ] || fail "$1 requires RUN"
            run=$2
            shift 2
            ;;
        --run-type)
            [ $# -ge 2 ] || fail "$1 requires TYPE"
            run_type=$2
            shift 2
            ;;
        --devices)
            parse_filter_values DEVICE_FILTER device "${@:2}"
            if [ "${PARSE_FILTER_MODE}" = "all" ]; then
                DEVICE_ALL=1
                DEVICE_FILTER=()
            else
                DEVICE_ALL=0
            fi
            shift "${PARSE_FILTER_CONSUMED}"
            ;;
        --fifos)
            parse_filter_values FIFO_FILTER fifo "${@:2}"
            if [ "${PARSE_FILTER_MODE}" = "all" ]; then
                FIFO_ALL=1
                FIFO_FILTER=()
            else
                FIFO_ALL=0
            fi
            shift "${PARSE_FILTER_CONSUMED}"
            ;;
        --allowed-spill-errors)
            [ $# -ge 2 ] || fail "$1 requires N"
            ALLOWED_SPILL_ERRORS=$2
            shift 2
            ;;
        --overwrite)
            OVERWRITE=1
            shift
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
case "${run_type}" in
    physics|testpulse)
        ;;
    *)
        fail "unsupported --run-type: ${run_type}"
        ;;
esac
[[ "${ALLOWED_SPILL_ERRORS}" =~ ^[0-9]+$ ]] || fail "--allowed-spill-errors must be a non-negative integer"

if [ ! -x "${DECODER}" ]; then
    fail "${DECODER} does not exist or is not executable"
fi

irpath="${actual_base}/${run_type}/${run}"
if [ ! -d "${irpath}" ]; then
   fail "${irpath} does not exist"
fi

orpath="${opath}/${run}"
mkdir -p "${orpath}"

device_count=0
for idpath in "${irpath}"/kc705* "${irpath}"/rdo*; do
    [ -d "${idpath}" ] || continue

    device=$(basename "${idpath}")
    if [ "${DEVICE_ALL}" -ne 1 ] && ! contains_value "${device}" "${DEVICE_FILTER[@]}"; then
        continue
    fi

    raw_files=("${idpath}"/raw/alcdaq.fifo_*.dat)
    if [ ${#raw_files[@]} -eq 0 ]; then
        echo " --- no raw files found for ${device} "
        continue
    fi

    odpath="${orpath}/${device}/decoded"
    mkdir -p "${odpath}"
    echo " --- decoding device ${device}: ${idpath} "
    device_count=$((device_count + 1))

    pids=()
    for fpath in "${raw_files[@]}"; do
        fname=$(basename "${fpath}")
        fifo_id=${fname#alcdaq.fifo_}; fifo_id=${fifo_id%.dat}
        [[ "${fifo_id}" =~ ^[0-9]+$ ]] || fail "could not extract numeric FIFO id from ${fname}"
        if [ "${FIFO_ALL}" -ne 1 ] && ! contains_value "${fifo_id}" "${FIFO_FILTER[@]}"; then
            continue
        fi

        output="${odpath}/${fname%.dat}.root"
        if [ "${OVERWRITE}" -ne 1 ] && [ -f "${output}" ]; then
            echo " --- decoded output exists, skipping: ${output}"
            continue
        fi

        "${DECODER}" \
            --input "${fpath}" \
            --output "${output}" \
            --allowed-spill-errors "${ALLOWED_SPILL_ERRORS}" &
        pids+=("$!")
    done

    if [ ${#pids[@]} -eq 0 ]; then
        echo " --- no decoder jobs started for ${device} "
        continue
    fi

    wait "${pids[@]}"
    echo " --- decoder jobs completed for ${device} "
done

if [ "${device_count}" -eq 0 ]; then
    echo " --- no devices decoded "
    exit 1
fi

echo " --- decoder workflow completed "
