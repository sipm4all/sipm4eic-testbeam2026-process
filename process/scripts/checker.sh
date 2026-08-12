#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

run_type="physics"
ipath="/data/2026-testbeam/process"
opath="/data/2026-testbeam/process"

CHECKER="${ROOT_DIR}/process/bin/checker"

run=""
DEVICE_FILTER=()
FIFO_FILTER=()
DEVICE_ALL=1
FIFO_ALL=1

usage()
{
    cat <<EOF_USAGE
usage:
  $0 --run RUN [--run-type TYPE] [--devices all|DEVICE ...] [--fifos all|FIFO ...]

required:
  --run, -r RUN          run name/directory

options:
  --run-type TYPE        accepted for compatibility with decoder.sh; checker.sh reads decoded files from /data/2026-testbeam/process/RUN
                          default: physics
  --devices DEVICE ...   device directory names to process, default: all
  --fifos FIFO ...       FIFO numbers to process, default: all
  --help, -h             show this help message

examples:
  $0 --run 12345
  $0 --run 12345 --run-type testpulse --devices kc705-200 rdo-{192..195} --fifos {0..16}
  $0 --run 12345 --devices all --fifos all
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


value_from_check()
{
    local file=$1
    local key=$2
    awk -v key="${key}" '$1 == key":" { print $2; found=1; exit } END { if (!found) print "" }' "${file}"
}

sum_key()
{
    local key=$1
    shift
    local total=0
    local file value
    for file in "$@"; do
        value=$(value_from_check "${file}" "${key}")
        [[ "${value}" =~ ^[0-9]+$ ]] || value=0
        total=$((total + value))
    done
    echo "${total}"
}

all_yes_key()
{
    local key=$1
    shift
    local file value
    for file in "$@"; do
        value=$(value_from_check "${file}" "${key}")
        if [ "${value}" != "yes" ]; then
            echo "no"
            return
        fi
    done
    echo "yes"
}

all_no_key()
{
    local key=$1
    shift
    local file value
    for file in "$@"; do
        value=$(value_from_check "${file}" "${key}")
        if [ "${value}" != "no" ]; then
            echo "yes"
            return
        fi
    done
    echo "no"
}

write_aggregate_check()
{
    local level=$1
    local name=$2
    local output=$3
    shift 3
    local files=("$@")

    local nfiles=${#files[@]}
    local entries start_spill end_spill alcor_hits trigger_tags unknown_words errors
    entries=$(sum_key entries "${files[@]}")
    start_spill=$(sum_key start_spill_type7 "${files[@]}")
    end_spill=$(sum_key end_spill_type15 "${files[@]}")
    alcor_hits=$(sum_key alcor_hits_type1 "${files[@]}")
    trigger_tags=$(sum_key trigger_tags_type9 "${files[@]}")
    unknown_words=$(sum_key unknown_words "${files[@]}")
    errors=$(sum_key errors "${files[@]}")

    local counters_ok open_ok balance_ok consistent
    counters_ok=$(all_yes_key spill_counter_consistent "${files[@]}")
    open_ok=$(all_no_key open_spill_at_eof "${files[@]}")
    balance_ok=$(all_yes_key spill_count_balance "${files[@]}")

    consistent="yes"
    if [ "${counters_ok}" != "yes" ] || [ "${open_ok}" != "no" ] || [ "${balance_ok}" != "yes" ] || [ "${unknown_words}" -ne 0 ] || [ "${errors}" -ne 0 ]; then
        consistent="no"
    fi

    {
        echo "level: ${level}"
        echo "name: ${name}"
        echo "files: ${nfiles}"
        echo "entries: ${entries}"
        echo "start_spill_type7: ${start_spill}"
        echo "end_spill_type15: ${end_spill}"
        echo "alcor_hits_type1: ${alcor_hits}"
        echo "trigger_tags_type9: ${trigger_tags}"
        echo "unknown_words: ${unknown_words}"
        echo "spill_counter_consistent: ${counters_ok}"
        echo "open_spill_at_eof: ${open_ok}"
        echo "spill_count_balance: ${balance_ok}"
        echo "consistent: ${consistent}"
        echo "errors: ${errors}"
        for file in "${files[@]}"; do
            echo "input_check: ${file}"
        done
        for file in "${files[@]}"; do
            local file_consistent file_errors
            file_consistent=$(value_from_check "${file}" consistent)
            file_errors=$(value_from_check "${file}" errors)
            [[ "${file_errors}" =~ ^[0-9]+$ ]] || file_errors=0
            if [ "${file_consistent}" = "no" ] || [ "${file_errors}" -ne 0 ]; then
                echo "problem_check: ${file}"
            fi
        done
    } > "${output}"
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

if [ ! -x "${CHECKER}" ]; then
    fail "${CHECKER} does not exist or is not executable"
fi

irpath="${ipath}/${run}"
if [ ! -d "${irpath}" ]; then
   fail "${irpath} does not exist; run decoder.sh first"
fi
orpath="${opath}/${run}"

device_checks=()
for device_path in "${irpath}"/kc705* "${irpath}"/rdo*; do
    [ -d "${device_path}" ] || continue

    device=$(basename "${device_path}")
    if [ "${DEVICE_ALL}" -ne 1 ] && ! contains_value "${device}" "${DEVICE_FILTER[@]}"; then
        continue
    fi

    idpath="${device_path}/decoded"
    odpath="${device_path}/decoded"
    echo " --- checking device ${device}: ${idpath} "

    if [ ! -d "${idpath}" ]; then
        echo " --- decoded directory not found for ${device}: ${idpath} "
        continue
    fi

    decoded_files=("${idpath}"/alcdaq.fifo_*.root)
    if [ ${#decoded_files[@]} -eq 0 ]; then
        echo " --- no decoded files found for ${device} "
        continue
    fi

    pids=()
    check_files=()
    for fpath in "${decoded_files[@]}"; do
        fname=$(basename "${fpath}")
        fifo_id=${fname#alcdaq.fifo_}; fifo_id=${fifo_id%.root}
        [[ "${fifo_id}" =~ ^[0-9]+$ ]] || fail "could not extract numeric FIFO id from ${fname}"
        if [ "${FIFO_ALL}" -ne 1 ] && ! contains_value "${fifo_id}" "${FIFO_FILTER[@]}"; then
            continue
        fi

        output="${odpath}/${fname%.root}.check"
        check_files+=("${output}")
        "${CHECKER}" --input "${fpath}" --output "${output}" &
        pids+=("$!")
    done

    if [ ${#pids[@]} -eq 0 ]; then
        echo " --- no checker jobs started for ${device} "
        continue
    fi

    wait "${pids[@]}"
    echo " --- per-FIFO checker jobs completed for ${device} "

    existing_check_files=()
    for check_file in "${check_files[@]}"; do
        [ -f "${check_file}" ] || fail "missing expected check file: ${check_file}"
        existing_check_files+=("${check_file}")
    done

    device_check="${odpath}/${device}.check"
    write_aggregate_check device "${device}" "${device_check}" "${existing_check_files[@]}"
    device_checks+=("${device_check}")
    echo " --- wrote device check ${device_check}"
done

if [ ${#device_checks[@]} -eq 0 ]; then
    echo " --- no device check files produced "
    exit 1
fi

run_check="${orpath}/${run}.check"
write_aggregate_check run "${run}" "${run_check}" "${device_checks[@]}"
echo " --- wrote run check ${run_check}"
echo " --- checker jobs completed "
