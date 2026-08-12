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
    [ -f "${file}" ] || fail "missing check file while aggregating: ${file}"
    awk -v key="${key}" '$1 == key":" { value=$2; gsub(/\r/, "", value); print value; found=1; exit } END { if (!found) print "" }' "${file}"
}

required_numeric_value()
{
    local file=$1
    local key=$2
    local value
    value=$(value_from_check "${file}" "${key}")
    if [[ "${value}" =~ ^[0-9]+$ ]]; then
        echo "${value}"
        return
    fi
    fail "check file ${file} does not contain numeric key '${key}'"
}

sum_key()
{
    local key=$1
    shift
    local total=0
    local file value
    for file in "$@"; do
        value=$(required_numeric_value "${file}" "${key}")
        total=$((total + value))
    done
    echo "${total}"
}

value_or_fallback()
{
    local file=$1
    local key=$2
    local fallback_key=$3
    local value
    value=$(value_from_check "${file}" "${key}")
    if [[ "${value}" =~ ^[0-9]+$ ]]; then
        echo "${value}"
        return
    fi
    value=$(value_from_check "${file}" "${fallback_key}")
    if [[ "${value}" =~ ^[0-9]+$ ]]; then
        echo "${value}"
        return
    fi
    fail "check file ${file} does not contain numeric key '${key}' or fallback key '${fallback_key}'"
}

min_spill_count_key()
{
    local key=$1
    case "${key}" in
        start_spill_type7)
            echo "min_start_spill_type7"
            ;;
        end_spill_type15)
            echo "min_end_spill_type15"
            ;;
        *)
            echo ""
            ;;
    esac
}

max_spill_count_key()
{
    local key=$1
    case "${key}" in
        start_spill_type7)
            echo "max_start_spill_type7"
            ;;
        end_spill_type15)
            echo "max_end_spill_type15"
            ;;
        *)
            echo ""
            ;;
    esac
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

uniform_spill_count_key()
{
    local key=$1
    shift
    local min_key max_key global_min global_max file file_min file_max
    min_key=$(min_spill_count_key "${key}")
    max_key=$(max_spill_count_key "${key}")
    global_min=""
    global_max=""
    for file in "$@"; do
        file_min=$(value_or_fallback "${file}" "${min_key}" "${key}")
        file_max=$(value_or_fallback "${file}" "${max_key}" "${key}")
        if [ -z "${global_min}" ] || [ "${file_min}" -lt "${global_min}" ]; then
            global_min=${file_min}
        fi
        if [ -z "${global_max}" ] || [ "${file_max}" -gt "${global_max}" ]; then
            global_max=${file_max}
        fi
    done
    if [ "${global_min:-0}" -eq "${global_max:-0}" ]; then
        echo "yes"
    else
        echo "no"
    fi
}

min_spill_count()
{
    local key=$1
    shift
    local min_key min="" file value
    min_key=$(min_spill_count_key "${key}")
    for file in "$@"; do
        value=$(value_or_fallback "${file}" "${min_key}" "${key}")
        if [ -z "${min}" ] || [ "${value}" -lt "${min}" ]; then
            min=${value}
        fi
    done
    echo "${min:-0}"
}

max_spill_count()
{
    local key=$1
    shift
    local max_key max="" file value
    max_key=$(max_spill_count_key "${key}")
    for file in "$@"; do
        value=$(value_or_fallback "${file}" "${max_key}" "${key}")
        if [ -z "${max}" ] || [ "${value}" -gt "${max}" ]; then
            max=${value}
        fi
    done
    echo "${max:-0}"
}

mode_spill_count()
{
    local key=$1
    shift
    local max_key file value candidate best_value best_count
    max_key=$(max_spill_count_key "${key}")
    declare -A counts=()
    for file in "$@"; do
        value=$(value_or_fallback "${file}" "${max_key}" "${key}")
        counts["${value}"]=$(( ${counts["${value}"]:-0} + 1 ))
    done

    best_value=0
    best_count=-1
    for candidate in "${!counts[@]}"; do
        if [ "${counts["${candidate}"]}" -gt "${best_count}" ] || { [ "${counts["${candidate}"]}" -eq "${best_count}" ] && [ "${candidate}" -gt "${best_value}" ]; }; then
            best_value=${candidate}
            best_count=${counts["${candidate}"]}
        fi
    done
    echo "${best_value}"
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
    alcor_hits=$(sum_key alcor_hits_type1 "${files[@]}")
    trigger_tags=$(sum_key trigger_tags_type9 "${files[@]}")
    unknown_words=$(sum_key unknown_words "${files[@]}")
    errors=$(sum_key errors "${files[@]}")

    local counters_ok open_ok balance_ok start_uniform end_uniform consistent
    local min_start max_start min_end max_end mode_start mode_end
    local discard_candidate_count would_pass_without_discard_candidates
    counters_ok=$(all_yes_key spill_counter_consistent "${files[@]}")
    open_ok=$(all_no_key open_spill_at_eof "${files[@]}")
    balance_ok=$(all_yes_key spill_count_balance "${files[@]}")
    start_uniform=$(uniform_spill_count_key start_spill_type7 "${files[@]}")
    end_uniform=$(uniform_spill_count_key end_spill_type15 "${files[@]}")
    min_start=$(min_spill_count start_spill_type7 "${files[@]}")
    max_start=$(max_spill_count start_spill_type7 "${files[@]}")
    min_end=$(min_spill_count end_spill_type15 "${files[@]}")
    max_end=$(max_spill_count end_spill_type15 "${files[@]}")
    mode_start=$(mode_spill_count start_spill_type7 "${files[@]}")
    mode_end=$(mode_spill_count end_spill_type15 "${files[@]}")
    start_spill=${mode_start}
    end_spill=${mode_end}

    discard_candidate_count=0
    would_pass_without_discard_candidates="yes"
    for file in "${files[@]}"; do
        local file_consistent file_errors file_unknown file_counters file_open file_balance
        local file_start_min file_start_max file_end_min file_end_max
        file_consistent=$(value_from_check "${file}" consistent)
        file_errors=$(required_numeric_value "${file}" errors)
        file_unknown=$(required_numeric_value "${file}" unknown_words)
        file_counters=$(value_from_check "${file}" spill_counter_consistent)
        file_open=$(value_from_check "${file}" open_spill_at_eof)
        file_balance=$(value_from_check "${file}" spill_count_balance)
        file_start_min=$(value_or_fallback "${file}" min_start_spill_type7 start_spill_type7)
        file_start_max=$(value_or_fallback "${file}" max_start_spill_type7 start_spill_type7)
        file_end_min=$(value_or_fallback "${file}" min_end_spill_type15 end_spill_type15)
        file_end_max=$(value_or_fallback "${file}" max_end_spill_type15 end_spill_type15)

        if [ "${file_start_min}" -ne "${mode_start}" ] || [ "${file_start_max}" -ne "${mode_start}" ] || [ "${file_end_min}" -ne "${mode_end}" ] || [ "${file_end_max}" -ne "${mode_end}" ]; then
            discard_candidate_count=$((discard_candidate_count + 1))
            continue
        fi

        if [ "${file_consistent}" = "no" ] || [ "${file_errors}" -ne 0 ] || [ "${file_unknown}" -ne 0 ] || [ "${file_counters}" != "yes" ] || [ "${file_open}" != "no" ] || [ "${file_balance}" != "yes" ]; then
            would_pass_without_discard_candidates="no"
        fi
    done
    if [ "${discard_candidate_count}" -ge "${nfiles}" ]; then
        would_pass_without_discard_candidates="no"
    fi

    consistent="yes"
    if [ "${counters_ok}" != "yes" ] || [ "${open_ok}" != "no" ] || [ "${balance_ok}" != "yes" ] || [ "${start_uniform}" != "yes" ] || [ "${end_uniform}" != "yes" ] || [ "${unknown_words}" -ne 0 ] || [ "${errors}" -ne 0 ]; then
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
        echo "spill_count_uniform_start: ${start_uniform}"
        echo "spill_count_uniform_end: ${end_uniform}"
        echo "min_start_spill_type7: ${min_start}"
        echo "max_start_spill_type7: ${max_start}"
        echo "min_end_spill_type15: ${min_end}"
        echo "max_end_spill_type15: ${max_end}"
        echo "consistent: ${consistent}"
        echo "errors: ${errors}"
        echo "discard_candidate_count: ${discard_candidate_count}"
        echo "would_pass_without_discard_candidates: ${would_pass_without_discard_candidates}"
        for file in "${files[@]}"; do
            echo "input_check: ${file}"
        done
        local reference_start reference_end file_start_min file_start_max file_end_min file_end_max
        reference_start="${mode_start}"
        reference_end="${mode_end}"
        for file in "${files[@]}"; do
            local file_consistent file_errors
            file_consistent=$(value_from_check "${file}" consistent)
            file_errors=$(value_from_check "${file}" errors)
            file_start_min=$(value_or_fallback "${file}" min_start_spill_type7 start_spill_type7)
            file_start_max=$(value_or_fallback "${file}" max_start_spill_type7 start_spill_type7)
            file_end_min=$(value_or_fallback "${file}" min_end_spill_type15 end_spill_type15)
            file_end_max=$(value_or_fallback "${file}" max_end_spill_type15 end_spill_type15)
            [[ "${file_errors}" =~ ^[0-9]+$ ]] || file_errors=0
            if [ "${file_consistent}" = "no" ] || [ "${file_errors}" -ne 0 ] || [ "${file_start_min}" -ne "${reference_start}" ] || [ "${file_start_max}" -ne "${reference_start}" ] || [ "${file_end_min}" -ne "${reference_end}" ] || [ "${file_end_max}" -ne "${reference_end}" ]; then
                echo "problem_check: ${file}"
                if [ "${file_start_min}" -ne "${reference_start}" ] || [ "${file_start_max}" -ne "${reference_start}" ] || [ "${file_end_min}" -ne "${reference_end}" ] || [ "${file_end_max}" -ne "${reference_end}" ]; then
                    echo "discard_candidate: ${file}"
                fi
                if [ "${file_consistent}" = "no" ]; then
                    echo "error: ${file}: input check is internally inconsistent"
                fi
                if [ "${file_errors}" -ne 0 ]; then
                    echo "error: ${file}: input check reports ${file_errors} errors"
                fi
                if [ "${file_start_min}" -ne "${reference_start}" ] || [ "${file_start_max}" -ne "${reference_start}" ]; then
                    echo "error: ${file}: START_SPILL count range ${file_start_min}..${file_start_max} differs from common count ${reference_start}"
                fi
                if [ "${file_end_min}" -ne "${reference_end}" ] || [ "${file_end_max}" -ne "${reference_end}" ]; then
                    echo "error: ${file}: END_SPILL count range ${file_end_min}..${file_end_max} differs from common count ${reference_end}"
                fi
            fi
        done
    } > "${output}"
}


fifo_id_from_check()
{
    local check_file=$1
    local name
    name=$(basename "${check_file}")
    name=${name%.check}
    name=${name#alcdaq.fifo_}
    [[ "${name}" =~ ^[0-9]+$ ]] || fail "could not extract FIFO id from check file ${check_file}"
    echo "${name}"
}

decoded_root_from_check()
{
    local check_file=$1
    local input
    input=$(value_from_check "${check_file}" input)
    if [ -n "${input}" ]; then
        echo "${input}"
    else
        echo "${check_file%.check}.root"
    fi
}

write_fifo_selection_lists()
{
    local run_name=$1
    local good_output=$2
    local bad_output=$3
    shift 3
    local device_checks=("$@")
    local tmp_good tmp_bad tmp_initial_good
    tmp_good="${good_output}.tmp"
    tmp_bad="${bad_output}.tmp"
    tmp_initial_good="${good_output}.initial.tmp"

    {
        echo "# run: ${run_name}"
        echo "# columns: device fifo decoded_root check_file"
    } > "${tmp_good}"
    {
        echo "# run: ${run_name}"
        echo "# columns: device fifo decoded_root check_file reason"
    } > "${tmp_bad}"
    : > "${tmp_initial_good}"

    declare -A selected_count_by_device=()
    local device_check device device_consistent device_repairable device_count
    local check_file fifo_id decoded_root
    for device_check in "${device_checks[@]}"; do
        device=$(value_from_check "${device_check}" name)
        [ -n "${device}" ] || device=$(basename "${device_check}" .check)
        device_consistent=$(value_from_check "${device_check}" consistent)
        device_repairable=$(value_from_check "${device_check}" would_pass_without_discard_candidates)
        device_count=$(required_numeric_value "${device_check}" start_spill_type7)

        declare -A discard_checks=()
        while read -r check_file; do
            [ -n "${check_file}" ] || continue
            discard_checks["${check_file}"]=1
        done < <(awk '$1 == "discard_candidate:" { print $2 }' "${device_check}")

        if [ "${device_consistent}" = "yes" ] || [ "${device_repairable}" = "yes" ]; then
            selected_count_by_device["${device}"]=${device_count}
            while read -r check_file; do
                [ -n "${check_file}" ] || continue
                fifo_id=$(fifo_id_from_check "${check_file}")
                decoded_root=$(decoded_root_from_check "${check_file}")
                if [ -n "${discard_checks["${check_file}"]+x}" ]; then
                    echo "${device} ${fifo_id} ${decoded_root} ${check_file} discard_candidate" >> "${tmp_bad}"
                else
                    echo "${device} ${fifo_id} ${device_count} ${decoded_root} ${check_file}" >> "${tmp_initial_good}"
                fi
            done < <(awk '$1 == "input_check:" { print $2 }' "${device_check}")
        else
            while read -r check_file; do
                [ -n "${check_file}" ] || continue
                fifo_id=$(fifo_id_from_check "${check_file}")
                decoded_root=$(decoded_root_from_check "${check_file}")
                echo "${device} ${fifo_id} ${decoded_root} ${check_file} device_not_repairable" >> "${tmp_bad}"
            done < <(awk '$1 == "input_check:" { print $2 }' "${device_check}")
        fi
        unset discard_checks
    done

    declare -A count_frequency=()
    local count run_count best_frequency candidate
    for device in "${!selected_count_by_device[@]}"; do
        count=${selected_count_by_device["${device}"]}
        count_frequency["${count}"]=$(( ${count_frequency["${count}"]:-0} + 1 ))
    done

    run_count=""
    best_frequency=-1
    for candidate in "${!count_frequency[@]}"; do
        if [ "${count_frequency["${candidate}"]}" -gt "${best_frequency}" ] || { [ "${count_frequency["${candidate}"]}" -eq "${best_frequency}" ] && { [ -z "${run_count}" ] || [ "${candidate}" -gt "${run_count}" ]; }; }; then
            run_count=${candidate}
            best_frequency=${count_frequency["${candidate}"]}
        fi
    done

    if [ -n "${run_count}" ]; then
        while read -r device fifo_id count decoded_root check_file; do
            [ -n "${device}" ] || continue
            if [ "${count}" -eq "${run_count}" ]; then
                echo "${device} ${fifo_id} ${decoded_root} ${check_file}" >> "${tmp_good}"
            else
                echo "${device} ${fifo_id} ${decoded_root} ${check_file} run_spill_count_outlier" >> "${tmp_bad}"
            fi
        done < "${tmp_initial_good}"
    else
        while read -r device fifo_id count decoded_root check_file; do
            [ -n "${device}" ] || continue
            echo "${device} ${fifo_id} ${decoded_root} ${check_file} no_repairable_run_selection" >> "${tmp_bad}"
        done < "${tmp_initial_good}"
    fi

    rm -f "${tmp_initial_good}"
    mv "${tmp_good}" "${good_output}"
    mv "${tmp_bad}" "${bad_output}"
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

good_fifo_list="${orpath}/${run}.good-fifos.list"
bad_fifo_list="${orpath}/${run}.bad-fifos.list"
write_fifo_selection_lists "${run}" "${good_fifo_list}" "${bad_fifo_list}" "${device_checks[@]}"
echo " --- wrote good FIFO list ${good_fifo_list}"
echo " --- wrote bad FIFO list ${bad_fifo_list}"
echo " --- checker jobs completed "
