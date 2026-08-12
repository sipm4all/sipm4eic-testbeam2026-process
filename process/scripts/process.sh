#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

run_type="physics"
ipath="/data/2026-testbeam/process"
opath="/data/2026-testbeam/process"

CALIBRATOR="${ROOT_DIR}/process/bin/calibrator"
SORTER="${ROOT_DIR}/process/bin/sorter"
APS="${ROOT_DIR}/process/bin/after-pulse-suppressor"
MERGER="${ROOT_DIR}/process/bin/merger"
TRIGGER="${ROOT_DIR}/process/bin/trigger"
HADD="hadd"

run=""
CALIBRATION_CONFIG=""
TRIGGER_CONFIGS=()
TRIGGER_TAGS=()
SORT_WINDOW=32768
APS_WINDOW=50
TRIGGER_WINDOW=256
DEVICE_FILTER=()
DEVICE_ALL=1
OVERWRITE=0
declare -A GOOD_FIFO=()
USE_GOOD_FIFO_LIST=0
GOOD_FIFO_LIST=""

WRITE_LOGS=0
CLEAN_DEVICE_SPILLS=1
CLEAN_MERGED_SPILLS=1
CLEAN_TRIGGERED_SPILLS=1

usage()
{
    cat <<EOF
usage:
  $0 --run RUN --calibration CALIBRATION.conf --trigger TRIGGER.conf TAG [options]

required:
  --run, -r RUN                  run name/directory
  --calibration, -c FILE         timing calibration configuration
  --trigger, -t FILE TAG         trigger configuration and output tag; may be repeated

options:
  --run-type TYPE                accepted for compatibility with decoder.sh; process.sh reads decoded files from /data/2026-testbeam/process/RUN
                                  default: physics
  --devices DEVICE ...           device directory names to process, default: all
  --good-fifo-list FILE          process only FIFOs listed by checker.sh
  --overwrite                    overwrite existing workflow outputs instead of skipping them
  --window, -w VALUE             trigger frame window, default: 256
  --help, -h                     show this help message

example:
  $0 --run 12345
  $0 --run 12345 --calibration process/config/calibration/calibration_example.conf --trigger process/config/trigger/trigger_range.conf range --trigger process/config/trigger/trigger_set.conf set --window 256
  $0 --run 20260618-183625 --run-type testpulse --devices rdo-192 --calibration calibration.conf --trigger trigger.conf calibcheck_rdo-192 --window 32
EOF
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

parse_device_filter_values()
{
    local values=()

    while [ $# -gt 0 ]; do
        case "$1" in
            --*|-*) break ;;
            *) values+=("$1"); shift ;;
        esac
    done

    [ ${#values[@]} -gt 0 ] || fail "--devices requires at least one value"

    PARSE_FILTER_MODE="list"
    PARSE_FILTER_CONSUMED=$((1 + ${#values[@]}))

    if [ ${#values[@]} -eq 1 ] && [ "${values[0]}" = "all" ]; then
        PARSE_FILTER_MODE="all"
        return
    fi

    for value in "${values[@]}"; do
        [ "${value}" != "all" ] || fail "'all' must be the only value for --devices"
        add_device_range_or_value DEVICE_FILTER "${value}"
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

load_good_fifo_list()
{
    local list=$1
    USE_GOOD_FIFO_LIST=0
    GOOD_FIFO=()

    if [ ! -f "${list}" ]; then
        fail "good FIFO list does not exist: ${list}"
    fi

    local device fifo decoded_root check_file
    while read -r device fifo decoded_root check_file; do
        [ -n "${device}" ] || continue
        case "${device}" in
            \#*) continue ;;
        esac
        [[ "${fifo}" =~ ^[0-9]+$ ]] || continue
        GOOD_FIFO["${device}:${fifo}"]=1
    done < "${list}"

    USE_GOOD_FIFO_LIST=1
    echo " --- using checker good-FIFO list ${list}"
}

good_fifo_allowed()
{
    local device=$1
    local fifo=$2
    if [ "${USE_GOOD_FIFO_LIST}" -ne 1 ]; then
        return 0
    fi
    [ -n "${GOOD_FIFO["${device}:${fifo}"]+x}" ]
}

cleanup_empty_device_dirs()
{
    if [ "${WRITE_LOGS}" -ne 0 ]; then
        return
    fi
    if [ "${CLEAN_DEVICE_SPILLS}" -ne 1 ]; then
        return
    fi

    for dpath in "${device_dirs[@]}"; do
        [ -d "${dpath}" ] || continue
        rmdir "${dpath}" 2>/dev/null || true
    done
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
            parse_device_filter_values "${@:2}"
            if [ "${PARSE_FILTER_MODE}" = "all" ]; then
                DEVICE_ALL=1
                DEVICE_FILTER=()
            else
                DEVICE_ALL=0
            fi
            shift "${PARSE_FILTER_CONSUMED}"
            ;;
        --good-fifo-list)
            [ $# -ge 2 ] || fail "$1 requires FILE"
            GOOD_FIFO_LIST=$2
            shift 2
            ;;
        --overwrite)
            OVERWRITE=1
            shift
            ;;
        --calibration|-c)
            [ $# -ge 2 ] || fail "$1 requires FILE"
            CALIBRATION_CONFIG=$2
            shift 2
            ;;
        --trigger|-t)
            [ $# -ge 3 ] || fail "$1 requires FILE TAG"
            TRIGGER_CONFIGS+=("$2")
            TRIGGER_TAGS+=("$3")
            shift 3
            ;;
        --window|-w)
            [ $# -ge 2 ] || fail "$1 requires VALUE"
            TRIGGER_WINDOW=$2
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
case "${run_type}" in
    physics|testpulse)
        ;;
    *)
        fail "unsupported --run-type: ${run_type}"
        ;;
esac
[ -n "${CALIBRATION_CONFIG}" ] || fail "missing --calibration"
[ ${#TRIGGER_CONFIGS[@]} -gt 0 ] || fail "at least one --trigger FILE TAG is required"
[[ "${TRIGGER_WINDOW}" =~ ^[0-9]+([.][0-9]+)?$ ]] || fail "--window must be a positive number"
awk "BEGIN { exit !(${TRIGGER_WINDOW} > 0) }" || fail "--window must be greater than zero"

if [ ! -f "${CALIBRATION_CONFIG}" ]; then
    fail "calibration config does not exist: ${CALIBRATION_CONFIG}"
fi

declare -A seen_trigger_tags=()
for i in "${!TRIGGER_CONFIGS[@]}"; do
    config=${TRIGGER_CONFIGS[$i]}
    tag=${TRIGGER_TAGS[$i]}
    [ -f "${config}" ] || fail "trigger config does not exist: ${config}"
    [ -n "${tag}" ] || fail "empty trigger tag for config ${config}"
    if [[ -n "${seen_trigger_tags[$tag]:-}" ]]; then
        fail "duplicate trigger tag: ${tag}"
    fi
    seen_trigger_tags[$tag]=1
done

irpath="${ipath}/${run}"
if [ ! -d "${irpath}" ]; then
   echo " ${irpath} does not exist; run decoder.sh first "
   exit 1
fi
orpath="${opath}/${run}"
if [ -n "${GOOD_FIFO_LIST}" ]; then
    load_good_fifo_list "${GOOD_FIFO_LIST}"
fi

merge_prefix="aps.sorted"
if [ "${DEVICE_ALL}" -ne 1 ]; then
    device_tag=$(IFS=_; echo "${DEVICE_FILTER[*]}")
    merge_prefix="aps.sorted.${device_tag}"
fi

if [ "${OVERWRITE}" -ne 1 ]; then
    final_outputs_done=1
    for tag in "${TRIGGER_TAGS[@]}"; do
        if [ ! -f "${orpath}/triggered.${tag}.root" ]; then
            final_outputs_done=0
            break
        fi
    done
    if [ "${final_outputs_done}" -eq 1 ]; then
        echo " --- final triggered outputs exist, nothing to do"
        exit 0
    fi
fi

### loop over device directories
merge_pids=()
processed_device_dirs=()
for device_path in "${irpath}"/kc705* "${irpath}"/rdo*; do

    [ -d "${device_path}" ] || continue

    device=$(basename "${device_path}")
    if [ "${DEVICE_ALL}" -ne 1 ] && ! contains_value "${device}" "${DEVICE_FILTER[@]}"; then
        continue
    fi
    idpath="${device_path}/decoded"
    echo " --- processing device ${device}: ${idpath} "

    if [ ! -d "${idpath}" ]; then
        echo " --- decoded directory not found for ${device}: ${idpath} "
        continue
    fi

    odpath="${device_path}/process"
    mkdir -p "${odpath}"

    device_spill_files=("${odpath}"/aps.sorted.spill_*.root)
    if [ "${OVERWRITE}" -ne 1 ] && [ ${#device_spill_files[@]} -gt 0 ]; then
        echo " --- device split-spill outputs exist, skipping device processing: ${odpath}/aps.sorted.spill_*.root"
        processed_device_dirs+=("${odpath}")
        continue
    fi

    decoded_files=("${idpath}"/alcdaq.fifo_*.root)
    if [ ${#decoded_files[@]} -eq 0 ]; then
        echo " --- no decoded files found for ${device} "
        continue
    fi

    ### loop over decoded files
    pids=()
    for fpath in "${decoded_files[@]}"; do

        fname=$(basename "${fpath}")
        fifo=${fname#alcdaq.}; fifo=${fifo%.root}
        fifo_id=${fname#alcdaq.fifo_}; fifo_id=${fifo_id%.root}
        [[ "${fifo_id}" =~ ^[0-9]+$ ]] || fail "could not extract numeric FIFO id from ${fname}"
        if ! good_fifo_allowed "${device}" "${fifo_id}"; then
            echo " --- skipping FIFO excluded by checker selection: ${device} fifo ${fifo_id}"
            continue
        fi

        ### calibrate, sort and AP suppress
        run_job "${odpath}/aps.sorted.${fifo}.log" bash -c '
            calibrator=$1
            sorter=$2
            aps=$3
            input=$4
            calibrated=$5
            sorted=$6
            output=$7
            calibration_config=$8
            sort_window=$9
            aps_window=${10}
            overwrite=${11}

            if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
                echo " --- output exists, skipping FIFO processing: ${output}"
                exit 0
            fi

            if [ "${overwrite}" -ne 1 ] && { [ -f "${calibrated}" ] || [ -f "${sorted}" ]; }; then
                echo " --- intermediate output exists, skipping FIFO processing: ${calibrated} ${sorted}"
                exit 0
            fi

            "${calibrator}" --input "${input}" --output "${calibrated}" --config "${calibration_config}"
            "${sorter}" --input "${calibrated}" --output "${sorted}" --window "${sort_window}"
            "${aps}" --input "${sorted}" --output "${output}" --window "${aps_window}"
            rm -f "${calibrated}" "${sorted}"
        ' _ "${CALIBRATOR}" "${SORTER}" "${APS}" \
            "${fpath}" \
            "${odpath}/calibrated.${fifo}.root" \
            "${odpath}/sorted.calibrated.${fifo}.root" \
            "${odpath}/aps.sorted.calibrated.${fifo}.root" \
            "${CALIBRATION_CONFIG}" \
            "${SORT_WINDOW}" \
            "${APS_WINDOW}" \
            "${OVERWRITE}" &
        pids+=($!)

    done

    if [ ${#pids[@]} -eq 0 ]; then
        echo " --- no FIFO processing jobs started for ${device} "
        continue
    fi
    wait "${pids[@]}"

    aps_files=()
    for fpath in "${decoded_files[@]}"; do
        fname=$(basename "${fpath}")
        fifo=${fname#alcdaq.}; fifo=${fifo%.root}
        fifo_id=${fname#alcdaq.fifo_}; fifo_id=${fifo_id%.root}
        [[ "${fifo_id}" =~ ^[0-9]+$ ]] || fail "could not extract numeric FIFO id from ${fname}"
        if ! good_fifo_allowed "${device}" "${fifo_id}"; then
            continue
        fi
        aps_file="${odpath}/aps.sorted.calibrated.${fifo}.root"
        if [ -f "${aps_file}" ]; then
            aps_files+=("${aps_file}")
        fi
    done
    if [ ${#aps_files[@]} -eq 0 ]; then
        echo " --- no selected after-pulse-suppressed files found for ${device} "
        continue
    fi

    processed_device_dirs+=("${odpath}")

    ### merge device, split by spill
    run_job "${odpath}/aps.sorted.log" bash -c '
        shopt -s nullglob
        merger=$1
        output=$2
        overwrite=$3
        shift 3

        existing=("${output%.root}".spill_*.root)
        if [ "${overwrite}" -ne 1 ] && [ ${#existing[@]} -gt 0 ]; then
            echo " --- device split-spill outputs exist, skipping device merge: ${output%.root}.spill_*.root"
            exit 0
        fi

        time -p "${merger}" --input "$@" --output "${output}" --split-spills
        rm -f "$@"
    ' _ "${MERGER}" "${odpath}/aps.sorted.root" "${OVERWRITE}" "${aps_files[@]}" &
    merge_pids+=($!)

done

echo " --- waiting for device merges to complete "
if [ ${#merge_pids[@]} -gt 0 ]; then
    wait "${merge_pids[@]}"
fi

### merge all, one parallel merge per spill
device_dirs=("${processed_device_dirs[@]}")

if [ ${#device_dirs[@]} -eq 0 ]; then
    echo " --- no device directories processed in this invocation "
    exit 1
fi

spill_files=()
for dpath in "${device_dirs[@]}"; do
    files=("${dpath}"/aps.sorted.spill_*.root)
    spill_files+=("${files[@]}")
done

if [ ${#spill_files[@]} -eq 0 ]; then
    echo " --- no device-level split-spill files found for devices processed in this invocation "
    exit 1
fi

spill_ids=()
while IFS= read -r spill_id; do
    spill_ids+=("${spill_id}")
done < <(
    for fpath in "${spill_files[@]}"; do
        fname=$(basename "${fpath}")
        spill_id=${fname#aps.sorted.spill_}
        spill_id=${spill_id%.root}
        echo "${spill_id}"
    done | sort -u
)

if [ ! -x "${TRIGGER}" ]; then
    echo " ${TRIGGER} does not exist or is not executable "
    exit 1
fi
if ! command -v "${HADD}" >/dev/null 2>&1; then
    echo " ${HADD} was not found "
    exit 1
fi
if [ ${#TRIGGER_CONFIGS[@]} -eq 0 ]; then
    echo " --- no trigger configurations defined "
    exit 1
fi
if [ ${#TRIGGER_CONFIGS[@]} -ne ${#TRIGGER_TAGS[@]} ]; then
    echo " --- TRIGGER_CONFIGS and TRIGGER_TAGS must have the same length "
    exit 1
fi
for config in "${TRIGGER_CONFIGS[@]}"; do
    if [ ! -f "${config}" ]; then
        echo " ${config} does not exist "
        exit 1
    fi
done

echo " --- final spill merges started "
final_pids=()
for spill_id in "${spill_ids[@]}"; do
    device_outputs=()
    for dpath in "${device_dirs[@]}"; do
        files=("${dpath}"/aps.sorted.spill_${spill_id}.root)
        device_outputs+=("${files[@]}")
    done
    if [ ${#device_outputs[@]} -eq 0 ]; then
        echo " --- no device-level files found for spill ${spill_id} "
        exit 1
    fi
    if [ ${#device_outputs[@]} -ne ${#device_dirs[@]} ]; then
        echo " --- spill ${spill_id} has ${#device_outputs[@]} files, expected ${#device_dirs[@]} "
        exit 1
    fi

    run_job "${orpath}/${merge_prefix}.spill_${spill_id}.log" bash -c '
        merger=$1
        output=$2
        clean_device_spills=$3
        overwrite=$4
        shift 4

        if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
            echo " --- merged spill exists, skipping final spill merge: ${output}"
            exit 0
        fi

        time -p "${merger}" --input "$@" --output "${output}"
        if [ "${clean_device_spills}" -eq 1 ]; then
            rm -f "$@"
        fi
    ' _ "${MERGER}" \
        "${orpath}/${merge_prefix}.spill_${spill_id}.root" \
        "${CLEAN_DEVICE_SPILLS}" \
        "${OVERWRITE}" \
        "${device_outputs[@]}" &
    final_pids+=($!)
done

wait "${final_pids[@]}"
cleanup_empty_device_dirs

echo " --- trigger jobs started "
for i in "${!TRIGGER_CONFIGS[@]}"; do
    tag=${TRIGGER_TAGS[$i]}
    config=${TRIGGER_CONFIGS[$i]}
    echo " --- trigger config ${tag}: ${config} "

    trigger_pids=()
    for spill_id in "${spill_ids[@]}"; do
        merged_spill="${orpath}/${merge_prefix}.spill_${spill_id}.root"
        if [ ! -f "${merged_spill}" ]; then
            echo " --- merged spill file not found: ${merged_spill} "
            exit 1
        fi

        triggered="${orpath}/triggered.${tag}.spill_${spill_id}.root"

        run_job "${orpath}/triggered.${tag}.spill_${spill_id}.log" bash -c '
            trigger=$1
            input=$2
            output=$3
            config=$4
            window=$5
            overwrite=$6

            if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
                echo " --- triggered spill exists, skipping trigger: ${output}"
                exit 0
            fi

            time -p "${trigger}" --input "${input}" \
                                --output "${output}" \
                                --config "${config}" \
                                --window "${window}"
        ' _ "${TRIGGER}" \
            "${merged_spill}" \
            "${triggered}" \
            "${config}" \
            "${TRIGGER_WINDOW}" \
            "${OVERWRITE}" &
        trigger_pids+=($!)
    done

    wait "${trigger_pids[@]}"

    triggered_files=()
    for spill_id in "${spill_ids[@]}"; do
        triggered="${orpath}/triggered.${tag}.spill_${spill_id}.root"
        if [ -f "${triggered}" ]; then
            triggered_files+=("${triggered}")
        fi
    done
    if [ ${#triggered_files[@]} -eq 0 ]; then
        echo " --- no triggered spill files found for tag ${tag} "
        exit 1
    fi
    if [ ${#triggered_files[@]} -ne ${#spill_ids[@]} ]; then
        echo " --- trigger tag ${tag} has ${#triggered_files[@]} files, expected ${#spill_ids[@]} "
        exit 1
    fi

    run_job "${orpath}/triggered.${tag}.log" bash -c '
        hadd=$1
        output=$2
        clean_triggered_spills=$3
        overwrite=$4
        shift 4

        if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
            echo " --- triggered output exists, skipping hadd: ${output}"
            exit 0
        fi

        time -p "${hadd}" -f "${output}" "$@"
        if [ "${clean_triggered_spills}" -eq 1 ]; then
            rm -f "$@"
        fi
    ' _ "${HADD}" "${orpath}/triggered.${tag}.root" "${CLEAN_TRIGGERED_SPILLS}" "${OVERWRITE}" "${triggered_files[@]}"
done

if [ "${CLEAN_MERGED_SPILLS}" -eq 1 ]; then
    rm -f "${orpath}"/${merge_prefix}.spill_*.root
fi
cleanup_empty_device_dirs
