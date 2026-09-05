#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

run_type="physics"
ipath="/data/2026-testbeam/process"
opath="/data/2026-testbeam/process"

CALIBRATOR="${ROOT_DIR}/process/bin/calibrator"
COORDINATOR="${ROOT_DIR}/process/bin/coordinator"
SORTER="${ROOT_DIR}/process/bin/sorter"
APS="${ROOT_DIR}/process/bin/after-pulse-suppressor"
MERGER="${ROOT_DIR}/process/bin/merger"
run=""
CALIBRATION_CONFIG=""
CLOCK_CONFIG=""
SORT_WINDOW=32768
APS_WINDOW=64
DEVICE_FILTER=()
DEVICE_ALL=1
OVERWRITE=0

WRITE_LOGS=0
CLEAN_DEVICE_SPILLS=1
CLEAN_MERGED_SPILLS=0

usage()
{
    cat <<EOF
usage:
  $0 --run RUN --calibration CALIBRATION.conf [options]

required:
  --run, -r RUN                  run name/directory
  --calibration, -c FILE         timing calibration configuration
  --clock FILE                   optional run-specific clock correction file

options:
  --run-type TYPE                accepted for compatibility with decoder.sh; process.sh reads decoded files from /data/2026-testbeam/process/RUN
                                  default: physics
  --devices DEVICE ...           device directory names to process, default: all
  --overwrite                    overwrite existing workflow outputs instead of skipping them
  --help, -h                     show this help message

example:
  $0 --run 12345 --calibration process/config/calibration/calibration_example.conf
  $0 --run 20260618-183625 --run-type testpulse --devices rdo-192 --calibration calibration.conf

After this workflow completes, run process/scripts/trigger.sh on the merged spill files.
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
    local _logfile=$1
    shift
    ( "$@" )
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
        --overwrite)
            OVERWRITE=1
            shift
            ;;
        --calibration|-c)
            [ $# -ge 2 ] || fail "$1 requires FILE"
            CALIBRATION_CONFIG=$2
            shift 2
            ;;
        --clock)
            [ $# -ge 2 ] || fail "$1 requires FILE"
            CLOCK_CONFIG=$2
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
if [ ! -f "${CALIBRATION_CONFIG}" ]; then
    fail "calibration config does not exist: ${CALIBRATION_CONFIG}"
fi
if [ -n "${CLOCK_CONFIG}" ] && [ ! -f "${CLOCK_CONFIG}" ]; then
    fail "clock correction file does not exist: ${CLOCK_CONFIG}"
fi

irpath="${ipath}/${run}"
if [ ! -d "${irpath}" ]; then
   echo " ${irpath} does not exist; run decoder.sh first "
   exit 1
fi
orpath="${opath}/${run}/process"
mkdir -p "${orpath}"

merge_prefix="aps.sorted"
if [ "${DEVICE_ALL}" -ne 1 ]; then
    device_tag=$(IFS=_; echo "${DEVICE_FILTER[*]}")
    merge_prefix="aps.sorted.${device_tag}"
fi

merged_spill_files=("${orpath}/${merge_prefix}".spill_*.root)
if [ "${OVERWRITE}" -ne 1 ] && [ ${#merged_spill_files[@]} -gt 0 ]; then
    echo " --- run-level merged spill outputs exist, skipping process workflow: ${orpath}/${merge_prefix}.spill_*.root"
    echo " --- pass --overwrite to regenerate them"
    exit 0
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

        ### calibrate, assign coordinates, sort and AP suppress
        run_job "${odpath}/aps.sorted.${fifo}.log" bash -c '
            set -euo pipefail

            calibrator=$1
            coordinator=$2
            sorter=$3
            aps=$4
            input=$5
            calibrated=$6
            coordinated=$7
            sorted=$8
            output=$9
            calibration_config=${10}
            sort_window=${11}
            aps_window=${12}
            overwrite=${13}
            clock_config=${14}
            run=${15}

            if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
                echo " --- output exists, skipping FIFO processing: ${output}"
                exit 0
            fi

            if [ "${overwrite}" -ne 1 ] && { [ -f "${calibrated}" ] || [ -f "${coordinated}" ] || [ -f "${sorted}" ]; }; then
                echo " --- intermediate output exists, skipping FIFO processing: ${calibrated} ${coordinated} ${sorted}"
                exit 0
            fi

            calibrator_args=(--input "${input}" --output "${calibrated}" --config "${calibration_config}")
            if [ -n "${clock_config}" ]; then
                calibrator_args+=(--clock "${clock_config}" --run "${run}")
            fi
            "${calibrator}" "${calibrator_args[@]}"
            "${coordinator}" --input "${calibrated}" --output "${coordinated}"
            "${sorter}" --input "${coordinated}" --output "${sorted}" --window "${sort_window}"
            "${aps}" --input "${sorted}" --output "${output}" --window "${aps_window}"
            rm -f "${calibrated}" "${coordinated}" "${sorted}"
        ' _ "${CALIBRATOR}" "${COORDINATOR}" "${SORTER}" "${APS}" \
            "${fpath}" \
            "${odpath}/calibrated.${fifo}.root" \
            "${odpath}/coordinated.calibrated.${fifo}.root" \
            "${odpath}/sorted.coordinated.calibrated.${fifo}.root" \
            "${odpath}/aps.sorted.coordinated.calibrated.${fifo}.root" \
            "${CALIBRATION_CONFIG}" \
            "${SORT_WINDOW}" \
            "${APS_WINDOW}" \
            "${OVERWRITE}" \
            "${CLOCK_CONFIG}" \
            "${run}" &
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
        aps_file="${odpath}/aps.sorted.coordinated.calibrated.${fifo}.root"
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
        set -euo pipefail
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

echo " --- final spill merges started "
final_pids=()
for spill_id in "${spill_ids[@]}"; do
    device_outputs=()
    for dpath in "${device_dirs[@]}"; do
        input_file="${dpath}/aps.sorted.spill_${spill_id}.root"
        if [ -f "${input_file}" ]; then
            device_outputs+=("${input_file}")
        else
            echo "WARNING: missing device spill, skipping input: ${input_file}"
        fi
    done
    if [ ${#device_outputs[@]} -eq 0 ]; then
        echo "WARNING: no device-level files found for spill ${spill_id}, skipping"
        continue
    fi
    run_job "${orpath}/${merge_prefix}.spill_${spill_id}.log" bash -c '
        set -euo pipefail

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

if [ "${CLEAN_MERGED_SPILLS}" -eq 1 ]; then
    rm -f "${orpath}"/${merge_prefix}.spill_*.root
fi
cleanup_empty_device_dirs

echo " --- processing completed"
echo " --- merged spill files: ${orpath}/${merge_prefix}.spill_*.root"
echo " --- run process/scripts/trigger.sh to produce triggered frame files"
