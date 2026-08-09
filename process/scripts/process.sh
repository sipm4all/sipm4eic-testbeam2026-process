#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

ipath="/data/2026-testbeam/actual/physics"
opath="/data/2026-testbeam/process"

CALIBRATOR="${ROOT_DIR}/process/bin/calibrator"
SORTER="${ROOT_DIR}/process/bin/sorter"
APS="${ROOT_DIR}/process/bin/after-pulse-suppressor"
MERGER="${ROOT_DIR}/process/bin/merger"
TRIGGER="${ROOT_DIR}/process/bin/trigger"
HADD="hadd"
CALIBRATION_CONFIG="${ROOT_DIR}/process/config/calibration/calibration_example.conf"
TRIGGER_CONFIGS=("${ROOT_DIR}/process/config/trigger/trigger_range.conf" "${ROOT_DIR}/process/config/trigger/trigger_set.conf")
TRIGGER_TAGS=("range" "set")
TRIGGER_WINDOW=256

WRITE_LOGS=0
CLEAN_DEVICE_SPILLS=1
CLEAN_MERGED_SPILLS=1
CLEAN_TRIGGERED_SPILLS=1

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

if [ $# -ne 1 ]; then
    echo " usage: $0 [run] "
    exit 1
fi
run=$1

irpath="${ipath}/${run}"
if [ ! -d "${irpath}" ]; then
   echo " ${irpath} does not exist "
   exit 1
fi
orpath="${opath}/${run}"
mkdir -p "${orpath}"

### loop over device directories
merge_pids=()
for idpath in "${irpath}"/kc705* "${irpath}"/rdo*; do

    [ -d "${idpath}" ] || continue

    device=$(basename "${idpath}")
    odpath="${orpath}/${device}"
    mkdir -p "${odpath}"
    echo " --- processing device ${device}: ${idpath} "

    decoded_files=("${idpath}"/decoded/alcdaq.fifo_*.root)
    if [ ${#decoded_files[@]} -eq 0 ]; then
        echo " --- no decoded files found for ${device} "
        continue
    fi

    ### loop over decoded files
    pids=()
    for fpath in "${decoded_files[@]}"; do

        fname=$(basename "${fpath}")
        fifo=${fname#alcdaq.}; fifo=${fifo%.root}

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

            "${calibrator}" --input "${input}" --output "${calibrated}" --config "${calibration_config}"
            "${sorter}" --input "${calibrated}" --output "${sorted}" --window 32768
            "${aps}" --input "${sorted}" --output "${output}" --window 50
            rm -f "${calibrated}" "${sorted}"
        ' _ "${CALIBRATOR}" "${SORTER}" "${APS}" \
            "${fpath}" \
            "${odpath}/calibrated.${fifo}.root" \
            "${odpath}/sorted.${fifo}.root" \
            "${odpath}/aps.sorted.${fifo}.root" \
            "${CALIBRATION_CONFIG}" &
        pids+=($!)

    done
    wait "${pids[@]}"

    aps_files=("${odpath}"/aps.sorted.fifo_*.root)
    if [ ${#aps_files[@]} -eq 0 ]; then
        echo " --- no after-pulse-suppressed files found for ${device} "
        continue
    fi

    ### merge device, split by spill
    run_job "${odpath}/aps.sorted.log" bash -c '
        merger=$1
        output=$2
        shift 2
        time -p "${merger}" --input "$@" --output "${output}" --split-spills
        rm -f "$@"
    ' _ "${MERGER}" "${odpath}/aps.sorted.root" "${aps_files[@]}" &
    merge_pids+=($!)

done

echo " --- waiting for device merges to complete "
if [ ${#merge_pids[@]} -gt 0 ]; then
    wait "${merge_pids[@]}"
fi

### merge all, one parallel merge per spill
device_dirs=()
for dpath in "${orpath}"/*; do
    [ -d "${dpath}" ] || continue
    files=("${dpath}"/aps.sorted.spill_*.root)
    if [ ${#files[@]} -gt 0 ]; then
        device_dirs+=("${dpath}")
    fi
done

if [ ${#device_dirs[@]} -eq 0 ]; then
    echo " --- no device directories with split-spill files found "
    exit 1
fi

spill_files=("${orpath}"/*/aps.sorted.spill_*.root)
if [ ${#spill_files[@]} -eq 0 ]; then
    echo " --- no device-level split-spill files found "
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
    device_outputs=("${orpath}"/*/aps.sorted.spill_${spill_id}.root)
    if [ ${#device_outputs[@]} -eq 0 ]; then
        echo " --- no device-level files found for spill ${spill_id} "
        exit 1
    fi
    if [ ${#device_outputs[@]} -ne ${#device_dirs[@]} ]; then
        echo " --- spill ${spill_id} has ${#device_outputs[@]} files, expected ${#device_dirs[@]} "
        exit 1
    fi

    run_job "${orpath}/aps.sorted.spill_${spill_id}.log" bash -c '
        merger=$1
        output=$2
        clean_device_spills=$3
        shift 3

        time -p "${merger}" --input "$@" --output "${output}"
        if [ "${clean_device_spills}" -eq 1 ]; then
            rm -f "$@"
        fi
    ' _ "${MERGER}" \
        "${orpath}/aps.sorted.spill_${spill_id}.root" \
        "${CLEAN_DEVICE_SPILLS}" \
        "${device_outputs[@]}" &
    final_pids+=($!)
done

wait "${final_pids[@]}"

echo " --- trigger jobs started "
for i in "${!TRIGGER_CONFIGS[@]}"; do
    tag=${TRIGGER_TAGS[$i]}
    config=${TRIGGER_CONFIGS[$i]}
    echo " --- trigger config ${tag}: ${config} "

    trigger_pids=()
    for spill_id in "${spill_ids[@]}"; do
        merged_spill="${orpath}/aps.sorted.spill_${spill_id}.root"
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

            time -p "${trigger}" --input "${input}" \
                                --output "${output}" \
                                --config "${config}" \
                                --window "${window}"
        ' _ "${TRIGGER}" \
            "${merged_spill}" \
            "${triggered}" \
            "${config}" \
            "${TRIGGER_WINDOW}" &
        trigger_pids+=($!)
    done

    wait "${trigger_pids[@]}"

    triggered_files=("${orpath}"/triggered.${tag}.spill_*.root)
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
        shift 3
        time -p "${hadd}" -f "${output}" "$@"
        if [ "${clean_triggered_spills}" -eq 1 ]; then
            rm -f "$@"
        fi
    ' _ "${HADD}" "${orpath}/triggered.${tag}.root" "${CLEAN_TRIGGERED_SPILLS}" "${triggered_files[@]}"
done

if [ "${CLEAN_MERGED_SPILLS}" -eq 1 ]; then
    rm -f "${orpath}"/aps.sorted.spill_*.root
fi
