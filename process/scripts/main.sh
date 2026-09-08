#!/usr/bin/env bash
set -euo pipefail

# Fixed installation and process-data locations used by this workflow.
PROCESS_DIR="/data/2026-testbeam/process"
SCRIPT_DIR="/data/2026-testbeam/process/sipm4eic-testbeam2026-process/process/scripts"
CONFIG_DIR="/data/2026-testbeam/process/sipm4eic-testbeam2026-process/process/config"

# User configuration.
CALIBRATION_CONFIG=""
CALIBRATION_MANIFEST="${CONFIG_DIR}/calibration/runs.conf"
CLOCK_CORRECTION_CONFIG=""
TRIGGER_CONFIG="${CONFIG_DIR}/trigger/timing.conf"
TRIGGER_TAG="timing"
FILTER_CONFIG="${CONFIG_DIR}/filter/recodata.conf"
FILTER_TAG="recodata"
WINDOW=32
DELTAT_MACRO="/data/2026-testbeam/process/sipm4eic-testbeam2026-process/macros/example/deltat.C"
RUN_TYPE="physics"
USE_GPU=1
OVERWRITE=0
STAGES=()


usage()
{
    cat <<EOF
usage:
  $0 RUN
  $0 RUNLIST

options:
  --step STAGE ...    run only these stages; may be repeated. Stages are:
                     decoder, checker, process, trigger, timing, filter, deltat, ring
                     With no --step, run the complete pipeline.
  --overwrite        pass --overwrite to selected stages
  --help, -h         show this help

The argument is interpreted as a run name when it is a directory or as a
newline-separated run list when it is a regular file.

Edit the configuration variables at the top of this script before running.
EOF
}

fail()
{
    echo "ERROR: $*" >&2
    usage >&2
    exit 1
}

resolve_calibration()
{
    local run=$1
    local manifest=${CALIBRATION_MANIFEST}
    [ -f "${manifest}" ] || fail "calibration manifest does not exist: ${manifest}"

    local row
    row=$(awk -v run="${run}" '
        /^[[:space:]]*#/ || NF == 0 { next }
        $1 == run { exact = $0; next }
        $1 == "*" { fallback = $0 }
        END { if (exact != "") print exact; else if (fallback != "") print fallback }
    ' "${manifest}")
    [ -n "${row}" ] || fail "no calibration entry for run ${run} in ${manifest}"

    local key calibration clock
    read -r key calibration clock _ <<< "${row}"
    [ -n "${calibration}" ] || fail "missing calibration filename for run ${run} in ${manifest}"
    CALIBRATION_CONFIG="${CONFIG_DIR}/calibration/${calibration}"
    CLOCK_CORRECTION_CONFIG=""
    if [ -n "${clock:-}" ] && [ "${clock}" != "-" ]; then
        CLOCK_CORRECTION_CONFIG="${CONFIG_DIR}/calibration/${clock}"
    fi
}

run_one()
{
    local run=$1
    echo " --- processing run: ${run}"
    resolve_calibration "${run}"

    local common=(--run "${run}" --run-type "${RUN_TYPE}")
    local overwrite=()
    local clock_option=()
    [ "${OVERWRITE}" -eq 1 ] && overwrite+=(--overwrite)
    [ -n "${CLOCK_CORRECTION_CONFIG}" ] && clock_option+=(--clock "${CLOCK_CORRECTION_CONFIG}")

    local do_decoder=0 do_checker=0 do_process=0 do_trigger=0
    local do_timing=0 do_ring=0 do_filter=0 do_deltat=0
    if [ ${#STAGES[@]} -eq 0 ]; then
        do_decoder=1; do_checker=1; do_process=1; do_trigger=1
        do_timing=1; do_filter=1
    else
        local stage
        for stage in "${STAGES[@]}"; do
            case "${stage}" in
                decoder) do_decoder=1 ;; checker) do_checker=1 ;;
                process) do_process=1 ;; trigger) do_trigger=1 ;;
                timing) do_timing=1 ;; ring|ring-finder) do_ring=1 ;;
                deltat) do_deltat=1 ;;
                filter) do_filter=1 ;;
                *) fail "unknown stage: ${stage}" ;;
            esac
        done
    fi

    [ "${do_decoder}" -eq 1 ] && "${SCRIPT_DIR}/decoder.sh" "${common[@]}" "${overwrite[@]}"
    [ "${do_checker}" -eq 1 ] && "${SCRIPT_DIR}/checker.sh" "${common[@]}"

    [ "${do_process}" -eq 1 ] && "${SCRIPT_DIR}/process.sh" "${common[@]}" \
        --calibration "${CALIBRATION_CONFIG}" \
        "${clock_option[@]}" \
        "${overwrite[@]}"

    [ "${do_trigger}" -eq 1 ] && "${SCRIPT_DIR}/trigger.sh" "${common[@]}" \
        --trigger "${TRIGGER_CONFIG}" "${TRIGGER_TAG}" \
        --window "${WINDOW}" "${overwrite[@]}"

    [ "${do_timing}" -eq 1 ] && "${SCRIPT_DIR}/timing.sh" "${common[@]}" \
        --trigger "${TRIGGER_TAG}" --parallel-spills --jobs 8 "${overwrite[@]}"

    [ "${do_timing}" -eq 1 ] && rm -f "${PROCESS_DIR}/${run}/trigger/triggered.${TRIGGER_TAG}.spill_"*.root

    local gpu=()
    [ "${USE_GPU}" -eq 1 ] && gpu+=(--gpu)

    [ "${do_filter}" -eq 1 ] && "${SCRIPT_DIR}/filter.sh" "${common[@]}" \
        --trigger "${TRIGGER_TAG}" \
        --filter "${FILTER_CONFIG}" "${FILTER_TAG}" "${gpu[@]}" "${overwrite[@]}"

    if [ "${do_deltat}" -eq 1 ]; then
        local analysis_dir="${PROCESS_DIR}/${run}/analysis"
        mkdir -p "${analysis_dir}"
        local deltat_output="${analysis_dir}/deltat.${TRIGGER_TAG}.root"
        if [ "${OVERWRITE}" -eq 1 ] || [ ! -f "${deltat_output}" ]; then
            local deltat_pids=() input spill output
            local timing_inputs=("${PROCESS_DIR}/${run}/trigger"/timing.${TRIGGER_TAG}.spill_*.root)
            [ ${#timing_inputs[@]} -gt 0 ] || fail "no timing spill files for deltat"
            for input in "${timing_inputs[@]}"; do
                spill=${input##*.spill_}; spill=${spill%.root}
                output="${analysis_dir}/deltat.${TRIGGER_TAG}.spill_${spill}.root"
                if [ "${OVERWRITE}" -eq 1 ] || [ ! -f "${output}" ]; then
                    root -l -b -q -e ".L ${DELTAT_MACRO}" \
                        -e "deltat(\"${input}\", std::make_shared<channel_target_t>(-1), std::make_shared<timing_reference_t>(\"T\"), {}, \"${output}\")" &
                    deltat_pids+=("$!")
                fi
            done
            local pid
            for pid in "${deltat_pids[@]}"; do wait "${pid}"; done
            local deltat_inputs=("${analysis_dir}"/deltat.${TRIGGER_TAG}.spill_*.root)
            [ ${#deltat_inputs[@]} -gt 0 ] || fail "deltat produced no spill files"
            hadd -f "${deltat_output}" "${deltat_inputs[@]}"
            rm -f -- "${deltat_inputs[@]}"
        else
            echo " --- deltat output exists, skipping: ${deltat_output}"
        fi
    fi

    [ "${do_ring}" -eq 1 ] && "${SCRIPT_DIR}/ring-finder.sh" "${common[@]}" \
        --trigger "${TRIGGER_TAG}" --parallel-spills --jobs 8 "${gpu[@]}" "${overwrite[@]}"

    [ "${do_ring}" -eq 1 ] && rm -f "${PROCESS_DIR}/${run}/trigger/timing.${TRIGGER_TAG}.spill_"*.root

}

while [ $# -gt 0 ]; do
    case "$1" in
        --step)
            [ $# -ge 2 ] || fail "--step requires STAGE"
            shift
            while [ $# -gt 0 ] && [[ "$1" != --* ]]; do
                case "$1" in
                    decoder|checker|process|trigger|timing|filter|deltat|ring|ring-finder)
                        STAGES+=("$1")
                        shift
                        ;;
                    *)
                        fail "unknown stage: $1"
                        ;;
                esac
            done
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
            [ -z "${run_argument:-}" ] || fail "only one RUN or RUNLIST is allowed"
            run_argument=$1
            shift
            ;;
    esac
done
if [ -z "${run_argument:-}" ]; then
    usage >&2
    exit 1
fi

if [ -f "${run_argument}" ]; then
    while IFS= read -r run || [ -n "${run}" ]; do
        [ -n "${run}" ] || continue
        [[ "${run}" == \#* ]] && continue
        run_one "${run}"
    done < "${run_argument}"
else
    run_one "${run_argument}"
fi
