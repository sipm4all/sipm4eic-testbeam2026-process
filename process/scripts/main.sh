#!/usr/bin/env bash
set -euo pipefail

# User configuration.
CONFIG_DIR="/data/2026-testbeam/process/sipm4eic-testbeam2026-process/process/config"
CALIBRATION_CONFIG="${CONFIG_DIR}/calibration/calibration.20260824.conf"
CLOCK_CORRECTIONS_CONFIG="${CONFIG_DIR}/calibration/clock-corrections.empty.conf"
TRIGGER_CONFIG="${CONFIG_DIR}/trigger/timing.conf"
TRIGGER_TAG="timing"
FILTER_CONFIG="${CONFIG_DIR}/filter/recodata.conf"
FILTER_TAG="recodata"
WINDOW=32
RUN_TYPE="physics"
USE_GPU=1
OVERWRITE=0
STAGES=()

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
SCRIPTS="${ROOT_DIR}/process/scripts"
DATA_ROOT="/data/2026-testbeam/process"

usage()
{
    cat <<EOF
usage:
  $0 RUN
  $0 RUNLIST

options:
  --step STAGE ...    run only these stages; may be repeated. Stages are:
                     decoder, checker, process, trigger, timing, ring, filter
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

run_one()
{
    local run=$1
    echo " --- processing run: ${run}"

    local common=(--run "${run}" --run-type "${RUN_TYPE}")
    local overwrite=()
    [ "${OVERWRITE}" -eq 1 ] && overwrite+=(--overwrite)

    local do_decoder=0 do_checker=0 do_process=0 do_trigger=0
    local do_timing=0 do_ring=0 do_filter=0
    if [ ${#STAGES[@]} -eq 0 ]; then
        do_decoder=1; do_checker=1; do_process=1; do_trigger=1
        do_timing=1; do_ring=1; do_filter=1
    else
        local stage
        for stage in "${STAGES[@]}"; do
            case "${stage}" in
                decoder) do_decoder=1 ;; checker) do_checker=1 ;;
                process) do_process=1 ;; trigger) do_trigger=1 ;;
                timing) do_timing=1 ;; ring|ring-finder) do_ring=1 ;;
                filter) do_filter=1 ;;
                *) fail "unknown stage: ${stage}" ;;
            esac
        done
    fi

    [ "${do_decoder}" -eq 1 ] && "${SCRIPTS}/decoder.sh" "${common[@]}" "${overwrite[@]}"
    [ "${do_checker}" -eq 1 ] && "${SCRIPTS}/checker.sh" "${common[@]}"

    [ "${do_process}" -eq 1 ] && "${SCRIPTS}/process.sh" "${common[@]}" \
        --calibration "${CALIBRATION_CONFIG}" \
        --clock "${CLOCK_CORRECTIONS_CONFIG}" \
        "${overwrite[@]}"

    [ "${do_trigger}" -eq 1 ] && "${SCRIPTS}/trigger.sh" "${common[@]}" \
        --trigger "${TRIGGER_CONFIG}" "${TRIGGER_TAG}" \
        --window "${WINDOW}" "${overwrite[@]}"

    [ "${do_timing}" -eq 1 ] && "${SCRIPTS}/timing.sh" "${common[@]}" \
        --trigger "${TRIGGER_TAG}" --parallel-spills --jobs 8 "${overwrite[@]}"

    [ "${do_timing}" -eq 1 ] && rm -f "${DATA_ROOT}/${run}/trigger/triggered.${TRIGGER_TAG}.spill_"*.root

    local gpu=()
    [ "${USE_GPU}" -eq 1 ] && gpu+=(--gpu)
    [ "${do_ring}" -eq 1 ] && "${SCRIPTS}/ring-finder.sh" "${common[@]}" \
        --trigger "${TRIGGER_TAG}" --parallel-spills --jobs 8 "${gpu[@]}" "${overwrite[@]}"

    [ "${do_ring}" -eq 1 ] && rm -f "${DATA_ROOT}/${run}/trigger/timing.${TRIGGER_TAG}.spill_"*.root

    [ "${do_filter}" -eq 1 ] && "${SCRIPTS}/filter.sh" "${common[@]}" \
        --trigger "${TRIGGER_TAG}" \
        --filter "${FILTER_CONFIG}" "${FILTER_TAG}" "${gpu[@]}" "${overwrite[@]}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --step)
            [ $# -ge 2 ] || fail "--step requires STAGE"
            shift
            while [ $# -gt 0 ] && [[ "$1" != --* ]]; do
                case "$1" in
                    decoder|checker|process|trigger|timing|ring|ring-finder|filter)
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
