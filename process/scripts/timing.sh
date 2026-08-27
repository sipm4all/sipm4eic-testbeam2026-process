#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

run_type="physics"
ipath="/data/2026-testbeam/process"
opath="/data/2026-testbeam/process"

TIMING="${ROOT_DIR}/process/bin/timing"

run=""
TRIGGER_TAGS=()
OVERWRITE=0
PARALLEL_SPILLS=0
JOBS=8

WRITE_LOGS=0
CLEAN_TIMING_SPILLS=0

usage()
{
    cat <<EOF
usage:
  $0 --run RUN --trigger TAG [options]

required:
  --run, -r RUN                  run name/directory
  --trigger, -t TAG              trigger output tag; may be repeated

options:
  --run-type TYPE                accepted for symmetry with the other workflows, default: physics
  --parallel-spills              run timing on triggered.<tag>.spill_*.root in parallel
  --jobs N                       maximum parallel spill timing jobs, default: 8
  --overwrite                    overwrite existing timing outputs instead of skipping them
  --clean-timing-spills          legacy option; split spill outputs are retained by this workflow
  --help, -h                     show this help message

examples:
  $0 --run 20260623-185238 --trigger fingers
  $0 --run 20260623-185238 --trigger fingers --parallel-spills --jobs 16
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

wait_for_jobs()
{
    local max_jobs=$1
    shift
    local -n pids_ref=$1

    while [ "${#pids_ref[@]}" -ge "${max_jobs}" ]; do
        local pid=${pids_ref[0]}
        wait "${pid}"
        pids_ref=("${pids_ref[@]:1}")
    done
}

wait_all_jobs()
{
    local -n pids_ref=$1
    local pid

    for pid in "${pids_ref[@]}"; do
        wait "${pid}"
    done
    pids_ref=()
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
        --trigger|-t)
            [ $# -ge 2 ] || fail "$1 requires TAG"
            TRIGGER_TAGS+=("$2")
            shift 2
            ;;
        --parallel-spills)
            PARALLEL_SPILLS=1
            shift
            ;;
        --jobs)
            [ $# -ge 2 ] || fail "$1 requires N"
            JOBS=$2
            shift 2
            ;;
        --overwrite)
            OVERWRITE=1
            shift
            ;;
        --clean-timing-spills)
            CLEAN_TIMING_SPILLS=1
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
[ ${#TRIGGER_TAGS[@]} -gt 0 ] || fail "at least one --trigger TAG is required"
[[ "${JOBS}" =~ ^[0-9]+$ ]] || fail "--jobs must be a positive integer"
[ "${JOBS}" -gt 0 ] || fail "--jobs must be greater than zero"

if [ ! -x "${TIMING}" ]; then
    fail "${TIMING} does not exist or is not executable"
fi

irpath="${ipath}/${run}/trigger"
if [ ! -d "${irpath}" ]; then
   fail "${irpath} does not exist; run trigger.sh first"
fi
orpath="${opath}/${run}/trigger"
mkdir -p "${orpath}"

declare -A seen_trigger_tags=()
for tag in "${TRIGGER_TAGS[@]}"; do
    [ -n "${tag}" ] || fail "empty trigger tag"
    if [[ -n "${seen_trigger_tags[$tag]:-}" ]]; then
        fail "duplicate trigger tag: ${tag}"
    fi
    seen_trigger_tags[$tag]=1
done

echo " --- timing workflow started"

for tag in "${TRIGGER_TAGS[@]}"; do
        output="${orpath}/timing.${tag}.root"

    if [ "${OVERWRITE}" -ne 1 ] && [ -f "${output}" ]; then
        echo " --- timing output exists, skipping tag ${tag}: ${output}"
        continue
    fi

    if [ "${PARALLEL_SPILLS}" -eq 0 ]; then
        input="${irpath}/triggered.${tag}.root"
        [ -f "${input}" ] || fail "triggered input does not exist: ${input}"

        echo " --- timing tag ${tag}: ${input}"
        run_job "${orpath}/timing.${tag}.log" \
            time -p "${TIMING}" --input "${input}" --output "${output}"
        continue
    fi

    spill_files=("${irpath}"/triggered.${tag}.spill_*.root)
    if [ ${#spill_files[@]} -eq 0 ]; then
        fail "no triggered split-spill files found: ${irpath}/triggered.${tag}.spill_*.root"
    fi

    echo " --- timing tag ${tag}: ${#spill_files[@]} spill files"
    timing_pids=()
    timing_spill_files=()

    for input in "${spill_files[@]}"; do
        fname=$(basename "${input}")
        spill_id=${fname#triggered.${tag}.spill_}
        spill_id=${spill_id%.root}
        spill_output="${orpath}/timing.${tag}.spill_${spill_id}.root"
        timing_spill_files+=("${spill_output}")

        wait_for_jobs "${JOBS}" timing_pids
        run_job "${orpath}/timing.${tag}.spill_${spill_id}.log" bash -c '
            timing=$1
            input=$2
            output=$3
            overwrite=$4

            if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
                echo " --- timing spill exists, skipping: ${output}"
                exit 0
            fi

            time -p "${timing}" --input "${input}" --output "${output}"
        ' _ "${TIMING}" "${input}" "${spill_output}" "${OVERWRITE}" &
        timing_pids+=($!)
    done

    wait_all_jobs timing_pids

    existing_timing_spills=()
    for spill_output in "${timing_spill_files[@]}"; do
        if [ -f "${spill_output}" ]; then
            existing_timing_spills+=("${spill_output}")
        fi
    done
    if [ ${#existing_timing_spills[@]} -ne ${#spill_files[@]} ]; then
        fail "timing tag ${tag} has ${#existing_timing_spills[@]} spill outputs, expected ${#spill_files[@]}"
    fi

    echo " --- timing tag ${tag}: kept ${#existing_timing_spills[@]} spill files"
done

echo " --- timing workflow completed"
