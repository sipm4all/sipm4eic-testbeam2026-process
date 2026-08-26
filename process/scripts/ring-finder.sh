#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
RING_FINDER="${ROOT_DIR}/process/bin/ring-finder-hough"
HADD="hadd"

run_type="physics"
ipath="/data/2026-testbeam/process"
opath="/data/2026-testbeam/process"
run=""
input_stage="timing"
trigger_tags=()
overwrite=0
parallel_spills=0
jobs=8
clean_ring_spills=0
use_gpu=0
ring_options=()

usage()
{
    cat <<EOF
usage:
  $0 --run RUN --trigger TAG [options]

required:
  --run, -r RUN                  run name/directory
  --trigger, -t TAG              input trigger tag; may be repeated

options:
  --run-type TYPE                default: physics
  --input-stage STAGE            trigger or timing, default: timing
  --gpu                          use the CUDA Hough backend
  --parallel-spills              process split-spill files in parallel and hadd
  --jobs N                       maximum parallel spill jobs, default: 8
  --overwrite                    overwrite existing ring outputs
  --clean-ring-spills            remove ring spill files after hadd, default: keep
  --min-inliers N                minimum ring inliers
  --max-rings N                  maximum rings per frame
  --max-shared-fraction VALUE    maximum shared fraction of the smaller ring
  --max-events N                 maximum frames to process; default: all
  --spatial-resolution VALUE     Hough spatial resolution in mm
  --time-resolution VALUE        Hough time resolution in native units
  --x0-step VALUE                Hough x step in mm
  --y0-step VALUE                Hough y step in mm
  --radius-step VALUE            Hough radius step in mm
  --t-step VALUE                 Hough time step in native units
  --min-e VALUE                  minimum ellipse eccentricity
  --max-e VALUE                  maximum ellipse eccentricity
  --e-step VALUE                 Hough eccentricity step
  --min-phi VALUE                minimum ellipse angle in radians
  --max-phi VALUE                maximum ellipse angle in radians
  --phi-step VALUE               Hough ellipse angle step in radians
  --ransac-tolerance VALUE       RANSAC spatial tolerance in mm
  --ransac-time-window VALUE     RANSAC/local-Hough time half-width
  --help, -h                     show this help message

examples:
  $0 --run 20260623-185238 --trigger timing --gpu
  $0 --run 20260623-185238 --trigger timing --parallel-spills --jobs 16 --gpu
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
    ( "$@" ) > "${logfile}" 2>&1
}

wait_for_jobs()
{
    local max_jobs=$1
    shift
    local -n pids_ref=$1
    while [ "${#pids_ref[@]}" -ge "${max_jobs}" ]; do
        wait "${pids_ref[0]}"
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

run=""
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
            trigger_tags+=("$2")
            shift 2
            ;;
        --input-stage)
            [ $# -ge 2 ] || fail "$1 requires STAGE"
            input_stage=$2
            shift 2
            ;;
        --gpu)
            use_gpu=1
            shift
            ;;
        --parallel-spills)
            parallel_spills=1
            shift
            ;;
        --jobs)
            [ $# -ge 2 ] || fail "$1 requires N"
            jobs=$2
            shift 2
            ;;
        --overwrite)
            overwrite=1
            shift
            ;;
        --clean-ring-spills)
            clean_ring_spills=1
            shift
            ;;
        --min-inliers|--max-rings|--max-shared-fraction|--max-events|--spatial-resolution|\
        --time-resolution|--x0-step|--y0-step|--radius-step|--t-step|\
        --min-e|--max-e|--e-step|--min-phi|--max-phi|--phi-step|\
        --ransac-tolerance|--ransac-time-window)
            [ $# -ge 2 ] || fail "$1 requires VALUE"
            ring_options+=("$1" "$2")
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
[ ${#trigger_tags[@]} -gt 0 ] || fail "at least one --trigger TAG is required"
case "${run_type}" in
    physics|testpulse) ;;
    *) fail "unsupported --run-type: ${run_type}" ;;
esac
case "${input_stage}" in
    trigger|timing) ;;
    *) fail "--input-stage must be trigger or timing" ;;
esac
[[ "${jobs}" =~ ^[0-9]+$ ]] || fail "--jobs must be a positive integer"
[ "${jobs}" -gt 0 ] || fail "--jobs must be greater than zero"
[ -x "${RING_FINDER}" ] || fail "${RING_FINDER} does not exist or is not executable"
command -v "${HADD}" >/dev/null 2>&1 || fail "${HADD} was not found"

declare -A seen_tags=()
for tag in "${trigger_tags[@]}"; do
    [ -n "${tag}" ] || fail "empty trigger tag"
    if [[ -n "${seen_tags[$tag]:-}" ]]; then
        fail "duplicate trigger tag: ${tag}"
    fi
    seen_tags[$tag]=1
done

irpath="${ipath}/${run}/trigger"
orpath="${opath}/${run}/trigger"
[ -d "${irpath}" ] || fail "${irpath} does not exist; run trigger.sh first"
mkdir -p "${orpath}"

echo " --- ring-finder workflow started"
for tag in "${trigger_tags[@]}"; do
    if [ "${input_stage}" = "timing" ]; then
        input_prefix="triggered.${tag}.timing"
        output_prefix="triggered.${tag}.timing.ring"
    else
        input_prefix="triggered.${tag}"
        output_prefix="triggered.${tag}.ring"
    fi
    output="${orpath}/${output_prefix}.root"

    if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
        echo " --- ring output exists, skipping tag ${tag}: ${output}"
        continue
    fi

    if [ "${parallel_spills}" -eq 0 ]; then
        input="${irpath}/${input_prefix}.root"
        [ -f "${input}" ] || fail "input does not exist: ${input}"
        run_job "${orpath}/${output_prefix}.log" bash -c '
            finder=$1
            input=$2
            output=$3
            overwrite=$4
            gpu=$5
            shift 5
            if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
                echo " --- ring output exists, skipping: ${output}"
                exit 0
            fi
            args=()
            [ "${gpu}" -eq 1 ] && args+=(--gpu)
            [ "${overwrite}" -eq 1 ] && args+=(--overwrite)
            time -p "${finder}" --input "${input}" --output "${output}" \
                "${args[@]}" "$@"
        ' _ "${RING_FINDER}" "${input}" "${output}" "${overwrite}" \
            "${use_gpu}" "${ring_options[@]}"
        continue
    fi

    spill_files=("${irpath}"/"${input_prefix}".spill_*.root)
    [ ${#spill_files[@]} -gt 0 ] || \
        fail "no split-spill files found: ${irpath}/${input_prefix}.spill_*.root"
    echo " --- ring tag ${tag}: ${#spill_files[@]} spill files"

    ring_pids=()
    ring_outputs=()
    for input in "${spill_files[@]}"; do
        fname=$(basename "${input}")
        spill_id=${fname#${input_prefix}.spill_}
        spill_id=${spill_id%.root}
        spill_output="${orpath}/${output_prefix}.spill_${spill_id}.root"
        ring_outputs+=("${spill_output}")
        wait_for_jobs "${jobs}" ring_pids
        run_job "${orpath}/${output_prefix}.spill_${spill_id}.log" bash -c '
            finder=$1
            input=$2
            output=$3
            overwrite=$4
            gpu=$5
            shift 5
            if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
                echo " --- ring spill exists, skipping: ${output}"
                exit 0
            fi
            args=()
            [ "${gpu}" -eq 1 ] && args+=(--gpu)
            [ "${overwrite}" -eq 1 ] && args+=(--overwrite)
            time -p "${finder}" --input "${input}" --output "${output}" \
                "${args[@]}" "$@"
        ' _ "${RING_FINDER}" "${input}" "${spill_output}" "${overwrite}" \
            "${use_gpu}" "${ring_options[@]}" &
        ring_pids+=( $! )
    done
    wait_all_jobs ring_pids

    existing_outputs=()
    for spill_output in "${ring_outputs[@]}"; do
        [ -f "${spill_output}" ] && existing_outputs+=("${spill_output}")
    done
    [ ${#existing_outputs[@]} -eq ${#spill_files[@]} ] || \
        fail "ring tag ${tag} has ${#existing_outputs[@]} outputs, expected ${#spill_files[@]}"

    run_job "${orpath}/${output_prefix}.log" bash -c '
        hadd=$1
        output=$2
        clean=$3
        overwrite=$4
        shift 4
        if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then
            echo " --- ring output exists, skipping hadd: ${output}"
            exit 0
        fi
        time -p "${hadd}" -f "${output}" "$@"
        if [ "${clean}" -eq 1 ]; then
            rm -f "$@"
        fi
    ' _ "${HADD}" "${output}" "${clean_ring_spills}" "${overwrite}" \
        "${existing_outputs[@]}"
done

echo " --- ring-finder workflow completed"
