#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
FILTER="${ROOT_DIR}/process/bin/filter"
HADD="hadd"
ipath="/data/2026-testbeam/process"; opath="/data/2026-testbeam/process"
run=""; trigger_tag=""; filter_tag=""; filter_config=""; overwrite=0; jobs=8
usage() { echo "usage: $0 --run RUN --trigger TAG --filter FILE TAG [--jobs N] [--overwrite]"; }
fail() { echo "ERROR: $*" >&2; usage >&2; exit 1; }
while [ $# -gt 0 ]; do
    case "$1" in
        --run|-r) [ $# -ge 2 ] || fail "$1 requires RUN"; run=$2; shift 2 ;;
        --run-type) [ $# -ge 2 ] || fail "$1 requires TYPE"; shift 2 ;;
        --trigger|-t) [ $# -ge 2 ] || fail "$1 requires TAG"; trigger_tag=$2; shift 2 ;;
        --filter|-f) [ $# -ge 3 ] || fail "$1 requires FILE TAG"; filter_config=$2; filter_tag=$3; shift 3 ;;
        --jobs) [ $# -ge 2 ] || fail "$1 requires N"; jobs=$2; shift 2 ;;
        --overwrite) overwrite=1; shift ;;
        --gpu) shift ;;
        --help|-h) usage; exit 0 ;;
        *) fail "unknown option: $1" ;;
    esac
done
[ -n "${run}" ] || fail "missing --run"; [ -n "${trigger_tag}" ] || fail "missing --trigger"; [ -n "${filter_tag}" ] || fail "missing filter tag"
[ -f "${filter_config}" ] || fail "filter config does not exist: ${filter_config}"
[ -x "${FILTER}" ] || fail "${FILTER} does not exist or is not executable"
[[ "${jobs}" =~ ^[0-9]+$ ]] && [ "${jobs}" -gt 0 ] || fail "--jobs must be positive"
input_dir="${ipath}/${run}/trigger"; output_dir="${opath}/${run}/trigger"; mkdir -p "${output_dir}"
inputs=("${input_dir}"/rings.${trigger_tag}.spill_*.root)
[ ${#inputs[@]} -gt 0 ] || fail "no ring spill files found: ${input_dir}/rings.${trigger_tag}.spill_*.root"
outputs=("${output_dir}"/filtered."${filter_tag}"."${trigger_tag}".spill_*.root)
merged_output="${output_dir}/filtered.${filter_tag}.${trigger_tag}.root"
if [ "${overwrite}" -ne 1 ] && [ -f "${merged_output}" ]; then
    echo " --- filtered output exists, skipping: ${merged_output}"
    exit 0
fi
pids=()
for input in "${inputs[@]}"; do
    spill=${input##*.spill_}; spill=${spill%.root}; output="${output_dir}/filtered.${filter_tag}.${trigger_tag}.spill_${spill}.root"
    if [ "${overwrite}" -ne 1 ] && [ -f "${output}" ]; then echo " --- filter output exists, skipping: ${output}"; continue; fi
    while [ ${#pids[@]} -ge "${jobs}" ]; do wait "${pids[0]}"; pids=("${pids[@]:1}"); done
    ( time -p "${FILTER}" --input "${input}" --output "${output}" --config "${filter_config}" ) > "${output%.root}.log" 2>&1 &
    pids+=("$!")
done
for pid in "${pids[@]}"; do wait "${pid}"; done
outputs=("${output_dir}"/filtered."${filter_tag}"."${trigger_tag}".spill_*.root)
[ ${#outputs[@]} -eq ${#inputs[@]} ] || fail "filtered spill count ${#outputs[@]} differs from input count ${#inputs[@]}"
if [ "${overwrite}" -eq 1 ] || [ ! -f "${merged_output}" ]; then
    command -v "${HADD}" >/dev/null 2>&1 || fail "${HADD} was not found"
    "${HADD}" -f "${merged_output}" "${outputs[@]}"
else
    echo " --- filtered output exists, skipping merge: ${merged_output}"
fi
echo " --- removing ${#outputs[@]} filtered spill files"
rm -f -- "${outputs[@]}"
remaining=("${output_dir}"/filtered."${filter_tag}"."${trigger_tag}".spill_*.root)
[ ${#remaining[@]} -eq 0 ] || fail "could not remove filtered spill files"
echo " --- filter workflow completed"
