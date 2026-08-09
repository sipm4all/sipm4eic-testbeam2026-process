#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

ipath="/data/2026-testbeam/process"
opath="/data/2026-testbeam/process"

TRIGGER="${ROOT_DIR}/process/bin/trigger"
CONFIG="${ROOT_DIR}/process/config/trigger/trigger_range.conf"
WINDOW=256

WRITE_LOGS=0

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

if [ ! -x "${TRIGGER}" ]; then
    echo " ${TRIGGER} does not exist or is not executable "
    exit 1
fi
if [ ! -f "${CONFIG}" ]; then
    echo " ${CONFIG} does not exist "
    exit 1
fi

spill_files=("${irpath}"/aps.sorted.spill_*.root)
if [ ${#spill_files[@]} -eq 0 ]; then
    echo " --- no split-spill files found "
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

echo " --- trigger jobs started "
pids=()
for spill_id in "${spill_ids[@]}"; do

    run_job "${orpath}/triggered.spill_${spill_id}.log" bash -c '
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
        "${irpath}/aps.sorted.spill_${spill_id}.root" \
        "${orpath}/triggered.spill_${spill_id}.root" \
        "${CONFIG}" \
        "${WINDOW}" &
    pids+=($!)
done

wait "${pids[@]}"
echo " --- trigger jobs completed "
