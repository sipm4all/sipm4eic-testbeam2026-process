#! /usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

run_type="physics"
ipath="/data/2026-testbeam/process"
opath="/data/2026-testbeam/process"

TRIGGER="${ROOT_DIR}/process/bin/trigger"

run=""
TRIGGER_CONFIGS=()
TRIGGER_TAGS=()
TRIGGER_WINDOW=256
DEVICE_FILTER=()
DEVICE_ALL=1
INPUT_PREFIX=""
OVERWRITE=0

WRITE_LOGS=0
CLEAN_TRIGGERED_SPILLS=0

usage()
{
    cat <<EOF
usage:
  $0 --run RUN --trigger TRIGGER.conf TAG [options]

required:
  --run, -r RUN                  run name/directory
  --trigger, -t FILE TAG         trigger configuration and output tag; may be repeated

options:
  --run-type TYPE                accepted for symmetry with process.sh, default: physics
  --devices DEVICE ...           device subset used by process.sh, default: all
  --input-prefix PREFIX          merged spill prefix, default derived from --devices
  --overwrite                    overwrite existing trigger outputs instead of skipping them
  --window, -w VALUE             trigger frame window, default: 256
  --clean-triggered-spills       legacy option; split spill outputs are retained by this workflow
  --help, -h                     show this help message

examples:
  $0 --run 12345 --trigger process/config/trigger/trigger_range.conf range --window 256
  $0 --run 20260618-183625 --run-type testpulse --devices rdo-192 --trigger trigger.conf calibcheck_rdo-192 --window 32
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
        --input-prefix)
            [ $# -ge 2 ] || fail "$1 requires PREFIX"
            INPUT_PREFIX=$2
            shift 2
            ;;
        --overwrite)
            OVERWRITE=1
            shift
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
        --clean-triggered-spills)
            CLEAN_TRIGGERED_SPILLS=1
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
[ ${#TRIGGER_CONFIGS[@]} -gt 0 ] || fail "at least one --trigger FILE TAG is required"
[[ "${TRIGGER_WINDOW}" =~ ^[0-9]+([.][0-9]+)?$ ]] || fail "--window must be a positive number"
awk "BEGIN { exit !(${TRIGGER_WINDOW} > 0) }" || fail "--window must be greater than zero"

if [ ! -x "${TRIGGER}" ]; then
    fail "${TRIGGER} does not exist or is not executable"
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

irpath="${ipath}/${run}/process"
if [ ! -d "${irpath}" ]; then
   fail "${irpath} does not exist; run process.sh first"
fi
orpath="${opath}/${run}/trigger"
mkdir -p "${orpath}"

merge_prefix="${INPUT_PREFIX}"
if [ -z "${merge_prefix}" ]; then
    merge_prefix="aps.sorted"
    if [ "${DEVICE_ALL}" -ne 1 ]; then
        device_tag=$(IFS=_; echo "${DEVICE_FILTER[*]}")
        merge_prefix="aps.sorted.${device_tag}"
    fi
fi

spill_files=("${irpath}"/${merge_prefix}.spill_*.root)
if [ ${#spill_files[@]} -eq 0 ]; then
    fail "no merged split-spill files found: ${irpath}/${merge_prefix}.spill_*.root"
fi

spill_ids=()
while IFS= read -r spill_id; do
    spill_ids+=("${spill_id}")
done < <(
    for fpath in "${spill_files[@]}"; do
        fname=$(basename "${fpath}")
        spill_id=${fname#${merge_prefix}.spill_}
        spill_id=${spill_id%.root}
        echo "${spill_id}"
    done | sort -u
)

echo " --- trigger input prefix: ${merge_prefix}"
echo " --- trigger jobs started"
for i in "${!TRIGGER_CONFIGS[@]}"; do
    tag=${TRIGGER_TAGS[$i]}
    config=${TRIGGER_CONFIGS[$i]}
    echo " --- trigger config ${tag}: ${config}"

    trigger_pids=()
    for spill_id in "${spill_ids[@]}"; do
        merged_spill="${irpath}/${merge_prefix}.spill_${spill_id}.root"
        if [ ! -f "${merged_spill}" ]; then
            fail "merged spill file not found: ${merged_spill}"
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
        fail "no triggered spill files found for tag ${tag}"
    fi
    if [ ${#triggered_files[@]} -ne ${#spill_ids[@]} ]; then
        fail "trigger tag ${tag} has ${#triggered_files[@]} files, expected ${#spill_ids[@]}"
    fi

    echo " --- trigger tag ${tag}: kept ${#triggered_files[@]} spill files"
done

echo " --- trigger jobs completed"
