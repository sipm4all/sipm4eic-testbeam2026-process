#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)
MODEL_DIR="${ROOT_DIR}/process/source/timing/python"
BUILD_DIR="${ROOT_DIR}/process/build"

NEVENTS=${1:-1000}
PYTHON=${PYTHON:-python3}

events=$(mktemp /tmp/timing_estimator_events.XXXXXX.csv)
reference=$(mktemp /tmp/timing_estimator_reference.XXXXXX.csv)
trap 'rm -f "${events}" "${reference}"' EXIT

"${PYTHON}" "${ROOT_DIR}/process/source/timing/test/make_validation_events.py" \
    --events "${NEVENTS}" \
    --output "${events}"

"${PYTHON}" "${ROOT_DIR}/process/source/timing/test/python_reference.py" \
    --model-dir "${MODEL_DIR}" \
    --input "${events}" \
    --output "${reference}"

"${BUILD_DIR}/validate_timing_estimator" "${events}" "${reference}"
