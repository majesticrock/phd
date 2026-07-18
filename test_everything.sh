#!/usr/bin/env bash
set -euo pipefail

# Default jobs: try nproc, then sysctl, else 1
DEFAULT_JOBS="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)"

USAGE="$(cat <<EOF
Usage: $(basename "$0") [--step <step>] [-j N | --jobs N] [--build-only]

Options:
  --step        Which step to run. One of: all, PhdUtility, Hubbard, Continuum, LatticeCUT
                Default: all
  -j, --jobs    Number of parallel jobs passed to make (default: CPU count)
  --build-only  Only build (make -jN). Do not run the 'test' target.
Examples:
  $(basename "$0") --step all
  $(basename "$0") --step Hubbard -j 8
  $(basename "$0") --step Continuum --build-only
EOF
)"

STEP="all"
JOBS="$DEFAULT_JOBS"
BUILD_ONLY="false"

# simple arg parsing
while [[ $# -gt 0 ]]; do
  case "$1" in
    --step)
      shift
      STEP="${1:-}"
      shift
      ;;
    -j|--jobs)
      shift
      JOBS="${1:-}"
      shift
      ;;
    --build-only)
      BUILD_ONLY="true"
      shift
      ;;
    -h|--help)
      echo "$USAGE"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "$USAGE" >&2
      exit 2
      ;;
  esac
done

# project directories (relative to repo root)
declare -A PROJECT_DIRS=(
  ["PhdUtility"]="PhdUtility"
  ["FermionCommute"]="cpp/FermionCommute"
  ["Hubbard"]="cpp/Hubbard"
  ["Continuum"]="cpp/ContinuumSystem"
  ["LatticeCUT"]="cpp/LatticeCUT"
)

check_dir() {
  local dir="$1"
  if [[ ! -d "$dir" ]]; then
    echo "Project directory not found: $dir" >&2
    exit 3
  fi
}

run_make_build() {
  local dir="$1"
  local jflag="-j${JOBS}"
  check_dir "$dir"
  echo
  echo "================================================================"
  echo " Building: make ${jflag}  (dir: ${dir})"
  echo "================================================================"
  (cd "$dir" && make "${jflag}")
}

run_make_test() {
  local dir="$1"
  check_dir "$dir"
  echo
  echo "================================================================"
  echo " Running tests: make test  (dir: ${dir})"
  echo "================================================================"
  # run 'make test' without -j to avoid mixing parallel build output with test harness output
  (cd "$dir" && make test)
}

# Combined convenience: build then (optionally) run test
build_and_maybe_test() {
  local dir="$1"
  run_make_build "$dir"
  if [[ "$BUILD_ONLY" == "false" ]]; then
    run_make_test "$dir"
  fi
}

case "$STEP" in
  all)
    # order preserved from your old script
    build_and_maybe_test "${PROJECT_DIRS[PhdUtility]}"

    # Reinstall PhdUtility so that the most recent version will be used
    (cd "${PROJECT_DIRS[PhdUtility]}" && make install)
    # Recompute the expectation values for the projects
    (cd "${PROJECT_DIRS[FermionCommute]}" && make "-j${JOBS}")

    # Run the project tests
    build_and_maybe_test "${PROJECT_DIRS[Hubbard]}"
    build_and_maybe_test "${PROJECT_DIRS[Continuum]}"
    build_and_maybe_test "${PROJECT_DIRS[LatticeCUT]}"
    ;;
  PhdUtility|Hubbard|Continuum|LatticeCUT)
    build_and_maybe_test "${PROJECT_DIRS[$STEP]}"
    ;;
  *)
    echo "Unknown step: $STEP" >&2
    echo "$USAGE" >&2
    exit 2
    ;;
esac

echo
if [[ "$BUILD_ONLY" == "true" ]]; then
  echo "All requested builds completed successfully."
else
  echo "All requested tests completed successfully."
fi