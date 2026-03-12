#!/usr/bin/env zsh
# =============================================================================
#  run_pipeline.zsh
#
#  phi SDME pipeline — two paths:
#
#  ── PATH 1: Closure test (no real data, no FastMC) ───────────────────────────
#    Generates events, extracts pseudo-data, fits, compares to truth params.
#
#    ./run_pipeline.zsh --mode closure
#    ./run_pipeline.zsh --mode closure --config config/schc.cfg
#    ./run_pipeline.zsh --mode closure --stat-err 0.03 --perturb 0.15
#
#  ── PATH 2: Real experimental data ───────────────────────────────────────────
#    Uses real CLAS12 data + MC acceptance from fastMC.jar.
#
#    # Step A: generate MC and run fastMC (do once):
#    ./run_pipeline.zsh --mode gen-fastmc
#
#    # Step B: fit real data against MC acceptance:
#    ./run_pipeline.zsh --mode real-data \
#        --data-xsec data/xsec.dat --data-bsa data/bsa.dat
#
#    # Or with moments:
#    ./run_pipeline.zsh --mode real-data \
#        --data-xsec data/xsec.dat --data-bsa data/bsa.dat \
#        --data-mom  data/moments.dat
#
#  ── Utility modes ────────────────────────────────────────────────────────────
#    ./run_pipeline.zsh --mode gen           # generate LUND only
#    ./run_pipeline.zsh --mode analysis      # analysis plots from existing LUND
#    ./run_pipeline.zsh --mode gen-analysis  # gen + analysis plots (truth-level)
#    ./run_pipeline.zsh --mode build         # build only
#
#  ── Options (all modes) ──────────────────────────────────────────────────────
#    --config    <cfg>     Physics parameter config (default: config/default.cfg)
#    --root      <path>    ROOT installation path (if not in PATH)
#    --jobs      <N>       Parallel make jobs (default: 4)
#    --skip-build          Skip cmake/make (use existing binaries)
#    --output-dir <dir>    Output directory (default: <project>/output)
#    --use-fastmc          Use fastMC.jar for analysis plots (--mode analysis)
#    --fastmc-jar <path>   Path to fastMC.jar (default: Gavalian's path on ifarm)
#    --stat-err  <frac>    Pseudo-data stat error fraction for closure (default: 0.05)
#    --perturb   <frac>    Starting offset for closure test (default: 0.10)
#    --data-xsec <file>    Real xsec data file (--mode real-data)
#    --data-bsa  <file>    Real BSA data file  (--mode real-data)
#    --data-mom  <file>    Real moments file   (--mode real-data, optional)
#
# =============================================================================

setopt ERR_EXIT        # exit on any error
setopt PIPE_FAIL       # fail if any pipe segment fails

# ── Colours ────────────────────────────────────────────────────────────────────
autoload -U colors && colors
info()  { print -P "%F{cyan}[INFO]%f  $*" }
ok()    { print -P "%F{green}[OK]%f    $*" }
warn()  { print -P "%F{yellow}[WARN]%f  $*" >&2 }
err()   { print -P "%F{red}[ERROR]%f $*" >&2; exit 1 }
sep()   { print -P "%F{white}$(printf '─%.0s' {1..60})%f" }

# ── Defaults ───────────────────────────────────────────────────────────────────
SCRIPT_DIR="${0:A:h}"
PROJECT_DIR="${SCRIPT_DIR}"
BUILD_DIR="${PROJECT_DIR}/build"
CONFIG="${PROJECT_DIR}/config/default.cfg"
MODE="closure"              # closure | real-data | gen | gen-fastmc | gen-analysis | analysis | build
ROOT_PATH=""
JOBS=4
SKIP_BUILD=false
OUTPUT_DIR="${PROJECT_DIR}/output"
LUND_FILE="${OUTPUT_DIR}/events.lund"
LUND_FASTMC_FILE="${OUTPUT_DIR}/events_fastmc.lund"
LOG_DIR="${OUTPUT_DIR}/logs"
USE_FASTMC=false
FASTMC_JAR="/home/gavalian/Software/fastMC/fastMC.jar"
DATA_XSEC=""
DATA_BSA=""
DATA_MOM=""
STAT_ERR="0.05"
PERTURB="0.10"
THREADS=0                   # 0 = let OpenMP decide (uses all available cores)

# ── Argument parsing ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)        MODE="$2";         shift 2 ;;
        --config)      CONFIG="$2";       shift 2 ;;
        --root)        ROOT_PATH="$2";    shift 2 ;;
        --jobs)        JOBS="$2";         shift 2 ;;
        --skip-build)  SKIP_BUILD=true;   shift   ;;
        --output-dir)  OUTPUT_DIR="$2";  shift 2 ;;
        --use-fastmc)  USE_FASTMC=true;   shift   ;;
        --fastmc-jar)  FASTMC_JAR="$2";  shift 2 ;;
        --data-xsec)   DATA_XSEC="$2";  shift 2 ;;
        --data-bsa)    DATA_BSA="$2";   shift 2 ;;
        --data-mom)    DATA_MOM="$2";   shift 2 ;;
        --stat-err)    STAT_ERR="$2";   shift 2 ;;
        --perturb)     PERTURB="$2";    shift 2 ;;
        --threads)     THREADS="$2";    shift 2 ;;
        --help|-h)
            sed -n '3,48p' "$0"
            exit 0 ;;
        *)  err "Unknown argument: $1" ;;
    esac
done

# Derived paths (after OUTPUT_DIR may have been overridden)
LUND_FILE="${OUTPUT_DIR}/events.lund"
LUND_FASTMC_FILE="${OUTPUT_DIR}/events_fastmc.lund"
LOG_DIR="${OUTPUT_DIR}/logs"

# ── Sanity checks ──────────────────────────────────────────────────────────────
[[ -f "${CONFIG}" ]]              || err "Config not found: ${CONFIG}"
[[ -d "${PROJECT_DIR}/include" ]] || err "Bad project dir: ${PROJECT_DIR}"

if [[ "${MODE}" == "real-data" ]]; then
    [[ -n "${DATA_XSEC}" ]] || err "--mode real-data requires --data-xsec <file>"
    [[ -n "${DATA_BSA}"  ]] || err "--mode real-data requires --data-bsa  <file>"
    [[ -f "${DATA_XSEC}" ]] || err "xsec file not found: ${DATA_XSEC}"
    [[ -f "${DATA_BSA}"  ]] || err "BSA  file not found: ${DATA_BSA}"
    [[ -z "${DATA_MOM}"  ]] || [[ -f "${DATA_MOM}" ]] || err "Moments file not found: ${DATA_MOM}"
fi

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

# ── Banner ─────────────────────────────────────────────────────────────────────
sep
print -P "%F{magenta}  phi SDME Pipeline  —  CLAS12 RGA phi electroproduction%f"
sep
info "Mode        : ${MODE}"
info "Config      : ${CONFIG}"

# ── Thread control ─────────────────────────────────────────────────────────────
if [[ "${THREADS}" -gt 0 ]]; then
    export OMP_NUM_THREADS="${THREADS}"
    info "Threads     : ${THREADS} (OMP_NUM_THREADS=${THREADS})"
else
    # Let OpenMP use all available cores; report how many that is
    local _ncpu
    _ncpu=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo "?")
    info "Threads     : auto (${_ncpu} cores available)"
fi
info "Output dir  : ${OUTPUT_DIR}"
if [[ "${MODE}" == "closure" ]]; then
    info "Stat error  : ${STAT_ERR}  |  Start perturb: ${PERTURB}"
elif [[ "${MODE}" == "real-data" ]]; then
    info "xsec data   : ${DATA_XSEC}"
    info "BSA  data   : ${DATA_BSA}"
    [[ -n "${DATA_MOM}" ]] && info "Moments     : ${DATA_MOM}"
    info "MC LUND     : ${LUND_FASTMC_FILE}  (fastmc=${USE_FASTMC})"
fi
sep

# =============================================================================
# STEP 0 — Source ROOT
# =============================================================================
step0_source_root() {
    sep
    info "STEP 0 — Locating ROOT"
    sep

    if [[ -n "${ROOT_PATH}" ]]; then
        local thisroot="${ROOT_PATH}/bin/thisroot.zsh"
        [[ -f "${thisroot}" ]] || thisroot="${ROOT_PATH}/bin/thisroot.sh"
        [[ -f "${thisroot}" ]] || err "Cannot find thisroot.sh under ${ROOT_PATH}"
        source "${thisroot}"
        ok "Sourced ROOT from ${ROOT_PATH}"
    elif command -v root-config &>/dev/null; then
        ok "ROOT already in PATH: $(root-config --version)"
    else
        # Try common installation locations
        local candidates=(
            "/opt/root/bin/thisroot.zsh"
            "/usr/local/root/bin/thisroot.zsh"
            "${HOME}/root/bin/thisroot.zsh"
            "/cvmfs/sft.cern.ch/lcg/releases/ROOT/latest/x86_64-centos7-gcc8-opt/bin/thisroot.zsh"
        )
        local found=false
        for c in "${candidates[@]}"; do
            if [[ -f "${c}" ]]; then
                source "${c}"
                ok "Auto-sourced ROOT from ${c}"
                found=true
                break
            fi
        done
        if ! ${found}; then
            err "ROOT not found. Either:\n  • Add ROOT to PATH before running this script\n  • Pass --root /path/to/root to this script"
        fi
    fi

    # Verify
    local rootver
    rootver=$(root-config --version 2>/dev/null) \
        || err "root-config not working after sourcing ROOT"
    info "ROOT version: ${rootver}"
}

# =============================================================================
# STEP 1 — Build
# =============================================================================
step1_build() {
    sep
    info "STEP 1 — Building"
    sep

    if ${SKIP_BUILD}; then
        warn "--skip-build set; skipping cmake/make"
        [[ -x "${BUILD_DIR}/phi_gen" ]] || err "phi_gen not found in ${BUILD_DIR}. Run without --skip-build first."
        [[ -x "${BUILD_DIR}/phi_fit" ]] || err "phi_fit not found in ${BUILD_DIR}."
        ok "Using existing binaries in ${BUILD_DIR}"
        return
    fi

    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"

    info "Running cmake..."
    cmake "${PROJECT_DIR}" \
          -DCMAKE_BUILD_TYPE=Release \
          2>&1 | tee "${LOG_DIR}/cmake.log" \
        || err "cmake failed — check ${LOG_DIR}/cmake.log"

    info "Building with make -j${JOBS}..."
    make -j"${JOBS}" 2>&1 | tee "${LOG_DIR}/make.log" \
        || err "make failed — check ${LOG_DIR}/make.log"

    [[ -x "./phi_gen" ]]      || err "phi_gen binary not produced"
    [[ -x "./phi_fit" ]]      || err "phi_fit binary not produced"
    [[ -x "./phi_analysis" ]] || err "phi_analysis binary not produced"

    ok "Build complete: phi_gen  phi_fit  phi_analysis"
    cd "${PROJECT_DIR}"
}

# =============================================================================
# STEP 2 — Generate LUND events
# =============================================================================
step2_generate() {
    sep
    info "STEP 2 — Generating LUND events"
    sep

    info "Config  : ${CONFIG}"
    info "Output  : ${LUND_FILE}"

    local gen_log="${LOG_DIR}/gen.log"

    # Run generator; progress → log, events → LUND file
    "${BUILD_DIR}/phi_gen" "${CONFIG}" \
        > "${LUND_FILE}" \
        2> "${gen_log}" \
        || err "phi_gen failed — check ${gen_log}"

    local n_events
    n_events=$(grep -c '^[[:space:]]*[0-9]' "${LUND_FILE}" 2>/dev/null || echo "?")
    ok "Generation complete"
    info "LUND file : ${LUND_FILE}"
    info "Gen log   : ${gen_log}"

    # Print last few lines of gen log (summary stats)
    sep
    info "Generator summary:"
    tail -5 "${gen_log}"
}

# =============================================================================
# STEP 2b — Analysis: dσ/dt' and BSA plots from LUND (no FastMC)
# =============================================================================
# =============================================================================
# STEP 2b — Analysis: dσ/dt' and BSA plots from LUND (optionally via FastMC)
# =============================================================================
step2b_analysis() {
    sep
    if ${USE_FASTMC}; then
        info "STEP 2b — Running LUND analysis WITH FastMC detector filter"
    else
        info "STEP 2b — Running LUND analysis (truth-level, no FastMC)"
    fi
    sep

    [[ -f "${LUND_FILE}" ]] || err "LUND file not found: ${LUND_FILE}. Run --mode gen first."

    local ana_input="${LUND_FILE}"
    local ana_log="${LOG_DIR}/analysis.log"
    local ana_extra_flags=""

    # ── Optional: run fastMC.jar on the LUND file first ──────────────────────
    if ${USE_FASTMC}; then
        info "FastMC jar  : ${FASTMC_JAR}"
        [[ -f "${FASTMC_JAR}" ]] || err "fastMC.jar not found: ${FASTMC_JAR}\n  → Run without --use-fastmc for truth-level analysis, or specify path with --fastmc-jar"

        which java &>/dev/null || err "java not found in PATH (run: module load clas12)"

        info "Running fastMC.jar on ${LUND_FILE}..."
        local fastmc_log="${LOG_DIR}/fastmc.log"
        java -jar "${FASTMC_JAR}" "${LUND_FILE}" "${LUND_FASTMC_FILE}" \
            2>&1 | tee "${fastmc_log}" \
            || err "fastMC.jar failed — check ${fastmc_log}"

        [[ -f "${LUND_FASTMC_FILE}" ]] || err "fastMC.jar did not produce output: ${LUND_FASTMC_FILE}"
        ok "fastMC.jar complete → ${LUND_FASTMC_FILE}"

        ana_input="${LUND_FASTMC_FILE}"
        ana_extra_flags="--fastmc"
    fi

    info "Input  : ${ana_input}"
    info "Output : ${OUTPUT_DIR}"
    info "Log    : ${ana_log}"
    ${USE_FASTMC} && info "Mode   : FastMC filtered (e'→ECAL, p→DC, K±→FTOF)"
    sep

    "${BUILD_DIR}/phi_analysis" "${ana_input}" \
        ${ana_extra_flags} \
        --output-dir "${OUTPUT_DIR}" \
        2>&1 | tee "${ana_log}" \
        || err "phi_analysis failed — check ${ana_log}"

    ok "Analysis complete"
    info "Cross-section plots : ${OUTPUT_DIR}/xsec/"
    info "BSA plots           : ${OUTPUT_DIR}/bsa/"
    info "ROOT file           : ${OUTPUT_DIR}/phi_analysis.root"
}

# =============================================================================
# STEP 3a — PATH 1: Closure test fit
# =============================================================================
step3_closure() {
    sep
    info "STEP 3 — PATH 1: Closure test fit"
    info "  LUND  : ${LUND_FILE}"
    info "  Config: ${CONFIG}"
    info "  stat_err=${STAT_ERR}  perturb=${PERTURB}"
    sep

    local fit_log="${LOG_DIR}/fit_closure.log"

    "${BUILD_DIR}/phi_fit" \
        --closure       "${LUND_FILE}" \
        --config        "${CONFIG}"   \
        --output-dir    "${OUTPUT_DIR}/closure" \
        --stat-err      "${STAT_ERR}" \
        --perturb       "${PERTURB}"  \
        2>&1 | tee "${fit_log}" \
        || err "phi_fit (closure) failed — check ${fit_log}"

    ok "Closure fit complete"
    info "Results → ${OUTPUT_DIR}/closure/"
    info "Plots   → ${OUTPUT_DIR}/closure/plots/"
    info "Log     → ${fit_log}"
    sep
    print -P "%F{yellow}  Check the pull table above: all |pull| < 1 = good closure.%f"
    print -P "%F{yellow}  Fit plots saved to: ${OUTPUT_DIR}/closure/plots/%f"
    sep
}

# =============================================================================
# STEP 3b — Run fastMC.jar on generated LUND (for real-data path)
# =============================================================================
step3b_run_fastmc() {
    sep
    info "STEP 3b — Running fastMC.jar"
    info "  Input  : ${LUND_FILE}"
    info "  Output : ${LUND_FASTMC_FILE}"
    info "  Jar    : ${FASTMC_JAR}"
    sep

    [[ -f "${LUND_FILE}"  ]] || err "LUND file not found: ${LUND_FILE}. Run --mode gen first."
    [[ -f "${FASTMC_JAR}" ]] || err "fastMC.jar not found: ${FASTMC_JAR}"
    which java &>/dev/null   || err "java not in PATH (run: module load clas12)"

    local fmc_log="${LOG_DIR}/fastmc.log"
    java -jar "${FASTMC_JAR}" "${LUND_FILE}" "${LUND_FASTMC_FILE}" \
        2>&1 | tee "${fmc_log}" \
        || err "fastMC.jar failed — check ${fmc_log}"

    [[ -f "${LUND_FASTMC_FILE}" ]] || err "fastMC output not found: ${LUND_FASTMC_FILE}"
    ok "fastMC.jar → ${LUND_FASTMC_FILE}"
}

# =============================================================================
# STEP 3c — PATH 2: Real data fit
# =============================================================================
step3_realdata() {
    sep
    info "STEP 3 — PATH 2: Real data fit"
    info "  xsec    : ${DATA_XSEC}"
    info "  BSA     : ${DATA_BSA}"
    [[ -n "${DATA_MOM}" ]] && info "  moments : ${DATA_MOM}"
    info "  MC LUND : ${LUND_FASTMC_FILE}  (fastmc=${USE_FASTMC})"
    info "  Config  : ${CONFIG}"
    sep

    local mc_lund_flag=""
    ${USE_FASTMC} && mc_lund_flag="--fastmc"

    local fit_log="${LOG_DIR}/fit_realdata.log"
    local fit_args=(
        --data       "${DATA_XSEC}" "${DATA_BSA}"
        --mc-lund    "${LUND_FASTMC_FILE}"
        --config     "${CONFIG}"
        --output-dir "${OUTPUT_DIR}/real"
    )
    ${USE_FASTMC} && fit_args+=(--fastmc)
    [[ -n "${DATA_MOM}" ]] && fit_args+=(--moments "${DATA_MOM}")

    "${BUILD_DIR}/phi_fit" "${fit_args[@]}" \
        2>&1 | tee "${fit_log}" \
        || err "phi_fit (real data) failed — check ${fit_log}"

    ok "Real data fit complete"
    info "Results → ${OUTPUT_DIR}/real/"
    info "Log     → ${fit_log}"
}

# =============================================================================
# STEP 4 — Summary
# =============================================================================
step4_summary() {
    sep
    print -P "%F{green}  Pipeline complete!%f"
    sep
    info "Output: ${OUTPUT_DIR}"
    ls -lh "${OUTPUT_DIR}" 2>/dev/null | awk 'NR>1{print "  " $0}'
    sep
    case "${MODE}" in
        closure)
            info "Next: inspect ${OUTPUT_DIR}/closure/fit_result_params.cfg"
            info "      compare to ${OUTPUT_DIR}/closure/fit_truth_params.cfg"
            info "      all |pull| < 1 → fit procedure validated"
            info "Plots: ${OUTPUT_DIR}/closure/plots/"
            info "       xsec_Q2_*.pdf   — dσ/dt fit + pull panel per Q² bin"
            info "       bsa_Q2_*.pdf    — A_LU  fit + pull panel per Q² bin"
            info "       param_summary.pdf — parameter pulls summary"
            ;;
        real-data|gen-fastmc)
            info "Next: inspect ${OUTPUT_DIR}/real/fit_result_params.cfg"
            info "      check chi2/ndf in ${LOG_DIR}/fit_realdata.log"
            info "Plots: ${OUTPUT_DIR}/real/plots/"
            ;;
        gen-analysis)
            info "Next: inspect ${OUTPUT_DIR}/xsec/ and ${OUTPUT_DIR}/bsa/"
            info "      run --mode closure to validate fit procedure"
            ;;
    esac
    sep
}

# =============================================================================
# Main — dispatch by mode
# =============================================================================
step0_source_root

case "${MODE}" in

    closure)
        # PATH 1: generate truth events → closure fit → verify param recovery
        step1_build
        step2_generate
        step3_closure
        step4_summary
        ;;

    real-data)
        # PATH 2: fit real data against existing fastMC-filtered MC LUND
        # Expects: --data-xsec, --data-bsa, and events_fastmc.lund already present
        step1_build
        step3_realdata
        step4_summary
        ;;

    gen-fastmc)
        # PATH 2 prep: generate MC events + run fastMC.jar
        # After this, run --mode real-data with your real data files
        step1_build
        step2_generate
        step3b_run_fastmc
        step4_summary
        ;;

    gen)
        # Generate LUND only
        step1_build
        step2_generate
        step4_summary
        ;;

    gen-analysis)
        # Generate + analysis plots (truth-level)
        step1_build
        step2_generate
        step2b_analysis
        step4_summary
        ;;

    analysis)
        # Analysis plots from existing LUND (truth or fastMC)
        step1_build
        step2b_analysis
        step4_summary
        ;;

    build)
        step1_build
        ok "Build complete."
        ;;

    *)
        err "Unknown mode '${MODE}'.\nValid: closure | real-data | gen-fastmc | gen | gen-analysis | analysis | build"
        ;;
esac
