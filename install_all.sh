#!/usr/bin/env bash
# =============================================================================
# NeHub — install all tools (or a chosen subset)
#
# The install commands mirror .github/workflows/ci.yml, which is tested on every
# pull request. Keep the two in sync (or let CI call this script).
#
# Usage:
#   bash install/install_all.sh                    # install everything
#   bash install/install_all.sh smcpp relate       # install only these
#   bash install/install_all.sh --force phlash     # remove and reinstall
#   bash install/install_all.sh --tools-dir /path/to/tools
#
# Tools:
#   conda envs : basic_Ne  bcftools  ts_tools  phlash  smcpp
#   standalone : gone  relate  stairway  easysfs
#   groups     : conda  standalone  all (default)
#
# Options:
#   --tools-dir DIR   where standalone tools are installed         [./tools]
#   --force           remove an existing env/folder and reinstall
#   -h, --help        show this help
#
# System requirements (install with your package manager / module system):
#   conda, git, unzip, cmake, make, a C++ compiler, Java (for Stairway Plot),
#   zlib/bz2/lzma/curl development headers (for building Relate)
#
# Output:
#   conf/tools.config — paths of the standalone tools. Either add
#   `includeConfig 'conf/tools.config'` to nextflow.config (after params.config),
#   or use the --flags printed at the end.
#
# Each tool is installed independently: a failure is reported in the summary
# but does not stop the remaining tools.
# =============================================================================
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="$REPO_DIR/conf/envs"
TOOLS_DIR="$REPO_DIR/tools"
TOOLS_CONFIG="$REPO_DIR/conf/tools.config"
FORCE=0

STAIRWAY_VERSION="2.1.3"
SMCPP_VERSION="1.15.4"

ALL_CONDA=(basic_Ne bcftools ts_tools phlash smcpp)
ALL_STANDALONE=(gone relate stairway easysfs)

# ---------------------------------------------------------------- arguments --
REQUESTED=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tools-dir)  mkdir -p "$2"; TOOLS_DIR="$(cd "$2" && pwd)"; shift 2 ;;
        --force)      FORCE=1; shift ;;
        -h|--help)    sed -n '2,36p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        conda)        REQUESTED+=("${ALL_CONDA[@]}"); shift ;;
        standalone)   REQUESTED+=("${ALL_STANDALONE[@]}"); shift ;;
        all)          REQUESTED+=("${ALL_CONDA[@]}" "${ALL_STANDALONE[@]}"); shift ;;
        -*)           echo "Unknown option: $1" >&2; exit 2 ;;
        *)            REQUESTED+=("$1"); shift ;;
    esac
done
[[ ${#REQUESTED[@]} -eq 0 ]] && REQUESTED=("${ALL_CONDA[@]}" "${ALL_STANDALONE[@]}")
mkdir -p "$TOOLS_DIR"

# ------------------------------------------------------------------ helpers --
log()  { printf '\n==> %s\n' "$*"; }
warn() { printf 'WARNING: %s\n' "$*" >&2; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "'$1' is required but not found in PATH"; }

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

# Returns 0 if the env should be created, 1 if it already exists and is kept.
fresh_env() {
    local name=$1
    if env_exists "$name"; then
        if [[ $FORCE -eq 1 ]]; then
            echo "Removing existing env '$name' (--force)"
            conda env remove -n "$name" -y >/dev/null
        else
            echo "Env '$name' already exists — skipping (use --force to reinstall)"
            return 1
        fi
    fi
    return 0
}

# Returns 0 if the folder should be (re)created, 1 if it exists and is kept.
fresh_dir() {
    local dir=$1
    if [[ -e "$dir" ]]; then
        if [[ $FORCE -eq 1 ]]; then
            rm -rf "$dir"
        else
            echo "$dir already exists — skipping (use --force to reinstall)"
            return 1
        fi
    fi
    return 0
}

# ------------------------------------------------------------ conda tools --
install_basic_Ne() {
    if fresh_env basic_Ne; then
        conda env create --name basic_Ne --file "$ENV_DIR/basic_Ne.yml"
    fi
}

install_bcftools() {
    if fresh_env bcftools; then
        conda env create --name bcftools --file "$ENV_DIR/bcftools.yml"
    fi
}

install_ts_tools() {
    if fresh_env ts_tools; then
        conda env create --name ts_tools --file "$ENV_DIR/ts_tools.yml"
    fi
    # tszip and msprime are not in ts_tools.yml but are needed by the ts_*.py
    # scripts and the test-data generator (no-op if already installed)
    conda install -n ts_tools -y -c conda-forge tszip msprime
}

install_phlash() {
    # dinopy from bioconda avoids building it from source (fails with Cython 3)
    if fresh_env phlash; then
    conda env create -n ts_tools --file "$ENV_DIR/phlash.yml"
    fi
}

install_smcpp() {
    # SMC++ only builds for Python <= 3.10; setuptools<81 keeps pkg_resources
    if fresh_env smcpp; then
        conda create -n smcpp --file "$ENV_DIR/smcpp.yml"
    fi
    conda run -n smcpp smc++ version
}

# ------------------------------------------------------- standalone tools --
install_gone() {
    need git
    local dest="$TOOLS_DIR/GONE"
    if fresh_dir "$dest"; then
        git clone --depth 1 https://github.com/esrud/GONE.git "$dest"
    fi
    chmod -R +x "$dest/Linux/PROGRAMMES/"
}

install_relate() {
    need git; need cmake; need make
    local dest="$TOOLS_DIR/relate"
    fresh_dir "$dest" || return 0

    local src; src="$(mktemp -d)/relate_src"
    git clone --depth 1 https://github.com/MyersGroup/relate.git "$src"
    (cd "$src" && rm -rf build && mkdir build && cd build && cmake .. && make -j"$(nproc 2>/dev/null || echo 2)")

    # copy binaries into <tools>/relate/bin — mkdir first so the layout is the
    # same whether or not the target folder already existed
    mkdir -p "$dest/bin"
    cp "$src"/bin/* "$dest/bin/"
    rm -rf "$(dirname "$src")"

    test -x "$dest/bin/RelateFileFormats" || die "Relate not installed at $dest/bin"
    echo "Relate built: $(ls "$dest/bin" | tr '\n' ' ')"
}

install_stairway() {
    need git; need unzip
    command -v java >/dev/null 2>&1 || warn "java not found — Stairway Plot needs Java at run time"
    local dest="$TOOLS_DIR/stairway_plot_es"
    fresh_dir "$dest" || return 0

    local tmp; tmp="$(mktemp -d)"
    git clone --depth 1 https://github.com/xiaoming-liu/stairway-plot-v2.git "$tmp/src"
    unzip -o -q "$tmp/src/stairway_plot_v${STAIRWAY_VERSION}.zip" -d "$tmp/src/stairway_src"
    mv "$tmp/src/stairway_src/stairway_plot_v${STAIRWAY_VERSION}/stairway_plot_es" "$dest"
    rm -rf "$tmp"
    echo "Stairway installed: $(ls "$dest" | tr '\n' ' ')"
}

install_easysfs() {
    need git
    local dest="$TOOLS_DIR/easySFS"
    if fresh_dir "$dest"; then
        git clone --depth 1 https://github.com/isaacovercast/easySFS.git "$dest"
    fi
    chmod +x "$dest/easySFS.py"
}

# --------------------------------------------------------------------- main --
for t in "${REQUESTED[@]}"; do
    if [[ " ${ALL_CONDA[*]} " == *" $t "* ]]; then
        need conda
        # shellcheck disable=SC1091
        source "$(conda info --base)/etc/profile.d/conda.sh"
        break
    fi
done

declare -A STATUS
for tool in "${REQUESTED[@]}"; do
    fn="install_${tool}"
    if ! declare -F "$fn" >/dev/null; then
        warn "unknown tool '$tool' (valid: ${ALL_CONDA[*]} ${ALL_STANDALONE[*]})"
        STATUS[$tool]="UNKNOWN"
        continue
    fi
    log "Installing $tool"
    # subshell + set -e: stop this tool at its first error, then carry on
    if ( set -e; "$fn" ); then STATUS[$tool]="OK"; else STATUS[$tool]="FAILED"; fi
done

# ---------------------------------------------------- write conf/tools.config --
# Same parameter names as ci_params.json in the CI workflow
declare -A PARAMS
[[ -d "$TOOLS_DIR/GONE" ]]              && PARAMS[GONE.installation_path]="$TOOLS_DIR/GONE"
[[ -d "$TOOLS_DIR/relate/bin" ]]        && PARAMS[Relate.Relate_path]="$TOOLS_DIR/relate"
[[ -f "$TOOLS_DIR/easySFS/easySFS.py" ]] && PARAMS[stairway.easysfs_path]="$TOOLS_DIR/easySFS/easySFS.py"
[[ -d "$TOOLS_DIR/stairway_plot_es" ]]  && PARAMS[stairway.stairway_path]="$TOOLS_DIR/stairway_plot_es"

if [[ ${#PARAMS[@]} -gt 0 ]]; then
    {
        echo "// Generated by install/install_all.sh on $(date '+%Y-%m-%d %H:%M')"
        echo "// Add to nextflow.config after params.config:  includeConfig 'conf/tools.config'"
        for key in $(printf '%s\n' "${!PARAMS[@]}" | sort); do
            echo "params.${key} = '${PARAMS[$key]}'"
        done
    } > "$TOOLS_CONFIG"
    log "Wrote $TOOLS_CONFIG"
    echo "Equivalent command-line flags:"
    for key in $(printf '%s\n' "${!PARAMS[@]}" | sort); do
        echo "    --${key} \"${PARAMS[$key]}\""
    done
fi

# ------------------------------------------------------------------- summary --
log "Summary"
rc=0
for tool in "${REQUESTED[@]}"; do
    printf '  %-10s %s\n' "$tool" "${STATUS[$tool]}"
    [[ ${STATUS[$tool]} != OK ]] && rc=1
done
exit $rc