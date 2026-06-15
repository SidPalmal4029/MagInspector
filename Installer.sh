#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# CONFIG
# =========================================================
ENV_NAME="MagInspector"
PYTHON_VER="3.12"
DEFAULT_THREADS=8
MAX_THREADS=16

# =========================================================
# HELP
# =========================================================
show_help() {
    cat << EOF
MAGReporter Installer

Usage:
  bash Installer.sh [options]

Options:
  -t, --threads INT    Number of threads to use (default: $DEFAULT_THREADS, max: $MAX_THREADS)
  -h, --help           Show this help message and exit

Examples:
  bash Installer.sh --threads 16
  bash Installer.sh -t 12
  THREADS=8 bash Installer.sh

Notes:
  - Threads are globally capped at $MAX_THREADS
  - Databases will be installed inside the conda environment
EOF
}


# =========================================================
# ARGUMENT PARSING
# =========================================================
THREADS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--threads)
            if [[ -n "${2:-}" && "$2" =~ ^[0-9]+$ ]]; then
                THREADS="$2"
                shift 2
            else
                echo "[ERROR] --threads requires a positive integer"
                exit 1
            fi
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            echo "Use --help for usage"
            exit 1
            ;;
    esac
done


# =========================================================
# THREAD NORMALIZATION (GLOBAL CAP)
# =========================================================
THREADS="${THREADS:-$DEFAULT_THREADS}"

THREADS=$(printf "%.0f" "$THREADS" 2>/dev/null || echo "$DEFAULT_THREADS")

if (( THREADS > MAX_THREADS )); then
    echo "[INFO] THREADS capped to $MAX_THREADS (requested: $THREADS)"
    THREADS=$MAX_THREADS
fi

(( THREADS < 1 )) && THREADS=1

echo "[INSTALL] Using $THREADS threads (max: $MAX_THREADS)"


# =========================================================
# INIT
# =========================================================
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "[INSTALL] Starting MAG QC setup"

# =========================================================
# CREATE ENV (IDEMPOTENT)
# =========================================================
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "[INSTALL] Environment exists → $ENV_NAME"
else
    echo "[INSTALL] Creating environment → $ENV_NAME"
    mamba create -y -n "$ENV_NAME" python="$PYTHON_VER"
fi

conda activate "$ENV_NAME"

# =========================================================
# CHANNELS (ENV-SCOPED)
# =========================================================
conda config --env --add channels defaults
conda config --env --add channels conda-forge
conda config --env --add channels bioconda

# =========================================================
# INSTALL TOOLS
# =========================================================
echo "[INSTALL] Installing tools"

mamba install -y \
    checkm2 \
    gtdbtk \
    coverm \
    barrnap \
    trnascan-se \
    bedtools \
    blast \
    bbmap \
    parallel \
    aria2 \
    seqkit

# =========================================================
# DIRECTORY STRUCTURE
# =========================================================
DB_ROOT="$CONDA_PREFIX/DB"
GTDB_DIR="$DB_ROOT/GTDBTK"
SILVA_DIR="$DB_ROOT/SILVA"
CHECKM2_DIR="$DB_ROOT/CHECKM2"

mkdir -p "$GTDB_DIR" "$SILVA_DIR" "$CHECKM2_DIR"

# =========================================================
# CHECKM2 DATABASE
# =========================================================
echo "[INSTALL] Setting up CheckM2 DB"

if [[ -f "$CHECKM2_DIR/database.dmnd" ]]; then
    echo "[INSTALL] CheckM2 DB exists → skipping"
else
    checkm2 database --download --path "$CHECKM2_DIR" || {
        echo "[ERROR] CheckM2 DB download failed"
        exit 1
    }
fi

# =========================================================
# GTDB-Tk DATABASE
# =========================================================
echo "[INSTALL] Setting up GTDB-Tk DB"

if [[ -n "$(ls -A "$GTDB_DIR" 2>/dev/null)" ]]; then
    echo "[INSTALL] GTDB already present → skipping"
else
    echo "[INSTALL] Downloading GTDB (~100GB)"

    download-db.sh -d "$GTDB_DIR" -t "$THREADS" || {
        echo "[ERROR] GTDB download failed"
        exit 1
    }
fi

# =========================================================
# SILVA DATABASE
# =========================================================
echo "[INSTALL] Setting up SILVA DB"

cd "$SILVA_DIR"

if [[ -f silva_nr99.nsq ]]; then
    echo "[INSTALL] SILVA DB exists → skipping"
else
    SILVA_URL="https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"

    DL_THREADS=$THREADS
    (( DL_THREADS < 2 )) && DL_THREADS=2

    echo "[INSTALL] Downloading SILVA"

    aria2c \
        -x "$DL_THREADS" \
        -s "$DL_THREADS" \
        -k 1M \
        --file-allocation=none \
        --continue=true \
        -o silva_nr99.fasta.gz \
        "$SILVA_URL"

    if [[ ! -s silva_nr99.fasta.gz ]]; then
        echo "[ERROR] SILVA download failed"
        exit 1
    fi

    echo "[INSTALL] Decompressing SILVA"
    gunzip -c silva_nr99.fasta.gz > silva_nr99.fasta

    echo "[INSTALL] Building BLAST DB"
    makeblastdb \
        -in silva_nr99.fasta \
        -dbtype nucl \
        -parse_seqids \
        -out silva_nr99
fi

cd - >/dev/null

# =========================================================
# CONDA ENV HOOKS (CRITICAL)
# =========================================================
echo "[INSTALL] Setting environment variables (env-local)"

ACTIVATE_DIR="$CONDA_PREFIX/etc/conda/activate.d"
DEACTIVATE_DIR="$CONDA_PREFIX/etc/conda/deactivate.d"

mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

cat << EOF > "$ACTIVATE_DIR/db_paths.sh"
export GTDBTK_DATA_PATH="$GTDB_DIR"
export SILVA_DB="$SILVA_DIR/silva_nr99"
export CHECKM2DB="$CHECKM2_DIR/database.dmnd"
EOF

cat << EOF > "$DEACTIVATE_DIR/db_paths.sh"
unset GTDBTK_DATA_PATH
unset SILVA_DB
unset CHECKM2DB
EOF

# =========================================================
# VALIDATION
# =========================================================
echo "[INSTALL] Validating setup"

[[ -f "$SILVA_DIR/silva_nr99.nsq" ]] || { echo "[ERROR] SILVA DB missing"; exit 1; }
[[ -f "$CHECKM2_DIR/database.dmnd" ]] || { echo "[ERROR] CheckM2 DB missing"; exit 1; }
[[ -d "$GTDB_DIR" ]] || { echo "[ERROR] GTDB missing"; exit 1; }

# =========================================================
# DONE
# =========================================================
echo ""
echo "[INSTALL] SUCCESS"
echo "----------------------------------------"
echo "Activate:"
echo "  conda activate $ENV_NAME"
echo ""
echo "Verify:"
echo "  bash verify.sh"
echo ""
echo "DBs (auto-set on activate):"
echo "  GTDBTK_DATA_PATH=$GTDB_DIR"
echo "  SILVA_DB=$SILVA_DIR/silva_nr99"
echo "  CHECKM2DB=$CHECKM2_DIR/database.dmnd"
