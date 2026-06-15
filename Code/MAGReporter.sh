#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# ARGUMENTS
# =========================================================
BINS=""
OUTDIR=""
THREADS=8
R1=""
R2=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --bins) BINS="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --reads1) R1="$2"; shift 2 ;;
        --reads2) R2="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown arg $1"; exit 1 ;;
    esac
done

[[ -z "$BINS" || -z "$OUTDIR" ]] && { echo "[ERROR] bins/outdir required"; exit 1; }

QC="$OUTDIR"
LOG="$QC/logs/mag_qc.log"

mkdir -p "$QC"/{logs,checkm2,gtdbtk,rrna,trna,silva,bbmap,coverage,tmp}

echo "[QC] START" | tee -a "$LOG"

# =========================================================
# ENV CHECK
# =========================================================

echo "[QC] checking environment" | tee -a "$LOG"

[[ -z "${GTDBTK_DATA_PATH:-}" ]] && { echo "[ERROR] GTDBTK_DATA_PATH not set"; exit 1; }

# =========================================================
# SILVA DB RESOLUTION
# =========================================================

# Default (only if not provided)
if [[ -z "${SILVA_DB:-}" ]]; then
    SILVA_DB="$CONDA_PREFIX/SILVA_DB/silva_nr99"
fi

# If directory provided → resolve DB prefix
if [[ -d "$SILVA_DB" ]]; then
    DB=$(find "$SILVA_DB" -name "*.nsq" | head -n1 | sed 's/.nsq$//')

    [[ -z "$DB" ]] && {
        echo "[ERROR] No BLAST DB found in $SILVA_DB" | tee -a "$LOG"
        exit 1
    }

    SILVA_DB="$DB"
fi

# Validate DB prefix
if [[ ! -f "${SILVA_DB}.nsq" ]]; then
    echo "[ERROR] Invalid SILVA DB prefix: $SILVA_DB" | tee -a "$LOG"
    exit 1
fi

echo "[QC] Using SILVA DB: $SILVA_DB" | tee -a "$LOG"


# =========================================================
# THREAD SPLIT
# =========================================================
if (( THREADS <= 8 )); then
    JOBS=$THREADS
elif (( THREADS <= 20 )); then
    JOBS=$(( THREADS / 2 ))
else
    JOBS=$(( THREADS - 4 ))
fi

# =========================================================
# FILE DISCOVERY (SAFE)
# =========================================================
FILES=()
while IFS= read -r f; do
    FILES+=("$f")
done < <(find "$BINS" -type f -name "*.fa" 2>/dev/null)

if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No .fa files found in $BINS" | tee -a "$LOG"
    exit 1
fi

echo "[QC] Found ${#FILES[@]} MAGs" | tee -a "$LOG"
echo "[QC] ENTERING PIPELINE" | tee -a "$LOG"

# =========================================================
# RUNNER
# =========================================================
run_stage() {
    stage="$1"; shift

    echo "[QC] ===== $stage =====" | tee -a "$LOG"
    "$@" 2>&1 | tee -a "$LOG" || { echo "[ERROR] $stage failed"; exit 1; }
    echo "[QC] $stage DONE" | tee -a "$LOG"
}

run_parallel() {
    stage="$1"
    mode="$2"

    echo "[QC] ===== $stage =====" | tee -a "$LOG"
    echo "[QC] $stage → ${#FILES[@]} MAGs | Jobs: $JOBS" | tee -a "$LOG"

    parallel -j "$JOBS" --bar --joblog "$QC/logs/${stage}.joblog" '
        f="{}"
        b=$(basename "$f" .fa)

        case "'"$mode"'" in

        barrnap)
            gff="'"$QC"'/rrna/${b}.gff"
            if [[ ! -f "$gff" ]]; then
                kingdom=$(awk -v b="$b" '\''$1==b{print $2}'\'' "'"$domain_map"'")

                if [[ "$kingdom" == "arc" ]]; then
                    barrnap --kingdom arc "$f" > "$gff"
                else
                    barrnap --kingdom bac "$f" > "$gff"
                fi
            fi
            ;;

        trna)
            out="'"$QC"'/trna/${b}.txt"
            [[ -f "$out" ]] || tRNAscan-SE "$f" -o "$out"
            ;;

        silva)
            gff="'"$QC"'/rrna/${b}.gff"
            [[ ! -s "$gff" ]] && exit 0

            for TYPE in 16S 23S 5S; do
                tmp="'"$QC"'/tmp/${b}_${TYPE}.fa"
                tmp_gff="'"$QC"'/tmp/${b}_${TYPE}.gff"
                out="'"$QC"'/silva/${b}_${TYPE}.tsv"

                awk -v t="$TYPE" '\''$3=="rRNA" && $9 ~ t'\'' "$gff" > "$tmp_gff"
                bedtools getfasta -fi "$f" -bed "$tmp_gff" -fo "$tmp"

                blastn -query "$tmp" -db "'"$SILVA_DB"'" \
                    -out "$out" \
                    -outfmt "6 qseqid stitle pident length evalue bitscore" \
                    -max_target_seqs 1
            done
            ;;

        bbmap)
            out="'"$QC"'/bbmap/${b}.stats"
            [[ -s "$out" ]] || stats.sh in="$f" out="$out"
            ;;

        esac
    ' ::: "${FILES[@]}" 2>&1 | tee -a "$LOG"

    echo "[QC] $stage DONE" | tee -a "$LOG"
}

# =========================================================
# CHECKM2
# =========================================================
if [[ ! -s "$QC/checkm2/quality_report.tsv" ]]; then
    run_stage "CheckM2" \
        checkm2 predict -x fa -i "$BINS" -o "$QC/checkm2" --threads "$THREADS"
fi

# =========================================================
# GTDB
# =========================================================
if [[ ! -s "$QC/gtdbtk/gtdbtk.bac120.summary.tsv" ]]; then
    run_stage "GTDB-Tk" \
        gtdbtk classify_wf \
        --genome_dir "$BINS" \
        --out_dir "$QC/gtdbtk" \
        --cpus "$THREADS" \
        --extension fa
fi

# =========================================================
# COVERAGE
# =========================================================
if [[ -n "$R1" && -n "$R2" ]]; then
    if [[ ! -s "$QC/coverage/coverage.tsv" ]]; then
        run_stage "Coverage" \
            coverm genome \
            --genome-fasta-directory "$BINS" \
            --genome-fasta-extension fa \
            --coupled "$R1" "$R2" \
            --threads "$THREADS" \
            --methods mean covered_fraction variance \
            -o "$QC/coverage/coverage.tsv"
    fi
fi

# =========================================================
# DOMAIN MAP
# =========================================================
domain_map="$QC/domain_map.tsv"
> "$domain_map"

if [[ -f "$QC/gtdbtk/gtdbtk.bac120.summary.tsv" ]]; then
    awk 'NR>1 {print $1"\tbac"}' "$QC/gtdbtk/gtdbtk.bac120.summary.tsv" >> "$domain_map"
    echo "[QC] Domain map built" | tee -a "$LOG"
fi

if [[ -f "$QC/gtdbtk/gtdbtk.ar53.summary.tsv" ]]; then
    awk 'NR>1 {print $1"\tarc"}' "$QC/gtdbtk/gtdbtk.ar53.summary.tsv" >> "$domain_map"
    echo "[QC] Domain map built" | tee -a "$LOG"
fi

# =========================================================
# FUNCTIONS (NO EXPORT NEEDED)
# =========================================================
run_barrnap() {
    f="$1"; b=$(basename "$f" .fa)
    gff="$QC/rrna/${b}.gff"

    [[ -f "$gff" ]] && return

    kingdom=$(awk -v b="$b" '$1==b{print $2}' "$domain_map")

    if [[ "$kingdom" == "arc" ]]; then
        barrnap --kingdom arc "$f" > "$gff"
    else
        barrnap --kingdom bac "$f" > "$gff"
    fi
}

run_trna() {
    f="$1"; b=$(basename "$f" .fa)
    out="$QC/trna/${b}.txt"
    [[ -f "$out" ]] && return
    tRNAscan-SE "$f" -o "$out"
}

run_silva() {
    f="$1"; b=$(basename "$f" .fa)
    gff="$QC/rrna/${b}.gff"
    [[ ! -s "$gff" ]] && return

    for TYPE in 16S 23S 5S; do
        tmp="$QC/tmp/${b}_${TYPE}.fa"
        tmp_gff="$QC/tmp/${b}_${TYPE}.gff"
        out="$QC/silva/${b}_${TYPE}.tsv"

        awk -v t="$TYPE" '$3=="rRNA" && $9 ~ t' "$gff" > "$tmp_gff"
        bedtools getfasta -fi "$f" -bed "$tmp_gff" -fo "$tmp"

        blastn -query "$tmp" -db "$SILVA_DB" \
            -out "$out" \
            -outfmt "6 qseqid stitle pident length evalue bitscore" \
            -max_target_seqs 1
    done
}

run_bbmap() {
    f="$1"; b=$(basename "$f" .fa)
    out="$QC/bbmap/${b}.stats"
    [[ -f "$out" ]] && return
    stats.sh in="$f" out="$out"
}

# =========================================================
# PARALLEL STAGES
# =========================================================
run_parallel "Barrnap" barrnap
run_parallel "tRNA" trna
run_parallel "SILVA" silva
run_parallel "BBMap" bbmap

# =========================================================
# FINAL AGGREGATION (HARDENED + CLEAN + SAFE)
# =========================================================
export FINAL_BINSET="$BINS" QC
python3 << 'PYCODE'
#!/usr/bin/env python3

import os, glob, sys, re
from collections import defaultdict

QC = os.environ.get("QC")
BINSET = os.environ.get("FINAL_BINSET")

if not QC or not BINSET:
    sys.exit("[ERROR] Missing QC or FINAL_BINSET")

# =========================================================
# HELPERS (STRICT + SAFE)
# =========================================================
def clean(x):
    return str(x).replace("\n","").replace("\t","").strip()

def r2(x):
    try: return f"{float(x):.2f}"
    except: return "0.00"

def i(x):
    try: return str(int(float(x)))
    except: return "0"

def mb(x):
    try: return f"{float(x)/1e6:.2f}"
    except: return "0.00"

def mag_class(c,t):
    try: c,t=float(c),float(t)
    except: return "LQ"
    if c>=90 and t<=5: return "HQ"
    if c>=50 and t<=10: return "MQ"
    return "LQ"

# =========================================================
# BBMAP PARSER (MATCHES YOUR REAL OUTPUT)
# =========================================================
def parse_bbmap(fp):

    stats = dict(contigs=0,size=0,n50=0,l50=0,l50_type="unknown",gc=0,gc_std=0)

    if not os.path.exists(fp):
        return stats

    with open(fp) as f:
        lines = f.readlines()

    # ---- GC ----
    for i,line in enumerate(lines):
        if line.startswith("A\tC\tG\tT"):
            try:
                vals = lines[i+1].split()
                stats["gc"] = float(vals[7])*100
                stats["gc_std"] = float(vals[8])*100
            except:
                pass
            break

    for line in lines:
        try:
            # Contigs
            if "scaffold total" in line:
                stats["contigs"] = int(line.split()[-1])

            # Genome size
            elif "scaffold sequence total" in line:
                m = re.search(r'([\d.]+)\s*Mb', line)
                if m:
                    stats["size"] = int(float(m.group(1))*1e6)

            # N/L50
            elif "scaffold N/L50" in line:
                try:
                    part = line.split(":", 1)[1].strip()

                    # Flexible split (handles spaces around "/")
                    pieces = re.split(r'\s*/\s*', part)

                    if len(pieces) != 2:
                        raise ValueError("Unexpected N/L50 format")

                    n50_raw, l50_raw = pieces

                    # ---- N50 ----
                    stats["n50"] = int(float(n50_raw.strip()))

                    # ---- L50 ----
                    l50_raw = l50_raw.strip()

                    # Extract numeric portion
                    match = re.search(r'([\d.]+)', l50_raw)
                    if not match:
                        raise ValueError("Invalid L50 value")

                    l50_val = float(match.group(1))

                    # Detect unit
                    if "kbp" in l50_raw.lower():
                        stats["l50_type"] = "length"
                        stats["l50"] = int(l50_val * 1000)

                    elif "bp" in l50_raw.lower():
                        stats["l50_type"] = "length"
                        stats["l50"] = int(l50_val)

                    else:
                        stats["l50_type"] = "count"
                        stats["l50"] = int(l50_val)

                except Exception:
                    print(f"[WARN] Failed to parse N/L50 line: {line}")
        except:
            continue

    return stats

# =========================================================
# CHECKM2
# =========================================================
checkm={}
fp=f"{QC}/checkm2/quality_report.tsv"

if os.path.exists(fp):
    for l in open(fp):
        if l.startswith("Name"): continue
        p=l.strip().split("\t")
        if len(p)>=3:
            try:
                checkm[p[0]]=(float(p[1]),float(p[2]),float(p[5]) if len(p)>5 else 0)
            except:
                continue

# =========================================================
# GTDB
# =========================================================
gtdb={}
for f in ["gtdbtk.bac120.summary.tsv","gtdbtk.ar53.summary.tsv"]:
    fp=f"{QC}/gtdbtk/{f}"
    if not os.path.exists(fp): continue

    for l in open(fp):
        if l.startswith("user_genome"): continue
        p=l.strip().split("\t")
        if len(p)>1:
            gtdb[p[0]]=p[1] if p[1] else "Unclassified"

# =========================================================
# COVERAGE (STRICT)
# =========================================================
coverage={}
fp=f"{QC}/coverage/coverage.tsv"

if os.path.exists(fp):
    for l in open(fp):
        if l.startswith("genome"): continue
        p=l.strip().split()

        if len(p)<4: continue

        try:
            coverage[p[0]]=(float(p[1]),float(p[2]),float(p[3]))
        except:
            continue

# =========================================================
# RNA PARSING (FIXED)
# =========================================================
def rrna(fp):
    c16=c23=c5=0
    if not os.path.exists(fp): return 0,0,0

    for l in open(fp):
        if l.startswith("#"): continue

        if "16S" in l: c16+=1
        elif "23S" in l: c23+=1
        elif "5S" in l: c5+=1

    return c16,c23,c5

def trna(fp):
    if not os.path.exists(fp): return 0
    return sum(1 for l in open(fp) if not l.startswith("#"))

# =========================================================
# OUTPUT
# =========================================================
bins=sorted(glob.glob(f"{BINSET}/*.fa"))
out_fp=f"{QC}/final_qc.tsv"

header=[
"MAG_ID","GTDB","Completeness","Contamination","MAG_Quality","MAG_Class",
"Mean_Coverage","Breadth","Coverage_Variance",
"Contigs","Genome_Size(Mb)","N50","L50_value","L50_unit",
"GC_Content","GC_std","Coding_Density",
"16S","23S","5S","tRNA","rRNA_total",
"rRNA_consensus","rRNA_status","rRNA_source"
]

with open(out_fp,"w") as out:
    out.write("\t".join(header)+"\n")

    for f in bins:
        b=os.path.basename(f).replace(".fa","")

        comp,cont,cds=checkm.get(b,(0,0,0))
        magq=comp-5*cont

        cov_mean,cov_br,cov_var=coverage.get(b,(0,0,0))

        bb=parse_bbmap(f"{QC}/bbmap/{b}.stats")

        r16,r23,r5=rrna(f"{QC}/rrna/{b}.gff")
        t=trna(f"{QC}/trna/{b}.txt")
        rt=r16+r23+r5

        row=[
            b,
            gtdb.get(b,"Unclassified"),
            r2(comp),r2(cont),r2(magq),mag_class(comp,cont),
            r2(cov_mean),r2(cov_br),r2(cov_var),
            i(bb["contigs"]),mb(bb["size"]),i(bb["n50"]),i(bb["l50"]),bb["l50_type"],
            r2(bb["gc"]),r2(bb["gc_std"]),r2(cds),
            i(r16),i(r23),i(r5),i(t),i(rt),
            "NA","missing","none"
        ]

        row=[clean(x) for x in row]

        if len(row)!=len(header):
            print(f"[WARN] Skipping {b} (bad row)")
            continue

        out.write("\t".join(row)+"\n")

print(f"[QC] final_qc.tsv written → {out_fp}")
PYCODE
echo "[QC] COMPLETE" | tee -a "$LOG"
