#!/usr/bin/env bash
# Strict, resumable SOAPdenovo2-based T-DNA workflow.
set -Eeuo pipefail
umask 002
export CONDA_SOLVER=classic
# Some SGE environments serialize shell functions incompletely. They are not
# needed by this workflow and are removed before subprocesses are launched.
for exported_func in $(env | sed -n 's/^BASH_FUNC_\([^=]*\)%%=.*/\1/p'); do unset -f "$exported_func" 2>/dev/null || true; done

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
: "${PROJECT:?Set PROJECT to the T-DNA analysis project directory}"
: "${SAMPLE:?Set SAMPLE to a filesystem-safe sample identifier}"
: "${SOAP_ROOT:?Set SOAP_ROOT to the runtime output directory}"
: "${BAM_PRIMARY:?Set BAM_PRIMARY to the primary coordinate-sorted BAM}"
: "${REF:?Set REF to the indexed plant reference FASTA}"
: "${TDNA:?Set TDNA to the supplied T-DNA FASTA}"
ROOT=$SOAP_ROOT
BAM_SECONDARY=${BAM_SECONDARY:-}
BAM=$BAM_PRIMARY
if [[ -n "$BAM_SECONDARY" ]]; then
  BAM="$ROOT/00.manifest/$SAMPLE.merged.rmdup.bam"
elif [[ "${NORMALIZE_BAM_HEADER:-0}" == 1 ]]; then
  BAM="$ROOT/00.manifest/$SAMPLE.sanitized.input.bam"
fi
SAMTOOLS=${SAMTOOLS:-samtools}
CONDA_BIN=${CONDA_BIN:-conda}
FALLBACK_TOOL_BIN=${FALLBACK_TOOL_BIN:-}
TOOL_ENV=${TOOL_ENV:-$ROOT/envs/assembly}
EXPECTED_TDNA_BP=${EXPECTED_TDNA_BP:-}
TDNA_SCOPE_NOTE=${TDNA_SCOPE_NOTE:-Not specified}
NORMALIZE_BAM_HEADER=${NORMALIZE_BAM_HEADER:-0}
THREADS=${NSLOTS:-24}
export PATH="$ROOT/envs/bin:$TOOL_ENV/bin:$PATH"
LOG="$ROOT/logs/pipeline.events.tsv"
PROGRESS="$ROOT/logs/progress.tsv"
HB="$ROOT/logs/heartbeat.tsv"
mkdir -p "$ROOT"/{00.manifest,01.qnames,02.fastq,03.assembly,04.alignments,05.junctions,06.validation,07.report,envs,logs,tmp,scripts}

now(){ date '+%F %T %z'; }
event(){ printf '%s\t%s\t%s\t%s\n' "$(now)" "${JOB_ID:-interactive}" "$1" "$2" >> "$LOG"; }
heartbeat(){ while :; do printf '%s\t%s\t%s\n' "$(now)" "${JOB_ID:-interactive}" "${CURRENT_STAGE:-startup}" >> "$HB"; sleep 60; done; }
heartbeat & HBPID=$!
trap 'rc=$?; kill "$HBPID" 2>/dev/null || true; event "${CURRENT_STAGE:-startup}" "FAILED rc=$rc line=$LINENO"; exit $rc' ERR INT TERM
trap 'kill "$HBPID" 2>/dev/null || true' EXIT

run_stage(){
  local stage=$1; shift
  CURRENT_STAGE=$stage; export CURRENT_STAGE
  if [[ -s "$ROOT/logs/$stage.done" ]]; then event "$stage" SKIP_DONE; return; fi
  local attempt
  attempt=$(( $(awk -F '\t' -v s="$stage" '$1==s{n++}END{print n+0}' "$ROOT/logs/attempts.tsv" 2>/dev/null || true) + 1 ))
  printf '%s\t%s\t%s\t%s\n' "$stage" "$attempt" "$(now)" "${JOB_ID:-interactive}" >> "$ROOT/logs/attempts.tsv"
  event "$stage" "START attempt=$attempt"
  local start end rc
  start=$(date +%s)
  set +e
  local timebin=""
  if [[ -x /usr/bin/time ]]; then
    timebin=/usr/bin/time
  elif [[ -x /bin/time ]]; then
    timebin=/bin/time
  fi
  if [[ -x "$timebin" ]]; then
    "$timebin" -v -o "$ROOT/logs/$stage.time.txt" bash "$SCRIPT_DIR/pipeline.sh" --run-stage "$stage" >"$ROOT/logs/$stage.stdout.log" 2>"$ROOT/logs/$stage.stderr.log"
  else
    bash "$SCRIPT_DIR/pipeline.sh" --run-stage "$stage" >"$ROOT/logs/$stage.stdout.log" 2>"$ROOT/logs/$stage.stderr.log"
    printf 'External GNU time unavailable on %s; wall time retained in progress.tsv.\n' "$(hostname)" > "$ROOT/logs/$stage.time.txt"
  fi
  rc=$?
  set -e
  end=$(date +%s)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$stage" "$(date -d "@$start" '+%F %T')" "$(date -d "@$end" '+%F %T')" "$((end-start))" "$rc" "${JOB_ID:-interactive}" >> "$PROGRESS"
  if (( rc != 0 )); then event "$stage" "FAILED rc=$rc"; return "$rc"; fi
  printf '%s\n' "$(now)" > "$ROOT/logs/$stage.done"
  event "$stage" "END wall_s=$((end-start))"
}

stage00(){
  [[ -s "$BAM_PRIMARY" && -s "$REF" && -s "$TDNA" ]] || { echo "Missing input" >&2; return 2; }
  "$SAMTOOLS" quickcheck -v "$BAM_PRIMARY"
  if [[ -n "$BAM_SECONDARY" ]]; then
    [[ -s "$BAM_SECONDARY" ]] || { echo "Missing secondary input" >&2; return 2; }
    "$SAMTOOLS" quickcheck -v "$BAM_SECONDARY"
    if [[ ! -s "$BAM" ]]; then
      "$SAMTOOLS" merge -@ "$THREADS" -f "$BAM.tmp" "$BAM_PRIMARY" "$BAM_SECONDARY"
      mv "$BAM.tmp" "$BAM"
      "$SAMTOOLS" index -@ "$THREADS" "$BAM"
    fi
  elif [[ "$NORMALIZE_BAM_HEADER" == 1 ]]; then
    # Normalize a workflow-owned copy when an input BAM header has trailing
    # whitespace; never alter or overwrite the source BAM.
    if [[ ! -s "$BAM" ]]; then
      "$SAMTOOLS" view --no-PG -H "$BAM_PRIMARY" | sed 's/[[:space:]]*$//' > "$ROOT/00.manifest/$SAMPLE.sanitized.header.sam"
      "$SAMTOOLS" reheader -P "$ROOT/00.manifest/$SAMPLE.sanitized.header.sam" "$BAM_PRIMARY" > "$BAM.tmp"
      mv "$BAM.tmp" "$BAM"
      "$SAMTOOLS" index -@ "$THREADS" "$BAM"
    fi
  fi
  "$SAMTOOLS" --version | head -2 > "$ROOT/00.manifest/samtools.version.txt"
  "$SAMTOOLS" quickcheck -v "$BAM"
  local tdna_bp
  tdna_bp=$(awk '!/^>/{gsub(/[[:space:]]/,"");n+=length}END{print n+0}' "$TDNA")
  if [[ -n "$EXPECTED_TDNA_BP" && "$tdna_bp" != "$EXPECTED_TDNA_BP" ]]; then
    echo "WARNING expected $EXPECTED_TDNA_BP bp, found $tdna_bp" >&2
  fi
  {
    printf 'field\tvalue\n'
    printf 'sample\t%s\ninput_bam_primary\t%s\ninput_bam_secondary\t%s\neffective_bam\t%s\nreference\t%s\ntdna_fasta\t%s\ntdna_length_bp\t%s\n' "$SAMPLE" "$BAM_PRIMARY" "${BAM_SECONDARY:-NONE}" "$BAM" "$REF" "$TDNA" "$tdna_bp"
    printf 'input_note\t%s\n' "${INPUT_NOTE:-rmdup BAM input}"
    printf 'expected_tdna_length_bp\t%s\n' "${EXPECTED_TDNA_BP:-NOT_SET}"
    printf 'tdna_scope\t%s\n' "$TDNA_SCOPE_NOTE"
    printf 'candidate_definition\tALL QNAMEs mapped to supplied T-DNA UNION ALL QNAMEs unmapped to plant reference, with every recoverable mate restored\n'
    printf 'homology_policy\trisk label only; never hard-delete candidates\n'
  } > "$ROOT/00.manifest/run_manifest.tsv"
  [[ -s "$PROGRESS" ]] || printf 'stage\tstart\tend\twall_seconds\texit_code\tjob_id\n' > "$PROGRESS"

  local bindir="$ROOT/envs/bin"; mkdir -p "$bindir"
  command -v python3 > "$ROOT/00.manifest/python3.path.txt"
  for tool in bwa minimap2 SOAPdenovo-63mer SOAPdenovo-127mer SOAPdenovo-Trans-63mer SOAPdenovo-Trans-127mer nucmer show-coords; do
    command -v "$tool" > "$ROOT/00.manifest/$tool.path.txt" 2>/dev/null || true
  done
  if ! command -v bwa >/dev/null && [[ -n "$FALLBACK_TOOL_BIN" && -x "$FALLBACK_TOOL_BIN/bwa" ]]; then ln -sf "$FALLBACK_TOOL_BIN/bwa" "$bindir/bwa"; fi
  if ! command -v minimap2 >/dev/null && [[ -n "$FALLBACK_TOOL_BIN" && -x "$FALLBACK_TOOL_BIN/minimap2" ]]; then ln -sf "$FALLBACK_TOOL_BIN/minimap2" "$bindir/minimap2"; fi
  export PATH="$bindir:$PATH"
  if ! command -v SOAPdenovo-63mer >/dev/null && ! command -v SOAPdenovo-Trans-63mer >/dev/null; then
    local conda
    conda=$(command -v "$CONDA_BIN" 2>/dev/null || true)
    [[ -n "$conda" && -x "$conda" ]] || { echo "No SOAPdenovo and conda unavailable" >&2; return 3; }
    "$conda" create -y --solver classic -p "$TOOL_ENV" --override-channels \
      -c https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda \
      -c https://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge \
      -c https://mirrors.ustc.edu.cn/anaconda/pkgs/main soapdenovo2 minimap2 bwa mummer4 || \
    "$conda" create -y --solver classic -p "$TOOL_ENV" --override-channels -c bioconda -c conda-forge soapdenovo2 minimap2 bwa mummer4 || \
    "$conda" create -y --solver classic -p "$TOOL_ENV" --override-channels -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge soapdenovo2 minimap2 bwa
  fi
  [[ -x "$TOOL_ENV/bin/SOAPdenovo-63mer" ]] && ln -sf "$TOOL_ENV/bin/SOAPdenovo-63mer" "$bindir/SOAPdenovo-63mer"
  [[ -x "$TOOL_ENV/bin/SOAPdenovo-Trans-63mer" ]] && ln -sf "$TOOL_ENV/bin/SOAPdenovo-Trans-63mer" "$bindir/SOAPdenovo-63mer"
  [[ -x "$TOOL_ENV/bin/minimap2" ]] && ln -sf "$TOOL_ENV/bin/minimap2" "$bindir/minimap2"
  [[ -x "$TOOL_ENV/bin/bwa" ]] && ln -sf "$TOOL_ENV/bin/bwa" "$bindir/bwa"
  for x in nucmer show-coords; do [[ -x "$TOOL_ENV/bin/$x" ]] && ln -sf "$TOOL_ENV/bin/$x" "$bindir/$x"; done
  command -v SOAPdenovo-63mer; command -v minimap2; command -v bwa
}

stage01(){
  export PATH="$ROOT/envs/bin:$TOOL_ENV/bin:$PATH"
  mkdir -p "$ROOT/01.qnames/vector_map_tmp"
  if [[ ! -s "$ROOT/01.qnames/tdna_idx.bwt" ]]; then bwa index -p "$ROOT/01.qnames/tdna_idx" "$TDNA"; fi
  IDX="$ROOT/01.qnames/tdna_idx"
  # Only QNAME membership is needed here. Streaming each primary record avoids
  # a costly whole-BAM name sort; strand is restored later from the original BAM.
  "$SAMTOOLS" view -@ "$THREADS" -F 2304 "$BAM" | \
    awk 'BEGIN{OFS="\n"}$10!="*"&&$11!="*"{q=$1; sub(/[[:space:]].*/,"",q); print "@"q,$10,"+",$11}' | \
    bwa mem -t "$THREADS" "$IDX" - | \
    "$SAMTOOLS" view -@ "$THREADS" -F 4 - | cut -f1 | \
    LC_ALL=C sort -u -S 24G -T "$ROOT/01.qnames/vector_map_tmp" > "$ROOT/01.qnames/tdna_mapped.qnames.txt.tmp"
  mv "$ROOT/01.qnames/tdna_mapped.qnames.txt.tmp" "$ROOT/01.qnames/tdna_mapped.qnames.txt"
  "$SAMTOOLS" view -@ "$THREADS" -f 4 -F 2304 "$BAM" | cut -f1 | \
    LC_ALL=C sort -u -S 8G -T "$ROOT/01.qnames/vector_map_tmp" > "$ROOT/01.qnames/genome_unmapped.qnames.txt.tmp"
  mv "$ROOT/01.qnames/genome_unmapped.qnames.txt.tmp" "$ROOT/01.qnames/genome_unmapped.qnames.txt"
  LC_ALL=C sort -u -S 24G -T "$ROOT/01.qnames/vector_map_tmp" \
    "$ROOT/01.qnames/tdna_mapped.qnames.txt" "$ROOT/01.qnames/genome_unmapped.qnames.txt" > "$ROOT/01.qnames/candidate.qnames.txt.tmp"
  mv "$ROOT/01.qnames/candidate.qnames.txt.tmp" "$ROOT/01.qnames/candidate.qnames.txt"
  {
    printf 'set\tunique_qnames\n'
    printf 'tdna_mapped\t%s\n' "$(wc -l < "$ROOT/01.qnames/tdna_mapped.qnames.txt")"
    printf 'plant_unmapped\t%s\n' "$(wc -l < "$ROOT/01.qnames/genome_unmapped.qnames.txt")"
    printf 'union\t%s\n' "$(wc -l < "$ROOT/01.qnames/candidate.qnames.txt")"
  } > "$ROOT/01.qnames/read_counts.tsv"
}

stage02(){
  local namebam="$ROOT/02.fastq/$SAMPLE.candidate.namesort.bam"
  "$SAMTOOLS" view -@ "$THREADS" -b -F 2304 -N "$ROOT/01.qnames/candidate.qnames.txt" "$BAM" | \
    "$SAMTOOLS" sort -@ "$THREADS" -n -m 2G -T "$ROOT/tmp/${SAMPLE}_namesort" -o "$namebam.tmp" -
  mv "$namebam.tmp" "$namebam"
  "$SAMTOOLS" fastq -@ "$THREADS" -n \
    -1 "$ROOT/02.fastq/$SAMPLE.candidates.R1.fastq.gz" \
    -2 "$ROOT/02.fastq/$SAMPLE.candidates.R2.fastq.gz" \
    -s "$ROOT/02.fastq/$SAMPLE.candidates.singleton.fastq.gz" \
    -0 "$ROOT/02.fastq/$SAMPLE.candidates.unclassified.fastq.gz" "$namebam"
  {
    printf 'file\treads\n'
    for f in "$ROOT/02.fastq"/*.fastq.gz; do printf '%s\t%s\n' "$(basename "$f")" "$(gzip -cd "$f" | awk 'END{print NR/4}')"; done
  } > "$ROOT/02.fastq/fastq_counts.tsv"
}

stage03(){
  export PATH="$ROOT/envs/bin:$TOOL_ENV/bin:$PATH"
  gzip -cd "$ROOT/02.fastq/$SAMPLE.candidates.singleton.fastq.gz" "$ROOT/02.fastq/$SAMPLE.candidates.unclassified.fastq.gz" | gzip -1 > "$ROOT/02.fastq/$SAMPLE.candidates.all_singletons.fastq.gz"
  local maxlen avgins
  # Consume each stream fully so pipefail does not treat an intentional early
  # downstream exit as a gzip/samtools SIGPIPE failure (exit 141).
  maxlen=$(gzip -cd "$ROOT/02.fastq/$SAMPLE.candidates.R1.fastq.gz" | awk 'NR%4==2&&length>m{m=length}END{print m+0}')
  local existing_stats="$PROJECT/$SAMPLE/01.qc/$SAMPLE.rmdup.stats.txt"
  if [[ -s "$existing_stats" ]]; then
    avgins=$(awk -F '\t' '$1=="SN" && $2=="insert size average:"{v=$3}END{if(v!="")printf "%.0f",v}' "$existing_stats")
  else
    avgins=$($SAMTOOLS view "$BAM" Chr1:1-10000000 | awk '$9!=0{x=$9<0?-$9:$9;if(x>=100&&x<=2000){s+=x;n++}}END{if(n)printf "%.0f",s/n}')
  fi
  [[ -n "$avgins" && "$avgins" -ge 100 ]] || avgins=350
  cat > "$ROOT/03.assembly/soap.config" <<EOF
max_rd_len=$maxlen
[LIB]
avg_ins=$avgins
reverse_seq=0
asm_flags=3
rank=1
pair_num_cutoff=3
map_len=32
q1=$ROOT/02.fastq/$SAMPLE.candidates.R1.fastq.gz
q2=$ROOT/02.fastq/$SAMPLE.candidates.R2.fastq.gz
q=$ROOT/02.fastq/$SAMPLE.candidates.all_singletons.fastq.gz
EOF
  local soap
  soap=$(command -v SOAPdenovo-63mer || command -v SOAPdenovo-Trans-63mer)
  "$soap" all -s "$ROOT/03.assembly/soap.config" -K 41 -p "$THREADS" -R -o "$ROOT/03.assembly/$SAMPLE.K41"
  [[ -s "$ROOT/03.assembly/$SAMPLE.K41.contig" ]] || return 4
  python3 "$SCRIPT_DIR/fasta_metrics.py" "$ROOT/03.assembly/$SAMPLE.K41.contig" "$ROOT/03.assembly/assembly_metrics.tsv"
}

stage04(){
  export PATH="$ROOT/envs/bin:$TOOL_ENV/bin:$PATH"
  local contigs="$ROOT/03.assembly/$SAMPLE.K41.contig"
  minimap2 -t "$THREADS" -x asm5 --secondary=yes "$REF" "$contigs" > "$ROOT/04.alignments/contigs_to_reference.paf"
  minimap2 -t "$THREADS" -x asm5 --secondary=yes "$TDNA" "$contigs" > "$ROOT/04.alignments/contigs_to_tdna.paf"
  minimap2 -t "$THREADS" -x asm5 --secondary=yes "$REF" "$TDNA" > "$ROOT/04.alignments/tdna_to_reference_homology.paf"
  if command -v nucmer >/dev/null && command -v show-coords >/dev/null; then
    if nucmer --maxmatch -t "$THREADS" -p "$ROOT/04.alignments/mummer.contigs_to_reference" "$REF" "$contigs" && \
       show-coords -THrd "$ROOT/04.alignments/mummer.contigs_to_reference.delta" > "$ROOT/04.alignments/mummer.contigs_to_reference.coords.tsv" && \
       nucmer --maxmatch -t "$THREADS" -p "$ROOT/04.alignments/mummer.contigs_to_tdna" "$TDNA" "$contigs" && \
       show-coords -THrd "$ROOT/04.alignments/mummer.contigs_to_tdna.delta" > "$ROOT/04.alignments/mummer.contigs_to_tdna.coords.tsv"; then
      printf 'MUMmer4 and minimap2 both completed; PAF is the normalized candidate-calling interface.\n' > "$ROOT/04.alignments/alignment_backend.txt"
    else
      printf 'MUMmer4 runtime failed; minimap2 PAF retained and used; MUMmer outputs/logs preserved for retry.\n' > "$ROOT/04.alignments/alignment_backend.txt"
    fi
  else
    printf 'minimap2 fallback used; rerun stage04 after installing MUMmer4 for orthogonal validation.\n' > "$ROOT/04.alignments/alignment_backend.txt"
  fi
}

stage05(){
  python3 "$SCRIPT_DIR/call_junctions.py" \
    --genome-paf "$ROOT/04.alignments/contigs_to_reference.paf" \
    --tdna-paf "$ROOT/04.alignments/contigs_to_tdna.paf" \
    --homology-paf "$ROOT/04.alignments/tdna_to_reference_homology.paf" \
    --sample "$SAMPLE" --output "$ROOT/05.junctions/$SAMPLE.junction_candidates.tsv"
  awk -F '\t' 'NR==1 || $1!=""' "$ROOT/05.junctions/$SAMPLE.junction_candidates.tsv" > "$ROOT/05.junctions/$SAMPLE.junction_candidates.checked.tsv"
}

stage06(){
  export PATH="$ROOT/envs/bin:$TOOL_ENV/bin:$PATH"
  local contigs="$ROOT/03.assembly/$SAMPLE.K41.contig"
  if [[ $(wc -l < "$ROOT/05.junctions/$SAMPLE.junction_candidates.tsv") -le 1 ]]; then
    python3 "$SCRIPT_DIR/validate_junctions.py" --samtools "$SAMTOOLS" --bam /dev/null \
      --candidates "$ROOT/05.junctions/$SAMPLE.junction_candidates.tsv" \
      --output "$ROOT/06.validation/$SAMPLE.junction_candidates.validated.tsv"
    printf 'No junction candidates; read-back mapping skipped.\n' > "$ROOT/06.validation/no_candidates.txt"
    return
  fi
  minimap2 -t "$THREADS" -ax sr "$contigs" "$ROOT/02.fastq/$SAMPLE.candidates.R1.fastq.gz" "$ROOT/02.fastq/$SAMPLE.candidates.R2.fastq.gz" | \
    "$SAMTOOLS" sort -@ "$THREADS" -m 2G -T "$ROOT/tmp/validate_pairs" -o "$ROOT/06.validation/pairs_to_contigs.bam" -
  minimap2 -t "$THREADS" -ax sr "$contigs" "$ROOT/02.fastq/$SAMPLE.candidates.all_singletons.fastq.gz" | \
    "$SAMTOOLS" sort -@ "$THREADS" -m 2G -T "$ROOT/tmp/validate_singletons" -o "$ROOT/06.validation/singletons_to_contigs.bam" -
  "$SAMTOOLS" merge -@ "$THREADS" -f -u - "$ROOT/06.validation/pairs_to_contigs.bam" "$ROOT/06.validation/singletons_to_contigs.bam" | \
    "$SAMTOOLS" sort -@ "$THREADS" -m 2G -T "$ROOT/tmp/validate_merge" -o "$ROOT/06.validation/candidates_to_contigs.bam.tmp" -
  mv "$ROOT/06.validation/candidates_to_contigs.bam.tmp" "$ROOT/06.validation/candidates_to_contigs.bam"
  "$SAMTOOLS" index -@ "$THREADS" "$ROOT/06.validation/candidates_to_contigs.bam"
  python3 "$SCRIPT_DIR/validate_junctions.py" --samtools "$SAMTOOLS" \
    --bam "$ROOT/06.validation/candidates_to_contigs.bam" \
    --candidates "$ROOT/05.junctions/$SAMPLE.junction_candidates.tsv" \
    --output "$ROOT/06.validation/$SAMPLE.junction_candidates.validated.tsv"
}

stage07(){
  cp "$ROOT/06.validation/$SAMPLE.junction_candidates.validated.tsv" "$ROOT/07.report/$SAMPLE.final_breakpoint_candidates.tsv"
  {
    echo "# $SAMPLE SOAPdenovo2 T-DNA analysis report"
    echo
    echo "Generated: $(now); SGE job: ${JOB_ID:-unknown}; host: $(hostname)"
    echo
    echo 'Candidate pool: every QNAME mapped to the supplied T-DNA plus every QNAME unmapped to the plant reference, deduplicated by QNAME, with recoverable mates restored.'
    echo
    echo "**Construct-reference scope:** $TDNA_SCOPE_NOTE"
    echo
    echo 'T-DNA-genome homology is retained as `TDNA_GENOME_HOMOLOGY`; it is never used as a hard deletion filter.'
    echo
    echo '## Counts'; sed 's/^/    /' "$ROOT/01.qnames/read_counts.tsv"; sed 's/^/    /' "$ROOT/02.fastq/fastq_counts.tsv"
    echo '## Assembly'; sed 's/^/    /' "$ROOT/03.assembly/assembly_metrics.tsv"
    echo '## Breakpoints'; awk 'END{print "validated candidate rows: " NR-1}' "$ROOT/07.report/$SAMPLE.final_breakpoint_candidates.tsv"
    echo
    echo 'Risk-labelled rows are retained for audit, but `RISK_FLAGGED` rows are not confirmed insertion breakpoints.'
    sed 's/^/    /' "$ROOT/07.report/$SAMPLE.final_breakpoint_candidates.tsv"
    if [[ -s "$ROOT/07.report/SGE_job_summary.tsv" ]]; then
      echo '## SGE jobs'; sed 's/^/    /' "$ROOT/07.report/SGE_job_summary.tsv"
    fi
    if [[ -s "$ROOT/07.report/failures_retries.tsv" ]]; then
      echo '## Failures and retries'; sed 's/^/    /' "$ROOT/07.report/failures_retries.tsv"
    fi
    echo '## Resources'; for f in "$ROOT"/logs/stage*.time.txt; do echo "### $(basename "$f")"; (grep -E 'Elapsed|User time|System time|Percent of CPU|Maximum resident|Exit status' "$f" || true) | sed 's/^/    /'; done
  } > "$ROOT/07.report/REPORT.md"
}

if [[ ${1:-} == --run-stage ]]; then
  [[ $# == 2 && $2 =~ ^stage0[0-7]$ ]] || { echo "Invalid stage selector" >&2; exit 64; }
  "$2"
  exit
fi

run_stage stage00 stage00
run_stage stage01 stage01
run_stage stage02 stage02
run_stage stage03 stage03
run_stage stage04 stage04
run_stage stage05 stage05
run_stage stage06 stage06
run_stage stage07 stage07
event pipeline COMPLETE
