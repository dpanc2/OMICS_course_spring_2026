#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

THREADS="${THREADS:-8}"
GENOME_FA="${GENOME_FA:-../day1_HiC_practice/data/reference/T2T_human.fna}"
CHROM_SIZES="${CHROM_SIZES:-../day1_HiC_practice/data/reference/chrom.sizes}"
GENOME_SIZE="${GENOME_SIZE:-hs}"

DEFAULT_SAMPLES=(
  MoPh7_H3K9me3
  MoPh7_H3K27Ac
  MoPh7_H2A119Ub
  MoPh7_CTCF
)

SAMPLES=("$@")
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  SAMPLES=("${DEFAULT_SAMPLES[@]}")
fi

mkdir -p results/alignments results/tracks results/macs results/logs

need_file() {
  if [[ ! -s "$1" ]]; then
    echo "Missing file: $1" >&2
    exit 1
  fi
}

ensure_bwa_index() {
  local suffix
  for suffix in amb ann bwt pac sa; do
    need_file "${GENOME_FA}.${suffix}"
  done
}

input_for_sample() {
  case "$1" in
    MoPh7_H3K9me3|MoPh7_H3K27Ac|MoPh7_H2A119Ub)
      printf '%s\n' "MoPh7_input"
      ;;
    MoPh7_CTCF)
      printf '%s\n' "MoPh7_CTCF_input"
      ;;
    *)
      echo "Unknown ChIP sample: $1" >&2
      exit 1
      ;;
  esac
}

macs_mode_for_sample() {
  case "$1" in
    MoPh7_H3K9me3|MoPh7_H2A119Ub)
      printf '%s\n' "broad"
      ;;
    *)
      printf '%s\n' "narrow"
      ;;
  esac
}

process_alignment() {
  local sample="$1"
  local raw_dir="data/raw/${sample}"
  local out_dir="results/alignments/${sample}"

  local r1="${raw_dir}/${sample}_R1.fastq.gz"
  local r2="${raw_dir}/${sample}_R2.fastq.gz"
  need_file "$r1"
  need_file "$r2"

  mkdir -p "$out_dir/FilteredReads" "$out_dir/logs"

  if [[ -s "${out_dir}/${sample}.q30.bam" && -s "${out_dir}/${sample}.q30.bam.bai" ]]; then
    echo "== ${sample}: q30 BAM exists, skip alignment"
  else
    echo "== ${sample}: trimming adapters"
    cutadapt \
      -j "$THREADS" \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      -m 70 \
      -o "${out_dir}/FilteredReads/R1.trim.fq.gz" \
      -p "${out_dir}/FilteredReads/R2.trim.fq.gz" \
      "$r1" "$r2" \
      --report=full > "${out_dir}/logs/${sample}.cutadapt_trim3.txt"

    cutadapt \
      --cores "$THREADS" \
      -g ^TGACTTatcatgtctgctcgaagcGTCGAGCTAAGTCGTATGGTATGCGCAATTACATTCAGCT \
      -g ^ATGGTATGCGCAATTACATTCAGCT \
      -a AGCTGAATGTAATTGCGCATACCAT \
      -A AGCTGAATGTAATTGCGCATACCAT \
      -e 0.2 \
      --overlap 10 \
      -o "${out_dir}/FilteredReads/R1.cln.fq.gz" \
      -p "${out_dir}/FilteredReads/R2.cln.fq.gz" \
      "${out_dir}/FilteredReads/R1.trim.fq.gz" \
      "${out_dir}/FilteredReads/R2.trim.fq.gz" \
      --report=full > "${out_dir}/logs/${sample}.cutadapt_filter5.txt"

    echo "== ${sample}: BWA alignment"
    bwa mem \
      -t "$THREADS" \
      "$GENOME_FA" \
      "${out_dir}/FilteredReads/R1.cln.fq.gz" \
      "${out_dir}/FilteredReads/R2.cln.fq.gz" \
      2> "${out_dir}/logs/${sample}.bwa.log" \
      | samtools view -b - \
      | samtools sort -@ "$THREADS" -o "${out_dir}/${sample}.sorted.bam"

    samtools index -@ "$THREADS" "${out_dir}/${sample}.sorted.bam"

    echo "== ${sample}: MAPQ >= 30"
    samtools view \
      -@ "$THREADS" \
      -b \
      -q 30 \
      -f 0x2 \
      "${out_dir}/${sample}.sorted.bam" \
      > "${out_dir}/${sample}.q30.bam"

    samtools index -@ "$THREADS" "${out_dir}/${sample}.q30.bam"
  fi

  echo "== ${sample}: bigWig"
  bedtools genomecov \
    -ibam "${out_dir}/${sample}.q30.bam" \
    -bg \
    > "results/tracks/${sample}.q30.bedGraph"

  sort -k1,1 -k2,2n \
    "results/tracks/${sample}.q30.bedGraph" \
    > "results/tracks/${sample}.q30.sorted.bedGraph"

  bedGraphToBigWig \
    "results/tracks/${sample}.q30.sorted.bedGraph" \
    "$CHROM_SIZES" \
    "results/tracks/${sample}.q30.bw"
}

call_macs3() {
  local chip="$1"
  local input="$2"
  local mode
  mode="$(macs_mode_for_sample "$chip")"

  local chip_bam="results/alignments/${chip}/${chip}.q30.bam"
  local input_bam="results/alignments/${input}/${input}.q30.bam"
  need_file "$chip_bam"
  need_file "$input_bam"

  local out_dir="results/macs/${chip}"
  mkdir -p "$out_dir"

  echo "== ${chip}: MACS3 callpeak (${mode})"
  if [[ "$mode" == "broad" ]]; then
    macs3 callpeak \
      -t "$chip_bam" \
      -c "$input_bam" \
      -f BAMPE \
      -g "$GENOME_SIZE" \
      -n "$chip" \
      --outdir "$out_dir" \
      --broad \
      --broad-cutoff 0.1 \
      -B \
      --SPMR
  else
    macs3 callpeak \
      -t "$chip_bam" \
      -c "$input_bam" \
      -f BAMPE \
      -g "$GENOME_SIZE" \
      -n "$chip" \
      --outdir "$out_dir" \
      -B \
      --SPMR
  fi

  local treat_bdg="${out_dir}/${chip}_treat_pileup.bdg"
  local control_bdg="${out_dir}/${chip}_control_lambda.bdg"
  need_file "$treat_bdg"
  need_file "$control_bdg"

  sort -k1,1 -k2,2n "$treat_bdg" > "${out_dir}/${chip}_treat_pileup.sorted.bdg"
  bedGraphToBigWig \
    "${out_dir}/${chip}_treat_pileup.sorted.bdg" \
    "$CHROM_SIZES" \
    "${out_dir}/${chip}.bw"

  echo "== ${chip}: MACS3 FE track"
  macs3 bdgcmp \
    -t "$treat_bdg" \
    -c "$control_bdg" \
    --outdir "$out_dir" \
    --o-prefix "${chip}_FE" \
    -m FE

  local fe_bdg="${out_dir}/${chip}_FE_FE.bdg"
  need_file "$fe_bdg"
  sort -k1,1 -k2,2n "$fe_bdg" > "${out_dir}/${chip}_FE.sorted.bdg"
  bedGraphToBigWig \
    "${out_dir}/${chip}_FE.sorted.bdg" \
    "$CHROM_SIZES" \
    "${out_dir}/${chip}_FE.bw"
}

need_file "$GENOME_FA"
need_file "$CHROM_SIZES"
ensure_bwa_index

for chip in "${SAMPLES[@]}"; do
  input="$(input_for_sample "$chip")"
  process_alignment "$chip"
  process_alignment "$input"
  call_macs3 "$chip" "$input"
done

echo "DONE"
