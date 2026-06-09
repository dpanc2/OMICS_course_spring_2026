# День 3. ChIP-seq: выравнивание, треки покрытий и MACS3

Нам нужно получить файлы, которые будем открывать в IGV:

- `q30.bam` — выравнивания после фильтрации по mapq30
- `q30.bw` — трек покрытия
- `MACS3 FE.bw` — fold enrichment относительно input.

На занятии руками разберем все шаги для чипа `H3K9me3`  
Для `H3K27Ac`, `H2A119Ub` и `CTCF` запустим готовый скрипт, который прогоняет такой же пайплайн.

## 1. Переходим в папку практики

```bash
cd ~/OMICS_course_spring_2026/day3_ChIPseq_practice
```

## 2. Активируем окружение

Используем то же окружение, что и в предыдущих практиках:

```bash
conda activate hic_practice
```

Если нужные пакеты еще не установлены:

```bash
mamba install -y -c conda-forge -c bioconda umi_tools setuptools
```

## 3. Скачиваем raw reads (идем на 10 шаг и запускаем выравнивание для H3K9me3)

Для разбора скачиваем `H3K9me3` и input

```bash
python3 scripts/download_raw_reads.py \
  --sample MoPh7_H3K9me3 \
  --sample MoPh7_input
```

После этого появятся файлы:

```text
data/raw/MoPh7_H3K9me3/MoPh7_H3K9me3_R1.fastq.gz
data/raw/MoPh7_H3K9me3/MoPh7_H3K9me3_R2.fastq.gz
data/raw/MoPh7_histone_input/MoPh7_input_R1.fastq.gz
data/raw/MoPh7_histone_input/MoPh7_input_R2.fastq.gz
```

Если нужно скачать для всех меток, используем скрипт без флагов:

```bash
python3 scripts/download_raw_reads.py
```

## 4. Быстрая проверка качества

Можно запустить, но необязателено :)

```bash
mkdir -p results/fastqc

fastqc \
  data/raw/MoPh7_H3K9me3/MoPh7_H3K9me3_R1.fastq.gz \
  data/raw/MoPh7_H3K9me3/MoPh7_H3K9me3_R2.fastq.gz \
  data/raw/MoPh7_histone_input/MoPh7_histone_input_R1.fastq.gz \
  data/raw/MoPh7_histone_input/MoPh7_histone_input_R2.fastq.gz \
  -o results/fastqc

multiqc results/fastqc -o results/fastqc
```

## 5. Проверяем BWA index

Используем T2T референс, который уже готовили в первый день:

```bash
GENOME_FA="../day1_HiC_practice/data/reference/T2T_human.fna"
CHROM_SIZES="../day1_HiC_practice/data/reference/chrom.sizes"

ls -lh "$GENOME_FA"
ls -lh "$CHROM_SIZES"
```

## 6. Обрезаем адаптеры

Сначала обработаем `H3K9me3`.

```bash
SAMPLE="MoPh7_H3K9me3"
THREADS=8

mkdir -p "results/alignments/${SAMPLE}/FilteredReads" \
         "results/alignments/${SAMPLE}/logs"

R1="data/raw/${SAMPLE}/${SAMPLE}_R1.fastq.gz"
R2="data/raw/${SAMPLE}/${SAMPLE}_R2.fastq.gz"
OUT="results/alignments/${SAMPLE}"
```

Обрезаем стандартные адаптеры Illumina:

```bash
cutadapt \
  -j "$THREADS" \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -m 70 \
  -o "${OUT}/FilteredReads/R1.cln.fq.gz" \
  -p "${OUT}/FilteredReads/R2.cln.fq.gz" \
  "$R1" "$R2" \
  --report=full > "${OUT}/logs/${SAMPLE}.cutadapt_trim3.txt"
```

## 7. Выравниваем риды на геном

```bash
BWA_INDEX="../day1_HiC_practice/data/reference/T2T_human.fna"

bwa mem \
  -t "$THREADS" \
  "$BWA_INDEX" \
  "${OUT}/FilteredReads/R1.cln.fq.gz" \
  "${OUT}/FilteredReads/R2.cln.fq.gz" \
  2> "${OUT}/logs/${SAMPLE}.bwa.log" \
  | samtools view -b - \
  | samtools sort -@ "$THREADS" -o "${OUT}/${SAMPLE}.sorted.bam"

samtools index -@ "$THREADS" "${OUT}/${SAMPLE}.sorted.bam"
```

## 8. Оставляем выравнивания с mapq >= 30

Оставляем только хорошо выровненные proper pairs:

```bash
samtools view \
  -@ "$THREADS" \
  -b \
  -q 30 \
  -f 0x2 \
  "${OUT}/${SAMPLE}.sorted.bam" \
  > "${OUT}/${SAMPLE}.q30.bam"

samtools index -@ "$THREADS" "${OUT}/${SAMPLE}.q30.bam"
```

Проверим число reads:

```bash
samtools flagstat "${OUT}/${SAMPLE}.q30.bam"
```

## 9. Делаем coverage bigWig

```bash
mkdir -p results/tracks

bedtools genomecov \
  -ibam "${OUT}/${SAMPLE}.q30.bam" \
  -bg \
  > "results/tracks/${SAMPLE}.q30.bedGraph"

sort -k1,1 -k2,2n \
  "results/tracks/${SAMPLE}.q30.bedGraph" \
  > "results/tracks/${SAMPLE}.q30.sorted.bedGraph"

bedGraphToBigWig \
  "results/tracks/${SAMPLE}.q30.sorted.bedGraph" \
  "$CHROM_SIZES" \
  "results/tracks/${SAMPLE}.q30.bw"
```

Файлы, которые можно открыть в IGV:

```text
results/alignments/MoPh7_H3K9me3/MoPh7_H3K9me3.q30.bam
results/tracks/MoPh7_H3K9me3.q30.bw
```

## 10. Повторяем alignment для input

MACS3 сравнивает ChIP с input, поэтому input должен пройти те же шаги:

```bash
bash scripts/run_all_chipseq_pipeline.sh MoPh7_H3K9me3
```

Эта команда автоматически обработает:

```text
MoPh7_H3K9me3
MoPh7_input
```

Если первые шаги для `H3K9me3` уже были сделаны руками, скрипт увидит готовый `q30.bam` и не будет заново запускать alignment для него.

## 11. MACS3 peak calling для H3K9me3

В скрипте MACS3 запустится автоматически, но полезно отдельно понимать, что происходит.

Для `H3K9me3` используем broad peak calling:

```bash
CHIP="MoPh7_H3K9me3"
INPUT="MoPh7_input"

CHIP_BAM="results/alignments/${CHIP}/${CHIP}.q30.bam"
INPUT_BAM="results/alignments/${INPUT}/${INPUT}.q30.bam"
MACS_OUT="results/macs/${CHIP}"

mkdir -p "$MACS_OUT"

macs3 callpeak \
  -t "$CHIP_BAM" \
  -c "$INPUT_BAM" \
  -g hs \
  -n "$CHIP" \
  --outdir "$MACS_OUT" \
  --broad \
  --broad-cutoff 0.1 \
  -B --SPMR
```

MACS3 создаст несколько файлов.  
Для IGV нам нужны:

```text
results/macs/MoPh7_H3K9me3/MoPh7_H3K9me3_treat_pileup.bdg
results/macs/MoPh7_H3K9me3/MoPh7_H3K9me3_control_lambda.bdg
results/macs/MoPh7_H3K9me3/MoPh7_H3K9me3_peaks.broadPeak
```

Теперь считаем fold enrichment относительно input:

```bash
macs3 bdgcmp \
  -t "${MACS_OUT}/${CHIP}_treat_pileup.bdg" \
  -c "${MACS_OUT}/${CHIP}_control_lambda.bdg" \
  --outdir "$MACS_OUT" \
  --o-prefix "${CHIP}_FE" \
  -m FE

sort -k1,1 -k2,2n \
  "${MACS_OUT}/${CHIP}_FE_FE.bdg" \
  > "${MACS_OUT}/${CHIP}_FE.sorted.bdg"

bedGraphToBigWig \
  "${MACS_OUT}/${CHIP}_FE.sorted.bdg" \
  "$CHROM_SIZES" \
  "${MACS_OUT}/${CHIP}_FE.bw"
```

Файлы для IGV:

```text
results/macs/MoPh7_H3K9me3/MoPh7_H3K9me3_FE.bw
results/macs/MoPh7_H3K9me3/MoPh7_H3K9me3_peaks.broadPeak
```

## 12. Запускаем остальные ChIP-seq samples

Скачиваем все данные:

```bash
python3 scripts/download_raw_reads.py
```

Запускаем полный pipeline:

```bash
bash scripts/run_all_chipseq_pipeline.sh
```

По умолчанию скрипт обработает:

```text
MoPh7_H3K9me3
MoPh7_H3K27Ac
MoPh7_H2A119Ub
MoPh7_CTCF
```

Для гистоновых меток будет использован общий input:

```text
MoPh7_input
```

Можно запустить только один sample:

```bash
bash scripts/run_all_chipseq_pipeline.sh MoPh7_H3K27Ac
```

## 13. Что должно получиться

Для каждого ChIP sample:

```text
results/alignments/<sample>/<sample>.q30.bam
results/alignments/<sample>/<sample>.q30.bam.bai
results/tracks/<sample>.q30.bw
results/macs/<sample>/<sample>.bw
results/macs/<sample>/<sample>_FE.bw
```

Peak files:

```text
results/macs/MoPh7_H3K9me3/MoPh7_H3K9me3_peaks.broadPeak
results/macs/MoPh7_H2A119Ub/MoPh7_H2A119Ub_peaks.broadPeak
results/macs/MoPh7_H3K27Ac/MoPh7_H3K27Ac_peaks.narrowPeak
results/macs/MoPh7_CTCF/MoPh7_CTCF_peaks.narrowPeak
```

## 14. Что смотреть в IGV

Откройте в IGV:

- `q30.bam` для проверки самих выравниваний
- `FE.bw` как fold enrichment относительно input
- `broadPeak` или `narrowPeak` как координаты пиков

Для `H3K9me3` ожидаем широкие домены  
Для `H3K27Ac` и `CTCF` пики обычно выглядят узкими

## 15. Задание

1. Откройте в IGV `H3K9me3.q30.bw`, `H3K9me3.bw` и `H3K9me3_FE.bw`
2. Найдите регион с сильным `H3K9me3` сигналом
3. Сравните этот регион с input coverage
4. Сравните форму сигнала `H3K9me3`, `H3K27Ac`, `H2A119Ub` и `CTCF`

