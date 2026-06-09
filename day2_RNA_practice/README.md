# День 2. RNA-seq и CAGE-seq: выравнивание и визуализация треков

В этой практике мы работаем с транскриптомными данными:

- **RNA-seq** — paired-end reads, будем выравнивать как обычный RNA-seq
- **CAGE-seq** — single-end reads

Основная цель занятия — получить:

- sorted bam
- mapq30 bam
- `bedGraph` coverage track
- `bigWig` coverage track

Для остальных образцов оставим готовый loop-скрипт, чтобы запустить его самостоятельно и сравнить треки между образцами

## 1. Переходим в папку практики

```bash
cd ~/OMICS_course_spring_2026/day2_RNA_practice
```

## 2. Активируем окружение

```bash
conda activate hic_practice

conda install -c conda-forge -c bioconda \
  star bwa samtools bedtools ucsc-bedgraphtobigwig requests tqdm
```

Здесь:

- `STAR` — splice-aware aligner для RNA-seq;
- `bwa` — обычный aligner для DNA-seq, используем только для сравнения;
- `samtools` — работа с BAM;
- `bedtools` — получение coverage;
- `bedGraphToBigWig` — конвертация `bedGraph` в `bigWig`.

## 3. Скачиваем raw reads

Для занятия достаточно скачать только `MoPh7`:

```bash
python3 scripts/download_raw_reads.py --sample MoPh7
```

После этого появятся файлы:

```text
data/raw/rnaseq/MoPh7_R1.fastq.gz
data/raw/rnaseq/MoPh7_R2.fastq.gz
data/raw/cage/MoPh7_R1.fastq.gz
```

Если нужно скачать все образцы:

```bash
python3 scripts/download_raw_reads.py
```

## 4. Быстрая проверка качества

Этот шаг можно запустить, но он необязателен

```bash
mkdir -p results/fastqc

fastqc \
  data/raw/rnaseq/MoPh7_R1.fastq.gz \
  data/raw/rnaseq/MoPh7_R2.fastq.gz \
  data/raw/cage/MoPh7_R1.fastq.gz \
  -o results/fastqc

multiqc results/fastqc -o results/fastqc
```

## 5. Готовим STAR index

STAR не выравнивает reads прямо на референс, сначала нужно один раз построить индекс

```bash
mkdir -p ../day1_HiC_practice/data/reference/star_index
GENOME_FA="../day1_HiC_practice/data/reference/T2T_human.fna"
```

Проверим, что есть .fna и `chrom.sizes`:

```bash
ls -lh "$GENOME_FA"
ls -lh ../day1_HiC_practice/data/reference/chrom.sizes
```

Строим индекс STAR:

```bash
STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir ../day1_HiC_practice/data/reference/star_index \
  --genomeFastaFiles ../day1_HiC_practice/data/reference/T2T_human.fna \
  --genomeSAindexNbases 13
```
Если индекс уже есть в `../day1_HiC_practice/data/reference/star_index`, этот шаг пропускаем.

## 6. Выравниваем RNA-seq для одного образца MoPh7

```bash
mkdir -p results/star/rnaseq results/logs

SAMPLE="MoPh7"
THREADS=8
STAR_INDEX="../day1_HiC_practice/data/reference/star_index"

R1="data/raw/rnaseq/${SAMPLE}_R1.fastq.gz"
R2="data/raw/rnaseq/${SAMPLE}_R2.fastq.gz"
PREFIX="results/star/rnaseq/${SAMPLE}_"

STAR \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$PREFIX" \
  --outSAMattrRGline "ID:${SAMPLE}_rnaseq" "SM:${SAMPLE}" "PL:ILLUMINA"
```

STAR создаст файл с длинным именем.
Переименуем его проще:

```bash
mv \
  "results/star/rnaseq/${SAMPLE}_Aligned.sortedByCoord.out.bam" \
  "results/star/rnaseq/${SAMPLE}.rnaseq.STAR.bam"
```

Индексируем BAM:

```bash
samtools index -@ "$THREADS" "results/star/rnaseq/${SAMPLE}.rnaseq.STAR.bam"
```

Сделаем отдельный BAM только с хорошо выровненными reads:

```bash
samtools view \
  -@ "$THREADS" \
  -b \
  -q 30 \
  "results/star/rnaseq/${SAMPLE}.rnaseq.STAR.bam" \
  > "results/star/rnaseq/${SAMPLE}.rnaseq.STAR.q30.bam"

samtools index -@ "$THREADS" "results/star/rnaseq/${SAMPLE}.rnaseq.STAR.q30.bam"
```

Посмотрим статистику STAR:

```bash
cat "results/star/rnaseq/${SAMPLE}_Log.final.out"
```

Полезные строки:

- `Number of input reads`
- `Uniquely mapped reads %`
- `% of reads mapped to multiple loci`
- `% of reads unmapped`
- `Number of splices`

## 7. Получаем RNA-seq bedGraph и bigWig

Для IGV удобно сделать coverage tracks.

```bash
mkdir -p results/tracks/rnaseq

BAM="results/star/rnaseq/${SAMPLE}.rnaseq.STAR.bam"
CHROM_SIZES="../day1_HiC_practice/data/reference/chrom.sizes"

bedtools genomecov \
  -ibam "$BAM" \
  -bg \
  -split \
  > "results/tracks/rnaseq/${SAMPLE}.rnaseq.STAR.bedGraph"
```

Опция `-split` важна для RNA-seq BAM, потому что reads могут быть выровнены через splice junctions.

Сортируем `bedGraph`:

```bash
sort -k1,1 -k2,2n \
  "results/tracks/rnaseq/${SAMPLE}.rnaseq.STAR.bedGraph" \
  > "results/tracks/rnaseq/${SAMPLE}.rnaseq.STAR.sorted.bedGraph"
```

Конвертируем в `bigWig`:

```bash
bedGraphToBigWig \
  "results/tracks/rnaseq/${SAMPLE}.rnaseq.STAR.sorted.bedGraph" \
  "$CHROM_SIZES" \
  "results/tracks/rnaseq/${SAMPLE}.rnaseq.STAR.bw"
```

## 8. Пример: выравниваем RNA-seq с помощью bwa

BWA не умеет разделять рид по сайту сплайсинга.
Мы используем его как пример: так будет видно, как меняется трек, если не учитывать сайты сплайсинга.

Выравниваем тот же RNA-seq образец `MoPh7`:

```bash
mkdir -p results/bwa/rnaseq

SAMPLE="MoPh7"
THREADS=8
BWA_INDEX="../day1_HiC_practice/data/reference/T2T_human.fna"

R1="data/raw/rnaseq/${SAMPLE}_R1.fastq.gz"
R2="data/raw/rnaseq/${SAMPLE}_R2.fastq.gz"

bwa mem \
  -t "$THREADS" \
  "$BWA_INDEX" \
  "$R1" "$R2" \
  | samtools sort \
      -@ "$THREADS" \
      -o "results/bwa/rnaseq/${SAMPLE}.rnaseq.BWA.bam"

samtools index -@ "$THREADS" "results/bwa/rnaseq/${SAMPLE}.rnaseq.BWA.bam"
```

Оставим только reads с `MAPQ >= 30`:

```bash
samtools view \
  -@ "$THREADS" \
  -b \
  -q 30 \
  "results/bwa/rnaseq/${SAMPLE}.rnaseq.BWA.bam" \
  > "results/bwa/rnaseq/${SAMPLE}.rnaseq.BWA.q30.bam"

samtools index -@ "$THREADS" "results/bwa/rnaseq/${SAMPLE}.rnaseq.BWA.q30.bam"
```

Делаем BWA coverage track для IGV:

```bash
mkdir -p results/tracks/bwa

BAM="results/bwa/rnaseq/${SAMPLE}.rnaseq.BWA.q30.bam"
CHROM_SIZES="../day1_HiC_practice/data/reference/chrom.sizes"

bedtools genomecov \
  -ibam "$BAM" \
  -bg \
  > "results/tracks/bwa/${SAMPLE}.rnaseq.BWA.q30.bedGraph"

sort -k1,1 -k2,2n \
  "results/tracks/bwa/${SAMPLE}.rnaseq.BWA.q30.bedGraph" \
  > "results/tracks/bwa/${SAMPLE}.rnaseq.BWA.q30.sorted.bedGraph"

bedGraphToBigWig \
  "results/tracks/bwa/${SAMPLE}.rnaseq.BWA.q30.sorted.bedGraph" \
  "$CHROM_SIZES" \
  "results/tracks/bwa/${SAMPLE}.rnaseq.BWA.q30.bw"
```

Теперь в IGV можно открыть рядом:

```text
results/tracks/rnaseq/MoPh7.rnaseq.STAR.bw
results/tracks/bwa/MoPh7.rnaseq.BWA.q30.bw
```

## 9. Выравниваем CAGE для MoPh7

CAGE reads single-end, поэтому передаем STAR один fastq

```bash
mkdir -p results/star/cage results/tracks/cage

SAMPLE="MoPh7"
R1="data/raw/cage/${SAMPLE}_R1.fastq.gz"
PREFIX="results/star/cage/${SAMPLE}_"

STAR \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$R1" \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$PREFIX" \
  --outSAMattrRGline "ID:${SAMPLE}_cage" "SM:${SAMPLE}" "PL:ILLUMINA"
```

Переименуем и проиндексируем BAM:

```bash
mv \
  "results/star/cage/${SAMPLE}_Aligned.sortedByCoord.out.bam" \
  "results/star/cage/${SAMPLE}.cage.STAR.bam"

samtools index -@ "$THREADS" "results/star/cage/${SAMPLE}.cage.STAR.bam"
```

Сделаем bam с `MAPQ >= 30`:

```bash
samtools view \
  -@ "$THREADS" \
  -b \
  -q 30 \
  "results/star/cage/${SAMPLE}.cage.STAR.bam" \
  > "results/star/cage/${SAMPLE}.cage.STAR.q30.bam"

samtools index -@ "$THREADS" "results/star/cage/${SAMPLE}.cage.STAR.q30.bam"
```

Делаем `bedGraph` и `bigWig`:

```bash
BAM="results/star/cage/${SAMPLE}.cage.STAR.bam"

bedtools genomecov \
  -ibam "$BAM" \
  -bg \
  -split \
  > "results/tracks/cage/${SAMPLE}.cage.STAR.bedGraph"

sort -k1,1 -k2,2n \
  "results/tracks/cage/${SAMPLE}.cage.STAR.bedGraph" \
  > "results/tracks/cage/${SAMPLE}.cage.STAR.sorted.bedGraph"

bedGraphToBigWig \
  "results/tracks/cage/${SAMPLE}.cage.STAR.sorted.bedGraph" \
  "$CHROM_SIZES" \
  "results/tracks/cage/${SAMPLE}.cage.STAR.bw"
```

## 11. Что должно получиться

Должны появиться:

```text
results/star/rnaseq/MoPh7.rnaseq.STAR.bam
results/star/rnaseq/MoPh7.rnaseq.STAR.bam.bai
results/star/rnaseq/MoPh7.rnaseq.STAR.q30.bam
results/star/rnaseq/MoPh7.rnaseq.STAR.q30.bam.bai
results/tracks/rnaseq/MoPh7.rnaseq.STAR.sorted.bedGraph
results/tracks/rnaseq/MoPh7.rnaseq.STAR.bw

results/bwa/rnaseq/MoPh7.rnaseq.BWA.bam
results/bwa/rnaseq/MoPh7.rnaseq.BWA.bam.bai
results/bwa/rnaseq/MoPh7.rnaseq.BWA.q30.bam
results/bwa/rnaseq/MoPh7.rnaseq.BWA.q30.bam.bai
results/tracks/bwa/MoPh7.rnaseq.BWA.q30.sorted.bedGraph
results/tracks/bwa/MoPh7.rnaseq.BWA.q30.bw

results/star/cage/MoPh7.cage.STAR.bam
results/star/cage/MoPh7.cage.STAR.bam.bai
results/star/cage/MoPh7.cage.STAR.q30.bam
results/star/cage/MoPh7.cage.STAR.q30.bam.bai
results/tracks/cage/MoPh7.cage.STAR.sorted.bedGraph
results/tracks/cage/MoPh7.cage.STAR.bw
```

Файлы `.bw` можно открыть в IGV рядом с Hi-C compartment tracks.

## 12. Запуск для всех образцов

На занятии полный прогон всех образцов может занять много времени.
Для самостоятельной работы есть готовый скрипт:

```bash
bash scripts/run_all_star_tracks.sh
```

Он проходит по образцам:

```text
MoPh7
MoPh11
MoPh14
MoPh15
```

Для каждого образца скрипт:

1. выравнивает RNA-seq paired-end reads;
2. выравнивает CAGE single-end reads;
3. индексирует BAM;
4. делает `MAPQ >= 30` BAM;
5. делает `bedGraph`;
6. делает `bigWig`.

Если BAM уже существует, alignment пропускается.
Если нужно пересчитать все заново:

```bash
FORCE=1 bash scripts/run_all_star_tracks.sh
```

## Задание

После самостоятельного прогона всех образцов:

1. Откройте в IGV `bigWig` треки RNA-seq и CAGE для `MoPh7`, `MoPh11`, `MoPh14`, `MoPh15`.
2. Найдите геномный регион, где сигнал отличается между образцами
3. Сравните RNA-seq и CAGE в этом регионе
4. Для `MoPh7` сравните STAR и BWA треки.
5. Подумайте, что лучше отражает общий уровень транскрипта, а что лучше показывает активность промотора.
