# Практика 1. Hi-C: подготовка ридов и получение `.hic`

Курс: методы анализа данных структурно-функциональной организации хроматина.

Этот репозиторий будет использоваться для практических занятий. В первой практике мы
разберем базовый путь от сырых paired-end Hi-C ридов до файла `.hic`, который можно
открыть в Juicebox и использовать дальше для анализа контактных карт.

## Что будет сделано

1. Подготовка ридов к выравниванию:
   - скачивание FASTQ;
   - первичный контроль качества через FastQC;
   - обрезка адаптеров и низкокачественных концов через cutadapt;
   - повторный FastQC после обрезки.
2. Получение Hi-C карты через Juicer:
   - установка Juicer tools;
   - подготовка структуры директорий;
   - запуск пайплайна Juicer на локальных данных;
   - проверка результата `.hic`.

Дальше на курсе мы будем работать с более глубокой готовой Hi-C картой и ноутбуками
на основе материалов из репозитория:
<https://github.com/dpanc2/BI_HiC_analysis>

## Требования

Практика рассчитана на локальный компьютер. Команды ниже предполагают macOS или Linux
и установленный `conda`/`mamba`.

Нужные инструменты:

- `fastqc`
- `cutadapt`
- `multiqc` опционально, но удобно для сводного отчета
- `bwa`
- `samtools`
- Java 8 или новее
- Juicer tools

Создадим отдельное окружение:

```bash
conda create -n hic_practice -c conda-forge -c bioconda \
  fastqc cutadapt multiqc bwa samtools openjdk=11 wget
conda activate hic_practice
```

Проверка:

```bash
fastqc --version
cutadapt --version
bwa 2>&1 | head
samtools --version
java -version
```

## Структура проекта

```text
chromatin-structure-function-course/
├── README.md
├── data/
│   ├── raw/
│   ├── trimmed/
│   ├── reference/
│   └── juicer/
├── results/
│   ├── fastqc_raw/
│   ├── fastqc_trimmed/
│   └── hic/
├── notebooks/
└── scripts/
```

Создать директории:

```bash
mkdir -p data/raw data/trimmed data/reference data/juicer
mkdir -p results/fastqc_raw results/fastqc_trimmed results/hic
mkdir -p notebooks scripts
```

## Шаг 1. Скачать сырые риды

В этой практике используются paired-end риды:

- `Copy of MoPh7_S85_L001_R1_001.fastq.gz`
- `Copy of MoPh7_S85_L001_R2_001.fastq.gz`

Скачаем их в `data/raw/` и сразу дадим короткие имена без пробелов:

```bash
wget -O data/raw/MoPh7_R1.fastq.gz \
  "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20MoPh7_S85_L001_R1_001.fastq.gz"

wget -O data/raw/MoPh7_R2.fastq.gz \
  "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20MoPh7_S85_L001_R2_001.fastq.gz"
```

Проверим, что файлы скачались:

```bash
ls -lh data/raw/
```

## Шаг 2. Первичный контроль качества FastQC

Запустим FastQC для обоих файлов:

```bash
fastqc \
  data/raw/MoPh7_R1.fastq.gz \
  data/raw/MoPh7_R2.fastq.gz \
  -o results/fastqc_raw
```

Опционально можно собрать общий отчет MultiQC:

```bash
multiqc results/fastqc_raw -o results/fastqc_raw
```

Что посмотреть в отчете:

- качество по позициям рида;
- наличие адаптеров;
- overrepresented sequences;
- распределение GC;
- длину ридов.

## Шаг 3. Обрезка адаптеров и низкокачественных концов cutadapt

Для paired-end данных запускаем `cutadapt` сразу на двух FASTQ.

Базовая команда:

```bash
cutadapt \
  -q 20 \
  -m 20 \
  -o data/trimmed/MoPh7_R1.trimmed.fastq.gz \
  -p data/trimmed/MoPh7_R2.trimmed.fastq.gz \
  data/raw/MoPh7_R1.fastq.gz \
  data/raw/MoPh7_R2.fastq.gz
```

Параметры:

- `-q 20` обрезает концы ридов с качеством ниже Q20;
- `-m 20` удаляет пары, где хотя бы один рид после обрезки стал короче 20 нуклеотидов;
- `-o` задает файл для первого рида;
- `-p` задает файл для второго рида.

Если FastQC показал конкретные адаптеры, их можно указать явно:

```bash
cutadapt \
  -a ADAPTER_FOR_R1 \
  -A ADAPTER_FOR_R2 \
  -q 20 \
  -m 20 \
  -o data/trimmed/MoPh7_R1.trimmed.fastq.gz \
  -p data/trimmed/MoPh7_R2.trimmed.fastq.gz \
  data/raw/MoPh7_R1.fastq.gz \
  data/raw/MoPh7_R2.fastq.gz
```

После обрезки снова запускаем FastQC:

```bash
fastqc \
  data/trimmed/MoPh7_R1.trimmed.fastq.gz \
  data/trimmed/MoPh7_R2.trimmed.fastq.gz \
  -o results/fastqc_trimmed

multiqc results/fastqc_trimmed -o results/fastqc_trimmed
```

## Шаг 4. Установка Juicer tools

Официальный репозиторий Juicer:
<https://github.com/aidenlab/juicer>

Для локальной практики нам нужен `juicer_tools.jar`.

```bash
mkdir -p tools
wget -O tools/juicer_tools.jar \
  https://github.com/aidenlab/Juicebox/releases/download/v2.20.00/juicer_tools.2.20.00.jar
```

Проверка:

```bash
java -jar tools/juicer_tools.jar
```

Если команда выводит справку Juicer tools, установка прошла успешно.

## Шаг 5. Подготовить геномный референс

Для полноценного запуска Juicer нужны:

- FASTA референсного генома;
- индекс `bwa`;
- файл размеров хромосом;
- файл сайтов рестрикции для выбранного фермента.

На занятии важно обсудить, что Hi-C библиотека зависит от протокола и фермента
рестрикции. Для Juicer это не техническая мелочь, а часть модели эксперимента.

Пример для небольшого учебного референса:

```bash
# FASTA кладем в data/reference/genome.fa
bwa index data/reference/genome.fa

samtools faidx data/reference/genome.fa
cut -f1,2 data/reference/genome.fa.fai > data/reference/chrom.sizes
```

Файл сайтов рестрикции обычно готовится скриптами из Juicer. В реальном анализе нужно
выбрать фермент, которым была приготовлена библиотека, например `MboI`, `DpnII`,
`HindIII`.

## Шаг 6. Подготовить структуру для Juicer

Juicer ожидает определенную структуру директорий. Для одного образца можно сделать так:

```bash
mkdir -p data/juicer/MoPh7/fastq

ln -s ../../trimmed/MoPh7_R1.trimmed.fastq.gz \
  data/juicer/MoPh7/fastq/MoPh7_R1.fastq.gz

ln -s ../../trimmed/MoPh7_R2.trimmed.fastq.gz \
  data/juicer/MoPh7/fastq/MoPh7_R2.fastq.gz
```

Если Juicer установлен из GitHub-репозитория:

```bash
git clone https://github.com/aidenlab/juicer.git tools/juicer
```

Основной shell-пайплайн находится в:

```text
tools/juicer/CPU/juicer.sh
```

## Шаг 7. Запуск Juicer

Пример команды для локального CPU-запуска:

```bash
bash tools/juicer/CPU/juicer.sh \
  -d data/juicer/MoPh7 \
  -z data/reference/genome.fa \
  -p data/reference/chrom.sizes \
  -y data/reference/restriction_sites.txt \
  -s MboI \
  -t 4
```

Параметры:

- `-d` директория эксперимента Juicer;
- `-z` FASTA референсного генома;
- `-p` размеры хромосом;
- `-y` файл сайтов рестрикции;
- `-s` фермент рестрикции;
- `-t` число потоков.

После успешного запуска Juicer итоговый `.hic` обычно появляется внутри директории
эксперимента, например:

```text
data/juicer/MoPh7/aligned/inter.hic
```

Скопируем результат в общую папку:

```bash
cp data/juicer/MoPh7/aligned/inter.hic results/hic/MoPh7.inter.hic
```

Проверка:

```bash
ls -lh results/hic/
java -jar tools/juicer_tools.jar dump observed NONE \
  results/hic/MoPh7.inter.hic 1 1 BP 1000000 | head
```

## Шаг 8. Установка Juicebox для визуализации

Juicebox нужен для интерактивного просмотра `.hic` карт.

Страница релизов:
<https://github.com/aidenlab/Juicebox/releases>

Для macOS или Linux можно скачать `.jar`:

```bash
wget -O tools/Juicebox.jar \
  https://github.com/aidenlab/Juicebox/releases/download/v2.20.00/Juicebox.jar
```

Запуск:

```bash
java -jar tools/Juicebox.jar
```

В Juicebox:

1. `File` -> `Open`;
2. выбрать `results/hic/MoPh7.inter.hic`;
3. выбрать хромосому или весь геном;
4. менять разрешение карты;
5. сравнить вид сырой и нормализованной матрицы, если нормализация доступна.

## Что сдать после этой части

Минимальный результат:

- HTML-отчеты FastQC до обрезки;
- HTML-отчеты FastQC после обрезки;
- короткий текстовый вывод `cutadapt`;
- файл `.hic` или описание, на каком шаге Juicer остановился;
- скриншот открытой карты в Juicebox.

## Дальше

В следующих блоках практики будут добавлены:

- работа с готовой глубокой Hi-C картой;
- чтение `.cool`/`.mcool`;
- анализ через `cooler` и `cooltools`;
- визуализация контактных карт в Jupyter Notebook;
- сопоставление Hi-C с RNA-seq и CAGE-seq.
