# 04.06. Hi-C: подготовка ридов и получение `.hic` файлов

Базовый анализ от сырых paired-end Hi-C ридов до файла `.hic`, который можно открыть в Juicebox и
использовать в дальнейшем для анализа карт контактов.

## План

1. Подготовка ридов к выравниванию:
   - скачивание fastq
   - первичный контроль качества через FastQC
   - обрезка адаптеров и низкокачественных хвостов через cutadapt
2. Получение Hi-C карты через Juicer:
   - установка Juicer tools
   - подготовка структуры директорий (важно для запуска Juicer)
   - запуск пайплайна Juicer
   - получение .hic

## Требования

Работаем локально на ноутбуке. Команды ниже предполагают macOS или Linux
и установленный `conda`/`mamba`

Нужные инструменты:

- `fastqc`
- `cutadapt`
- `bwa`
- `samtools`
- `hic2cool`
- `cooler`
- `cooltools`
- `requests`
- `tqdm`
- Java 8 или новее
- Juicer tools

Создадим отдельное окружение:

```bash
conda create -n hic_practice -c conda-forge -c bioconda \
  fastqc cutadapt bwa samtools hic2cool cooler cooltools requests tqdm openjdk=11 wget
conda activate hic_practice
```

Проверка:

```bash
fastqc --version
cutadapt --version
bwa 2>&1 | head -3
samtools --version
python3 -m hic2cool --help
cooler --version
python3 -c "import requests, tqdm"
java -version
```

## Шаг 1. Скачиваем сырые риды

В практике 1го дня используются paired-end риды

- `Copy of MoPh7_S85_L001_R1_001.fastq.gz`
- `Copy of MoPh7_S85_L001_R2_001.fastq.gz`

Скачаем их в `data/raw/`

```bash
wget -O data/raw/MoPh7_R1.fastq.gz \
  "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20MoPh7_S85_L001_R1_001.fastq.gz"

wget -O data/raw/MoPh7_R2.fastq.gz \
  "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20MoPh7_S85_L001_R2_001.fastq.gz"
```

Проверим, что файлы скачались

```bash
ls -lh data/raw/
```

## Шаг 2. FastQC

Запустим FastQC для обоих файлов:

```bash
fastqc \
  data/raw/MoPh7_R1.fastq.gz \
  data/raw/MoPh7_R2.fastq.gz \
  -o results/fastqc_raw
```

Что посмотреть в отчете:

- качество по позициям рида
- наличие адаптеров
- overrepresented sequences
- распределение GC
- длину ридов

## Шаг 3. Обрезка адаптеров и низкокачественных концов cutadapt

Для paired-end данных запускаем `cutadapt` сразу на двух FASTQ.

Базовая команда:

```bash
cutadapt \
  -q 20 \
  -m 70 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -o data/trimmed/MoPh7_R1.trimmed.fastq.gz \
  -p data/trimmed/MoPh7_R2.trimmed.fastq.gz \
  data/raw/MoPh7_R1.fastq.gz \
  data/raw/MoPh7_R2.fastq.gz
```

Параметры:

- `-q 20` обрезает концы ридов с качеством ниже Q20;
- `-m 20` удаляет пары, где хотя бы один рид после обрезки стал короче 70 нуклеотидов;
- `-a` обрезать адаптер Illumina;
- `-o` задает файл для первого рида;
- `-p` задает файл для второго рида.


## Шаг 4. Установка Juicer

Официальный репозиторий Juicer:
<https://github.com/aidenlab/juicer>

Скачаем Juicer целиком в папку `tools/`:

```bash
mkdir -p tools
git clone https://github.com/aidenlab/juicer.git tools/juicer
```

Проверка:

```bash
ls tools/juicer/CPU/juicer.sh
```

Если команда показывает путь к `juicer.sh`, установка прошла успешно.

## Шаг 5. Референсный геном

Для запуска Juicer нам нужны

- fasta референсного генома
- индекс `bwa`
- файл размеров хромосом `chrom.sizes`
- файл сайтов рестрикции для выбранного фермента

В качестве референса берем T2T геном человека:

```bash
# кладем в data/reference/
wget -O data/references/T2T_human.fna https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz 

bwa index data/references/T2T_human.fna

samtools faidx data/references/T2T_human.fna
cut -f1,2 data/references/T2T_human.fna.fai > data/references/chrom.sizes
```

## Шаг 6. Подготовка к запуску Juicer

Juicer ожидает определенную структуру директорий. Для одного образца можно сделать так:

```bash
mkdir -p data/juicer/MoPh7/fastq

ln -s ../../trimmed/MoPh7_R1.trimmed.fastq.gz \
  data/juicer/MoPh7/fastq/MoPh7_R1.fastq.gz

ln -s ../../trimmed/MoPh7_R2.trimmed.fastq.gz \
  data/juicer/MoPh7/fastq/MoPh7_R2.fastq.gz
```

Файл сайтов рестрикции обычно готовится скриптами из Juicer. В реальном анализе нужно
выбрать фермент, которым была приготовлена библиотека, например `MboI`, `DpnII`,
`HindIII`. В нашем случае мы используем фермент `DpnII`

Сгенерируем файл позиций сайтов рестрикции для T2T референса:

```bash
python tools/juicer/misc/generate_site_positions.py \
  DpnII \
  T2T_human \
  data/references/T2T_human.fna

mv T2T_human_DpnII.txt data/references/restriction_sites_DpnII.txt
```

Скрипт принимает три аргумента: название фермента, короткое имя генома и FASTA-файл.
На выходе он создает файл формата `<genome>_<enzyme>.txt`, который передается в
Juicer через параметр `-y`.

## Шаг 7. Запуск Juicer

Пример команды для локального CPU-запуска:

```bash
bash tools/juicer/CPU/juicer.sh \
  -d data/juicer/MoPh7 \
  -z data/references/T2T_human.fna \
  -p data/references/chrom.sizes \
  -y data/references/restriction_sites_DpnII.txt \
  -s DpnII \
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

Сохраним его в общей папке результатов:

```bash
mkdir -p results/hic
cp data/juicer/MoPh7/aligned/inter.hic results/hic/MoPh7.inter.hic
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

Запуск

```bash
java -jar tools/Juicebox.jar
```

В Juicebox:

1. `File` -> `Open`;
2. выбрать `results/hic/MoPh7.inter.hic`;
3. выбрать хромосому или весь геном;
4. менять разрешение карты;
5. сравнить вид сырой и нормализованной матрицы, если нормализация доступна.

Проверить список разрешений внутри `.mcool` можно так:

```bash
cooler ls data/Control_inter_30.mcool
```


# Задание первого дня

Напишите bash-pipeline для обработки нескольких Hi-C образцов от сырых paired-end
ридов до `.hic` файлов

### Входные данные

Общая папка с сырыми ридами:

```text
https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/
```

Для самостоятельной обработки используйте три дополнительных образца:

| Образец | R1 | R2 |
| --- | --- | --- |
| `MoPh11` | `Copy of MoPh11_S86_L001_R1_001.fastq.gz` | `Copy of MoPh11_S86_L001_R2_001.fastq.gz` |
| `MoPh14` | `Copy of MoPh14_S87_L001_R1_001.fastq.gz` | `Copy of MoPh14_S87_L001_R2_001.fastq.gz` |
| `MoPh15` | `Copy of MoPh15_S88_L001_R1_001.fastq.gz` | `Copy of MoPh15_S88_L001_R2_001.fastq.gz` |

Образец `MoPh7`, разобранный в практике, можно использовать как пример и контроль
для структуры директорий и команд.

### Что должен делать пайплайн

Для каждого образца:

1. Скачать paired-end FASTQ файлы в `data/raw/`.
2. Запустить `FastQC` для сырых ридов.
3. Обрезать адаптеры и низкокачественные хвосты с помощью `cutadapt`.
4. Подготовить директорию `data/juicer/<sample>/fastq/`.
5. Запустить `Juicer` и получить файл `inter_30.hic`.
6. Скопировать итоговую карту в `results/hic/<sample>.inter_30.hic`.

### Ожидаемый результат

В конце работы пайплайна должны появиться `30.hic` файлы для четырех образцов:

```text
results/hic/MoPh7.inter.hic
results/hic/MoPh11.inter.hic
results/hic/MoPh14.inter.hic
results/hic/MoPh15.inter.hic
```

### Финальный вопрос

После получения карт сравните четыре образца между собой в Juicebox:

- отличаются ли карты визуально
- есть ли крупные перестройки
- на каких хромосомах они заметны
- какие дополнительные проверки нужны, чтобы убедиться, что это не артефакт
  покрытия, качества ридов или выравнивания
