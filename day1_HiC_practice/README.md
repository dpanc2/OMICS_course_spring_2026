# 04.06. Hi-C: подготовка ридов и получение `.hic` файлов

Базовый анализ от сырых paired-end Hi-C ридов до файла `.hic`, который можно открыть
в Juicebox и использовать в дальнейшем для анализа карт контактов.

## План

1. Подготовка референсного генома:
   - скачивание FASTA;
   - распаковка `.fna.gz`;
   - индексация `bwa`;
   - создание `chrom.sizes`;
   - создание файла сайтов рестрикции для Juicer.
2. Подготовка ридов к выравниванию:
   - скачивание FASTQ;
   - первичный контроль качества через FastQC;
   - обрезка адаптеров и низкокачественных хвостов через cutadapt.
3. Получение Hi-C карты через Juicer:
   - установка Juicer;
   - подготовка структуры директорий;
   - запуск пайплайна Juicer;
   - получение `.hic`.

## Требования

Работаем локально на ноутбуке. Команды ниже предполагают macOS или Linux
и установленный `conda`/`mamba`.

Нужные инструменты:

- `fastqc`
- `cutadapt`
- `bwa`
- `samtools`
- `wget`
- `gzip`
- Java 8 или новее
- Juicer
- Juicebox

Создадим отдельное окружение:

```bash
conda create -n hic_practice -c conda-forge -c bioconda \
  fastqc cutadapt bwa samtools openjdk=11 wget
conda activate hic_practice
```

Проверка:

```bash
fastqc --version
cutadapt --version
bwa 2>&1 | head -3
samtools --version
java -version
```

Создадим рабочие папки:

```bash
mkdir -p data/raw data/trimmed data/reference data/juicer
mkdir -p results/fastqc_raw results/cutadapt results/hic
mkdir -p tools
```

## Шаг 1. Референсный геном

Начнем с референса, потому что скачивание и индексация генома занимают больше всего
времени.

Для запуска Juicer нам нужны:

- FASTA референсного генома;
- индекс `bwa`;
- файл размеров хромосом `chrom.sizes`;
- файл сайтов рестрикции для выбранного фермента.

В качестве референса берем T2T геном человека. На сервере NCBI файл лежит в gzip
формате, поэтому скачиваем его как `.fna.gz`, а затем распаковываем рабочую копию
`.fna`.

```bash
wget --no-check-certificate \
  -O data/reference/T2T_human.fna.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

gzip -dkf data/reference/T2T_human.fna.gz
```

Проверим, что появились оба файла:

```bash
ls -lh data/reference/T2T_human.fna.gz data/reference/T2T_human.fna
```

Индексируем распакованный FASTA:

```bash
bwa index data/reference/T2T_human.fna

samtools faidx data/reference/T2T_human.fna
cut -f1,2 data/reference/T2T_human.fna.fai > data/reference/chrom.sizes
```

## Шаг 2. Установка Juicer

Официальный репозиторий Juicer:
[https://github.com/aidenlab/juicer](https://github.com/aidenlab/juicer)

Скачаем Juicer целиком в папку `tools/`:

```bash
if [ ! -d tools/juicer ]; then
  git clone https://github.com/aidenlab/juicer.git tools/juicer
fi
```

Проверка:

```bash
ls tools/juicer/CPU/juicer.sh
ls tools/juicer/misc/generate_site_positions.py
```

Если команды показывают пути к `juicer.sh` и `generate_site_positions.py`, установка
прошла успешно.

## Шаг 3. Файл сайтов рестрикции

Файл сайтов рестрикции готовится скриптом из Juicer. В реальном анализе нужно выбрать
фермент, которым была приготовлена библиотека, например `MboI`, `DpnII`, `HindIII`.
В нашем случае используем фермент `DpnII`.

Для скрипта нужен распакованный FASTA, поэтому используем `data/reference/T2T_human.fna`,
а не `.fna.gz`.

```bash
python3 tools/juicer/misc/generate_site_positions.py \
  DpnII \
  T2T_human \
  data/reference/T2T_human.fna

mv T2T_human_DpnII.txt data/reference/restriction_sites_DpnII.txt
```

Скрипт принимает три аргумента: название фермента, короткое имя генома и FASTA-файл.
На выходе он создает файл формата `<genome>_<enzyme>.txt`, который передается в
Juicer через параметр `-y`.

## Шаг 4. Скачиваем сырые риды

В практике первого дня используются paired-end риды:

- `Copy of MoPh7_S85_L001_R1_001.fastq.gz`
- `Copy of MoPh7_S85_L001_R2_001.fastq.gz`

Скачаем их в `data/raw/`. Для этого сервера добавляем `--no-check-certificate`, чтобы
не упасть на ошибке проверки сертификата.

```bash
wget --no-check-certificate \
  -O data/raw/MoPh7_R1.fastq.gz \
  "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20MoPh7_S85_L001_R1_001.fastq.gz"

wget --no-check-certificate \
  -O data/raw/MoPh7_R2.fastq.gz \
  "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20MoPh7_S85_L001_R2_001.fastq.gz"
```

Проверим, что файлы скачались:

```bash
ls -lh data/raw/
```

## Шаг 5. FastQC

Запустим FastQC для обоих файлов:

```bash
fastqc \
  data/raw/MoPh7_R1.fastq.gz \
  data/raw/MoPh7_R2.fastq.gz \
  -o results/fastqc_raw
```

Что посмотреть в отчете:

- качество по позициям рида;
- наличие адаптеров;
- overrepresented sequences;
- распределение GC;
- длину ридов.

## Шаг 6. Обрезка адаптеров и низкокачественных концов cutadapt

Для paired-end данных запускаем `cutadapt` сразу на двух FASTQ. Отчет `cutadapt`
сохраним в отдельный log-файл.

```bash
cutadapt \
  -q 20 \
  -m 70 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -o data/trimmed/MoPh7_R1.trimmed.fastq.gz \
  -p data/trimmed/MoPh7_R2.trimmed.fastq.gz \
  data/raw/MoPh7_R1.fastq.gz \
  data/raw/MoPh7_R2.fastq.gz \
  > results/cutadapt/MoPh7.cutadapt.log 2>&1
```

Параметры:

- `-q 20` обрезает концы ридов с качеством ниже Q20;
- `-m 70` удаляет пары, где хотя бы один рид после обрезки стал короче 70 нуклеотидов;
- `-a` обрезает адаптер Illumina;
- `-o` задает файл для первого рида;
- `-p` задает файл для второго рида.

Проверим log:

```bash
less results/cutadapt/MoPh7.cutadapt.log
```

## Шаг 7. Подготовка к запуску Juicer

Juicer ожидает определенную структуру директорий. Для одного образца можно сделать так:

```bash
mkdir -p data/juicer/MoPh7/fastq

ln -s "$(pwd)/data/trimmed/MoPh7_R1.trimmed.fastq.gz" \
  data/juicer/MoPh7/fastq/MoPh7_R1.fastq.gz

ln -s "$(pwd)/data/trimmed/MoPh7_R2.trimmed.fastq.gz" \
  data/juicer/MoPh7/fastq/MoPh7_R2.fastq.gz
```

Проверим, что ссылки ведут на существующие файлы:

```bash
ls -lh data/juicer/MoPh7/fastq/
```

## Шаг 8. Запуск Juicer

Пример команды для локального CPU-запуска:

```bash
bash tools/juicer/CPU/juicer.sh \
  -d data/juicer/MoPh7 \
  -z data/reference/T2T_human.fna \
  -p data/reference/chrom.sizes \
  -y data/reference/restriction_sites_DpnII.txt \
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
эксперимента:

```text
data/juicer/MoPh7/aligned/inter.hic
```

Сохраним его в общей папке результатов:

```bash
cp data/juicer/MoPh7/aligned/inter.hic results/hic/MoPh7.inter.hic
```

## Шаг 9. Установка Juicebox для визуализации

Juicebox нужен для интерактивного просмотра `.hic` карт.

Страница релизов:
[https://github.com/aidenlab/Juicebox/releases](https://github.com/aidenlab/Juicebox/releases)

Для macOS или Linux можно скачать `.jar`:

```bash
wget --no-check-certificate \
  -O tools/Juicebox.jar \
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

## Задание первого дня

Напишите bash-pipeline для обработки нескольких Hi-C образцов от сырых paired-end
ридов до `.hic` файлов.

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
4. Сохранить log-файл `cutadapt` в `results/cutadapt/`.
5. Подготовить директорию `data/juicer/<sample>/fastq/`.
6. Запустить `Juicer` и получить файл `inter.hic`.
7. Скопировать итоговую карту в `results/hic/<sample>.inter.hic`.

### Ожидаемый результат

В конце работы пайплайна должны появиться `.hic` файлы для четырех образцов:

```text
results/hic/MoPh7.inter.hic
results/hic/MoPh11.inter.hic
results/hic/MoPh14.inter.hic
results/hic/MoPh15.inter.hic
```

### Финальный вопрос

После получения карт сравните четыре образца между собой в Juicebox:

- отличаются ли карты визуально;
- есть ли крупные перестройки;
- на каких хромосомах они заметны;
- какие дополнительные проверки нужны, чтобы убедиться, что это не артефакт
  покрытия, качества ридов или выравнивания.
