# День 4. WGBS/Bismark: анализ метилирования и пересечения с эпигеномными треками

В этой практике мы не запускаем Bismark, это слишком долго. Работаем с готовым output Bismark для `MoPh7`

Цели практики:

- прочитать `*.bismark.cov.gz`
- посчитать coverage, beta-value и M-value
- построить QC-графики
- сделать `bedGraph` и `bigWig` треки для IGV
- сравнить метилирование с H3K27ac и H3K9me3 signal из ChIP-seq практики

## 1. Переходим в папку практики

```bash
cd ~/OMICS_course_spring_2026/day4_WGBS_practice
```

## 2. Окружение

Используем то же окружение курса:

```bash
conda activate hic_practice
```

Если пакетов не хватает:

```bash
mamba install -y -c conda-forge -c bioconda \
  pandas numpy matplotlib seaborn scipy statsmodels tqdm \
  pybigwig pysam jupyterlab
```

Или через pip для базовых notebooks:

```bash
pip install -r requirements.txt
```

## 3. Входные данные

Для практики используем только основной файл Bismark coverage:

```text
MoPh7_1_bismark_bt2_pe.bismark.cov.gz
```

Не копируем уже дополнительно обработанные файлы вроде `*_avg_beta.bedGraph`, `*_CpG.bed`, `*_high_meth.bed`.

## 4. Скачиваем подготовленные файлы

Запускаем скрипт

```bash
python3 scripts/download_bismark_data.py
```

После этого в `data/bismark/` появятся все файлы со страницы:

```text
data/bismark/MoPh7_1_bismark_bt2_PE_report.txt
data/bismark/MoPh7_1_bismark_bt2_pe.CpG_report.txt
data/bismark/MoPh7_1_bismark_bt2_pe.M-bias.txt
data/bismark/MoPh7_1_bismark_bt2_pe.bedGraph.gz
data/bismark/MoPh7_1_bismark_bt2_pe.bismark.cov.gz
data/bismark/MoPh7_1_bismark_bt2_pe.cytosine_context_summary.txt
data/bismark/MoPh7_1_bismark_bt2_pe_splitting_report.txt
```

Основной файл для работы:

```text
data/bismark/MoPh7_1_bismark_bt2_pe.bismark.cov.gz
```

Его не переименовываем: во всех notebooks используем это имя.

Формат `bismark.cov.gz`:

```text
chrom  start  end  methylation_percentage  count_methylated  count_unmethylated
```

## 5. Что лежит в output Bismark

В папке `data/bismark/` лежат не только координаты метилирования, но и отчеты Bismark. Их полезно посмотреть перед анализом

### Основной отчет выравнивания

```bash
less data/bismark/MoPh7_1_bismark_bt2_PE_report.txt
```

В этом файле общая статистика: сколько ридов было выровнено, сколько пар прошло выоавнивание, сколько прочтений оказалось уникально выровнено

### CpG report

```bash
less data/bismark/MoPh7_1_bismark_bt2_pe.CpG_report.txt
```

Формат колонок

```text
chromosome  position  strand  methylated_count  unmethylated_count  context  trinucleotide
```

### M-bias report

```bash
less data/bismark/MoPh7_1_bismark_bt2_pe.M-bias.txt
```

Формат основных колонок:

```text
position  count_methylated  count_unmethylated  methylation_percent  coverage
```

M-bias показывает, меняется ли оценка метилирования вдоль рида. Это важно для QC: иногда первые или последние позиции рида дают смещенную оценку

### Cytosine context summary

```bash
less data/bismark/MoPh7_1_bismark_bt2_pe.cytosine_context_summary.txt
```

Формат колонок:

```text
upstream  C-context  full_context  count_methylated  count_unmethylated  percent_methylation
```

Этот файл помогает посмотреть метилирование в разных контекстах: CpG, CHG, CHH

### bedGraph от Bismark

```bash
zless data/bismark/MoPh7_1_bismark_bt2_pe.bedGraph.gz
```

Это готовый genome track с процентом метилирования. Мы сделаем свои `bedGraph` и `bigWig` треки из `bismark.cov.gz`

## 6. H3K27ac и H3K9me3 треки

Для сравнения используем ChIP-seq tracks из 3го дня

```text
../day3_ChIPseq_practice/results/macs/MoPh7_H3K27Ac/MoPh7_H3K27Ac_FE.bw
../day3_ChIPseq_practice/results/macs/MoPh7_H3K9me3/MoPh7_H3K9me3_FE.bw
```

## 7. Порядок запуска notebooks

Запускаем по порядку

```text
notebooks/bismark_qc_and_methylation_values.ipynb
notebooks/create_bigwig_tracks_for_igv.ipynb
notebooks/integrate_methylation_with_h3k27ac.ipynb
```

## 8. Результаты

После запуска notebooks появятся:

```text
results/tables/MoPh7_cpg_methylation_values.tsv.gz
results/tracks/MoPh7_beta_methylation.bedGraph
results/tracks/MoPh7_beta_methylation.bw
results/tracks/MoPh7_m_value.bedGraph
results/tracks/MoPh7_m_value.bw
results/tracks/MoPh7_coverage.bedGraph
results/tracks/MoPh7_coverage.bw
results/tracks/T2T_gc_content_100bp.bw
results/tracks/T2T_cpg_obs_exp_100bp.bw
results/tables/MoPh7_methylation_chipseq_signals.tsv.gz
results/tables/MoPh7_methylation_chipseq_correlations.tsv
```

## 9. Что открыть в IGV

```text
results/tracks/MoPh7_beta_methylation.bw
results/tracks/MoPh7_coverage.bw
results/tracks/T2T_gc_content_100bp.bw
results/tracks/T2T_cpg_obs_exp_100bp.bw
../day3_ChIPseq_practice/results/macs/MoPh7_H3K27Ac/MoPh7_H3K27Ac_FE.bw
../day3_ChIPseq_practice/results/macs/MoPh7_H3K9me3/MoPh7_H3K9me3_FE.bw
```

## 10. Идеи дополнительных заданий

- Посчитать среднее метилирование промоторов
- Сравнить метилирование промоторов/энхансеров/CpG островков/фона
- Построить methylation profile вокруг TSS
- Добавить несколько samples и сравнить профиль метилирования между ними

## 11. Вопросы

- Почему CpG с низким coverage лучше удалить?
- Чем beta-value отличается от M-value?
- Почему H3K27ac может быть связан с низким метилированием в регуляторных областях?

