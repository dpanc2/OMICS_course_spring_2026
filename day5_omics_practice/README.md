# День 5. Простая интеграция omics-данных

Сегодня работаем с готовыми файлами, которые получили в предыдущие дни: bam, bw, bed peaks, compartments и boundaries

## Что будем делать

1. Читать RNAseq bam через `pysam` и bw через `pyBigWig`
2. Сравнивать ChIPseq пики разных гистоновых меток с A/B компартментами
3. Смотреть, связаны ли CTCF пики с границами доменов
4. Строить метапрофили ChIPseq вокруг заданного геномного региона

## Окружение

Используем то же окружение курса:

```bash
conda activate hic_practice
```

Если не устанавливали пакеты

```bash
pip install -r requirements.txt
```

## Ноутбуки

Запускаем по порядку:

```text
notebooks/1_pysam_pybigwig_intro.ipynb
notebooks/2_histone_marks_vs_compartments.ipynb
notebooks/3_ctcf_vs_domain_boundaries.ipynb
```

## Откуда берутся данные

Используют файлы из предыдущих практик

- `../day2_RNA_practice/results/` — RNAseq bam и bw треки
- `../day2_HiC_practice/results/compartments/` — A/B компартменты
- `../day3_ChIPseq_practice/results/` — ChIPseq пики и bw треки

Большие файлы не храним в этой папке. Если хочется работать с копиями, положите их в `data/` по описанию в `data/README.md`.

## Базовый минимум, который нужно с собой унести

- как открыть bam и посмотреть риды
- что такое MAPQ, CIGAR, duplicate, secondary и supplementary reads
- как посчитать покрытие по bam
- как читать bw и доставать сигналы из региона
- как пересекать bed файлы(таблицы)
- как сравнить ChIPseq пики с A/B компартментами
- как оценить связь CTCF пики с границами доменов

