# Data for day5 omics practice

Большие файлы не коммитим в git. Ноутбуки по умолчанию читают результаты из предыдущих дней курса.

Если нужно положить файлы локально в `day5_omics_practice/data/`, используйте такую структуру:

```text
data/bam/             RNA-seq или ChIP-seq BAM + BAI
data/bigwig/          bigWig tracks
data/chipseq/         ChIP-seq peaks: narrowPeak/broadPeak/BED
data/compartments/    BED/bedGraph с A/B compartments
data/boundaries/      BED с boundaries/domain borders
data/annotation/      TSS/genes BED, если есть
```

Минимальный формат compartments:

```text
chrom  start  end  compartment  score
```

Минимальный формат boundaries:

```text
chrom  start  end  boundary_id  score
```

Если готового boundaries BED нет, notebook 03 умеет построить простые boundary-like points из смены A/B compartment.
