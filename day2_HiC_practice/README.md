# День 2. Анализ готовых Hi-C карт

Во второй день мы работаем не с сырыми ридами, а с уже собранными Hi-C картами

Цель практики:

- скачать готовые `.hic` карты
- конвертировать их в multi-resolution cooler формат `.mcool`
- сбалансировать карты
- посмотреть контактные матрицы
- посчитать зависимость контактов от геномного расстояния
- посчитать insulation score и границы доменов
- посчитать A/B-компартменты и saddle plots
- сохранить E1-компартментный трек в `bedGraph` и `bigWig`

## 1. Переходим в папку практики

```bash
cd ~/OMICS_course_spring_2026/day2_HiC_practice
```

## 2. Активируем окружение

Можно использовать окружении с первого дня, главное, чтобы были установлены python-библиотеки:

- `cooler`
- `cooltools`
- `hic2cool`
- `bioframe`
- `pyBigWig`
- `pysam`
- `jupyter`
- `matplotlib`
- `seaborn`
- `pandas`
- `numpy`

```bash
conda activate hic_practice
conda install cooler cooltools hic2cool bioframe pybigwig pysam
```

## 3. Скачиваем и подготавливаем карты

Запускаем скрипт из папкт второго дня (скачиваем новые данные)

```bash
python3 scripts/download_data.py
```

Скрипт скачивает четыре `.hic` карты

- `MoPh7_enr_v2.hic`
- `MoPh11_enr_v2.hic`
- `MoPh14_enr_v2.hic`
- `MoPh15_enr_v2.hic`

Для каждой карты скрипт делает

1. скачивание в папку `data/`;
2. конвертацию `.hic` -> `.mcool`;
3. балансировку нужных разрешений: `10 kb`, `100 kb`, `1 Mb`.

Если файлы уже существуют, скрипт не скачивает их заново.
Если разрешение уже сбалансировано, скрипт тоже пропускает этот шаг.

## 4. Проверим результат

Посмотрим, какие разрешения есть в `.mcool` файле:

```bash
python3 -m cooler ls data/MoPh7_enr_v2.mcool
```

Проверим одно разрешение:

```bash
python3 -m cooler info data/MoPh7_enr_v2.mcool::/resolutions/100000
```

Если в таблице bins есть колонка `weight`, значит карта сбалансирована.

## 5. Открываем ноутбуки

Порядок работы:

1. `visualization.ipynb` Простая визуализация и отрисовка карт
2. `contacts_vs_distance.ipynb` Считаем, как средняя частота контактов зависит от геномного расстояния
3. `insulation_and_boundaries.ipynb` Инсуляция
4. `compartments_and_saddles.ipynb` А/В компартментализация

В конце ноутбук `compartments_and_saddles.ipynb` сохраняет E1-треки для всех образцов

```text
results/compartments/*_E1_res100000.bedGraph
results/compartments/*_E1_res100000.bw
```

## 6. Можно проверить `.hic` в Juicebox

Оригинальные `.hic` файлы остаются в папке `data/`
