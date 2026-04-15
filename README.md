# Clonogenic Scan Pipeline

Строгий scan-only пайплайн для анализа клоногенного теста на TIFF-сканах планшетов.

Этот репозиторий предназначен только для сканов. Фото планшетов, JPEG/PNG со смартфона, низкое разрешение и смешанные наборы данных здесь не поддерживаются сознательно. Если передать не-скан, Java-анализатор завершится ошибкой.

## Что делает проект

1. Находит лунки на TIFF-скане по заданной пользователем схеме.
2. Сегментирует окрашенные колонии внутри каждой лунки.
3. Пытается разделять слипшиеся пятна только там, где это оправдано.
4. Считает колонии по каждой лунке и строит overlay с контурами и номерами.
5. Собирает средние по тем группам, которые пользователь задал в manifest/config.
6. При необходимости строит CDI-таблицы и heatmap.

## Что лежит в репозитории

```text
clonogenic_scan_pipeline/
├── README.md
├── requirements.txt
├── .gitignore
├── data/
│   └── .gitkeep
├── colab/
│   └── Clonogenic_Scan_Pipeline_Colab.ipynb
├── examples/
│   ├── settings.properties
│   ├── well_layout.csv
│   └── plate_manifest.csv
├── results/
│   └── .gitkeep
├── scripts/
│   ├── run_clonogenic_java.sh
│   ├── run_material_pipeline.sh
│   ├── build_cdi_table.py
│   └── heatmap.py
├── src/
│   └── java/
│       ├── ClonogenicAnalyzer.java
│       └── WellMaskRefiner.java
└── templates/
    └── CDI UiO-66 (Zr).xlsx
```

## Зависимости

### Java

- Java 17+ рекомендуется.
- Fiji / ImageJ с доступом к `ij-*.jar`.

Скрипт `scripts/run_clonogenic_java.sh` ищет jar так:

1. `FIJI_IJ_JAR`
2. `IJ_JAR`
3. `/home/maksegr/Applications/Fiji.app/jars/ij-1.54p.jar`

Перед запуском на новой машине задайте:

```bash
export FIJI_IJ_JAR=/absolute/path/to/Fiji.app/jars/ij-1.54p.jar
```

### Python

- Python 3.10+ рекомендуется.
- Пакеты из `requirements.txt`:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Формат входных файлов

Поддерживаются только `.tif` / `.tiff`.

Есть два режима работы.

### Режим 1. Современный и рекомендуемый

Используйте три внешних файла:

- `settings.properties` — основные параметры анализа;
- `well_layout.csv` — геометрия лунок;
- `plate_manifest.csv` — метаданные по каждой лунке каждого изображения.

Это основной режим для других пользователей, других клеточных линий, других доз, других концентраций, другого числа повторностей и другой раскладки лунок.

Примеры лежат в:

- `examples/settings.properties`
- `examples/well_layout.csv`
- `examples/plate_manifest.csv`

### Режим 2. Наследуемый

Если `manifest` не передан, анализатор может читать часть метаданных из имени файла. Имя должно иметь вид:

```text
<dose>_<unit>_<left-column concentration>_<right-column concentration>[_<capture-index>].tif
```

Примеры:

- `0_Gy_0,3_0,15.tif`
- `2_Gy_0,0_0,075.tif`
- `4_Gy_0_0,75_2.tif`

### Что задает пользователь в `settings.properties`

Минимально рекомендуемые поля:

- `dataset_mode=hf|zr|auto`
- `min_colony_area=...`
- `summary_group_fields=concentration`

Через этот файл удобно менять:

- минимальный размер колонии;
- основные пороги сегментации;
- поля, по которым надо усреднять лунки внутри одной пластины.

Если вы включаете CDI-режим, в `summary_group_fields` должен участвовать `concentration`, потому что CDI-скрипт ожидает именно это поле.

### Что задает пользователь в `well_layout.csv`

CSV задает стартовые центры лунок.

Обязательные колонки:

- `well_index`
- `x_fraction`
- `y_fraction`
- `radius_fraction`

Смысл:

- `x_fraction` и `y_fraction` — координаты центра лунки в долях ширины и высоты рабочей области;
- `radius_fraction` — стартовый радиус в долях `min(width, height)`.

Это позволяет работать не только с классической схемой 6 лунок.

### Что задает пользователь в `plate_manifest.csv`

Это таблица с одной строкой на одну лунку одного изображения.

Обязательные колонки:

- `image_name`
- `well_index`

Допустимы как CSV с запятыми, так и CSV с точкой с запятой. Для русской локали и имен файлов с запятыми практичнее использовать `;`.

Практически обязательные для нормальной биологической интерпретации:

- `dose_gy`
- `dose_unit`
- `replicate`
- `condition_name`
- `condition_value`
- `cell_line`

Допустимо добавлять любые собственные поля:

- `sample_name`
- `drug_name`
- `plate_id`
- `batch`
- `operator`
- любые другие колонки

Они будут перенесены в `well_counts.csv`, `segments.csv` и могут участвовать в `group_summary.csv`.

### Наследуемый режим и старое правило Hf / Zr

В CSV исторически используются метки `top` и `bottom`, но для сканов это фактически:

- `well 1-3` = левая колонка лунок
- `well 4-6` = правая колонка лунок

Для Hf действует базовое правило:

- левая колонка = `top_concentration`
- правая колонка = `bottom_concentration`

Для Zr действует исключение:

- у всех пластин из папки `Clon ZR Egor` концентрации инвертируются справа налево
- исключение: контроль `0_Gy_0_0,75.tif` остается в обычной ориентации

Это правило уже зашито в `ClonogenicAnalyzer.java`. Не дублируйте его вручную в Excel.

### Явный выбор набора данных

Чтобы не зависеть от имени папки, анализатор принимает:

```bash
--dataset auto|hf|zr
```

Рекомендуемое правило:

- локально и в Google Colab всегда передавать `--dataset hf` или `--dataset zr` явно;
- `auto` оставлять только для старых запусков, где структура папок сохранена как раньше.

## Быстрый запуск

### 1. Универсальный count-only запуск

Этот вариант подходит для большинства новых пользователей.

```bash
scripts/run_material_pipeline.sh \
  --input-dir /absolute/path/to/tiff_scans \
  --output-dir results/generic_run \
  --config examples/settings.properties \
  --layout examples/well_layout.csv \
  --manifest examples/plate_manifest.csv \
  --skip-cdi
```

### 2. Только Java-анализ

```bash
scripts/run_clonogenic_java.sh --output-dir results/hf --dataset hf /absolute/path/to/*.tif
```

### 3. Полный пайплайн с CDI

Этот режим имеет смысл только для экспериментов, которые действительно укладываются в матрицу доза × условие.

#### Hf

```bash
scripts/run_material_pipeline.sh \
  --input-dir /absolute/path/to/Clon_Egor \
  --output-dir results/hf \
  --material "UiO-66 (Hf)" \
  --config examples/settings.properties \
  --dataset hf
```

#### Zr

```bash
scripts/run_material_pipeline.sh \
  --input-dir /absolute/path/to/Clon\ ZR\ Egor \
  --output-dir results/zr \
  --material "UiO-66 (Zr)" \
  --config examples/settings.properties \
  --dataset zr
```

## Google Colab

Для пользователей без локальной настройки среды подготовлен ноутбук:

- `colab/Clonogenic_Scan_Pipeline_Colab.ipynb`

Как использовать:

1. Залить этот репозиторий на GitHub.
2. Открыть ноутбук `colab/Clonogenic_Scan_Pipeline_Colab.ipynb` в Google Colab.
3. В первой ячейке вставить URL своего GitHub-репозитория.
4. Выполнить ячейки сверху вниз:
   - установка зависимостей;
   - upload TIFF;
   - upload или редактирование `well_layout.csv` и `plate_manifest.csv`;
   - выбор основных параметров вроде `DATASET`, `MIN_COLONY_AREA`, `RUN_CDI`;
   - запуск анализа;
   - скачивание ZIP с результатами.

Ноутбук специально использует тот же shell-пайплайн и тот же Java-код, что и локальный запуск. Это важно: локально и в Colab должны считаться одинаковые результаты при одинаковых входных TIFF.

## Что появится в `results/<run>`

- `*/well_counts.csv` — counts по лункам
- `*/segments.csv` — все сегменты с морфометрией
- `*/group_summary.csv` — средние по группе внутри одной пластины
- `*/plate_overlay.png` — исходное изображение с контурами и номерами
- `*/well_XX_overlay.png` — overlay по отдельной лунке
- `duplicate_well_means.csv` — среднее по дубликатам одной и той же пластины
- `duplicate_group_means.csv` — среднее по группе после объединения дубликатов
- `run_summary.json` — сводка по полному запуску
- `cdi_*_long.csv` — длинная таблица CDI, если CDI включен
- `cdi_*_verification.csv` — плоская таблица для проверки, если CDI включен
- `cdi_*_verification.xlsx` — таблица с форматированием шаблона, если CDI включен
- `cdi_*_heatmap.png` — финальная heatmap, если CDI включен

## Где крутить параметры

Для обычной работы:

- `examples/settings.properties`
- `examples/well_layout.csv`
- `examples/plate_manifest.csv`

Для продвинутой доработки:

- `src/java/ClonogenicAnalyzer.java`

Если стандартных полей `settings.properties` уже не хватает, тогда имеет смысл менять `TUNE_*` в Java-коде.

## Правила доработки проекта

1. Не ломать scan-only режим.
   Если нужен анализ фото, делайте отдельную ветку или отдельный репозиторий.

2. Не менять silently правило ориентации Zr.
   Любое изменение раскладки должно сопровождаться комментарием в README и отдельной проверкой на `0_Gy_0_0,75.tif`.

3. Не ломать режим `manifest + layout`.
   Для новых пользователей это основной интерфейс.

4. Не редактировать CDI вручную.
   Источник истины — `group_summary.csv` и `summary.json`.

5. Любую новую эвристику добавлять так, чтобы она была параметризована через `TUNE_*`, а не зашита напрямую.

6. Перед commit обязательно проверять:
   - хотя бы один Hf прогон
   - хотя бы один Zr прогон
   - что `cdi_*_long.csv` и `cdi_*_heatmap.png` пересчитаны заново

## Минимальная проверка после изменений

```bash
scripts/run_clonogenic_java.sh \
  --output-dir results/smoke_test \
  --config examples/settings.properties \
  --layout examples/well_layout.csv \
  --manifest examples/plate_manifest.csv \
  /absolute/path/to/0_Gy_0,3_0,15.tif
```

Если этот запуск не собирается или не создает `plate_overlay.png`, изменения нельзя считать готовыми.

## Что не хранить в git

Не коммитьте:

- исходные TIFF-сканы
- `.rar` архивы
- локальные `results/*`
- `build/`
- `.venv/`
- временные probe-папки
