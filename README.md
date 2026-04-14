# Clonogenic Scan Pipeline

Строгий scan-only пайплайн для анализа клоногенного теста на TIFF-сканах 6-луночных планшетов.

Этот репозиторий предназначен только для сканов. Фото планшетов, JPEG/PNG со смартфона, низкое разрешение и смешанные наборы данных здесь не поддерживаются сознательно. Если передать не-скан, Java-анализатор завершится ошибкой.

## Что делает проект

1. Находит 6 лунок на TIFF-скане.
2. Сегментирует окрашенные колонии внутри каждой лунки.
3. Пытается разделять слипшиеся пятна только там, где это оправдано.
4. Считает колонии по каждой лунке и строит overlay с контурами и номерами.
5. Собирает средние по повторностям.
6. Строит CDI-таблицы и heatmap.

## Что лежит в репозитории

```text
clonogenic_scan_pipeline/
├── README.md
├── requirements.txt
├── .gitignore
├── data/
│   └── .gitkeep
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

Имя файла должно иметь вид:

```text
<dose>_<unit>_<left-column concentration>_<right-column concentration>[_<capture-index>].tif
```

Примеры:

- `0_Gy_0,3_0,15.tif`
- `2_Gy_0,0_0,075.tif`
- `4_Gy_0_0,75_2.tif`

### Важное замечание по раскладке лунок

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

## Быстрый запуск

### 1. Только Java-анализ

```bash
scripts/run_clonogenic_java.sh --output-dir results/hf /absolute/path/to/*.tif
```

### 2. Полный пайплайн по материалу

#### Hf

```bash
scripts/run_material_pipeline.sh \
  --input-dir /absolute/path/to/Clon_Egor \
  --output-dir results/hf \
  --material "UiO-66 (Hf)"
```

#### Zr

```bash
scripts/run_material_pipeline.sh \
  --input-dir /absolute/path/to/Clon\ ZR\ Egor \
  --output-dir results/zr \
  --material "UiO-66 (Zr)"
```

## Что появится в `results/<material>`

- `*/well_counts.csv` — counts по лункам
- `*/segments.csv` — все сегменты с морфометрией
- `*/group_summary.csv` — средние по группе внутри одной пластины
- `*/plate_overlay.png` — исходное изображение с контурами и номерами
- `*/well_XX_overlay.png` — overlay по отдельной лунке
- `duplicate_well_means.csv` — среднее по дубликатам одной и той же пластины
- `duplicate_group_means.csv` — среднее по группе после объединения дубликатов
- `run_summary.json` — сводка по полному запуску
- `cdi_4t1_long.csv` — длинная таблица CDI
- `cdi_4t1_verification.csv` — плоская таблица для проверки
- `cdi_4t1_verification.xlsx` — таблица с форматированием шаблона
- `cdi_4t1_heatmap.png` — финальная heatmap

## Где крутить параметры сегментации

Все ручки лежат в верхнем блоке `TUNE_*` в:

- `src/java/ClonogenicAnalyzer.java`

Менять нужно только параметры в этом блоке. Не разбрасывайте магические числа по коду.

Критические группы параметров:

- `TUNE_MIN_COLONY_AREA`
- `TUNE_MIN_AREA_SCALE_SCAN`
- `TUNE_SCAN_STRONG_THRESHOLD_FACTOR`
- `TUNE_SCAN_RESCUE_*`
- `TUNE_SPLIT_*`

## Правила доработки проекта

1. Не ломать scan-only режим.
   Если нужен анализ фото, делайте отдельную ветку или отдельный репозиторий.

2. Не менять silently правило ориентации Zr.
   Любое изменение раскладки должно сопровождаться комментарием в README и отдельной проверкой на `0_Gy_0_0,75.tif`.

3. Не менять формат имен входных файлов без обновления `parseMetadata`.

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
