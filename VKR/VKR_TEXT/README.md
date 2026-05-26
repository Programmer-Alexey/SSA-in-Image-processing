# VKR LaTeX project

Главный файл проекта: `Report.tex`.

Структура:

- `Report.tex` -- корневой файл, подключает титульник, главы, заключение и библиографию.
- `preamble.tex` -- пакеты, математические операторы, окружения и общие настройки.
- `title/` -- титульные листы.
- `chapters/introduction/` -- введение.
- `chapters/chapter1_algorithms/` -- глава 1, описание модели и алгоритмов.
- `chapters/chapter2_cssa_preprocessing/` -- глава 2, теоретическое обоснование CSSA-предобработки.
- `chapters/chapter3_experiments/` -- глава 3, численные эксперименты.
- `chapters/chapter3_experiments/tables/` -- таблицы для главы 3.
- `chapters/conclusion/` -- заключение.
- `assets/images/` -- изображения, сгруппированные по главам.
- `additional/` -- файлы стиля и библиография.
- `scripts/` -- вспомогательные скрипты и R Markdown-файлы для генерации результатов.

## Сборка отчета

```powershell
pdflatex -interaction=nonstopmode -halt-on-error Report.tex
bibtex Report
pdflatex -interaction=nonstopmode -halt-on-error Report.tex
pdflatex -interaction=nonstopmode -halt-on-error Report.tex
```

## Воспроизведение таблиц и рисунков

Все команды ниже запускаются из папки `VKR/VKR_TEXT`.

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' scripts\make_report_figures.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' scripts\make_discretization_step_tables.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' scripts\make_chapter3_tables.R
```

Назначение файлов:

- `scripts/00_common.R` -- общие функции, пути проекта, CSSA/ESPRIT/HT-утилиты.
- `scripts/make_report_figures.R` -- генерирует рисунки для глав 1 и 2 в `assets/images/`.
- `scripts/make_discretization_step_tables.R` -- генерирует таблицу вклада дискретизации Хафа.
- `scripts/make_chapter3_tables.R` -- генерирует `.tex`-таблицы главы 3 и сохраняет сырые выборки в `chapters/chapter3_experiments/tables/data/`.
- `scripts/make_chapter3_frequency_figure.R` -- дополнительный скрипт для эксперимента по локализации максимума в строке.

По умолчанию `make_chapter3_tables.R` использует `n=100` реализаций и записывает только те `.tex`-таблицы,
которые подключены в главе 3. Сырые выборки и сводная таблица сохраняются в CSV. Если нужно вывести все
комбинации уровней шума и порогов, перед запуском можно задать `$env:VKR_TABLE_SCOPE='all'`.

Для быстрой проверки можно временно задать меньшее число повторов:

```powershell
$env:VKR_N_REP='10'
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' scripts\make_chapter3_tables.R
Remove-Item Env:\VKR_N_REP
```

В демонстрационных рисунках для медианного фильтра используется модифицированный line-preserving вариант, чтобы тонкая прямая была видна на иллюстрации. В численных таблицах главы 3 используется обычный медианный фильтр с окном `3x3`.
