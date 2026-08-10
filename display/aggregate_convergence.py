#!/usr/bin/env python3
"""Поячеечная свёртка лестницы чекпойнтов (см. config_converge.yaml).

Выгрузка лестницы — сотни миллионов строк по 100 начальных условий на ячейку,
и целиком в parquet она не помещается. Но для вопроса «с какого t_trans карта
перестаёт меняться» сами траектории не нужны: нужна статистика ансамбля в
ячейке. Здесь из выгрузки получается по одной строке на
(тип связи, eps, t_trans, delta1, delta2) — это в сто раз компактнее и читается
в ноутбуке целиком.

Что считается и почему именно это:

* cap01/cap02/cap12 — доля начальных условий, где пара захвачена, то есть
  наклон разности фаз ниже порога theta = 2*pi/window. cap_any — доля, где
  захвачена хотя бы одна пара; именно так «в захвате» определено в статье.
* L, A, P — медиана по ансамблю, а не среднее: в мультистабильной ячейке
  распределение многомодальное, и среднее размазывается между режимами.
* P_spread = max - min по ансамблю — мера мультистабильности: если все
  начальные условия сходятся к одному режиму, разброс нулевой.

Через DuckDB, а не через поток на питоне: связка zstdcat | grep | парсер
разбирает CSV в один поток и упирается в него на 22 гигабайтах, тогда как
DuckDB распаковывает, разбирает и группирует на всех ядрах сразу и не держит
выгрузку в памяти целиком.

  ./aggregate_convergence.py <выход.csv> <кусок.csv.zst> [кусок2.csv.zst ...]
"""
import math
import os
import sys

import duckdb

# Порог захвата. Он задан длиной окна наблюдения: наклон меньше 2*pi/window
# неотличим от нуля, потому что за окно разность фаз не успевает набрать и
# одного оборота. Окно во всех чекпойнтах одно и то же, поэтому порог общий.
WINDOW = 1360.0
THETA = 2 * math.pi / WINDOW

# Ключ ячейки. Расстройки сравниваются как числа с плавающей точкой, но это
# безопасно: во всей выгрузке они получены одним и тем же выражением
# d_min + i*d_step, поэтому совпадают побитово.
KEYS = "coupling_type, eps, t_trans, delta1, delta2"

QUERY = """
COPY (
  SELECT
    {keys},
    count(*) AS n,
    avg(CASE WHEN s01 < {th} THEN 1.0 ELSE 0.0 END) AS cap01,
    avg(CASE WHEN s02 < {th} THEN 1.0 ELSE 0.0 END) AS cap02,
    avg(CASE WHEN s12 < {th} THEN 1.0 ELSE 0.0 END) AS cap12,
    avg(CASE WHEN least(s01, s02, s12) < {th} THEN 1.0 ELSE 0.0 END) AS cap_any,
    median(L) AS L_med,
    median(A) AS A_med,
    median(P) AS P_med,
    min(P) AS P_min,
    max(P) AS P_max,
    max(P) - min(P) AS P_spread,
    median(s01) AS s01_med,
    median(s02) AS s02_med,
    median(s12) AS s12_med
  FROM read_csv({src}, comment='#')
  GROUP BY {keys}
  ORDER BY {keys}
) TO '{dst}' (FORMAT CSV, HEADER)
"""


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    dst, sources = sys.argv[1], sys.argv[2:]
    print(f"порог захвата theta = {THETA:.4e} (окно {WINDOW:g})", flush=True)

    con = duckdb.connect()
    # Временные файлы — рядом с выгрузкой: в корне счётной машины меньше двух
    # гигабайт, а сортировка под медиану может вылиться на диск.
    con.execute(f"SET temp_directory='{os.path.dirname(os.path.abspath(dst))}'")
    src = "[" + ", ".join(f"'{s}'" for s in sources) + "]"
    con.execute(QUERY.format(keys=KEYS, th=repr(THETA), src=src, dst=dst))

    n = con.execute(f"SELECT count(*) FROM read_csv('{dst}')").fetchone()[0]
    size = os.path.getsize(dst) / 2 ** 20
    print(f"итого {n:,} строк, {size:.1f} МБ -> {dst}")


if __name__ == "__main__":
    main()
