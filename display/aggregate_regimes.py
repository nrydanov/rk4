#!/usr/bin/env python3
"""Поячеечная статистика режимов синхронизации.

Дополняет aggregate_convergence.py: там доли захвата по парам по отдельности, а
режим ячейки — это комбинация всех трёх пар сразу, и из отдельных долей она не
восстанавливается. Здесь для каждой ячейки считается распределение начальных
условий по режимам и из него — доминирующий режим, его доля и мультистабильность.

Код режима собирается битами из захвата пар: 4*(0-1) + 2*(0-2) + 1*(1-2).
Динамически различимых режимов пять, а не восемь: наклоны связаны тождеством
s01 = s02 - s12, поэтому если две пары захвачены, то захвачена и третья. Порог
это свойство ломает — |s12| <= |s01| + |s02| < 2*theta, и третий наклон может
попасть в [theta, 2*theta), — так что коды с двумя битами всё же появляются.
Все они сводятся в класс 9 «полный частотный захват», и обязательно до взятия
моды, иначе один режим дробится на несколько кодов и мода уходит к другому.

  ./aggregate_regimes.py <выход.csv> <источник.parquet|кусок.csv.zst ...>
"""
import math
import os
import sys

import duckdb

WINDOW = 1360.0
THETA = 2 * math.pi / WINDOW

KEYS = "coupling_type, eps, delta1, delta2"

QUERY = """
WITH b AS (
  SELECT {keys},
         (CASE WHEN s01 < {th} THEN 1 ELSE 0 END) AS b01,
         (CASE WHEN s02 < {th} THEN 1 ELSE 0 END) AS b02,
         (CASE WHEN s12 < {th} THEN 1 ELSE 0 END) AS b12
  FROM {src} {where}
),
r AS (
  -- Два и более захваченных бита сводятся в один класс 9 «полный частотный
  -- захват»: комбинация «ровно две пары из трёх» динамически невозможна
  -- (s01 = s02 - s12 тождественно), а порог её порождает, когда третий наклон
  -- попадает в [theta, 2*theta). Свести надо ДО взятия моды, иначе один режим
  -- искусственно дробится на несколько кодов и мода уходит к другому классу.
  SELECT {keys},
         CASE WHEN b01 + b02 + b12 >= 2 THEN 9
              ELSE 4 * b01 + 2 * b02 + b12 END AS code
  FROM b
),
c AS (
  SELECT {keys}, code, count(*) AS n FROM r GROUP BY {keys}, code
)
SELECT {keys},
       sum(n) AS n_total,
       arg_max(code, n) AS code_dom,
       max(n) AS n_dom,
       count(*) AS n_codes,
       sum(CASE WHEN code = 9 THEN n ELSE 0 END) AS n_full,
       sum(CASE WHEN code = 0 THEN n ELSE 0 END) AS n_none
FROM c GROUP BY {keys}
"""


def parts_of(path):
    """Части, на которые бьётся проход: у parquet — пары (тип связи, eps)."""
    if path.endswith(".parquet"):
        con = duckdb.connect()
        pairs = con.execute(
            f"SELECT DISTINCT coupling_type, eps FROM read_parquet('{path}') "
            "ORDER BY coupling_type, eps").fetchall()
        return [(f"read_parquet('{path}')",
                 f"WHERE coupling_type = {ct} AND eps = {eps!r}",
                 f"тип {ct}, eps {eps:g}") for ct, eps in pairs]
    return [(f"read_csv('{path}', comment='#')", "", os.path.basename(path))]


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    dst, sources = sys.argv[1], sys.argv[2:]
    print(f"порог захвата theta = {THETA:.4e}", flush=True)

    con = duckdb.connect()
    con.execute(f"SET temp_directory='{os.path.dirname(os.path.abspath(dst))}'")
    parts = [p for s in sources for p in parts_of(s)]
    for i, (src, where, label) in enumerate(parts, 1):
        select = QUERY.format(keys=KEYS, th=repr(THETA), src=src, where=where)
        con.execute(f"CREATE TABLE reg AS {select}" if i == 1
                    else f"INSERT INTO reg {select}")
        done = con.execute("SELECT count(*) FROM reg").fetchone()[0]
        print(f"[{i}/{len(parts)}] {label}: всего {done:,} ячеек", flush=True)

    # Мультистабильность — доля начальных условий вне доминирующего режима,
    # ровно как определено в статье
    con.execute(f"""COPY (
      SELECT {KEYS}, n_total, code_dom, n_codes,
             1.0 - n_dom::DOUBLE / n_total AS multistab,
             n_full::DOUBLE / n_total AS full_frac,
             n_none::DOUBLE / n_total AS none_frac
      FROM reg ORDER BY {KEYS}) TO '{dst}' (FORMAT CSV, HEADER)""")
    n = con.execute("SELECT count(*) FROM reg").fetchone()[0]
    print(f"итого {n:,} строк, {os.path.getsize(dst)/2**20:.1f} МБ -> {dst}")


if __name__ == "__main__":
    main()
