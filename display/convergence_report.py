#!/usr/bin/env python3
"""Насколько карта меняется от чекпойнта к чекпойнту и где перестаёт.

На вход идёт свёртка из aggregate_convergence.py. Для каждой панели (тип связи,
eps) и каждой величины карта на чекпойнте сравнивается с картой на предыдущем
чекпойнте: берётся медиана и 99-й процентиль модуля разности по ячейкам, а для
доли захвата — ещё и число ячеек, сменивших ответ на вопрос «захвачено».

Медиана показывает, сдвинулась ли карта целиком, 99-й процентиль — остались ли
отдельные ячейки, которые всё ещё едут. Максимум не берётся: одна ячейка на
границе языка может дрожать сколь угодно долго, и по ней судить о карте нельзя.

Сходимость объявляется не по одному пересечению порога, а по плато: величина
считается сошедшейся с того чекпойнта, начиная с которого изменение остаётся
ниже допуска на всех последующих. Допуски: для наклонов — theta/10, потому что
точнее наклон и не измеряется; для долей и для P, L, A — 1% размаха.

  ./convergence_report.py <свёртка.csv> [выход.csv]
"""
import os
import sys

import duckdb

THETA = 2 * 3.141592653589793 / 1360.0

# Величина -> допуск, ниже которого изменение считается неразличимым
METRICS = {
    "cap_any": 0.01,      # доля ячеек, 1%
    "P_med": 0.01,
    "P_spread": 0.01,
    "L_med": 0.01 * 3.6,  # L нормирована на N^2*A^2 = 3.6
    "A_med": 0.01,
    "s01_med": THETA / 10,
    "s02_med": THETA / 10,
    "s12_med": THETA / 10,
}

DIFF = """
WITH ladder AS (
  SELECT *, lag({m}) OVER (
      PARTITION BY coupling_type, eps, delta1, delta2 ORDER BY t_trans
    ) AS prev
  FROM agg
)
SELECT coupling_type, eps, t_trans,
       median(abs({m} - prev)) AS d_med,
       quantile_cont(abs({m} - prev), 0.99) AS d_p99
FROM ladder WHERE prev IS NOT NULL
GROUP BY coupling_type, eps, t_trans
ORDER BY coupling_type, eps, t_trans
"""


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else None

    con = duckdb.connect()
    con.execute(f"CREATE TABLE agg AS SELECT * FROM read_csv('{src}')")
    checkpoints = [r[0] for r in
                   con.execute("SELECT DISTINCT t_trans FROM agg "
                               "ORDER BY t_trans").fetchall()]
    print("чекпойнты:", ", ".join(f"{t:g}" for t in checkpoints))

    rows = []
    for metric, tol in METRICS.items():
        d = con.execute(DIFF.format(m=metric)).fetchall()
        for ct, eps, t, d_med, d_p99 in d:
            rows.append((metric, tol, int(ct), eps, t, d_med, d_p99,
                         d_p99 <= tol))

    # Панель считается сошедшейся с того чекпойнта, начиная с которого допуск
    # держится и дальше: одиночное пересечение порога сходимостью не является.
    print(f"\n{'величина':<10}{'связь':>6}{'eps':>7}  сошлось с t_trans")
    verdict = []
    for metric in METRICS:
        for ct in sorted({r[2] for r in rows if r[0] == metric}):
            for eps in sorted({r[3] for r in rows if r[0] == metric
                               and r[2] == ct}):
                seq = [r for r in rows if r[0] == metric and r[2] == ct
                       and r[3] == eps]
                seq.sort(key=lambda r: r[4])
                first = None
                for i, r in enumerate(seq):
                    if all(s[7] for s in seq[i:]):
                        first = r[4]
                        break
                verdict.append((metric, ct, eps, first))
                mark = f"{first:g}" if first is not None else "не сошлось"
                print(f"{metric:<10}{ct:>6}{eps:>7g}  {mark}")

    if dst:
        con.execute("CREATE TABLE diffs (metric VARCHAR, tol DOUBLE, "
                    "coupling_type INT, eps DOUBLE, t_trans DOUBLE, "
                    "d_med DOUBLE, d_p99 DOUBLE, within BOOLEAN)")
        con.executemany("INSERT INTO diffs VALUES (?,?,?,?,?,?,?,?)", rows)
        con.execute(f"COPY diffs TO '{dst}' (FORMAT CSV, HEADER)")
        print(f"\nразности по чекпойнтам -> {dst} "
              f"({os.path.getsize(dst) / 2**10:.0f} КБ)")


if __name__ == "__main__":
    main()
