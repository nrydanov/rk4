#!/usr/bin/env python3
"""Поячеечная свёртка лестницы чекпойнтов (см. config_converge.yaml).

Выгрузка лестницы — сотни миллионов строк по 100 начальных условий на ячейку,
и целиком в parquet она не помещается. Но для вопроса «с какого t_trans карта
перестаёт меняться» сами траектории не нужны: нужна статистика ансамбля в
ячейке. Здесь из каждого куска получается по одной строке на
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

Только numpy и pyarrow: на счётной машине pandas нет, а ставить его туда ради
одной свёртки незачем. Группировка сделана сортировкой по целочисленному ключу,
а медиана — разворотом в матрицу (ячейки, начальные условия), поэтому свёртка
идёт векторно, без прохода по группам в питоне.

  ./aggregate_convergence.py <выход.csv> <кусок.csv.zst> [кусок2.csv.zst ...]
"""
import math
import os
import subprocess
import sys

import numpy as np
import pyarrow as pa
import pyarrow.csv as pacsv

BLOCK = 256 << 20  # порция чтения, байт

# Порог захвата. Он задан длиной окна наблюдения: наклон меньше 2*pi/window
# неотличим от нуля, потому что за окно разность фаз не успевает набрать и
# одного оборота. Окно во всех чекпойнтах одно и то же, поэтому порог общий.
WINDOW = 1360.0
THETA = 2 * math.pi / WINDOW

# Читаются только те колонки, которые нужны свёртке. Начальные и конечные
# состояния — это 12 колонок из 23, и без них поток вдвое легче.
COLS = ["delta1", "delta2", "eps", "coupling_type", "t_trans",
        "L", "A", "P", "s01", "s02", "s12"]

OUT_COLS = ["coupling_type", "eps", "t_trans", "delta1", "delta2", "n",
            "cap01", "cap02", "cap12", "cap_any",
            "L_med", "A_med", "P_med", "P_min", "P_max", "P_spread",
            "s01_med", "s02_med", "s12_med"]


def read_chunk(path):
    """Кусок целиком в таблицу; zstd распаковывается на лету."""
    cat = "zstdcat" if path.endswith(".zst") else "cat"
    pipe = subprocess.Popen(f"{cat} {path!r} | grep -v '^#'",
                            shell=True, stdout=subprocess.PIPE)
    reader = pacsv.open_csv(
        pipe.stdout,
        read_options=pacsv.ReadOptions(block_size=BLOCK),
        convert_options=pacsv.ConvertOptions(include_columns=COLS))
    try:
        table = pa.Table.from_batches(list(reader))
    finally:
        pipe.wait()
    return {c: table.column(c).to_numpy(zero_copy_only=False) for c in COLS}


def aggregate(col):
    """Свёртка ансамбля в одну строку на ячейку и чекпойнт."""
    # Расстройки, чекпойнты и типы связи — сетка, а не непрерывные величины,
    # поэтому группировка идёт по их индексам. Это заодно снимает вопрос о
    # сравнении чисел с плавающей точкой на равенство.
    axes = {}
    idx = {}
    for name in ("coupling_type", "eps", "t_trans", "delta1", "delta2"):
        axes[name] = np.unique(col[name])
        idx[name] = np.searchsorted(axes[name], col[name])

    key = idx["coupling_type"]
    for name in ("eps", "t_trans", "delta1", "delta2"):
        key = key * len(axes[name]) + idx[name]

    order = np.argsort(key, kind="stable")
    key_sorted = key[order]
    starts = np.flatnonzero(np.r_[True, key_sorted[1:] != key_sorted[:-1]])
    sizes = np.diff(np.r_[starts, len(key_sorted)])
    n_ens = sizes[0]
    if not np.all(sizes == n_ens):
        sys.exit(f"ансамбли разного размера: от {sizes.min()} до {sizes.max()}; "
                 "свёртка рассчитана на одинаковые")

    # Все группы одного размера, поэтому отсортированный столбец разворачивается
    # в матрицу (ячейки, начальные условия) и сводится вдоль второй оси разом.
    def by_cell(name):
        return col[name][order].reshape(-1, n_ens)

    out = {}
    keys_at = order[starts]
    for name in ("coupling_type", "eps", "t_trans", "delta1", "delta2"):
        out[name] = col[name][keys_at]
    out["n"] = np.full(len(starts), n_ens, dtype=np.int64)

    # Наклоны бинарник уже пишет по модулю
    caps = {}
    for pair in ("01", "02", "12"):
        caps[pair] = by_cell("s" + pair) < THETA
        out["cap" + pair] = caps[pair].mean(axis=1)
    out["cap_any"] = (caps["01"] | caps["02"] | caps["12"]).mean(axis=1)

    for name in ("L", "A", "P"):
        out[name + "_med"] = np.median(by_cell(name), axis=1)
    P = by_cell("P")
    out["P_min"] = P.min(axis=1)
    out["P_max"] = P.max(axis=1)
    out["P_spread"] = out["P_max"] - out["P_min"]
    for pair in ("01", "02", "12"):
        out["s" + pair + "_med"] = np.median(by_cell("s" + pair), axis=1)
    return out


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    dst, sources = sys.argv[1], sys.argv[2:]
    print(f"порог захвата theta = {THETA:.4e} (окно {WINDOW:g})", flush=True)

    with open(dst, "w") as out:
        out.write(",".join(OUT_COLS) + "\n")
        total = 0
        for i, src in enumerate(sources, 1):
            col = read_chunk(src)
            n_raw = len(col["P"])
            agg = aggregate(col)
            del col
            rows = np.column_stack([agg[c] for c in OUT_COLS])
            np.savetxt(out, rows, delimiter=",", fmt="%.9g")
            total += len(rows)
            print(f"[{i}/{len(sources)}] {os.path.basename(src)}: "
                  f"{n_raw:,} строк -> {len(rows):,} ячеек", flush=True)

    size = os.path.getsize(dst) / 2 ** 20
    print(f"итого {total:,} строк, {size:.1f} МБ -> {dst}")


if __name__ == "__main__":
    main()
