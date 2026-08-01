#!/usr/bin/env python3
"""CSV (в том числе сжатый zstd) -> parquet порциями.

Память не растёт с размером входа: файл читается потоком и пишется
построчными группами, поэтому конвертируются и стогигабайтные выгрузки.

Шапка provenance (строки, начинающиеся на #) вырезается до парсера:
pyarrow, в отличие от pandas, комментарии в CSV не понимает и падает на
первой же такой строке.

Принимает несколько источников — это нужно для выгрузок, посчитанных по
кускам (см. simulation/experiments/multistability/run_chunked.py): все они
собираются в один parquet. Схема берётся из первого куска, остальные должны
ей соответствовать.

  ./to_parquet.py <выход.parquet> <вход.csv[.zst]> [вход2.csv[.zst] ...]
"""
import os
import subprocess
import sys

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

BLOCK = 256 << 20  # порция чтения, байт


def rows_of(path):
    """Поток батчей из одного файла; zstd распаковывается на лету."""
    cat = "zstdcat" if path.endswith(".zst") else "cat"
    pipe = subprocess.Popen(f"{cat} {path!r} | grep -v '^#'",
                            shell=True, stdout=subprocess.PIPE)
    reader = pacsv.open_csv(pipe.stdout,
                            read_options=pacsv.ReadOptions(block_size=BLOCK))
    try:
        for batch in reader:
            yield batch
    finally:
        pipe.wait()


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    dst, sources = sys.argv[1], sys.argv[2:]

    writer, total = None, 0
    for i, src in enumerate(sources, 1):
        n = 0
        for batch in rows_of(src):
            if writer is None:
                writer = pq.ParquetWriter(dst, batch.schema, compression="zstd")
            writer.write_table(pa.Table.from_batches([batch]))
            n += batch.num_rows
        total += n
        print(f"[{i}/{len(sources)}] {os.path.basename(src)}: {n:,} строк", flush=True)
    if writer:
        writer.close()
    size = os.path.getsize(dst) / 2 ** 30 if os.path.exists(dst) else 0
    print(f"итого {total:,} строк -> {dst} ({size:.2f} ГБ)")


if __name__ == "__main__":
    main()
