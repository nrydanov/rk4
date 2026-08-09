#!/usr/bin/env python3
"""Прогон по кускам с возобновлением.

Один кусок = одна группа связи (инерционная или диссипативная) при одном
значении eps. Начальные условия засеиваются как seed_seq{seed, i1, i2} и от
eps с типом связи не зависят, поэтому нарезка даёт побитово тот же результат,
что и цельный прогон (проверено сравнением).

Готовый кусок пишется во временный .part и переименовывается только после
успешного завершения, поэтому оборванный кусок не будет принят за готовый.
Повторный запуск пропускает готовые куски и продолжает с первого недостающего.

  ./run_chunked.py <config.yaml> <выходной каталог> [путь к бинарнику]
"""
import os, subprocess, sys, time
import yaml

cfg_path = sys.argv[1]
out_dir = sys.argv[2]
binary = sys.argv[3] if len(sys.argv) > 3 else \
    os.path.join(os.path.dirname(__file__), "..", "..", "build", "vdp_multistability")

with open(cfg_path) as f:
    base = yaml.safe_load(f)
sim = base["sim"]
eps_i = sim.get("epsilons_inertial", sim.get("epsilons", []))
eps_d = sim.get("epsilons_dissipative", sim.get("epsilons", []))

chunks = [("inert", e) for e in eps_i] + [("diss", e) for e in eps_d]
os.makedirs(out_dir, exist_ok=True)
cfg_dir = os.path.join(out_dir, "configs")
os.makedirs(cfg_dir, exist_ok=True)

print(f"кусков всего: {len(chunks)}  ({len(eps_i)} инерционных + {len(eps_d)} диссипативных)")
done_before = sum(1 for g, e in chunks
                  if os.path.exists(os.path.join(out_dir, f"{g}-eps{e:g}.csv.zst")))
if done_before:
    print(f"уже готово: {done_before}, они будут пропущены")

t_all = time.time()
for k, (group, eps) in enumerate(chunks, 1):
    name = f"{group}-eps{eps:g}"
    final = os.path.join(out_dir, f"{name}.csv.zst")
    if os.path.exists(final):
        print(f"[{k}/{len(chunks)}] {name}: уже готов, пропуск")
        continue

    chunk = {**base, "sim": {**sim,
                             "epsilons_inertial": [eps] if group == "inert" else [],
                             "epsilons_dissipative": [eps] if group == "diss" else []}}
    chunk["sim"].pop("epsilons", None)
    cfg_chunk = os.path.join(cfg_dir, f"{name}.yaml")
    with open(cfg_chunk, "w") as f:
        yaml.safe_dump(chunk, f, allow_unicode=True, sort_keys=False)

    part = final + ".part"
    print(f"[{k}/{len(chunks)}] {name}: счёт...", flush=True)
    t0 = time.time()
    rc = subprocess.call(
        # -12, а не -19: на 63 ядрах -19 выдаёт 4.4 МБ/с и становится узким
        # местом (часы на прогон) ради выигрыша в 9% против -12.
        f'"{binary}" "{cfg_chunk}" -o /dev/stdout | zstd -T0 -12 -q -o "{part}" -f',
        shell=True)
    if rc != 0:
        print(f"[{k}/{len(chunks)}] {name}: ОШИБКА, код {rc}; .part оставлен для разбора")
        sys.exit(rc)
    os.replace(part, final)
    sz = os.path.getsize(final) / 2**30
    print(f"[{k}/{len(chunks)}] {name}: готов за {(time.time()-t0)/60:.1f} мин, {sz:.2f} ГБ",
          flush=True)

print(f"все куски готовы за {(time.time()-t_all)/3600:.2f} ч -> {out_dir}")
