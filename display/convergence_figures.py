#!/usr/bin/env python3
"""Страница с картами по чекпойнтам и кривыми сходимости.

На вход идут свёртка из aggregate_convergence.py и разности из
convergence_report.py. На выходе — самодостаточный HTML: карты вкладываются
попиксельно (151x151, ровно сетка расчёта) и масштабируются браузером без
интерполяции, поэтому страница остаётся лёгкой и ничего не приукрашивает.

  ./convergence_figures.py <свёртка.csv> <разности.csv> <выход.html>
"""
import base64
import io
import sys

import duckdb
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap

# Последовательная шкала — один тон от светлого к тёмному (см. руководство по
# графике). Радуги быть не должно: величина непрерывная, а не категориальная.
BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#256abf",
        "#184f95", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("blue", BLUE)
# Шаги той же шкалы под линии: eps — величина упорядоченная, не набор меток
LINE_STEPS = ["#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#256abf",
              "#184f95", "#0d366b"]

NAMES = {0: "Инерционная", 1: "Инерционная норм.",
         2: "Диссипативная", 3: "Диссипативная норм."}


def grid(con, ct, eps, t, col="L_med"):
    """Карта величины по сетке расстроек для одной панели и чекпойнта."""
    rows = con.execute(
        f"SELECT delta1, delta2, {col} FROM agg WHERE coupling_type=? "
        "AND eps=? AND t_trans=? ORDER BY delta1, delta2",
        [ct, eps, t]).fetchall()
    d1 = sorted({r[0] for r in rows})
    d2 = sorted({r[1] for r in rows})
    m = np.full((len(d1), len(d2)), np.nan)
    i1 = {v: i for i, v in enumerate(d1)}
    i2 = {v: i for i, v in enumerate(d2)}
    for a, b, v in rows:
        m[i1[a], i2[b]] = v
    return m


def png(arr, vmin, vmax, cmap=CMAP):
    """Массив как PNG в data-URI, один элемент сетки — один пиксель."""
    buf = io.BytesIO()
    mpimg.imsave(buf, arr.T[::-1], vmin=vmin, vmax=vmax, cmap=cmap,
                 format="png")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def filmstrip(con, ct, eps, checkpoints, col="L_med"):
    """Ряд карт по чекпойнтам плюс ряд разностей с предыдущим."""
    maps = [grid(con, ct, eps, t, col) for t in checkpoints]
    vmin = min(np.nanmin(m) for m in maps)
    vmax = max(np.nanmax(m) for m in maps)
    diffs = [np.abs(maps[i] - maps[i - 1]) for i in range(1, len(maps))]
    dmax = max(np.nanquantile(d, 0.999) for d in diffs)
    return ([png(m, vmin, vmax) for m in maps],
            [png(d, 0, dmax) for d in diffs],
            vmin, vmax, dmax)


def svg_curves(con, metric, tol):
    """Малые кратные: убывание разности карт по лестнице, по типам связи."""
    W, H, PAD_L, PAD_B, PAD_T = 250, 165, 46, 28, 22
    checkpoints = [r[0] for r in con.execute(
        "SELECT DISTINCT t_trans FROM d WHERE t_trans>120 ORDER BY t_trans"
    ).fetchall()]
    lo, hi = 1e-6, 1e-1
    def y(v):
        v = max(v, lo)
        f = (np.log10(v) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
        return PAD_T + (1 - f) * (H - PAD_T - PAD_B)
    def x(i):
        return PAD_L + i * (W - PAD_L - 12) / (len(checkpoints) - 1)

    out = []
    for ct in (0, 1, 2, 3):
        eps_list = [r[0] for r in con.execute(
            "SELECT DISTINCT eps FROM d WHERE coupling_type=? ORDER BY eps",
            [ct]).fetchall()]
        s = [f'<svg viewBox="0 0 {W} {H}" role="img" '
             f'aria-label="{NAMES[ct]}">']
        s.append(f'<text x="{PAD_L}" y="13" class="ft">{NAMES[ct]}</text>')
        for p in (-6, -5, -4, -3, -2, -1):
            yy = y(10.0 ** p)
            s.append(f'<line x1="{PAD_L}" y1="{yy:.1f}" x2="{W-12}" '
                     f'y2="{yy:.1f}" class="grid"/>')
            s.append(f'<text x="{PAD_L-6}" y="{yy+3:.1f}" class="ax" '
                     f'text-anchor="end">1e{p}</text>')
        yt = y(tol)
        s.append(f'<line x1="{PAD_L}" y1="{yt:.1f}" x2="{W-12}" y2="{yt:.1f}" '
                 f'class="tol"/>')
        for i, t in enumerate(checkpoints):
            s.append(f'<text x="{x(i):.1f}" y="{H-10}" class="ax" '
                     f'text-anchor="middle">{t:g}</text>')
        for k, eps in enumerate(eps_list):
            pts = []
            for i, t in enumerate(checkpoints):
                v = con.execute(
                    "SELECT d_med FROM d WHERE metric=? AND coupling_type=? "
                    "AND eps=? AND t_trans=?", [metric, ct, eps, t]).fetchone()
                if v:
                    pts.append((x(i), y(v[0]), t, v[0]))
            color = LINE_STEPS[k % len(LINE_STEPS)]
            d = " ".join(f"{'M' if j == 0 else 'L'}{px:.1f},{py:.1f}"
                         for j, (px, py, _, _) in enumerate(pts))
            s.append(f'<path d="{d}" fill="none" stroke="{color}" '
                     f'stroke-width="2" stroke-linejoin="round"/>')
            for px, py, t, v in pts:
                s.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="4" '
                         f'fill="{color}" stroke="var(--surface-1)" '
                         f'stroke-width="2"><title>eps={eps:g}, '
                         f't_trans={t:g}: {v:.2g}</title></circle>')
            if pts:
                s.append(f'<text x="{pts[-1][0]+6:.1f}" y="{pts[-1][1]+3:.1f}" '
                         f'class="lbl">{eps:g}</text>')
        s.append("</svg>")
        out.append("".join(s))
    return out


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    agg_path, diff_path, dst = sys.argv[1], sys.argv[2], sys.argv[3]
    con = duckdb.connect()
    con.execute(f"CREATE TABLE agg AS SELECT * FROM read_csv('{agg_path}')")
    con.execute(f"CREATE TABLE d AS SELECT * FROM read_csv('{diff_path}')")
    checkpoints = [r[0] for r in con.execute(
        "SELECT DISTINCT t_trans FROM agg ORDER BY t_trans").fetchall()]
    print("чекпойнты:", checkpoints, flush=True)

    panels = [(2, 0.005, "самая медленная панель"),
              (0, 0.30, "самая быстрая панель")]
    strips = []
    for ct, eps, note in panels:
        maps, diffs, vmin, vmax, dmax = filmstrip(con, ct, eps, checkpoints)
        strips.append(dict(ct=ct, eps=eps, note=note, maps=maps, diffs=diffs,
                           vmin=vmin, vmax=vmax, dmax=dmax))
        print(f"карты готовы: {NAMES[ct]} eps={eps:g}", flush=True)

    curves = svg_curves(con, "L_med", 0.036)
    html = build_html(checkpoints, strips, curves, con)
    with open(dst, "w") as f:
        f.write(html)
    print(f"страница -> {dst} ({len(html)/2**20:.1f} МБ)")


CSS = """
:root {
  --surface-1:#fcfcfb; --plane:#f9f9f7; --ink:#0b0b0b; --ink-2:#52514e;
  --muted:#898781; --grid:#e1e0d9; --rule:#c3c2b7; --accent:#2a78d6;
  --warn:#d03b3b;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --surface-1:#1a1a19; --plane:#0d0d0d; --ink:#fff; --ink-2:#c3c2b7;
    --muted:#898781; --grid:#2c2c2a; --rule:#383835; --accent:#3987e5;
    --warn:#e66767;
  }
}
:root[data-theme="dark"] {
  --surface-1:#1a1a19; --plane:#0d0d0d; --ink:#fff; --ink-2:#c3c2b7;
  --muted:#898781; --grid:#2c2c2a; --rule:#383835; --accent:#3987e5;
  --warn:#e66767;
}
body { background:var(--plane); color:var(--ink); margin:0;
  font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }
.wrap { max-width:1080px; margin:0 auto; padding:40px 24px 80px; }
h1 { font-size:30px; line-height:1.25; margin:0 0 8px; letter-spacing:-.02em; }
h2 { font-size:20px; margin:48px 0 6px; letter-spacing:-.01em; }
p { color:var(--ink-2); margin:8px 0 0; max-width:68ch; }
.lede { font-size:18px; color:var(--ink-2); }
.card { background:var(--surface-1); border:1px solid var(--grid);
  border-radius:12px; padding:20px; margin-top:16px; }
.answer { border-left:3px solid var(--accent); }
.answer b { color:var(--ink); }
.small { font-size:14px; color:var(--muted); }
.facets { display:grid; grid-template-columns:repeat(2,minmax(0,1fr));
  gap:8px; }
@media (max-width:720px) { .facets { grid-template-columns:1fr; } }
svg { width:100%; height:auto; display:block; }
.ft { fill:var(--ink); font-size:11px; font-weight:600; }
.ax { fill:var(--muted); font-size:9px; }
.lbl { fill:var(--ink-2); font-size:9px; }
.grid { stroke:var(--grid); stroke-width:1; }
.tol { stroke:var(--warn); stroke-width:1.5; stroke-dasharray:4 3; }
.strip { overflow-x:auto; }
.row { display:flex; gap:10px; min-width:min-content; }
.cell { flex:0 0 150px; }
.cell img { width:150px; height:150px; image-rendering:pixelated;
  border:1px solid var(--grid); border-radius:4px; display:block; }
.cell .cap { font-size:12px; color:var(--ink-2); margin-top:6px;
  text-align:center; }
.ramp { height:8px; border-radius:4px; margin:10px 0 4px;
  background:linear-gradient(90deg,#cde2fb,#9ec5f4,#6da7ec,#3987e5,#2a78d6,
  #256abf,#184f95,#0d366b); }
table { border-collapse:collapse; width:100%; font-size:14px; margin-top:8px; }
th, td { text-align:right; padding:6px 10px;
  border-bottom:1px solid var(--grid); }
th:first-child, td:first-child { text-align:left; }
th { color:var(--muted); font-weight:600; font-size:12px;
  text-transform:uppercase; letter-spacing:.04em; }
td.ok { color:var(--accent); } td.no { color:var(--warn); }
"""


def build_html(checkpoints, strips, curves, con):
    cps = "".join(f"<th>{t:g}</th>" for t in checkpoints)

    # Таблица: медиана разности соседних карт по величинам
    metrics = [r[0] for r in con.execute(
        "SELECT DISTINCT metric FROM d ORDER BY metric").fetchall()]

    def table(col):
        trs = []
        for m in metrics:
            tol = con.execute("SELECT max(tol) FROM d WHERE metric=?",
                              [m]).fetchone()[0]
            cells = ""
            for t in checkpoints[1:]:
                v = con.execute(f"SELECT median({col}) FROM d WHERE metric=? "
                                "AND t_trans=?", [m, t]).fetchone()[0]
                cls = "ok" if v is not None and v <= tol else "no"
                cells += (f'<td class="{cls}">{v:.2g}</td>' if v
                          else "<td>—</td>")
            trs.append(f"<tr><td>{m}</td><td>{tol:.2g}</td>{cells}</tr>")
        head = ("<thead><tr><th>величина</th><th>допуск</th>"
                + "".join(f"<th>{t:g}</th>" for t in checkpoints[1:])
                + "</tr></thead>")
        return f"<table>{head}<tbody>{''.join(trs)}</tbody></table>"

    tbl_med, tbl_p99 = table("d_med"), table("d_p99")

    strip_html = ""
    for s in strips:
        maps = "".join(
            f'<div class="cell"><img src="{u}" alt="карта при t_trans={t:g}">'
            f'<div class="cap">{t:g}</div></div>'
            for u, t in zip(s["maps"], checkpoints))
        diffs = "".join(
            f'<div class="cell"><img src="{u}" alt="разность карт">'
            f'<div class="cap">{a:g}&thinsp;&rarr;&thinsp;{b:g}</div></div>'
            for u, a, b in zip(s["diffs"], checkpoints, checkpoints[1:]))
        strip_html += f"""
<h2>{NAMES[s['ct']]}, eps = {s['eps']:g} — {s['note']}</h2>
<div class="card">
  <p class="small">Карта когерентности L по расстройкам, медиана по 100
  начальным условиям. Каждый пиксель — ячейка сетки, δ₁ по горизонтали,
  δ₂ по вертикали, диапазон ±0.15. Шкала общая для всего ряда:
  {s['vmin']:.2f} … {s['vmax']:.2f}.</p>
  <div class="ramp"></div>
  <div class="strip"><div class="row">{maps}</div></div>
  <p class="small" style="margin-top:18px">Модуль разности с предыдущим
  чекпойнтом, шкала 0 … {s['dmax']:.3f}. Если бы карта сходилась, ряд бы
  темнел справа налево и гас.</p>
  <div class="strip"><div class="row">{diffs}</div></div>
</div>"""

    facets = "".join(f"<div>{c}</div>" for c in curves)

    return f"""<title>Сходимость карт по времени переходного процесса</title>
<style>{CSS}</style>
<div class="wrap">
<h1>Сходимость карт по времени переходного процесса</h1>
<p class="lede">Лестница из шести значений t_trans при неизменном окне
наблюдения 1360. Сетка 151×151, 100 начальных условий на ячейку, 12 значений
силы связи, четыре типа связи — 328 млн траекторий, свёрнутых в 3.3 млн ячеек.</p>

<div class="card answer">
<p><b>Основная часть карты перестаёт меняться при t_trans ≈ 640.</b> До этого
разность соседних карт падает (медиана по ячейкам для L: 0.0031 → 0.0021 →
0.0017), после — выходит на полку и дальше не убывает, сколько транзиент ни
удлиняй.</p>
<p><b>Полка не равна нулю, и удлинением транзиента она не убирается.</b>
Остаточные 0.0017 держатся одинаково на переходах 480→640, 640→860 и
860→1720, то есть не зависят ни от длины транзиента, ни от того, насколько
разнесены сравниваемые окна. Значит это не недосчитанный переходный
процесс: его вклад убывал бы с ростом t_trans. Чем именно задана полка —
конечной длиной окна наблюдения или медленным блужданием самих траекторий —
по этим данным не различить, нужен отдельный опыт (сравнить два
непересекающихся окна при одном и том же t_trans).</p>
<p><b>Хвост ячеек не сходится вовсе.</b> 99-й процентиль разности не убывает по
лестнице ни для одной величины. Это ячейки на границах языка захвата и в
мультистабильных областях, где ансамбль распределён по нескольким режимам;
там и не должно быть сходимости отдельной ячейки — сходится статистика.</p>
</div>

<h2>Как убывает разность соседних карт</h2>
<p>Медиана модуля разности по ячейкам для карты когерентности L, по типам
связи; линия — значение eps, светлее значит слабее связь. Пунктир — допуск
в 1% размаха L. Шкала логарифмическая.</p>
<div class="card"><div class="facets">{facets}</div></div>
<p class="small">Убывание кончается на третьем-четвёртом чекпойнте у всех
четырёх типов связи. Слабая связь (светлые линии) стоит выше сильной: время
фазовой релаксации обратно пропорционально силе связи, и при eps = 0.005 оно
на два порядка больше, чем при eps = 0.3.</p>

{strip_html}

<h2>Все величины разом</h2>
<p>Разность соседних карт, сведённая по 24 панелям. Синим — уложилось в
допуск, красным — нет.</p>
<div class="card">
<p class="small"><b>Медиана по ячейкам.</b> Уложилась в допуск везде и сразу,
с самого первого чекпойнта: типичная ячейка карты не меняется вовсе. Ниже
видно, что убывать эти числа перестают после t_trans ≈ 640.</p>
{tbl_med}
<p class="small" style="margin-top:22px"><b>99-й процентиль по ячейкам.</b>
Тот же расчёт по хвосту. Здесь допуск не выдержан нигде и ни на одном
чекпойнте, и по лестнице эти числа не убывают. Именно из-за этого хвоста
формальная проверка на плато объявляет несошедшимися 15–21 панель из 24.</p>
{tbl_p99}
</div>

<h2>Оговорка, без которой цифры читаются неверно</h2>
<div class="card">
<p>Окна соседних чекпойнтов перекрываются, и по-разному: [480, 1840] и
[640, 2000] делят 88% отсчётов, [640, 2000] и [860, 2220] — 84%, а последняя
пара, [860, 2220] и [1720, 3080], только 37%. Разность двух сильно
перекрытых окон занижена, поэтому столбцы таблицы не вполне сопоставимы
между собой.</p>
<p>Читается это так. Медианы по лестнице стоят на месте, хотя перекрытие
падает более чем вдвое, — если бы полку задавал шум конечного окна, на
последнем переходе она бы заметно выросла. А 99-й процентиль как раз растёт
(у доли захвата с 0.04 до 0.065), то есть хвост упрямых ячеек к перекрытию
чувствителен. Это ещё один довод, что основная часть карты и её хвост живут
по разным законам.</p>
<p class="small">Данные: коммит edb0ef1, конфиг config_converge.yaml, порог
захвата θ = 2π/1360 = 4.62e-3.</p>
</div>
</div>"""


if __name__ == "__main__":
    main()
