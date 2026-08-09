#include "config_array.hpp"
#include "coupling.hpp"
#include "forces.hpp"
#include "goals.hpp"
#include "phase_diff_slopes.hpp"
#include "provenance.hpp"
#include "vdp_ensemble.hpp"
#include <CLI11.hpp>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <optional>
#include <random>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

template <int N> struct Config {
  std::array<double, N> mus;             // параметр нелинейного затухания (μ)
  std::array<std::array<int, N>, N> adj; // матрица смежности
  double dt;                             // шаг интегрирования
  double d_min;                          // минимальная расстройка частоты
  double d_step;                         // шаг сетки расстроек
  // Времена переходного процесса, при которых снимаются метрики. Траектория
  // интегрируется один раз до самого позднего окна, а более ранние чекпойнты
  // читаются по дороге: прогон с меньшим t_trans — это та же самая траектория,
  // просто окно взято раньше. Отдельные запуски дали бы те же числа втридорога.
  std::vector<double> t_trans_list;
  // Длина окна наблюдения, общая для всех чекпойнтов. Фиксирована намеренно:
  // от неё зависит порог theta = 2*pi/window и шумовой пол оценки наклонов,
  // поэтому менять её вместе с t_trans значит менять две вещи разом и потерять
  // возможность судить о сходимости по транзиенту.
  double window;
  // Значения силы связи задаются отдельно для инерционных и диссипативных
  // вариантов: у диссипативной связи полуширина языка захвата примерно равна
  // eps, у инерционной вчетверо меньше, а порог гашения (mu - eps*k_n = 0)
  // существует только у диссипативной. Общий список одинаково плохо
  // обслуживал бы оба.
  std::vector<double> eps_inertial;
  std::vector<double> eps_dissipative;
  int grid_size;                         // число точек по каждой оси (delta)
  int n_ic;                              // число случайных начальных условий
  double ic_range;                       // начальные условия из [-ic_range, ic_range]
  uint32_t seed;                         // seed генератора начальных условий
};

struct Result {
  double delta1, delta2, eps; // параметры системы
  int coupling_type; // тип связи
  double t_trans;    // чекпойнт, с которого началось окно наблюдения
  double x0, y0, x1, y1, x2, y2;       // начальные условия
  double xf0, yf0, xf1, yf1, xf2, yf2; // конечное состояние
  double L, A, P; // значения рассматриваемых целевых функций
  double s01, s02, s12;
};

void print_progress(size_t done, size_t total,
                    std::chrono::steady_clock::time_point start,
                    const std::string &label) {
  double frac = (double)done / total;
  int bar_width = 35;
  int filled = (int)(frac * bar_width);
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                     std::chrono::steady_clock::now() - start)
                     .count();
  long eta = done > 0 ? (long)(elapsed * (double)(total - done) / done) : 0;

  std::cerr << "\r" << std::left << std::setw(16) << label << " [";
  for (int i = 0; i < bar_width; ++i)
    std::cerr << (i < filled ? '=' : (i == filled ? '>' : ' '));
  std::cerr << "] " << std::fixed << std::setprecision(1) << (frac * 100)
            << "% (" << done << "/" << total << ")";
  if (done > 0 && done < total)
    std::cerr << " ETA: " << eta / 60 << "m" << eta % 60 << "s";
  std::cerr << "   " << std::flush;
}

void write_result(std::ostream &out, const Result &r) {
  out << r.delta1 << "," << r.delta2 << "," << r.eps << "," << r.coupling_type
      << "," << r.t_trans
      << "," << r.x0 << "," << r.y0 << "," << r.x1 << "," << r.y1 << "," << r.x2
      << "," << r.y2 << "," << r.xf0 << "," << r.yf0 << "," << r.xf1 << ","
      << r.yf1 << "," << r.xf2 << "," << r.yf2 << "," << r.L << "," << r.A
      << "," << r.P << "," << r.s01 << "," << r.s02 << "," << r.s12 << "\n";
}

// Один полный прогон симуляции для конкретного типа связи
// Sink намеренно шаблонный параметр-функтор, чтобы код работал не только
// с std::ostream
template <int N, class CouplingFunc, class Sink>
void sweep(int coupling_type_id, const std::string &label, const Config<N> &cfg,
           const std::vector<double> &epsilons, Sink &&sink) {
  const size_t n_tasks =
      (size_t)cfg.grid_size * cfg.grid_size * epsilons.size() * cfg.n_ic;
  const int n_cp = (int)cfg.t_trans_list.size();
  const long n_win = std::lround(cfg.window / cfg.dt);
  // Окно чекпойнта c — это выборки с номерами (cp_start[c], cp_start[c]+n_win].
  // Нумерация с единицы, потому что метрика снимается после шага, а не до него.
  std::vector<long> cp_start(n_cp);
  long n_steps = 0, first_open = std::numeric_limits<long>::max();
  for (int c = 0; c < n_cp; ++c) {
    cp_start[c] = std::lround(cfg.t_trans_list[c] / cfg.dt);
    n_steps = std::max(n_steps, cp_start[c] + n_win);
    first_open = std::min(first_open, cp_start[c]);
  }
  std::vector<Result> results(n_tasks * n_cp);
  // Счетчик текущего прогресса
  std::atomic<size_t> progress{0};
  // Чтобы снизить падение производительности из-за постоянного вывода состояния,
  // мы делаем это дискретно раз в промежуток времени
  const size_t report_every = std::max((size_t)1, n_tasks / 50);
  auto start = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(static) collapse(3)
  for (size_t e = 0; e < epsilons.size(); ++e) {
    for (int i1 = 0; i1 < cfg.grid_size; ++i1) {
      for (int i2 = 0; i2 < cfg.grid_size; ++i2) {
        const size_t base = (e * (size_t)cfg.grid_size * cfg.grid_size +
                             i1 * cfg.grid_size + i2) *
                            cfg.n_ic;
        // Расстройка первого осциллятора
        const double delta1 = cfg.d_min + i1 * cfg.d_step;
        // Расстройка второго осциллятора
        const double delta2 = cfg.d_min + i2 * cfg.d_step;

        // TODO(nrydanov): это не работает для N отличных от трех, нужно придумать как генерализировать
        std::array<double, N> freqs{1.0, 1.0 + delta1, 1.0 + delta2};
        std::array<double, N> eps_coupling;
        eps_coupling.fill(epsilons[e]);

        // NOTE(nrydanov): я не уверен как лучше подбирать seed, чтобы
        // эксперимент был статистически корректным, сейчас для каждой
        // итерации используется свой, отдельным набор начальных условий
        std::seed_seq seq{cfg.seed, static_cast<uint32_t>(i1),
                          static_cast<uint32_t>(i2)};
        std::mt19937 rng(seq);
        std::uniform_real_distribution<double> dist(-cfg.ic_range,
                                                    cfg.ic_range);

        using Solver =
            vdp_ensemble::VdPEnsembleSolver<N, forces::NoopForce, CouplingFunc>;
        thread_local std::optional<Solver> solver;
        // По накопителю наклонов на чекпойнт: окна перекрываются, поэтому
        // одновременно открытых бывает несколько. Всё остальное — состояние
        // шага, оно общее.
        thread_local std::vector<phase_diff_slopes<double>> pds;
        thread_local std::optional<phase_t<double, int>> phase_goal;
        thread_local std::vector<double> L_acc, A_acc, P_acc;
        thread_local std::vector<std::array<double, 2 * N>> yf_cp;
        ampl_t<double, int> ampl_goal(2, N);
        if ((int)pds.size() != n_cp) {
          pds.clear();
          pds.reserve(n_cp);
          for (int c = 0; c < n_cp; ++c)
            pds.emplace_back(N, 2, 0, 1);
          L_acc.resize(n_cp);
          A_acc.resize(n_cp);
          P_acc.resize(n_cp);
          yf_cp.resize(n_cp);
        }

        std::array<double, 2 * N> y0;

        // Итетируемся по начальным условиям
        for (int ic = 0; ic < cfg.n_ic; ++ic) {
          // Генерируем начальные условия и кладем их в y0
          for (auto &v : y0)
            v = dist(rng);

          // Для экономии вычислений при первой итерации инициализируем solver,
          // а в дальнейшем переиспользуем его повторно
          if (!solver)
            solver.emplace(y0, freqs, cfg.mus, eps_coupling, cfg.adj,
                           forces::NoopForce{});
          else
            solver->reset(y0, freqs, eps_coupling);

          if (!phase_goal)
            phase_goal.emplace(2, N);
          for (int c = 0; c < n_cp; ++c) {
            pds[c].reset();
            L_acc[c] = A_acc[c] = P_acc[c] = 0.0;
          }

          for (long k = 1; k <= n_steps; ++k) {
            solver->step(cfg.dt);
            // До первого чекпойнта считать нечего — идём вхолостую
            if (k <= first_open)
              continue;
            const double t = k * cfg.dt;
            const auto &state = solver->getState();
            // Метрики шага не зависят от чекпойнта, поэтому считаются один раз
            // и раздаются всем открытым окнам
            std::array<double, N> raw_phase;
            for (int i = 0; i < N; ++i)
              raw_phase[i] = std::atan2(state[2 * i + 1], state[2 * i]);
            const double L_k = goals::coherence(state.data(), N);
            const double A_k = ampl_goal(state.data());
            const double P_k = phase_goal->from_raw_phases(raw_phase.data());
            for (int c = 0; c < n_cp; ++c) {
              if (k <= cp_start[c] || k > cp_start[c] + n_win)
                continue;
              pds[c].push_raw(t, raw_phase.data());
              L_acc[c] += L_k;
              A_acc[c] += A_k;
              P_acc[c] += P_k;
              // Конечное состояние — то, на котором окно закрылось
              if (k == cp_start[c] + n_win)
                yf_cp[c] = state;
            }
          }

          for (int c = 0; c < n_cp; ++c) {
            const auto &yf = yf_cp[c];
            results[(base + ic) * n_cp + c] = {delta1,
                              delta2,
                              epsilons[e],
                              coupling_type_id,
                              cfg.t_trans_list[c],
                              y0[0],
                              y0[1],
                              y0[2],
                              y0[3],
                              y0[4],
                              y0[5],
                              yf[0],
                              yf[1],
                              yf[2],
                              yf[3],
                              yf[4],
                              yf[5],
                              2.0 * L_acc[c] / n_win,
                              A_acc[c] / n_win,
                              P_acc[c] / n_win,
                              pds[c](0, 1),
                              pds[c](0, 2),
                              pds[c](1, 2)};
          }

          size_t done = progress.fetch_add(1, std::memory_order_relaxed) + 1;
          if (done % report_every == 0)
            print_progress(done, n_tasks, start, label);
        }
      }
    }
  }

  print_progress(n_tasks, n_tasks, start, label);
  std::cerr << "\n";
  for (const auto &r : results)
    sink(r);
}

template <int N>
int run(const YAML::Node &yaml, const std::string &config_path,
        const std::string &output_path) {
  const auto &s = yaml["sim"];
  double d_min = s["d_min"].as<double>();
  double d_max = s["d_max"].as<double>();
  double d_step = s["d_step"].as<double>();

  // Ключи t_trans_list и window задают лестницу чекпойнтов при неизменном окне
  // наблюдения; при их отсутствии работает прежняя пара t_transition и T, то
  // есть один чекпойнт — так продолжают считаться старые конфиги.
  std::vector<double> t_trans_list;
  double window;
  if (s["t_trans_list"]) {
    t_trans_list = s["t_trans_list"].as<std::vector<double>>();
    window = s["window"].as<double>();
  } else {
    t_trans_list = {s["t_transition"].as<double>()};
    window = s["T"].as<double>() - t_trans_list[0];
  }
  if (window <= 0.0) {
    std::cerr << "Observation window must be positive (got " << window << ")\n";
    return 1;
  }

  Config<N> cfg{
      .mus = to_array<N>(yaml["mus"].as<std::vector<double>>()),
      .adj = to_adj<N>(yaml["adj"].as<std::vector<std::vector<int>>>()),
      .dt = s["dt"].as<double>(),
      .d_min = d_min,
      .d_step = d_step,
      .t_trans_list = t_trans_list,
      .window = window,
      // Ключ epsilons задаёт общий список для всех типов связи (так устроены
      // прежние конфиги); epsilons_inertial и epsilons_dissipative задают их
      // раздельно и имеют приоритет.
      .eps_inertial = s[s["epsilons_inertial"] ? "epsilons_inertial" : "epsilons"]
                          .as<std::vector<double>>(),
      .eps_dissipative =
          s[s["epsilons_dissipative"] ? "epsilons_dissipative" : "epsilons"]
              .as<std::vector<double>>(),
      .grid_size =
          static_cast<int>(std::lround(std::abs(d_max - d_min) / d_step) + 1),
      .n_ic = s["n_ic"].as<int>(),
      .ic_range = s["ic_range"].as<double>(),
      .seed = s["seed"].as<uint32_t>(),
  };

  std::cerr << "Grid: " << cfg.grid_size << "x" << cfg.grid_size
            << "  epsilons: " << cfg.eps_inertial.size() << " inert / "
            << cfg.eps_dissipative.size() << " diss"
            << "  n_ic: " << cfg.n_ic
            << "  OpenMP threads: " << omp_get_max_threads() << "\n";
  std::cerr << "Window: " << cfg.window << "  checkpoints t_trans:";
  for (double t : cfg.t_trans_list)
    std::cerr << " " << t;
  std::cerr << "\n";

  std::ofstream out(output_path);
  if (!out.is_open()) {
    std::cerr << "Failed to open output file\n";
    return 1;
  }
  // Записываем параметры конфига в заголовок CSV для воспроизводимости
  provenance::write_header(out, "vdp_multistability", config_path, yaml);
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << "delta1,delta2,eps,coupling_type,t_trans,x0,y0,x1,y1,x2,y2,"
         "xf0,yf0,xf1,yf1,xf2,yf2,L,A,P,s01,s02,s12\n";

  auto sink = [&out](const Result &r) { write_result(out, r); };
  sweep<N, coupling::Inertial>(0, "Inertial", cfg, cfg.eps_inertial, sink);
  sweep<N, coupling::InertialNorm>(1, "InertialNorm", cfg, cfg.eps_inertial, sink);
  sweep<N, coupling::Dissipative>(2, "Dissipative", cfg, cfg.eps_dissipative, sink);
  sweep<N, coupling::DissipativeNorm>(3, "DissipativeNorm", cfg, cfg.eps_dissipative,
                                      sink);
  return 0;
}

int main(int argc, char **argv) {
  CLI::App app{"Van der Pol Multistability Experiment"};
  std::string config_path, output_path;
  app.add_option("config", config_path)->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", output_path)
      ->check(CLI::NonexistentPath | CLI::ExistingPath);
  CLI11_PARSE(app, argc, argv);

  YAML::Node yaml = YAML::LoadFile(config_path);
  const int N = yaml["N"].as<int>();
  switch (N) {
  case 3:
    return run<3>(yaml, config_path, output_path);
  default:
    std::cerr << "multistability supports only N=3 (got " << N << ")\n";
    return 1;
  }
}
