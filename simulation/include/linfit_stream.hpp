/*
  26.08.17 Аппроксимация данных прямой по методу наименьших квадратов
  (linear regression)

  Сначала накапливаем через push. Потом загружаем коэффиценты при
  помощи load

  Вариант линейной регрессии без B (то есть без intersept) не реализован

  В процессе счёта Dl может получаться отрицательным. Это значит, что
  S-коэффициенты накоплены с большими ошибками.

  25.03.26 Заменил вызовы my_assert на runtime_error, добавил функцию
  kv чтобы убрать ссылку на мою библиотеку
*/

#ifndef linfit_stream_hpp
#define linfit_stream_hpp

#include <cmath>
#include <stdexcept>

template<bool WITH_B, bool NEED_ERR, class Real, class LongReal=long double>
struct linfit_stream {

  LongReal Sx;
  LongReal Sxx;
  LongReal Sy;
  LongReal Syy;
  LongReal Sxy;
  std::size_t N;

  Real A,B;      // Коэффициенты линейной интерполяции y=A*x+B
  Real dA,dB;    // Погрешности A и B
  Real D;        // Сумма квадратов расстояний sum((A*x_i+B - y_i)^2)

  linfit_stream(void) {
    static_assert(WITH_B,"Not implemented without B");
    reset();
  }

  void reset(void) {Sx=Sxx=Sy=Syy=Sxy=0.0; N=0; A=B=D=dA=dB=0;}

  void push(Real xx, Real yy) {
    ++N;
    Sx+=xx; Sxx+=kv(xx);
    Sy+=yy; Sxy+=xx*yy;
    if (NEED_ERR) Syy+=kv(yy);
  }

  void load(void) {
    if (N == 0) throw std::runtime_error("N == 0");
    auto SxN = Sx / N;
    auto SyN = Sy / N;
    auto Vx = Sxx - SxN * Sx;
    if (Vx == 0.0) throw std::runtime_error("Vx == 0.0");
    auto Al = (Sxy - SxN * Sy) / Vx;
    auto Bl = (Sxx * SyN - SxN * Sxy) / Vx;
    if (NEED_ERR) {
      auto N2 = N * N;
      auto SxxN = Sxx / N2;
      auto SyyN = Syy / N2;
      auto SxyN = Sxy / N2;
      auto Dl = kv(Al) * SxxN + kv(Bl)/N + SyyN +
        2 * (Al * Bl * (SxN/N) - (Al * SxyN + Bl * (SyN/N)));
      if (Dl > -1e-10 && Dl < 0.0) Dl = 0.0; // прощаем небольшие погрешности
      // В принипе может вылезти отрицательное число. Это значит, что Sx, Sy, Sxx, Syy, Sxy накоплены с ошибками
      if (Dl < 0) throw std::runtime_error("Dl < 0");
      Dl *= N2;
      dA = static_cast<Real>(sqrt(Dl / ((N - 2) * Vx)));
      dB = static_cast<Real>(dA * sqrt(Sxx / N));
      D = static_cast<Real>(Dl);
    }
    A = static_cast<Real>(Al);
    B = static_cast<Real>(Bl);
  }

  Real kv(Real xx) {return xx*xx;}

  // void load(void) {
  //   LongReal Vx = Sxx*N-kv(Sx); my_assert(Vx!=0.0);
  //   A=static_cast<Real>((Sxy*N-Sx*Sy)/Vx);
  //   B=static_cast<Real>((Sxx*Sy-Sx*Sxy)/Vx);
  //   if (NEED_ERR) {
  //     D=kv(A)*Sxx + N*kv(B) + Syy + 2*(A*B*Sx - A*Sxy - B*Sy);
  //     dA=sqrt(N*D/((N-2)*Vx));
  //     dB=dA*sqrt(Sxx/N);
  //   }
  // } 
  
};

#endif

      // if (std::isnan(dA)) {
      //   std::cout << std::endl
      //             << "A=" << A
      //             << " B=" << B
      //             << " D=" << D
      //             << " dA=" << dA
      //             << " dB=" << dB
      //             << " Vx=" << Vx
      //             << " N=" << N
      //             << " Sx=" << Sx
      //             << " Sxx=" << Sxx
      //             << " Sy=" << Sy
      //             << " Syy=" << Syy
      //             << " Sxy=" << Sxy
      //             << std::endl;
      // }
