#define _USE_MATH_DEFINES
#include "../header/update.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/struct.h"

__device__ void Txx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t) {
  int i, j, k;
  int Txximax = sr.Txx.x, Txxjmax = sr.Txx.y, Txxkmax = sr.Txx.z;
  if (ip.mode == E_SINE) {
    ip.Txx[ip.in.x][ip.in.y][ip.in.z] = 0;
  } else if (ip.mode == E_RCOS) {
    if (t < 1. / ip.freq / dif.dt) {
      ip.Txx[ip.in.x][ip.in.y][ip.in.z] = 0;///* 8.e3 * 0.5 * */(1. - cos(2. * M_PI * ip.freq * (double)t * dif.dt)) / 2.;
    } else {
      ip.Txx[ip.in.x][ip.in.y][ip.in.z] = 0.;
    }
  } else {
    ip.Txx[ip.in.x][ip.in.y][ip.in.z] = 0.;
  }
#pragma omp parallel for private(i, j, k)
  for (k = 1; k <= Txxkmax - 1; k++) {
    for (j = 1; j <= Txxjmax - 1; j++) {
      for (i = 1; i <= Txximax - 1; i++) {
        aft->sa.Txxx[i][j][k] = (2. - ma.zetadx[i][j][k] * dif.dt) / (2. + ma.zetadx[i][j][k] * dif.dt) * bef->sa.Txxx[i][j][k]
         + 2. * (ma.c11[i][j][k] * dif.dt + ma.xi11[i][j][k]) / (2. + ma.zetadx[i][j][k] * dif.dt) * (aft->va.Vx[i][j][k] - aft->va.Vx[i - 1][j][k]) / dif.dx 
          - 2. * ma.xi11[i][j][k] / (2. + ma.zetadx[i][j][k] * dif.dt) * (bef->va.Vx[i][j][k] - bef->va.Vx[i - 1][j][k]) / dif.dx;

        aft->sa.Txxy[i][j][k] = (2. - ma.zetady[i][j][k] * dif.dt) / (2. + ma.zetady[i][j][k] * dif.dt) * bef->sa.Txxy[i][j][k]
         + 2. * (ma.ramda[i][j][k] * dif.dt + ma.khi[i][j][k]) / (2. + ma.zetady[i][j][k] * dif.dt) * (aft->va.Vy[i][j][k] - aft->va.Vy[i][j - 1][k]) / dif.dy
          - 2. * ma.khi[i][j][k] / (2. + ma.zetady[i][j][k] * dif.dt) * (bef->va.Vy[i][j][k] - bef->va.Vy[i][j - 1][k]) / dif.dy;

        aft->sa.Txxz[i][j][k] = (2. - ma.zetadz[i][j][k] * dif.dt) / (2. + ma.zetadz[i][j][k] * dif.dt) * bef->sa.Txxz[i][j][k]
         + 2. * (ma.ramda[i][j][k] * dif.dt + ma.khi[i][j][k]) / (2. + ma.zetadz[i][j][k] * dif.dt) * (aft->va.Vz[i][j][k] - aft->va.Vz[i][j][k - 1]) / dif.dz
          - 2. * ma.khi[i][j][k] / (2. + ma.zetadz[i][j][k] * dif.dt) * (bef->va.Vz[i][j][k] - bef->va.Vz[i][j][k - 1]) / dif.dz;
      }
    }
  }
  // // ここから
#pragma omp parallel for private(i, j)
  for (j = 0; j <= Txxjmax; j++) {
    for (i = 0; i <= Txximax; i++) {
      aft->sa.Txxx[i][j][0] = 0.;
      aft->sa.Txxx[i][j][Txxkmax] = 0.;
      aft->sa.Txxy[i][j][0] = 0.;
      aft->sa.Txxy[i][j][Txxkmax] = 0.;
      aft->sa.Txxz[i][j][0] = 0.;
      aft->sa.Txxz[i][j][Txxkmax] = 0.;
    }
  }
#pragma omp parallel for private(j, k)
  for (k = 1; k <= Txxkmax - 1; k++) {
    for (j = 0; j <= Txxjmax; j++) {
      aft->sa.Txxx[0][j][k] = 0.;
      aft->sa.Txxx[Txximax][j][k] = 0.;
      aft->sa.Txxy[0][j][k] = 0.;
      aft->sa.Txxy[Txximax][j][k] = 0.;
      aft->sa.Txxz[0][j][k] = 0.;
      aft->sa.Txxz[Txximax][j][k] = 0.;
    }
  }
#pragma omp parallel for private(i, k)
  for (k = 1; k <= Txxkmax - 1; k++) {
    for (i = 1; i <= Txximax - 1; i++) {
      aft->sa.Txxx[i][0][k] = 0.;
      aft->sa.Txxx[i][Txxjmax][k] = 0.;
      aft->sa.Txxy[i][0][k] = 0.;
      aft->sa.Txxy[i][Txxjmax][k] = 0.;
      aft->sa.Txxz[i][0][k] = 0.;
      aft->sa.Txxz[i][Txxjmax][k] = 0.;
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Txxkmax; k++) {
    for (j = 0; j <= Txxjmax; j++) {
      for (i = 0; i <= Txximax; i++) {
        aft->sa.Txx[i][j][k] = aft->sa.Txxx[i][j][k] + aft->sa.Txxy[i][j][k] + aft->sa.Txxz[i][j][k] + ip.Txx[i][j][k];
      }
    }
  }
}

__device__ void Tyy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t) {
  int i, j, k;
  int Tyyimax = sr.Tyy.x, Tyyjmax = sr.Tyy.y, Tyykmax = sr.Tyy.z;
  if (ip.mode == E_SINE) {
    ip.Tyy[ip.in.x][ip.in.y][ip.in.z] = 0;
  } else if (ip.mode == E_RCOS) {
    if (t < 1. / ip.freq / dif.dt) {
      ip.Tyy[ip.in.x][ip.in.y][ip.in.z] = 0;///* 8.e3 * 0.5 * */(1. - cos(2. * M_PI * ip.freq * (double)t * dif.dt)) / 2.;
    } else {
      ip.Tyy[ip.in.x][ip.in.y][ip.in.z] = 0.;
    }
  } else {
    ip.Tyy[ip.in.x][ip.in.y][ip.in.z] = 0.;
  }
#pragma omp parallel for private(i, j, k)
  for (k = 1; k <= Tyykmax - 1; k++) {
    for (j = 1; j <= Tyyjmax - 1; j++) {
      for (i = 1; i <= Tyyimax - 1; i++) {
        aft->sa.Tyyx[i][j][k] = (2. - ma.zetadx[i][j][k] * dif.dt) / (2. + ma.zetadx[i][j][k] * dif.dt) * bef->sa.Tyyx[i][j][k]
         + 2. * (ma.ramda[i][j][k] * dif.dt + ma.khi[i][j][k]) / (2. + ma.zetadx[i][j][k] * dif.dt) * (aft->va.Vx[i][j][k] - aft->va.Vx[i - 1][j][k]) / dif.dx
          - 2. * ma.khi[i][j][k] / (2. + ma.zetadx[i][j][k] * dif.dt) * (bef->va.Vx[i][j][k] - bef->va.Vx[i - 1][j][k]) / dif.dx;

        aft->sa.Tyyy[i][j][k] = (2. - ma.zetady[i][j][k] * dif.dt) / (2. + ma.zetady[i][j][k] * dif.dt) * bef->sa.Tyyy[i][j][k]
         + 2. * (ma.c11[i][j][k] * dif.dt + ma.xi11[i][j][k]) / (2. + ma.zetady[i][j][k] * dif.dt) * (aft->va.Vy[i][j][k] - aft->va.Vy[i][j - 1][k]) / dif.dy
          - 2. * ma.xi11[i][j][k] / (2. + ma.zetady[i][j][k] * dif.dt) * (bef->va.Vy[i][j][k] - bef->va.Vy[i][j - 1][k]) / dif.dy;

        aft->sa.Tyyz[i][j][k] = (2. - ma.zetadz[i][j][k] * dif.dt) / (2. + ma.zetadz[i][j][k] * dif.dt) * bef->sa.Tyyz[i][j][k]
         + 2. * (ma.ramda[i][j][k] * dif.dt + ma.khi[i][j][k]) / (2. + ma.zetadz[i][j][k] * dif.dt) * (aft->va.Vz[i][j][k] - aft->va.Vz[i][j][k - 1]) / dif.dz
          - 2. * ma.khi[i][j][k] / (2. + ma.zetadz[i][j][k] * dif.dt)  * (bef->va.Vz[i][j][k] - bef->va.Vz[i][j][k - 1]) / dif.dz;
      }
    }
  }
#pragma omp parallel for private(i, j)
  for (j = 0; j <= Tyyjmax; j++) {
    for (i = 0; i <= Tyyimax; i++) {
      aft->sa.Tyyx[i][j][0] = 0.;
      aft->sa.Tyyx[i][j][Tyykmax] = 0.;
      aft->sa.Tyyy[i][j][0] = 0.;
      aft->sa.Tyyy[i][j][Tyykmax] = 0.;
      aft->sa.Tyyz[i][j][0] = 0.;
      aft->sa.Tyyz[i][j][Tyykmax] = 0.;
    }
  }
#pragma omp parallel for private(j, k)
  for (k = 1; k <= Tyykmax - 1; k++) {
    for (j = 0; j <= Tyyjmax; j++) {
      aft->sa.Tyyx[0][j][k] = 0.;
      aft->sa.Tyyx[Tyyimax][j][k] = 0.;
      aft->sa.Tyyy[0][j][k] = 0.;
      aft->sa.Tyyy[Tyyimax][j][k] = 0.;
      aft->sa.Tyyz[0][j][k] = 0.;
      aft->sa.Tyyz[Tyyimax][j][k] = 0.;
    }
  }
#pragma omp parallel for private(i, k)
  for (k = 1; k <= Tyykmax - 1; k++) {
    for (i = 1; i <= Tyyimax - 1; i++) {
      aft->sa.Tyyx[i][0][k] = 0.;
      aft->sa.Tyyx[i][Tyyjmax][k] = 0.;
      aft->sa.Tyyy[i][0][k] = 0.;
      aft->sa.Tyyy[i][Tyyjmax][k] = 0.;
      aft->sa.Tyyz[i][0][k] = 0.;
      aft->sa.Tyyz[i][Tyyjmax][k] = 0.;
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Tyykmax; k++) {
    for (j = 0; j <= Tyyjmax; j++) {
      for (i = 0; i <= Tyyimax; i++) {
        aft->sa.Tyy[i][j][k] = aft->sa.Tyyx[i][j][k] + aft->sa.Tyyy[i][j][k] + aft->sa.Tyyz[i][j][k] + ip.Tyy[i][j][k];
      }
    }
  }
}

__device__ void Tzz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t) {
  int i, j, k;
  int Tzzimax = sr.Tzz.x, Tzzjmax = sr.Tzz.y, Tzzkmax = sr.Tzz.z;
  if (ip.mode == E_SINE) {
    ip.Tzz[ip.in.x][ip.in.y][ip.in.z] = (-1) * sin(2. * M_PI * ip.freq * (double)t * dif.dt) / 2.;
  } else if (ip.mode == E_RCOS) {
    if (t < 1. / ip.freq / dif.dt) {
      ip.Tzz[ip.in.x][ip.in.y][ip.in.z] = 8.e3 * 0.5 * (-1) * (1. - cos(2. * M_PI * ip.freq * (double)t * dif.dt)) / 2.;
    } else {
      ip.Tzz[ip.in.x][ip.in.y][ip.in.z] = 0.;
    }
  } else {
    ip.Tzz[ip.in.x][ip.in.y][ip.in.z] = 0.;
  }

  // Tzzの更新式
#pragma omp parallel for private(i, j, k)
  for (k = 1; k <= Tzzkmax - 1; k++) {
    for (j = 1; j <= Tzzjmax - 1; j++) {
      for (i = 1; i <= Tzzimax - 1; i++) {
        aft->sa.Tzzx[i][j][k] = (2. - ma.zetadx[i][j][k] * dif.dt) / (2. + ma.zetadx[i][j][k] * dif.dt) * bef->sa.Tzzx[i][j][k]
         + 2. * (ma.ramda[i][j][k] * dif.dt + ma.khi[i][j][k]) / (2. + ma.zetadx[i][j][k] * dif.dt) * (aft->va.Vx[i][j][k] - aft->va.Vx[i - 1][j][k]) / dif.dx
          - 2. * ma.khi[i][j][k] / (2. + ma.zetadx[i][j][k] * dif.dt) * (bef->va.Vx[i][j][k] - bef->va.Vx[i - 1][j][k]) / dif.dx;

        aft->sa.Tzzy[i][j][k] = (2. - ma.zetady[i][j][k] * dif.dt) / (2. + ma.zetady[i][j][k] * dif.dt) * bef->sa.Tzzy[i][j][k]
         + 2. * (ma.ramda[i][j][k] * dif.dt + ma.khi[i][j][k]) / (2. + ma.zetady[i][j][k] * dif.dt) * (aft->va.Vy[i][j][k] - aft->va.Vy[i][j - 1][k]) / dif.dy
          - 2. * ma.khi[i][j][k] / (2. + ma.zetady[i][j][k] * dif.dt)  * (bef->va.Vy[i][j][k] - bef->va.Vy[i][j - 1][k]) / dif.dy;

        aft->sa.Tzzz[i][j][k] = (2. - ma.zetadz[i][j][k] * dif.dt) / (2. + ma.zetadz[i][j][k] * dif.dt) * bef->sa.Tzzz[i][j][k]
         + 2. * (ma.c11[i][j][k] * dif.dt + ma.xi11[i][j][k]) / (2. + ma.zetadz[i][j][k] * dif.dt) * (aft->va.Vz[i][j][k] - aft->va.Vz[i][j][k - 1]) / dif.dz
          - 2. * ma.xi11[i][j][k] / (2. + ma.zetadz[i][j][k] * dif.dt) * (bef->va.Vz[i][j][k] - bef->va.Vz[i][j][k - 1]) / dif.dz;
      }
    }
  }
#pragma omp parallel for private(i, j)
  for (j = 0; j <= Tzzjmax; j++) {
    for (i = 0; i <= Tzzimax; i++) {
      aft->sa.Tzzx[i][j][0] = 0.;
      aft->sa.Tzzx[i][j][Tzzkmax] = 0.;
      aft->sa.Tzzy[i][j][0] = 0.;
      aft->sa.Tzzy[i][j][Tzzkmax] = 0.;
      aft->sa.Tzzz[i][j][0] = 0.;
      aft->sa.Tzzz[i][j][Tzzkmax] = 0.;
    }
  }
#pragma omp parallel for private(j, k)
  for (k = 1; k <= Tzzkmax - 1; k++) {
    for (j = 0; j <= Tzzjmax; j++) {
      aft->sa.Tzzx[0][j][k] = 0.;
      aft->sa.Tzzx[Tzzimax][j][k] = 0.;
      aft->sa.Tzzy[0][j][k] = 0.;
      aft->sa.Tzzy[Tzzimax][j][k] = 0.;
      aft->sa.Tzzz[0][j][k] = 0.;
      aft->sa.Tzzz[Tzzimax][j][k] = 0.;
    }
  }
#pragma omp parallel for private(i, k)
  for (k = 1; k <= Tzzkmax - 1; k++) {
    for (i = 1; i <= Tzzimax - 1; i++) {
      aft->sa.Tzzx[i][0][k] = 0.;
      aft->sa.Tzzx[i][Tzzjmax][k] = 0.;
      aft->sa.Tzzy[i][0][k] = 0.;
      aft->sa.Tzzy[i][Tzzjmax][k] = 0.;
      aft->sa.Tzzz[i][0][k] = 0.;
      aft->sa.Tzzz[i][Tzzjmax][k] = 0.;
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Tzzkmax; k++) {
    for (j = 0; j <= Tzzjmax; j++) {
      for (i = 0; i <= Tzzimax; i++) {
        aft->sa.Tzz[i][j][k] = aft->sa.Tzzx[i][j][k] + aft->sa.Tzzy[i][j][k] + aft->sa.Tzzz[i][j][k] + ip.Tzz[i][j][k];
      }
    }
  }
}
//垂直応力計算
__global__ void Sig(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t) {
  Txx(aft, bef, ma, dif, sr, ip, t);
  Tyy(aft, bef, ma, dif, sr, ip, t);
  Tzz(aft, bef, ma, dif, sr, ip, t);
}

void Txy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr) {
  int i, j, k;
  int Txyimax = tr.Txy.x, Txyjmax = tr.Txy.y, Txykmax = tr.Txy.z;
  double Hzetadx, Hzetady, Hzetadz, Hmu, Hgamma;
#pragma omp parallel for private(i, j, k)
  for (k = 1; k <= Txykmax - 1; k++) {
    for (j = 0; j <= Txyjmax; j++) {
      for (i = 0; i <= Txyimax; i++) {
        //PML:減衰係数,計算領域:摩擦定数
        Hzetadx = 4. * pow((1. / ma.zetadx[i + 1][j + 1][k]) + (1. / ma.zetadx[i][j + 1][k]) + (1. / ma.zetadx[i + 1][j][k]) + (1. / ma.zetadx[i][j][k]), -1.);
        //PML:減衰係数,計算領域:摩擦定数
        Hzetady = 4. * pow((1. / ma.zetady[i + 1][j + 1][k]) + (1. / ma.zetady[i][j + 1][k]) + (1. / ma.zetady[i + 1][j][k]) + (1. / ma.zetady[i][j][k]), -1.);
        //第2ラメ，横弾性係数(剛性率)
        Hmu     = 4. * pow((1. /     ma.mu[i + 1][j + 1][k]) + (1. /     ma.mu[i][j + 1][k]) + (1. /     ma.mu[i + 1][j][k]) + (1. /     ma.mu[i][j][k]), -1.);
        //第１粘性定数
        Hgamma  = 4. * pow((1. /  ma.gamma[i + 1][j + 1][k]) + (1. /  ma.gamma[i][j + 1][k]) + (1. /  ma.gamma[i + 1][j][k]) + (1. /  ma.gamma[i][j][k]), -1.);
        aft->ta.Txyx[i][j][k] = (2. - Hzetadx * dif.dt) / (2. + Hzetadx * dif.dt) * bef->ta.Txyx[i][j][k]
         + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetadx * dif.dt) * (aft->va.Vy[i + 1][j][k] - aft->va.Vy[i][j][k]) / dif.dx
          - 2. * Hgamma / (2. + Hzetadx * dif.dt) * (bef->va.Vy[i + 1][j][k] - bef->va.Vy[i][j][k]) / dif.dx;

        aft->ta.Txyy[i][j][k] = (2. - Hzetady * dif.dt) / (2. + Hzetady * dif.dt) * bef->ta.Txyy[i][j][k]
         + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetady * dif.dt) * (aft->va.Vx[i][j + 1][k] - aft->va.Vx[i][j][k]) / dif.dy
          - 2. * Hgamma / (2. + Hzetady * dif.dt) * (bef->va.Vx[i][j + 1][k] - bef->va.Vx[i][j][k]) / dif.dy;
      }
    }
  }
#pragma omp parallel for private(i, j)
  for (j = 0; j <= Txyjmax; j++) {
    for (i = 0; i <= Txyimax; i++) {
      aft->ta.Txyx[i][j][0] = 0.;
      aft->ta.Txyx[i][j][Txykmax] = 0.;
      aft->ta.Txyy[i][j][0] = 0.;
      aft->ta.Txyy[i][j][Txykmax] = 0.;
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Txykmax; k++) {
    for (j = 0; j <= Txyjmax; j++) {
      for (i = 0; i <= Txyimax; i++) {
        aft->ta.Txy[i][j][k] = aft->ta.Txyx[i][j][k] + aft->ta.Txyy[i][j][k];
      }
    }
  }
}

void Tyz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr) {
  int i, j, k;
  int Tyzimax = tr.Tyz.x, Tyzjmax = tr.Tyz.y, Tyzkmax = tr.Tyz.z;
  double Hzetadx, Hzetady, Hzetadz, Hmu, Hgamma;
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Tyzkmax; k++) {
    for (j = 0; j <= Tyzjmax; j++) {
      for (i = 1; i <= Tyzimax - 1; i++) {
        Hzetady = 4. * pow((1. / ma.zetady[i][j + 1][k + 1]) + (1. / ma.zetady[i][j + 1][k]) + (1. / ma.zetady[i][j][k + 1]) + (1. / ma.zetady[i][j][k]), -1.);
        Hzetadz = 4. * pow((1. / ma.zetadz[i][j + 1][k + 1]) + (1. / ma.zetadz[i][j + 1][k]) + (1. / ma.zetadz[i][j][k + 1]) + (1. / ma.zetadz[i][j][k]), -1.);
        Hmu     = 4. * pow((1. /     ma.mu[i][j + 1][k + 1]) + (1. /     ma.mu[i][j + 1][k]) + (1. /     ma.mu[i][j][k + 1]) + (1. /     ma.mu[i][j][k]), -1.);
        Hgamma  = 4. * pow((1. /  ma.gamma[i][j + 1][k + 1]) + (1. /  ma.gamma[i][j + 1][k]) + (1. /  ma.gamma[i][j][k + 1]) + (1. /  ma.gamma[i][j][k]), -1.);
        aft->ta.Tyzy[i][j][k] = (2. - Hzetady * dif.dt) / (2. + Hzetady * dif.dt) * bef->ta.Tyzy[i][j][k]
         + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetady * dif.dt) * (aft->va.Vz[i][j + 1][k] - aft->va.Vz[i][j][k]) / dif.dy
          - 2. * Hgamma / (2. + Hzetady * dif.dt) * (bef->va.Vz[i][j + 1][k] - bef->va.Vz[i][j][k]) / dif.dy;

        aft->ta.Tyzz[i][j][k] = (2. - Hzetadz * dif.dt) / (2. + Hzetadz * dif.dt) * bef->ta.Tyzz[i][j][k]
         + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetadz * dif.dt) * (aft->va.Vy[i][j][k + 1] - aft->va.Vy[i][j][k]) / dif.dz 
          - 2. * Hgamma / (2. + Hzetadz * dif.dt) * (bef->va.Vy[i][j][k + 1] - bef->va.Vy[i][j][k]) / dif.dz;
      }
    }
  }
#pragma omp parallel for private(j, k)
  for (k = 0; k <= Tyzkmax; k++) {
    for (j = 0; j <= Tyzjmax; j++) {
      aft->ta.Tyzy[0][j][k] = 0.;
      aft->ta.Tyzy[Tyzimax][j][k] = 0.;
      aft->ta.Tyzz[0][j][k] = 0.;
      aft->ta.Tyzz[Tyzimax][j][k] = 0.;
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Tyzkmax; k++) {
    for (j = 0; j <= Tyzjmax; j++) {
      for (i = 0; i <= Tyzimax; i++) {
        aft->ta.Tyz[i][j][k] = aft->ta.Tyzy[i][j][k] + aft->ta.Tyzz[i][j][k];
      }
    }
  }
}

void Tzx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr) {
  int i, j, k;
  int Tzximax = tr.Tzx.x, Tzxjmax = tr.Tzx.y, Tzxkmax = tr.Tzx.z;
  double Hzetadx, Hzetady, Hzetadz, Hmu, Hgamma;
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Tzxkmax; k++) {
    for (j = 1; j <= Tzxjmax - 1; j++) {
      for (i = 0; i <= Tzximax; i++) {
        Hzetadx = 4. * pow((1. / ma.zetadx[i + 1][j][k + 1]) + (1. / ma.zetadx[i + 1][j][k]) + (1. / ma.zetadx[i][j][k + 1]) + (1. / ma.zetadx[i][j][k]), -1.);
        Hzetadz = 4. * pow((1. / ma.zetadz[i + 1][j][k + 1]) + (1. / ma.zetadz[i + 1][j][k]) + (1. / ma.zetadz[i][j][k + 1]) + (1. / ma.zetadz[i][j][k]), -1.);
        Hmu     = 4. * pow((1. /     ma.mu[i + 1][j][k + 1]) + (1. /     ma.mu[i + 1][j][k]) + (1. /     ma.mu[i][j][k + 1]) + (1. /     ma.mu[i][j][k]), -1.);
        Hgamma  = 4. * pow((1. /  ma.gamma[i + 1][j][k + 1]) + (1. /  ma.gamma[i + 1][j][k]) + (1. /  ma.gamma[i][j][k + 1]) + (1. /  ma.gamma[i][j][k]), -1.);
        aft->ta.Tzxz[i][j][k] = (2. - Hzetadz * dif.dt) / (2. + Hzetadz * dif.dt) * bef->ta.Tzxz[i][j][k]
         + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetadz * dif.dt) * (aft->va.Vx[i][j][k + 1] - aft->va.Vx[i][j][k]) / dif.dz
          - 2. * Hgamma / (2. + Hzetadz * dif.dt) * (bef->va.Vx[i][j][k + 1] - bef->va.Vx[i][j][k]) / dif.dz;

        aft->ta.Tzxx[i][j][k] = (2. - Hzetadx * dif.dt) / (2. + Hzetadx * dif.dt) * bef->ta.Tzxx[i][j][k]
         + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetadx * dif.dt) * (aft->va.Vz[i + 1][j][k] - aft->va.Vz[i][j][k]) / dif.dx
          - 2. * Hgamma / (2. + Hzetadx * dif.dt) * (bef->va.Vz[i + 1][j][k] - bef->va.Vz[i][j][k]) / dif.dx;
      }
    }
  }
#pragma omp parallel for private(i, k)
  for (k = 0; k <= Tzxkmax; k++) {
    for (i = 0; i <= Tzximax; i++) {
      aft->ta.Tzxx[i][0][k] = 0.;
      aft->ta.Tzxx[i][Tzxjmax][k] = 0.;
      aft->ta.Tzxz[i][0][k] = 0.;
      aft->ta.Tzxz[i][Tzxjmax][k] = 0.;
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Tzxkmax; k++) {
    for (j = 0; j <= Tzxjmax; j++) {
      for (i = 0; i <= Tzximax; i++) {
        aft->ta.Tzx[i][j][k] = aft->ta.Tzxx[i][j][k] + aft->ta.Tzxz[i][j][k];
      }
    }
  }
}
//せん断応力計算
void Tau(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr) {
  Txy(aft, bef, ma, dif, tr);
  Tyz(aft, bef, ma, dif, tr);
  Tzx(aft, bef, ma, dif, tr);
}

void Vx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr) {
  int i, j, k;
  int Vximax = vr.Vx.x, Vxjmax = vr.Vx.y, Vxkmax = vr.Vx.z;
  double Azetaxx, Azetaxy, Azetaxz, Arho;
#pragma omp parallel for private(i, j, k)
  for (k = 1; k <= Vxkmax - 1; k++) {
    for (j = 1; j <= Vxjmax - 1; j++) {
      for (i = 0; i <= Vximax; i++) {
        Azetaxx = (ma.zetaxx[i + 1][j][k] + ma.zetaxx[i][j][k]) / 2.;
        Azetaxy = (ma.zetaxy[i + 1][j][k] + ma.zetaxy[i][j][k]) / 2.;
        Azetaxz = (ma.zetaxz[i + 1][j][k] + ma.zetaxz[i][j][k]) / 2.;
        Arho    = (   ma.rho[i + 1][j][k] +    ma.rho[i][j][k]) / 2.;
        aft->va.Vxx[i][j][k] = (2. * Arho - Azetaxx * dif.dt) / (2. * Arho + Azetaxx * dif.dt) * bef->va.Vxx[i][j][k]
         + 2. * dif.dt / (2. * Arho + Azetaxx * dif.dt) * (bef->sa.Txx[i + 1][j][k] - bef->sa.Txx[i][j][k]) / dif.dx;

        aft->va.Vxy[i][j][k] = (2. * Arho - Azetaxy * dif.dt) / (2. * Arho + Azetaxy * dif.dt) * bef->va.Vxy[i][j][k]
         + 2. * dif.dt / (2. * Arho + Azetaxy * dif.dt) * (bef->ta.Txy[i][j][k] - bef->ta.Txy[i][j - 1][k]) / dif.dy;

        aft->va.Vxz[i][j][k] = (2. * Arho - Azetaxz * dif.dt) / (2. * Arho + Azetaxz * dif.dt) * bef->va.Vxz[i][j][k]
         + 2. * dif.dt / (2. * Arho + Azetaxz * dif.dt) * (bef->ta.Tzx[i][j][k] - bef->ta.Tzx[i][j][k - 1]) / dif.dz;
      }
    }
  }
#pragma omp parallel for private(i, k)
  for (k = 0; k <= Vxkmax; k++) {
    for (i = 0; i <= Vximax; i++) {
      aft->va.Vxx[i][0][k] = 0.;
      aft->va.Vxx[i][Vxjmax][k] = 0.;
      aft->va.Vxy[i][0][k] = 0.;
      aft->va.Vxy[i][Vxjmax][k] = 0.;
      aft->va.Vxz[i][0][k] = 0.;
      aft->va.Vxz[i][Vxjmax][k] = 0.;
    }
  }
#pragma omp parallel for private(i, j)
  for (j = 1; j <= Vxjmax - 1; j++) {
    for (i = 0; i <= Vximax; i++) {
      aft->va.Vxx[i][j][0] = 0.;
      aft->va.Vxx[i][j][Vxkmax] = 0.;
      aft->va.Vxy[i][j][0] = 0.;
      aft->va.Vxy[i][j][Vxkmax] = 0.;
      aft->va.Vxz[i][j][0] = 0.;
      aft->va.Vxz[i][j][Vxkmax] = 0.;
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Vxkmax; k++) {
    for (j = 0; j <= Vxjmax; j++) {
      for (i = 0; i <= Vximax; i++) {
        aft->va.Vx[i][j][k] = aft->va.Vxx[i][j][k] + aft->va.Vxy[i][j][k] + aft->va.Vxz[i][j][k];
      }
    }
  }
}

void Vy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr) {
  int i, j, k;
  int Vyimax = vr.Vy.x, Vyjmax = vr.Vy.y, Vykmax = vr.Vy.z;
  double Azetayx, Azetayy, Azetayz, Arho;
#pragma omp parallel for private(i, j, k)
  for (k = 1; k <= Vykmax - 1; k++) {
    for (j = 0; j <= Vyjmax; j++) {
      for (i = 1; i <= Vyimax - 1; i++) {
        Azetayx = (ma.zetayx[i][j + 1][k] + ma.zetayx[i][j][k]) / 2.;
        Azetayy = (ma.zetayy[i][j + 1][k] + ma.zetayy[i][j][k]) / 2.;
        Azetayz = (ma.zetayz[i][j + 1][k] + ma.zetayz[i][j][k]) / 2.;
        Arho    = (   ma.rho[i][j + 1][k] +    ma.rho[i][j][k]) / 2.;
        aft->va.Vyx[i][j][k] = (2. * Arho - Azetayx * dif.dt) / (2. * Arho + Azetayx * dif.dt) * bef->va.Vyx[i][j][k]
         + 2. * dif.dt / (2. * Arho + Azetayx * dif.dt) * (bef->ta.Txy[i][j][k] - bef->ta.Txy[i - 1][j][k]) / dif.dx;
        aft->va.Vyy[i][j][k] = (2. * Arho - Azetayy * dif.dt) / (2. * Arho + Azetayy * dif.dt) * bef->va.Vyy[i][j][k]
         + 2. * dif.dt / (2. * Arho + Azetayy * dif.dt) * (bef->sa.Tyy[i][j + 1][k] - bef->sa.Tyy[i][j][k]) / dif.dy;
        aft->va.Vyz[i][j][k] = (2. * Arho - Azetayz * dif.dt) / (2. * Arho + Azetayz * dif.dt) * bef->va.Vyz[i][j][k]
         + 2. * dif.dt / (2. * Arho + Azetayz * dif.dt) * (bef->ta.Tyz[i][j][k] - bef->ta.Tyz[i][j][k - 1]) / dif.dz;
      }
    }
  }
#pragma omp parallel for private(j, k)
  for (k = 0; k <= Vykmax; k++) {
    for (j = 0; j <= Vyjmax; j++) {
      aft->va.Vyx[0][j][k] = 0.;
      aft->va.Vyx[Vyimax][j][k] = 0.;
      aft->va.Vyy[0][j][k] = 0.;
      aft->va.Vyy[Vyimax][j][k] = 0.;
      aft->va.Vyz[0][j][k] = 0.;
      aft->va.Vyz[Vyimax][j][k] = 0.;
    }
  }
#pragma omp parallel for private(i, j)
  for (j = 0; j <= Vyjmax; j++) {
    for (i = 1; i <= Vyimax - 1; i++) {
      aft->va.Vyx[i][j][0] = 0.;
      aft->va.Vyx[i][j][Vykmax] = 0.;
      aft->va.Vyy[i][j][0] = 0.;
      aft->va.Vyy[i][j][Vykmax] = 0.;
      aft->va.Vyz[i][j][0] = 0.;
      aft->va.Vyz[i][j][Vykmax] = 0.;
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Vykmax; k++) {
    for (j = 0; j <= Vyjmax; j++) {
      for (i = 0; i <= Vyimax; i++) {
        aft->va.Vy[i][j][k] = aft->va.Vyx[i][j][k] + aft->va.Vyy[i][j][k] + aft->va.Vyz[i][j][k];
      }
    }
  }
}

void Vz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr) {
  int i, j, k;
  int Vzimax = vr.Vz.x, Vzjmax = vr.Vz.y, Vzkmax = vr.Vz.z;
  double Azetazx, Azetazy, Azetazz, Arho;
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Vzkmax; k++) {
    for (j = 1; j <= Vzjmax - 1; j++) {
      for (i = 1; i <= Vzimax - 1; i++) {
        Azetazx = (ma.zetazx[i][j][k + 1] + ma.zetazx[i][j][k]) / 2.;
        Azetazy = (ma.zetazy[i][j][k + 1] + ma.zetazy[i][j][k]) / 2.;
        Azetazz = (ma.zetazz[i][j][k + 1] + ma.zetazz[i][j][k]) / 2.;
        Arho    = (   ma.rho[i][j][k + 1] +    ma.rho[i][j][k]) / 2.;
        aft->va.Vzx[i][j][k] = (2. * Arho - Azetazx * dif.dt) / (2. * Arho + Azetazx * dif.dt) * bef->va.Vzx[i][j][k]
         + 2. * dif.dt / (2. * Arho + Azetazx * dif.dt) * (bef->ta.Tzx[i][j][k] - bef->ta.Tzx[i - 1][j][k]) / dif.dx;
        aft->va.Vzy[i][j][k] = (2. * Arho - Azetazy * dif.dt) / (2. * Arho + Azetazy * dif.dt) * bef->va.Vzy[i][j][k]
         + 2. * dif.dt / (2. * Arho + Azetazy * dif.dt) * (bef->ta.Tyz[i][j][k] - bef->ta.Tyz[i][j - 1][k]) / dif.dy;
        aft->va.Vzz[i][j][k] = (2. * Arho - Azetazz * dif.dt) / (2. * Arho + Azetazz * dif.dt) * bef->va.Vzz[i][j][k]
         + 2. * dif.dt / (2. * Arho + Azetazz * dif.dt) * (bef->sa.Tzz[i][j][k + 1] - bef->sa.Tzz[i][j][k]) / dif.dz;
      }
    }
  }
#pragma omp parallel for private(j, k)
  for (k = 0; k <= Vzkmax; k++) {
    for (j = 0; j <= Vzjmax; j++) {
      aft->va.Vzx[0][j][k] = 0.;
      aft->va.Vzx[Vzimax][j][k] = 0.;
      aft->va.Vzy[0][j][k] = 0.;
      aft->va.Vzy[Vzimax][j][k] = 0.;
      aft->va.Vzz[0][j][k] = 0.;
      aft->va.Vzz[Vzimax][j][k] = 0.;
    }
  }
#pragma omp parallel for private(i, k)
  for (k = 0; k <= Vzkmax; k++) {
    for (i = 1; i <= Vzimax - 1; i++) {
      aft->va.Vzx[i][0][k] = 0.;
      aft->va.Vzx[i][Vzjmax][k] = 0.;
      aft->va.Vzy[i][0][k] = 0.;
      aft->va.Vzy[i][Vzjmax][k] = 0.;
      aft->va.Vzz[i][0][k] = 0.;
      aft->va.Vzz[i][Vzjmax][k] = 0.;
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Vzkmax; k++) {
    for (j = 0; j <= Vzjmax; j++) {
      for (i = 0; i <= Vzimax; i++) {
        aft->va.Vz[i][j][k] = aft->va.Vzx[i][j][k] + aft->va.Vzy[i][j][k] + aft->va.Vzz[i][j][k];
      }
    }
  }
}
//粒子速度計算
void Vel(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr) {
  Vx(aft, bef, ma, dif, vr);
  Vy(aft, bef, ma, dif, vr);
  Vz(aft, bef, ma, dif, vr);
}

void Acc(Coord_acc *A,BefAft *aft, BefAft *bef, Diff dif, Coord out, Coord center){
  A->x = ((aft->va.Vx[out.x - 1][out.y][out.z] - bef->va.Vx[out.x - 1][out.y][out.z]) / dif.dt  + (aft->va.Vx[out.x][out.y][out.z] - bef->va.Vx[out.x][out.y][out.z]) / dif.dt) / 2;
  A->y = ((aft->va.Vy[out.x][out.y - 1][out.z] - bef->va.Vy[out.x][out.y - 1][out.z]) / dif.dt  + (aft->va.Vy[out.x][out.y][out.z] - bef->va.Vy[out.x][out.y][out.z]) / dif.dt) / 2;
  A->z = ((aft->va.Vz[out.x][out.y][out.z - 1] - bef->va.Vz[out.x][out.y][out.z - 1]) / dif.dt  + (aft->va.Vz[out.x][out.y][out.z] - bef->va.Vz[out.x][out.y][out.z]) / dif.dt) / 2;

}
//更新
void swapBefAft(BefAft *aft, BefAft *bef, Range ran) {
  int i, j, k;
  int Txximax = ran.sr.Txx.x, Txxjmax = ran.sr.Txx.y, Txxkmax = ran.sr.Txx.z;
  int Tyyimax = ran.sr.Tyy.x, Tyyjmax = ran.sr.Tyy.y, Tyykmax = ran.sr.Tyy.z;
  int Tzzimax = ran.sr.Tzz.x, Tzzjmax = ran.sr.Tzz.y, Tzzkmax = ran.sr.Tzz.z;
  int Txyimax = ran.tr.Txy.x, Txyjmax = ran.tr.Txy.y, Txykmax = ran.tr.Txy.z;
  int Tyzimax = ran.tr.Tyz.x, Tyzjmax = ran.tr.Tyz.y, Tyzkmax = ran.tr.Tyz.z;
  int Tzximax = ran.tr.Tzx.x, Tzxjmax = ran.tr.Tzx.y, Tzxkmax = ran.tr.Tzx.z;
  int Vximax = ran.vr.Vx.x, Vxjmax = ran.vr.Vx.y, Vxkmax = ran.vr.Vx.z;
  int Vyimax = ran.vr.Vy.x, Vyjmax = ran.vr.Vy.y, Vykmax = ran.vr.Vy.z;
  int Vzimax = ran.vr.Vz.x, Vzjmax = ran.vr.Vz.y, Vzkmax = ran.vr.Vz.z;
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Txxkmax; k++) {
    for (j = 0; j <= Txxjmax; j++) {
      for (i = 0; i <= Txximax; i++) {
        bef->sa.Txx[i][j][k] = aft->sa.Txx[i][j][k];
        bef->sa.Txxx[i][j][k] = aft->sa.Txxx[i][j][k];
        bef->sa.Txxy[i][j][k] = aft->sa.Txxy[i][j][k];
        bef->sa.Txxz[i][j][k] = aft->sa.Txxz[i][j][k];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Tyykmax; k++) {
    for (j = 0; j <= Tyyjmax; j++) {
      for (i = 0; i <= Tyyimax; i++) {
        bef->sa.Tyy[i][j][k] = aft->sa.Tyy[i][j][k];
        bef->sa.Tyyx[i][j][k] = aft->sa.Tyyx[i][j][k];
        bef->sa.Tyyy[i][j][k] = aft->sa.Tyyy[i][j][k];
        bef->sa.Tyyz[i][j][k] = aft->sa.Tyyz[i][j][k];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Tzzkmax; k++) {
    for (j = 0; j <= Tzzjmax; j++) {
      for (i = 0; i <= Tzzimax; i++) {
        bef->sa.Tzz[i][j][k] = aft->sa.Tzz[i][j][k];
        bef->sa.Tzzx[i][j][k] = aft->sa.Tzzx[i][j][k];
        bef->sa.Tzzy[i][j][k] = aft->sa.Tzzy[i][j][k];
        bef->sa.Tzzz[i][j][k] = aft->sa.Tzzz[i][j][k];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (i = 0; i <= Txyimax; i++) {
    for (j = 0; j <= Txyjmax; j++) {
      for (k = 0; k <= Txykmax; k++) {
        bef->ta.Txy[i][j][k] = aft->ta.Txy[i][j][k];
        bef->ta.Txyx[i][j][k] = aft->ta.Txyx[i][j][k];
        bef->ta.Txyy[i][j][k] = aft->ta.Txyy[i][j][k];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (i = 0; i <= Tyzimax; i++) {
    for (j = 0; j <= Tyzjmax; j++) {
      for (k = 0; k <= Tyzkmax; k++) {
        bef->ta.Tyz[i][j][k] = aft->ta.Tyz[i][j][k];
        bef->ta.Tyzy[i][j][k] = aft->ta.Tyzy[i][j][k];
        bef->ta.Tyzz[i][j][k] = aft->ta.Tyzz[i][j][k];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (i = 0; i <= Tzximax; i++) {
    for (j = 0; j <= Tzxjmax; j++) {
      for (k = 0; k <= Tzxkmax; k++) {
        bef->ta.Tzx[i][j][k] = aft->ta.Tzx[i][j][k];
        bef->ta.Tzxz[i][j][k] = aft->ta.Tzxz[i][j][k];
        bef->ta.Tzxx[i][j][k] = aft->ta.Tzxx[i][j][k];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Vxkmax; k++) {
    for (j = 0; j <= Vxjmax; j++) {
      for (i = 0; i <= Vximax; i++) {
        bef->va.Vx[i][j][k] = aft->va.Vx[i][j][k];
        bef->va.Vxx[i][j][k] = aft->va.Vxx[i][j][k];
        bef->va.Vxy[i][j][k] = aft->va.Vxy[i][j][k];
        bef->va.Vxz[i][j][k] = aft->va.Vxz[i][j][k];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Vykmax; k++) {
    for (j = 0; j <= Vyjmax; j++) {
      for (i = 0; i <= Vyimax; i++) {
        bef->va.Vy[i][j][k] = aft->va.Vy[i][j][k];
        bef->va.Vyx[i][j][k] = aft->va.Vyx[i][j][k];
        bef->va.Vyy[i][j][k] = aft->va.Vyy[i][j][k];
        bef->va.Vyz[i][j][k] = aft->va.Vyz[i][j][k];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k <= Vzkmax; k++) {
    for (j = 0; j <= Vzjmax; j++) {
      for (i = 0; i <= Vzimax; i++) {
        bef->va.Vz[i][j][k] = aft->va.Vz[i][j][k];
        bef->va.Vzx[i][j][k] = aft->va.Vzx[i][j][k];
        bef->va.Vzy[i][j][k] = aft->va.Vzy[i][j][k];
        bef->va.Vzz[i][j][k] = aft->va.Vzz[i][j][k];
      }
    }
  }
}