#define _USE_MATH_DEFINES
#include "../header/init.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/struct.h"

#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))

void initMedium(Medium *med) {
  for (int mednum = 0; mednum < E_M_END; mednum++) {
    switch (mednum) {
      case E_AIR:
        med[mednum].rho = 1.205;  // 密度/////////////////////////////////////////////////////
        med[mednum].K = 1.422e5;  // 体積弾性率
        med[mednum].E = 0.;       // ヤング率
        med[mednum].G = 0.;       // 剛性率/////////////////////////////////////////////////////
        med[mednum].nu = 0.;     // ポアソン比
        med[mednum].ramda = med[mednum].K - 2. / 3. * med[mednum].G;  // 第1ラメ定数/////////////////////////////////////////////////////
        med[mednum].zeta = 5.;    //セル間摩擦係数/////////////////////////////////////////////////////
        med[mednum].gamma = 1.8e-5;   //第1粘性/////////////////////////////////////////////////////
        med[mednum].khi = 0.;         //第2粘性/////////////////////////////////////////////////////
        med[mednum].eta = 0.;
        med[mednum].omega = 0.;
        break;
      case E_CON:
        med[mednum].rho = 2400.;  // 密度/////////////////////////////////////////////////////
        med[mednum].E = 2.4e10;   // ヤング率
        med[mednum].nu = 0.2;     // ポアソン比
        med[mednum].G = med[mednum].E / 2. / (1. + med[mednum].nu);  // 剛性率////////////////////////////////////第2ラメ
        med[mednum].K = med[mednum].E / 3. / (1. - 2. * med[mednum].nu);  // 体積弾性率
        med[mednum].ramda = med[mednum].E * med[mednum].nu / (1. + med[mednum].nu) / (1. - 2. * med[mednum].nu);  // 第1ラメ定数//////////////
        med[mednum].zeta = 2.5e4;//セル間摩擦係数////////////////////////////////////////////////////////////////////////////////////
        med[mednum].eta = 0.005;//粘性定数算出係数(損失係数)
        med[mednum].omega = 2. * M_PI * 32.;//粘性定数算出係数(角周波数)
        med[mednum].gamma = med[mednum].eta * med[mednum].G / med[mednum].omega;//第1粘性定数//////////////////////////////////////
        med[mednum].khi = med[mednum].eta * med[mednum].ramda / med[mednum].omega;//第2粘性定数/////////////////////////////////////
        break;
      default:
        break;
    }
  }
}

void initCoord(Coord *co, int x, int y, int z) {
  co->x = x;
  co->y = y;
  co->z = z;
}

//差分間隔
void initDiff(Diff *dif, Medium *med) {
  dif->dx = 0.005;
  dif->dy = 0.005;
  dif->dz = 0.005;
  double tmp;
  
  for(int i = E_AIR; i < E_M_END - 1; i++){
    tmp = MAX(sqrt((med[i].K + 4. / 3. * med[i].G) / med[i].rho),tmp);
  }
  printf("v = %lf\n", tmp);
  dif->dt = dif->dx / tmp / 10000.;
}

void initPml(Pml *pml, Medium *med, Diff dif) {
  pml->ta = 4.;
  pml->fm = 3.574e4;
  double R = 1.e-20;
  double tmp,tmp_v;//max
  initCoord(&pml->pl1, 32, 32, 32);
  initCoord(&pml->pl2, 32, 32, 32);
  //計算領域内最高速度
  for(int i = E_AIR; i < E_M_END - 1; i++){
    tmp_v = MAX(sqrt((med[i].K + 4. / 3. * med[i].G) / med[i].rho),tmp);
  }
  //減衰係数最大値(PML層)
  for (int i = E_AIR + 1; i < E_M_END; i++) {
    tmp = tmp_v * (pml->ta + 1) / (2. * (double)pml->pl1.x * dif.dx) * log(1/R);
    pml->fm = MAX(tmp, pml->fm);
  }
}

void initConcrete(Object *con, Medium med, Pml pml, int spx, int spy, int spz, int ranx, int rany, int ranz) {
  con->med = med;
  // 0スタートだから-1
  // 改善できたらしたい
  spx = spx + pml.pl1.x - 1, spy = spy + pml.pl1.y - 1, spz = spz + pml.pl1.z - 1;//ok
  initCoord(&con->sp, spx, spy, spz);
  initCoord(&con->range, ranx, rany, ranz);
}

void initRange(Range *ran, int x, int y, int z, Pml pml) {
  x = x - 1, y = y - 1, z = z - 1;//ok
  initCoord(&ran->sr.Txx, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->sr.Tyy, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->sr.Tzz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->tr.Txy, x + pml.pl1.x + pml.pl2.x - 1, y + pml.pl1.y + pml.pl2.y - 1, z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->tr.Tyz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y - 1, z + pml.pl1.z + pml.pl2.z - 1);
  initCoord(&ran->tr.Tzx, x + pml.pl1.x + pml.pl2.x - 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z - 1);
  initCoord(&ran->vr.Vx , x + pml.pl1.x + pml.pl2.x - 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->vr.Vy , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y - 1, z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->vr.Vz , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z - 1);
} 

void initDeviceSigArr(SigArr *sa, SigRan sr) {
  int i, j;

  // 1,2次元目のメモリ確保では，実際にアクセスするのはhost
  // 3次元目の配列にアクセスするのはdevice
  // このことからmallc,cudaMallocを使い分けている

  // cudaMalloc((void ***)&sa->Txx , (sr.Txx.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&sa->Txxx, (sr.Txx.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&sa->Txxy, (sr.Txx.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&sa->Txxz, (sr.Txx.x + 1) * sizeof(double **));
  sa->Txx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  sa->Txxx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  sa->Txxy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  sa->Txxz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  
  for (i = 0; i < sr.Txx.x + 1; i++) {
    // cudaMalloc((void ***)&sa->Txx[i] , (sr.Txx.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&sa->Txxx[i], (sr.Txx.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&sa->Txxy[i], (sr.Txx.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&sa->Txxz[i], (sr.Txx.y + 1) * sizeof(double *));
    sa->Txx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    sa->Txxx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    sa->Txxy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    sa->Txxz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
  }
  for(i = 0; i < sr.Txx.x + 1; i++){
    // 3次元目のポインタもメモリ割り当て
    for (j = 0; j < sr.Txx.y + 1; j++) {
      cudaMalloc((void ***)&sa->Txx[i][j] , (sr.Txx.z + 1) * sizeof(double));
      cudaMalloc((void ***)&sa->Txxx[i][j], (sr.Txx.z + 1) * sizeof(double));
      cudaMalloc((void ***)&sa->Txxy[i][j], (sr.Txx.z + 1) * sizeof(double));
      cudaMalloc((void ***)&sa->Txxz[i][j], (sr.Txx.z + 1) * sizeof(double));
    }
  }
    
  // cudaMalloc((void ***)&sa->Tyy , (sr.Tyy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&sa->Tyyx, (sr.Tyy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&sa->Tyyy, (sr.Tyy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&sa->Tyyz, (sr.Tyy.x + 1) * sizeof(double **));
  sa->Tyy = (double ***)malloc(sizeof(double **) * (sr.Tyy.x + 1));
  sa->Tyyx = (double ***)malloc(sizeof(double **) * (sr.Tyy.x + 1));
  sa->Tyyy = (double ***)malloc(sizeof(double **) * (sr.Tyy.x + 1));
  sa->Tyyz = (double ***)malloc(sizeof(double **) * (sr.Tyy.x + 1));
  for (i = 0; i < sr.Tyy.x + 1; i++) {
    // cudaMalloc((void ***)&sa->Tyy[i] , (sr.Tyy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&sa->Tyyx[i], (sr.Tyy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&sa->Tyyy[i], (sr.Tyy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&sa->Tyyz[i], (sr.Tyy.y + 1) * sizeof(double *));
    sa->Tyy[i] = (double **)malloc(sizeof(double *) * (sr.Tyy.y + 1));
    sa->Tyyx[i] = (double **)malloc(sizeof(double *) * (sr.Tyy.y + 1));
    sa->Tyyy[i] = (double **)malloc(sizeof(double *) * (sr.Tyy.y + 1));
    sa->Tyyz[i] = (double **)malloc(sizeof(double *) * (sr.Tyy.y + 1));
  }
  for(i = 0; i < sr.Tyy.x + 1; i++){
    // 3次元目のポインタもメモリ割り当て
    for (j = 0; j < sr.Tyy.y + 1; j++) {
      cudaMalloc((void ***)&sa->Tyy[i][j] , (sr.Tyy.z + 1) * sizeof(double));
      cudaMalloc((void ***)&sa->Tyyx[i][j], (sr.Tyy.z + 1) * sizeof(double));
      cudaMalloc((void ***)&sa->Tyyy[i][j], (sr.Tyy.z + 1) * sizeof(double));
      cudaMalloc((void ***)&sa->Tyyz[i][j], (sr.Tyy.z + 1) * sizeof(double));
    }
  }
    
  // cudaMalloc((void ***)&sa->Tzz , (sr.Tzz.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&sa->Tzzx, (sr.Tzz.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&sa->Tzzy, (sr.Tzz.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&sa->Tzzz, (sr.Tzz.x + 1) * sizeof(double **));
  sa->Tzz = (double ***)malloc(sizeof(double **) * (sr.Tzz.x + 1));
  sa->Tzzx = (double ***)malloc(sizeof(double **) * (sr.Tzz.x + 1));
  sa->Tzzy = (double ***)malloc(sizeof(double **) * (sr.Tzz.x + 1));
  sa->Tzzz = (double ***)malloc(sizeof(double **) * (sr.Tzz.x + 1));
  for (i = 0; i < sr.Tzz.x + 1; i++) {
    // cudaMalloc((void ***)&sa->Tzz[i] , (sr.Tzz.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&sa->Tzzx[i], (sr.Tzz.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&sa->Tzzy[i], (sr.Tzz.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&sa->Tzzz[i], (sr.Tzz.y + 1) * sizeof(double *));
    sa->Tzz[i] = (double **)malloc(sizeof(double *) * (sr.Tzz.y + 1));
    sa->Tzzx[i] = (double **)malloc(sizeof(double *) * (sr.Tzz.y + 1));
    sa->Tzzy[i] = (double **)malloc(sizeof(double *) * (sr.Tzz.y + 1));
    sa->Tzzz[i] = (double **)malloc(sizeof(double *) * (sr.Tzz.y + 1));
  }
  for(i = 0; i < sr.Tzz.x + 1; i++){
    // 3次元目のポインタもメモリ割り当て
    for (j = 0; j < sr.Tzz.y + 1; j++) {
      cudaMalloc((void ***)&sa->Tzz[i][j] , (sr.Tzz.z + 1) * sizeof(double));
      cudaMalloc((void ***)&sa->Tzzx[i][j], (sr.Tzz.z + 1) * sizeof(double));
      cudaMalloc((void ***)&sa->Tzzy[i][j], (sr.Tzz.z + 1) * sizeof(double));
      cudaMalloc((void ***)&sa->Tzzz[i][j], (sr.Tzz.z + 1) * sizeof(double));
    }
  }
}

void initHostSigArr(SigArr *sa, SigRan sr) {
  int i, j;
  sa->Txx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  sa->Txxx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  sa->Txxy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  sa->Txxz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  for (i = 0; i <= sr.Txx.x; i++) {
    sa->Txx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    sa->Txxx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    sa->Txxy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    sa->Txxz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
  }
  for (i = 0; i <= sr.Txx.x; i++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      sa->Txx[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      sa->Txxx[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      sa->Txxy[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      sa->Txxz[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
    }
  }
  sa->Tyy = (double ***)malloc(sizeof(double **) * (sr.Tyy.x + 1));
  sa->Tyyx = (double ***)malloc(sizeof(double **) * (sr.Tyy.x + 1));
  sa->Tyyy = (double ***)malloc(sizeof(double **) * (sr.Tyy.x + 1));
  sa->Tyyz = (double ***)malloc(sizeof(double **) * (sr.Tyy.x + 1));
  for (i = 0; i <= sr.Tyy.x; i++) {
    sa->Tyy[i] = (double **)malloc(sizeof(double *) * (sr.Tyy.y + 1));
    sa->Tyyx[i] = (double **)malloc(sizeof(double *) * (sr.Tyy.y + 1));
    sa->Tyyy[i] = (double **)malloc(sizeof(double *) * (sr.Tyy.y + 1));
    sa->Tyyz[i] = (double **)malloc(sizeof(double *) * (sr.Tyy.y + 1));
  }
  for (i = 0; i <= sr.Tyy.x; i++) {
    for (j = 0; j <= sr.Tyy.y; j++) {
      sa->Tyy[i][j] = (double *)malloc(sizeof(double) * (sr.Tyy.z + 1));
      sa->Tyyx[i][j] = (double *)malloc(sizeof(double) * (sr.Tyy.z + 1));
      sa->Tyyy[i][j] = (double *)malloc(sizeof(double) * (sr.Tyy.z + 1));
      sa->Tyyz[i][j] = (double *)malloc(sizeof(double) * (sr.Tyy.z + 1));
    }
  }
  sa->Tzz = (double ***)malloc(sizeof(double **) * (sr.Tzz.x + 1));
  sa->Tzzx = (double ***)malloc(sizeof(double **) * (sr.Tzz.x + 1));
  sa->Tzzy = (double ***)malloc(sizeof(double **) * (sr.Tzz.x + 1));
  sa->Tzzz = (double ***)malloc(sizeof(double **) * (sr.Tzz.x + 1));
  for (i = 0; i <= sr.Tzz.x; i++) {
    sa->Tzz[i] = (double **)malloc(sizeof(double *) * (sr.Tzz.y + 1));
    sa->Tzzx[i] = (double **)malloc(sizeof(double *) * (sr.Tzz.y + 1));
    sa->Tzzy[i] = (double **)malloc(sizeof(double *) * (sr.Tzz.y + 1));
    sa->Tzzz[i] = (double **)malloc(sizeof(double *) * (sr.Tzz.y + 1));
  }
  for (i = 0; i <= sr.Tzz.x; i++) {
    for (j = 0; j <= sr.Tzz.y; j++) {
      sa->Tzz[i][j] = (double *)malloc(sizeof(double) * (sr.Tzz.z + 1));
      sa->Tzzx[i][j] = (double *)malloc(sizeof(double) * (sr.Tzz.z + 1));
      sa->Tzzy[i][j] = (double *)malloc(sizeof(double) * (sr.Tzz.z + 1));
      sa->Tzzz[i][j] = (double *)malloc(sizeof(double) * (sr.Tzz.z + 1));
    }
  }
}

void initDeviceTauArr(TauArr *ta, TauRan tr) {
  int i, j;
  // 1次元目のポインタをデバイスメモリに割り当て (X方向)
  // cudaMalloc((void ***)ta->Txy , (tr.Txy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)ta->Txyx, (tr.Txy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)ta->Txyy, (tr.Txy.x + 1) * sizeof(double **));
  ta->Txy = (double ***)malloc(sizeof(double **) * (tr.Txy.x + 1));
  ta->Txyx = (double ***)malloc(sizeof(double **) * (tr.Txy.x + 1));
  ta->Txyy = (double ***)malloc(sizeof(double **) * (tr.Txy.x + 1));
  // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
  for (i = 0; i <= tr.Txy.x; i++) {
    // cudaMalloc((void ***)ta->Txy, (tr.Txy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)ta->Txyx, (tr.Txy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)ta->Txyy, (tr.Txy.y + 1) * sizeof(double *));
    ta->Txy[i] = (double **)malloc(sizeof(double *) * (tr.Txy.y + 1));
    ta->Txyx[i] = (double **)malloc(sizeof(double *) * (tr.Txy.y + 1));
    ta->Txyy[i] = (double **)malloc(sizeof(double *) * (tr.Txy.y + 1));
    // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
    for (j = 0; j <= tr.Txy.y; j++) {
      cudaMalloc((void ***)ta->Txy, (tr.Txy.z + 1) * sizeof(double));
      cudaMalloc((void ***)ta->Txyx, (tr.Txy.z + 1) * sizeof(double));
      cudaMalloc((void ***)ta->Txyy, (tr.Txy.z + 1) * sizeof(double));
    }
  }
    
  // 1次元目のポインタをデバイスメモリに割り当て (X方向)
  // cudaMalloc((void ***)ta->Tyz , (tr.Txy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)ta->Tyzy, (tr.Txy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)ta->Tyzz, (tr.Txy.x + 1) * sizeof(double **));
  ta->Tyz = (double ***)malloc(sizeof(double **) * (tr.Tyz.x + 1));
  ta->Tyzy = (double ***)malloc(sizeof(double **) * (tr.Tyz.x + 1));
  ta->Tyzz = (double ***)malloc(sizeof(double **) * (tr.Tyz.x + 1));
  // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
  for (int i = 0; i <= tr.Txy.x; i++) {
    // cudaMalloc((void ***)ta->Tyz , (tr.Txy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)ta->Tyzy, (tr.Txy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)ta->Tyzz, (tr.Txy.y + 1) * sizeof(double *));
    ta->Tyz[i] = (double **)malloc(sizeof(double *) * (tr.Tyz.y + 1));
    ta->Tyzy[i] = (double **)malloc(sizeof(double *) * (tr.Tyz.y + 1));
    ta->Tyzz[i] = (double **)malloc(sizeof(double *) * (tr.Tyz.y + 1));
    // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
    for (int j = 0; j <= tr.Txy.y; j++) {
      cudaMalloc((void ***)ta->Tyz , (tr.Txy.z + 1) * sizeof(double));
      cudaMalloc((void ***)ta->Tyzy, (tr.Txy.z + 1) * sizeof(double));
      cudaMalloc((void ***)ta->Tyzz, (tr.Txy.z + 1) * sizeof(double));
    }
  }

  // 1次元目のポインタをデバイスメモリに割り当て (X方向)
  // cudaMalloc((void ***)ta->Tzx , (tr.Txy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)ta->Tzxz, (tr.Txy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)ta->Tzxx, (tr.Txy.x + 1) * sizeof(double **));
  ta->Tzx = (double ***)malloc(sizeof(double **) * (tr.Tzx.x + 1));
  ta->Tzxz = (double ***)malloc(sizeof(double **) * (tr.Tzx.x + 1));
  ta->Tzxx = (double ***)malloc(sizeof(double **) * (tr.Tzx.x + 1));
  // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
  for (int i = 0; i <= tr.Txy.x; i++) {
    // cudaMalloc((void ***)ta->Tzx, (tr.Txy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)ta->Tzxz, (tr.Txy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)ta->Tzxx, (tr.Txy.y + 1) * sizeof(double *));
    ta->Tzx[i] = (double **)malloc(sizeof(double *) * (tr.Tzx.y + 1));
    ta->Tzxz[i] = (double **)malloc(sizeof(double *) * (tr.Tzx.y + 1));
    ta->Tzxx[i] = (double **)malloc(sizeof(double *) * (tr.Tzx.y + 1));
    // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
    for (int j = 0; j <= tr.Txy.y; j++) {
      cudaMalloc((void ***)ta->Tzx, (tr.Txy.z + 1) * sizeof(double));
      cudaMalloc((void ***)ta->Tzxz, (tr.Txy.z + 1) * sizeof(double));
      cudaMalloc((void ***)ta->Tzxx, (tr.Txy.z + 1) * sizeof(double));
    }
  }
}

void initHostTauArr(TauArr *ta, TauRan tr) {
  int i, j;
  ta->Txy = (double ***)malloc(sizeof(double **) * (tr.Txy.x + 1));
  ta->Txyx = (double ***)malloc(sizeof(double **) * (tr.Txy.x + 1));
  ta->Txyy = (double ***)malloc(sizeof(double **) * (tr.Txy.x + 1));
  for (i = 0; i <= tr.Txy.x; i++) {
    ta->Txy[i] = (double **)malloc(sizeof(double *) * (tr.Txy.y + 1));
    ta->Txyx[i] = (double **)malloc(sizeof(double *) * (tr.Txy.y + 1));
    ta->Txyy[i] = (double **)malloc(sizeof(double *) * (tr.Txy.y + 1));
  }
  for (i = 0; i <= tr.Txy.x; i++) {
    for (j = 0; j <= tr.Txy.y; j++) {
      ta->Txy[i][j] = (double *)malloc(sizeof(double) * (tr.Txy.z + 1));
      ta->Txyx[i][j] = (double *)malloc(sizeof(double) * (tr.Txy.z + 1));
      ta->Txyy[i][j] = (double *)malloc(sizeof(double) * (tr.Txy.z + 1));
    }
  }
  ta->Tyz = (double ***)malloc(sizeof(double **) * (tr.Tyz.x + 1));
  ta->Tyzy = (double ***)malloc(sizeof(double **) * (tr.Tyz.x + 1));
  ta->Tyzz = (double ***)malloc(sizeof(double **) * (tr.Tyz.x + 1));
  for (i = 0; i <= tr.Tyz.x; i++) {
    ta->Tyz[i] = (double **)malloc(sizeof(double *) * (tr.Tyz.y + 1));
    ta->Tyzy[i] = (double **)malloc(sizeof(double *) * (tr.Tyz.y + 1));
    ta->Tyzz[i] = (double **)malloc(sizeof(double *) * (tr.Tyz.y + 1));
  }
  for (i = 0; i <= tr.Tyz.x; i++) {
    for (j = 0; j <= tr.Tyz.y; j++) {
      ta->Tyz[i][j] = (double *)malloc(sizeof(double) * (tr.Tyz.z + 1));
      ta->Tyzy[i][j] = (double *)malloc(sizeof(double) * (tr.Tyz.z + 1));
      ta->Tyzz[i][j] = (double *)malloc(sizeof(double) * (tr.Tyz.z + 1));
    }
  }
  ta->Tzx = (double ***)malloc(sizeof(double **) * (tr.Tzx.x + 1));
  ta->Tzxz = (double ***)malloc(sizeof(double **) * (tr.Tzx.x + 1));
  ta->Tzxx = (double ***)malloc(sizeof(double **) * (tr.Tzx.x + 1));
  for (i = 0; i <= tr.Tzx.x; i++) {
    ta->Tzx[i] = (double **)malloc(sizeof(double *) * (tr.Tzx.y + 1));
    ta->Tzxz[i] = (double **)malloc(sizeof(double *) * (tr.Tzx.y + 1));
    ta->Tzxx[i] = (double **)malloc(sizeof(double *) * (tr.Tzx.y + 1));
  }
  for (i = 0; i <= tr.Tzx.x; i++) {
    for (j = 0; j <= tr.Tzx.y; j++) {
      ta->Tzx[i][j] = (double *)malloc(sizeof(double) * (tr.Tzx.z + 1));
      ta->Tzxz[i][j] = (double *)malloc(sizeof(double) * (tr.Tzx.z + 1));
      ta->Tzxx[i][j] = (double *)malloc(sizeof(double) * (tr.Tzx.z + 1));
    }
  }
}

void initDeviceVelArr(VelArr *va, VelRan vr) {
  int i, j;

  // 1次元目のポインタをデバイスメモリに割り当て (X方向)
  // cudaMalloc((void ***)&va->Vx , (vr.Vx.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&va->Vxx, (vr.Vx.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&va->Vxy, (vr.Vx.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&va->Vxz, (vr.Vx.x + 1) * sizeof(double **));
  va->Vx = (double ***)malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxx = (double ***)malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxy = (double ***)malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxz = (double ***)malloc(sizeof(double **) * (vr.Vx.x + 1));
  // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
  for (i = 0; i <= vr.Vx.x; i++) {
    // cudaMalloc((void ***)&va->Vx[i] , (vr.Vx.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&va->Vxx[i], (vr.Vx.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&va->Vxy[i], (vr.Vx.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&va->Vxz[i], (vr.Vx.y + 1) * sizeof(double *));
    va->Vx[i] = (double **)malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxx[i] = (double **)malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxy[i] = (double **)malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxz[i] = (double **)malloc(sizeof(double *) * (vr.Vx.y + 1));
    
    // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
    for (j = 0; j <= vr.Vx.y; j++) {
    cudaMalloc((void ***)&va->Vx[i][j] , (vr.Vx.z + 1) * sizeof(double));
    cudaMalloc((void ***)&va->Vxx[i][j], (vr.Vx.z + 1) * sizeof(double));
    cudaMalloc((void ***)&va->Vxy[i][j], (vr.Vx.z + 1) * sizeof(double));
    cudaMalloc((void ***)&va->Vxz[i][j], (vr.Vx.z + 1) * sizeof(double));
    }
  }

  // 1次元目のポインタをデバイスメモリに割り当て (X方向)
  // cudaMalloc((void ***)&va->Vy , (vr.Vy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&va->Vyx, (vr.Vy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&va->Vyy, (vr.Vy.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&va->Vyz, (vr.Vy.x + 1) * sizeof(double **));
  va->Vy = (double ***)malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyx = (double ***)malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyy = (double ***)malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyz = (double ***)malloc(sizeof(double **) * (vr.Vy.x + 1));
  // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
  for (i = 0; i <= vr.Vy.x; i++) {
    // cudaMalloc((void ***)&va->Vy[i] , (vr.Vy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&va->Vyx[i], (vr.Vy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&va->Vyy[i], (vr.Vy.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&va->Vyz[i], (vr.Vy.y + 1) * sizeof(double *));
    va->Vy[i] = (double **)malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyx[i] = (double **)malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyy[i] = (double **)malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyz[i] = (double **)malloc(sizeof(double *) * (vr.Vy.y + 1));
    // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
    for (j = 0; j <= vr.Vy.y; j++) {
    cudaMalloc((void ***)&va->Vy[i][j] , (vr.Vy.z + 1) * sizeof(double));
    cudaMalloc((void ***)&va->Vyx[i][j], (vr.Vy.z + 1) * sizeof(double));
    cudaMalloc((void ***)&va->Vyy[i][j], (vr.Vy.z + 1) * sizeof(double));
    cudaMalloc((void ***)&va->Vyz[i][j], (vr.Vy.z + 1) * sizeof(double));
    }
  }
  // 1次元目のポインタをデバイスメモリに割り当て (X方向)
  // cudaMalloc((void ***)&va->Vz , (vr.Vz.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&va->Vzx, (vr.Vz.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&va->Vzy, (vr.Vz.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&va->Vzz, (vr.Vz.x + 1) * sizeof(double **));
  va->Vz = (double ***)malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzx = (double ***)malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzy = (double ***)malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzz = (double ***)malloc(sizeof(double **) * (vr.Vz.x + 1));
  // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
  for (i = 0; i <= vr.Vz.x; i++) {
    // cudaMalloc((void ***)&va->Vz[i] , (vr.Vz.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&va->Vzx[i], (vr.Vz.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&va->Vzy[i], (vr.Vz.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&va->Vzz[i], (vr.Vz.y + 1) * sizeof(double *));
    va->Vz[i] = (double **)malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzx[i] = (double **)malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzy[i] = (double **)malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzz[i] = (double **)malloc(sizeof(double *) * (vr.Vz.y + 1));
    // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
    for (j = 0; j <= vr.Vz.y; j++) {
    cudaMalloc((void ***)&va->Vz[i][j] , (vr.Vz.z + 1) * sizeof(double));
    cudaMalloc((void ***)&va->Vzx[i][j], (vr.Vz.z + 1) * sizeof(double));
    cudaMalloc((void ***)&va->Vzy[i][j], (vr.Vz.z + 1) * sizeof(double));
    cudaMalloc((void ***)&va->Vzz[i][j], (vr.Vz.z + 1) * sizeof(double));
    }
  }
}

void initHostVelArr(VelArr *va, VelRan vr) {
  int i, j;
  va->Vx = (double ***)malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxx = (double ***)malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxy = (double ***)malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxz = (double ***)malloc(sizeof(double **) * (vr.Vx.x + 1));
  for (i = 0; i <= vr.Vx.x; i++) {
    va->Vx[i] = (double **)malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxx[i] = (double **)malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxy[i] = (double **)malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxz[i] = (double **)malloc(sizeof(double *) * (vr.Vx.y + 1));
  }
  for (i = 0; i <= vr.Vx.x; i++) {
    for (j = 0; j <= vr.Vx.y; j++) {
      va->Vx[i][j] = (double *)malloc(sizeof(double) * (vr.Vx.z + 1));
      va->Vxx[i][j] = (double *)malloc(sizeof(double) * (vr.Vx.z + 1));
      va->Vxy[i][j] = (double *)malloc(sizeof(double) * (vr.Vx.z + 1));
      va->Vxz[i][j] = (double *)malloc(sizeof(double) * (vr.Vx.z + 1));
    }
  }
  va->Vy = (double ***)malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyx = (double ***)malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyy = (double ***)malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyz = (double ***)malloc(sizeof(double **) * (vr.Vy.x + 1));
  for (i = 0; i <= vr.Vy.x; i++) {
    va->Vy[i] = (double **)malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyx[i] = (double **)malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyy[i] = (double **)malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyz[i] = (double **)malloc(sizeof(double *) * (vr.Vy.y + 1));
  }
  for (i = 0; i <= vr.Vy.x; i++) {
    for (j = 0; j <= vr.Vy.y; j++) {
      va->Vy[i][j] = (double *)malloc(sizeof(double) * (vr.Vy.z + 1));
      va->Vyx[i][j] = (double *)malloc(sizeof(double) * (vr.Vy.z + 1));
      va->Vyy[i][j] = (double *)malloc(sizeof(double) * (vr.Vy.z + 1));
      va->Vyz[i][j] = (double *)malloc(sizeof(double) * (vr.Vy.z + 1));
    }
  }
  va->Vz = (double ***)malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzx = (double ***)malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzy = (double ***)malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzz = (double ***)malloc(sizeof(double **) * (vr.Vz.x + 1));
  for (i = 0; i <= vr.Vz.x; i++) {
    va->Vz[i] = (double **)malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzx[i] = (double **)malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzy[i] = (double **)malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzz[i] = (double **)malloc(sizeof(double *) * (vr.Vz.y + 1));
  }
  for (i = 0; i <= vr.Vz.x; i++) {
    for (j = 0; j <= vr.Vz.y; j++) {
      va->Vz[i][j] = (double *)malloc(sizeof(double) * (vr.Vz.z + 1));
      va->Vzx[i][j] = (double *)malloc(sizeof(double) * (vr.Vz.z + 1));
      va->Vzy[i][j] = (double *)malloc(sizeof(double) * (vr.Vz.z + 1));
      va->Vzz[i][j] = (double *)malloc(sizeof(double) * (vr.Vz.z + 1));
    }
  }
}

void initDeviceBefAft(BefAft *ba, Range ran) {
  initDeviceSigArr(&ba->sa, ran.sr);
  initDeviceTauArr(&ba->ta, ran.tr);
  initDeviceVelArr(&ba->va, ran.vr);
}

void initHostBefAft(BefAft *ba, Range ran) {
  initHostSigArr(&ba->sa, ran.sr);
  initHostTauArr(&ba->ta, ran.tr);
  initHostVelArr(&ba->va, ran.vr);
}

void initDeviceInpalse(Inpaluse *ip, SigRan sr, Pml pml, int mode, int x, int y, int z, double freq) {
  ip->freq = freq;
  ip->mode = mode;
  // int i, j, k;
  // initCoord(&ip->in, x + pml.pl1.x - 1, y + pml.pl1.y - 1, z + pml.pl1.z - 1);//ok
  // 1次元目のポインタをデバイスメモリに割り当て (X方向)
  // cudaMalloc((void ***)&ip->Txx, (sr.Txx.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&ip->Tyy, (sr.Txx.x + 1) * sizeof(double **));
  // cudaMalloc((void ***)&ip->Tzz, (sr.Txx.x + 1) * sizeof(double **));
  ip->Txx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ip->Tyy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ip->Tzz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
  for (int i = 0; i <= sr.Txx.x; i++) {
    // cudaMalloc((void ***)&ip->Txx[i], (sr.Txx.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&ip->Tyy[i], (sr.Txx.y + 1) * sizeof(double *));
    // cudaMalloc((void ***)&ip->Tzz[i], (sr.Txx.y + 1) * sizeof(double *));
    ip->Txx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ip->Tyy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ip->Tzz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
    for (int j = 0; j <= sr.Txx.y; j++) {
      cudaMalloc((void ***)&ip->Txx[i][j], (sr.Txx.z + 1) * sizeof(double));
      cudaMalloc((void ***)&ip->Tyy[i][j], (sr.Txx.z + 1) * sizeof(double));
      cudaMalloc((void ***)&ip->Tzz[i][j], (sr.Txx.z + 1) * sizeof(double));
    }
  }
  // for (k = 0; k <= sr.Txx.z; k++) {
  //   for (j = 0; j <= sr.Txx.y; j++) {
  //     for (i = 0; i <= sr.Txx.x; i++) {
  //       ip->Txx[i][j][k] = 0.;
  //       ip->Tyy[i][j][k] = 0.;
  //       ip->Tzz[i][j][k] = 0.;
  //     }
  //   }
  // }  
}

void initHostInpalse(Inpaluse *ip, SigRan sr, Pml pml, int mode, int x, int y, int z, double freq) {
  ip->freq = freq;
  ip->mode = mode;
  int i, j, k;
  // initCoord(&ip->in, x + pml.pl1.x - 1, y + pml.pl1.y - 1, z + pml.pl1.z - 1);//ok
  ip->Txx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ip->Tyy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ip->Tzz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  for (i = 0; i <= sr.Txx.x; i++) {
    ip->Txx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ip->Tyy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ip->Tzz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
  }
  for (i = 0; i <= sr.Txx.x; i++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      ip->Txx[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ip->Tyy[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ip->Tzz[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
    }
  }
  for (k = 0; k <= sr.Txx.z; k++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      for (i = 0; i <= sr.Txx.x; i++) {
        ip->Txx[i][j][k] = 0.;
        ip->Tyy[i][j][k] = 0.;
        ip->Tzz[i][j][k] = 0.;
      }
    }
  }  
}

void initDeviceMedArr(MedArr *ma, SigRan sr) {
  int i, j;
  // cudaMalloc((void ***)&ma->ramda, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->mu, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->c11, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->rho, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetaxx, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetaxy, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetaxz, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetayx, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetayy, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetayz, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetazx, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetazy, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetazz, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->gamma, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->khi, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->xi11, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetadx, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetady, sizeof(double **) * (sr.Txx.x + 1));
  // cudaMalloc((void ***)&ma->zetadz, sizeof(double **) * (sr.Txx.x + 1));
  ma->ramda = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->mu = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->c11 = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->rho = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->gamma = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->khi = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->xi11 = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetadx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetady = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetadz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  for (i = 0; i <= sr.Txx.x; i++) {
    // cudaMalloc((void ***)&ma->ramda[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->mu[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->c11[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->rho[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetaxx[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetaxy[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetaxz[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetayx[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetayy[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetayz[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetazx[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetazy[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetazz[i], sizeof(double **) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->gamma[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->khi[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->xi11[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetadx[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetady[i], sizeof(double *) * (sr.Txx.x + 1));
    // cudaMalloc((void ***)&ma->zetadz[i], sizeof(double *) * (sr.Txx.x + 1));
    ma->ramda[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->mu[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->c11[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->rho[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->gamma[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->khi[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->xi11[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetadx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetady[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetadz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    for (j = 0; j <= sr.Txx.y; j++) {
      cudaMalloc((void **)&ma->ramda[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->mu[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->c11[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->rho[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetaxx[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetaxy[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetaxz[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetayx[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetayy[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetayz[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetazx[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetazy[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetazz[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->gamma[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->khi[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->xi11[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetadx[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetady[i][j], sizeof(double) * (sr.Txx.x + 1));
      cudaMalloc((void ***)&ma->zetadz[i][j], sizeof(double) * (sr.Txx.x + 1));
    }
  }
}

void initHostMedArr(MedArr *ma, SigRan sr) {
  int i, j;
  ma->ramda = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->mu = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->c11 = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->rho = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazy = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->gamma = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->khi = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->xi11 = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetadx = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetady = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetadz = (double ***)malloc(sizeof(double **) * (sr.Txx.x + 1));
  for (i = 0; i <= sr.Txx.x; i++) {
    ma->ramda[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->mu[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->c11[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->rho[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazy[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->gamma[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->khi[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->xi11[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetadx[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetady[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetadz[i] = (double **)malloc(sizeof(double *) * (sr.Txx.y + 1));
  }
  for (i = 0; i <= sr.Txx.x; i++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      ma->ramda[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->mu[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->c11[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->rho[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetaxx[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetaxy[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetaxz[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetayx[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetayy[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetayz[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetazx[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetazy[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetazz[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->gamma[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->khi[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->xi11[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetadx[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetady[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetadz[i][j] = (double *)malloc(sizeof(double) * (sr.Txx.z + 1));
    }
  }
}

void initrandom(Coord con_size, Coord *clack, int ratio) {
  if(con_size.x < 3 || con_size.y < 3 || con_size.z < 3) {
    printf("Cannot place defects.\n");
    return;
  }
  int count = 0;
  // コンクリートセル数
  int max_Patern = con_size.x * con_size.y * con_size.z;
  // 内部欠陥パターン数
  int max_ClackPatern = (con_size.x - 2) * (con_size.y - 2) * (con_size.z - 2);
  // 割合による欠陥数
  int clack_count = max_Patern * ratio / 100;
  if(clack_count > max_ClackPatern){
    printf("The number of internal defects is insufficient.\n");
    return;
  }
  
  // 乱数の種を初期化
  srand(time(NULL));

  while (count < clack_count) {
    // 新しい乱数の組み合わせを生成
    int rand1 = rand() % (con_size.x - 2) + 1;
    int rand2 = rand() % (con_size.y - 2) + 1;
    int rand3 = rand() % (con_size.z - 2) + 1;

    // 重複がないかチェック
    int is_unique = 1;
    for (int i = 0; i < count; i++) {
      if (clack[i].x == rand1 && clack[i].y == rand2 && clack[i].z == rand3) {
        is_unique = 0;
        break;
      }
    }

    // 重複がなければ保存
    if (is_unique) {
      clack[count].x = rand1;
      clack[count].y = rand2;
      clack[count].z = rand3;
      count++;
    }
  }
}

void initClack(Object *clack, Medium med, Pml *pml, int spx, int spy, int spz, int ranx, int rany, int ranz) {
  clack->med = med;
  spx = spx + pml->pl1.x - 1, spy = spy + pml->pl1.y - 1, spz = spz + pml->pl1.z - 1;//ok
  initCoord(&clack->sp, spx, spy, spz);
  initCoord(&clack->range, ranx, rany, ranz);
}