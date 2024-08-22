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
  spx = spx + pml.pl1.x - 1, spy = spy + pml.pl1.y - 1, spz = spz + pml.pl1.z - 1;//ok
  initCoord(&con->sp, spx, spy, spz);
  initCoord(&con->range, ranx, rany, ranz);
}
/////////steeeeel
void initSteel(Object *steel, Medium med, Pml *pml, int spx, int spy, int spz, int ranx, int rany, int ranz) {
  steel->med = med;
  spx = spx + pml->pl1.x - 1, spy = spy + pml->pl1.y - 1, spz = spz + pml->pl1.z - 1;//ok
  initCoord(&steel->sp, spx, spy, spz);
  initCoord(&steel->range, ranx, rany, ranz);
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

void initSigArr(SigArr *sa, SigRan sr) {
    int i, j, k;
    cudaMalloc((void ***)&sa->Txx , (sr.Txx.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&sa->Txxx, (sr.Txx.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&sa->Txxy, (sr.Txx.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&sa->Txxz, (sr.Txx.x + 1) * sizeof(double **));
    
    for (int i = 0; i < sr.Txx.x + 1; ++i) {
        cudaMalloc((void ***)&sa->Txx[i] , (sr.Txx.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&sa->Txxx[i], (sr.Txx.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&sa->Txxy[i], (sr.Txx.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&sa->Txxz[i], (sr.Txx.y + 1) * sizeof(double *));

        cudaMemcpy(&(sa->Txx[i]), sa->Txx, sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(sa->Txxx[i]), sa->Txxx, sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(sa->Txxy[i]), sa->Txxy, sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(sa->Txxz[i]), sa->Txxz, sizeof(double **), cudaMemcpyHostToDevice);

        // 3次元目のポインタもメモリ割り当て
        for (int j = 0; j < sr.Txx.y + 1; ++j) {
            cudaMalloc((void ***)&sa->Txx[i][j] , (sr.Txx.z + 1) * sizeof(double));
            cudaMalloc((void ***)&sa->Txxx[i][j], (sr.Txx.z + 1) * sizeof(double));
            cudaMalloc((void ***)&sa->Txxy[i][j], (sr.Txx.z + 1) * sizeof(double));
            cudaMalloc((void ***)&sa->Txxz[i][j], (sr.Txx.z + 1) * sizeof(double));

            cudaMemcpy(&(sa->Txx[i][j]), sa->Txx, sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(sa->Txxx[i][j]), sa->Txxx, sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(sa->Txxy[i][j]), sa->Txxy, sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(sa->Txxz[i][j]), sa->Txxz, sizeof(double *), cudaMemcpyHostToDevice);
        }
    }
    
    cudaMalloc((void ***)&sa->Tyy , (sr.Tyy.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&sa->Tyyx, (sr.Tyy.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&sa->Tyyy, (sr.Tyy.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&sa->Tyyz, (sr.Tyy.x + 1) * sizeof(double **));

    for (int i = 0; i < sr.Tyy.x + 1; ++i) {
        cudaMalloc((void ***)&sa->Tyy[i] , (sr.Tyy.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&sa->Tyyx[i], (sr.Tyy.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&sa->Tyyy[i], (sr.Tyy.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&sa->Tyyz[i], (sr.Tyy.y + 1) * sizeof(double *));
        // 3次元目のポインタもメモリ割り当て
        for (int j = 0; j < sr.Tyy.y + 1; ++j) {
            cudaMalloc((void ***)&sa->Tyy[i][j] , (sr.Tyy.z + 1) * sizeof(double));
            cudaMalloc((void ***)&sa->Tyyx[i][j], (sr.Tyy.z + 1) * sizeof(double));
            cudaMalloc((void ***)&sa->Tyyy[i][j], (sr.Tyy.z + 1) * sizeof(double));
            cudaMalloc((void ***)&sa->Tyyz[i][j], (sr.Tyy.z + 1) * sizeof(double));
        }
    }
    
    cudaMalloc((void ***)&sa->Tzz , (sr.Tzz.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&sa->Tzzx, (sr.Tzz.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&sa->Tzzy, (sr.Tzz.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&sa->Tzzz, (sr.Tzz.x + 1) * sizeof(double **));

    for (int i = 0; i < sr.Tzz.x + 1; ++i) {
        cudaMalloc((void ***)&sa->Tzz[i] , (sr.Tzz.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&sa->Tzzx[i], (sr.Tzz.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&sa->Tzzy[i], (sr.Tzz.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&sa->Tzzz[i], (sr.Tzz.y + 1) * sizeof(double *));
        // 3次元目のポインタもメモリ割り当て
        for (int j = 0; j < sr.Tzz.y + 1; ++j) {
            cudaMalloc((void ***)&sa->Tzz[i][j] , (sr.Tzz.z + 1) * sizeof(double));
            cudaMalloc((void ***)&sa->Tzzx[i][j], (sr.Tzz.z + 1) * sizeof(double));
            cudaMalloc((void ***)&sa->Tzzy[i][j], (sr.Tzz.z + 1) * sizeof(double));
            cudaMalloc((void ***)&sa->Tzzz[i][j], (sr.Tzz.z + 1) * sizeof(double));
        }
    }
}

void initTauArr(TauArr *ta, TauRan tr) {
    int i, j, k;
    // 1次元目のポインタをデバイスメモリに割り当て (X方向)
    cudaMalloc((void ***)ta->Txy , (tr.Txy.x + 1) * sizeof(double **));
    cudaMalloc((void ***)ta->Txyx, (tr.Txy.x + 1) * sizeof(double **));
    cudaMalloc((void ***)ta->Txyy, (tr.Txy.x + 1) * sizeof(double **));

    // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
    for (int i = 0; i <= tr.Txy.x; i++) {
        cudaMalloc((void ***)ta->Txy, (tr.Txy.y + 1) * sizeof(double *));
        cudaMalloc((void ***)ta->Txyx, (tr.Txy.y + 1) * sizeof(double *));
        cudaMalloc((void ***)ta->Txyy, (tr.Txy.y + 1) * sizeof(double *));
        
        cudaMemcpy(&(ta->Txy[i]), ta->Txy, sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(ta->Txyx[i]), ta->Txyx, sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(ta->Txyy[i]), ta->Txyy, sizeof(double **), cudaMemcpyHostToDevice);
        
        // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
        for (int j = 0; j <= tr.Txy.y; j++) {
            cudaMalloc((void ***)ta->Txy, (tr.Txy.z + 1) * sizeof(double));
            cudaMalloc((void ***)ta->Txyx, (tr.Txy.z + 1) * sizeof(double));
            cudaMalloc((void ***)ta->Txyy, (tr.Txy.z + 1) * sizeof(double));
            
            cudaMemcpy(&(ta->Txy[i][j]), ta->Txy, sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(ta->Txyx[i][j]), ta->Txyx, sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(ta->Txyy[i][j]), ta->Txyy, sizeof(double *), cudaMemcpyHostToDevice);
        }
    }
    
    // 1次元目のポインタをデバイスメモリに割り当て (X方向)
    cudaMalloc((void ***)ta->Tyz , (tr.Txy.x + 1) * sizeof(double **));
    cudaMalloc((void ***)ta->Tyzy, (tr.Txy.x + 1) * sizeof(double **));
    cudaMalloc((void ***)ta->Tyzz, (tr.Txy.x + 1) * sizeof(double **));

    // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
    for (int i = 0; i <= tr.Txy.x; i++) {
        cudaMalloc((void ***)ta->Tyz , (tr.Txy.y + 1) * sizeof(double *));
        cudaMalloc((void ***)ta->Tyzy, (tr.Txy.y + 1) * sizeof(double *));
        cudaMalloc((void ***)ta->Tyzz, (tr.Txy.y + 1) * sizeof(double *));
        
        cudaMemcpy(&(ta->Tyz[i]) , ta->Tyz , sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(ta->Tyzy[i]), ta->Tyzy, sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(ta->Tyzz[i]), ta->Tyzz, sizeof(double **), cudaMemcpyHostToDevice);
        
        // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
        for (int j = 0; j <= tr.Txy.y; j++) {
            cudaMalloc((void ***)ta->Tyz , (tr.Txy.z + 1) * sizeof(double));
            cudaMalloc((void ***)ta->Tyzy, (tr.Txy.z + 1) * sizeof(double));
            cudaMalloc((void ***)ta->Tyzz, (tr.Txy.z + 1) * sizeof(double));
            
            cudaMemcpy(&(ta->Tyz[i][j]) , ta->Tyz , sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(ta->Tyzy[i][j]), ta->Tyzy, sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(ta->Tyzz[i][j]), ta->Tyzz, sizeof(double *), cudaMemcpyHostToDevice);
        }
    }

    // 1次元目のポインタをデバイスメモリに割り当て (X方向)
    cudaMalloc((void ***)ta->Tzx , (tr.Txy.x + 1) * sizeof(double **));
    cudaMalloc((void ***)ta->Tzxz, (tr.Txy.x + 1) * sizeof(double **));
    cudaMalloc((void ***)ta->Tzxx, (tr.Txy.x + 1) * sizeof(double **));

    // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
    for (int i = 0; i <= tr.Txy.x; i++) {
        cudaMalloc((void ***)ta->Tzx, (tr.Txy.y + 1) * sizeof(double *));
        cudaMalloc((void ***)ta->Tzxz, (tr.Txy.y + 1) * sizeof(double *));
        cudaMalloc((void ***)ta->Tzxx, (tr.Txy.y + 1) * sizeof(double *));
        
        cudaMemcpy(&(ta->Tzx[i]), ta->Tzx, sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(ta->Tzxz[i]), ta->Tzxz, sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(ta->Tzxx[i]), ta->Tzxx, sizeof(double **), cudaMemcpyHostToDevice);
        
        // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
        for (int j = 0; j <= tr.Txy.y; j++) {
            cudaMalloc((void ***)ta->Tzx, (tr.Txy.z + 1) * sizeof(double));
            cudaMalloc((void ***)ta->Tzxz, (tr.Txy.z + 1) * sizeof(double));
            cudaMalloc((void ***)ta->Tzxx, (tr.Txy.z + 1) * sizeof(double));
            
            cudaMemcpy(&(ta->Tzx[i][j]), ta->Tzx, sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(ta->Tzxz[i][j]), ta->Tzxz, sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(ta->Tzxx[i][j]), ta->Tzxx, sizeof(double *), cudaMemcpyHostToDevice);
        }
    }
}

void initVelArr(VelArr *va, VelRan vr) {
    int i, j, k;

    // 1次元目のポインタをデバイスメモリに割り当て (X方向)
    cudaMalloc((void ***)&va->Vx , (vr.Vx.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&va->Vxx, (vr.Vx.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&va->Vxy, (vr.Vx.x + 1) * sizeof(double **));
    cudaMalloc((void ***)&va->Vxz, (vr.Vx.x + 1) * sizeof(double **));

    // 2次元目のポインタをデバイスメモリに割り当て (Y方向)
    for (int i = 0; i <= vr.Vx.x; i++) {
        cudaMalloc((void ***)&va->Vx[i], (vr.Vx.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&va->Vxx[i], (vr.Vx.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&va->Vxy[i], (vr.Vx.y + 1) * sizeof(double *));
        cudaMalloc((void ***)&va->Vxz[i], (vr.Vx.y + 1) * sizeof(double *));
        
        cudaMemcpy(&( va->Vx[i]),  &va->Vx[i], sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(va->Vxx[i]), &va->Vxx[i], sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(va->Vxy[i]), &va->Vxy[i], sizeof(double **), cudaMemcpyHostToDevice);
        cudaMemcpy(&(va->Vxz[i]), &va->Vxz[i], sizeof(double **), cudaMemcpyHostToDevice);
        
        // 3次元目のポインタをデバイスメモリに割り当て (Z方向)
        for (int j = 0; j <= vr.Vx.y; j++) {
            cudaMalloc((void ***)&va->Vx[i][j] , (vr.Vx.z + 1) * sizeof(double));
            cudaMalloc((void ***)&va->Vxx[i][j], (vr.Vx.z + 1) * sizeof(double));
            cudaMalloc((void ***)&va->Vxy[i][j], (vr.Vx.z + 1) * sizeof(double));
            cudaMalloc((void ***)&va->Vxz[i][j], (vr.Vx.z + 1) * sizeof(double));
            
            cudaMemcpy(&(va->Vx[i][j]) , &va->Vx[i][j] , sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(va->Vxx[i][j]), &va->Vxx[i][j], sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(va->Vxy[i][j]), &va->Vxy[i][j], sizeof(double *), cudaMemcpyHostToDevice);
            cudaMemcpy(&(va->Vxz[i][j]), &va->Vxz[i][j], sizeof(double *), cudaMemcpyHostToDevice);
        }
    }

  va->Vx = malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxx = malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxy = malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxz = malloc(sizeof(double **) * (vr.Vx.x + 1));
  for (i = 0; i <= vr.Vx.x; i++) {
    va->Vx[i] = malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxx[i] = malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxy[i] = malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxz[i] = malloc(sizeof(double *) * (vr.Vx.y + 1));
  }
  for (i = 0; i <= vr.Vx.x; i++) {
    for (j = 0; j <= vr.Vx.y; j++) {
      va->Vx[i][j] = malloc(sizeof(double) * (vr.Vx.z + 1));
      va->Vxx[i][j] = malloc(sizeof(double) * (vr.Vx.z + 1));
      va->Vxy[i][j] = malloc(sizeof(double) * (vr.Vx.z + 1));
      va->Vxz[i][j] = malloc(sizeof(double) * (vr.Vx.z + 1));
    }
  }
  va->Vy = malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyx = malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyy = malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyz = malloc(sizeof(double **) * (vr.Vy.x + 1));
  for (i = 0; i <= vr.Vy.x; i++) {
    va->Vy[i] = malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyx[i] = malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyy[i] = malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyz[i] = malloc(sizeof(double *) * (vr.Vy.y + 1));
  }
  for (i = 0; i <= vr.Vy.x; i++) {
    for (j = 0; j <= vr.Vy.y; j++) {
      va->Vy[i][j] = malloc(sizeof(double) * (vr.Vy.z + 1));
      va->Vyx[i][j] = malloc(sizeof(double) * (vr.Vy.z + 1));
      va->Vyy[i][j] = malloc(sizeof(double) * (vr.Vy.z + 1));
      va->Vyz[i][j] = malloc(sizeof(double) * (vr.Vy.z + 1));
    }
  }
  va->Vz = malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzx = malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzy = malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzz = malloc(sizeof(double **) * (vr.Vz.x + 1));
  for (i = 0; i <= vr.Vz.x; i++) {
    va->Vz[i] = malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzx[i] = malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzy[i] = malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzz[i] = malloc(sizeof(double *) * (vr.Vz.y + 1));
  }
  for (i = 0; i <= vr.Vz.x; i++) {
    for (j = 0; j <= vr.Vz.y; j++) {
      va->Vz[i][j] = malloc(sizeof(double) * (vr.Vz.z + 1));
      va->Vzx[i][j] = malloc(sizeof(double) * (vr.Vz.z + 1));
      va->Vzy[i][j] = malloc(sizeof(double) * (vr.Vz.z + 1));
      va->Vzz[i][j] = malloc(sizeof(double) * (vr.Vz.z + 1));
    }
  }
}

void initBefAft(BefAft *ba, Range ran) {
  initSigArr(&ba->sa, ran.sr);
  initTauArr(&ba->ta, ran.tr);
  initVelArr(&ba->va, ran.vr);
}

void initInpalse(Inpaluse *ip, SigRan sr, Pml pml, int mode, int x, int y, int z, double freq) {
  ip->freq = freq;
  ip->mode = mode;
  int i, j, k;
  // initCoord(&ip->in, x + pml.pl1.x - 1, y + pml.pl1.y - 1, z + pml.pl1.z - 1);//ok
  ip->Txx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ip->Tyy = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ip->Tzz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  for (i = 0; i <= sr.Txx.x; i++) {
    ip->Txx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ip->Tyy[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ip->Tzz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
  }
  for (i = 0; i <= sr.Txx.x; i++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      ip->Txx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ip->Tyy[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ip->Tzz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
    }
  }
  for (k = 0; k <= sr.Txx.z; k++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      for (i = 0; i <= sr.Txx.x; i++) {
        ip->Txx[i][j][k] = 0.;
        ip->Tzz[i][j][k] = 0.;
      }
    }
  }  
}

void initMedArr(MedArr *ma, SigRan sr) {
  int i, j, k;
  ma->ramda = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->mu = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->c11 = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->rho = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxy = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayy = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazy = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->gamma = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->khi = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->xi11 = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetadx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetady = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetadz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  for (i = 0; i <= sr.Txx.x; i++) {
    ma->ramda[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->mu[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->c11[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->rho[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxy[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayy[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazy[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->gamma[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->khi[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->xi11[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetadx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetady[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetadz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
  }
  for (i = 0; i <= sr.Txx.x; i++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      ma->ramda[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->mu[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->c11[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->rho[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetaxx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetaxy[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetaxz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetayx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetayy[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetayz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetazx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetazy[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetazz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->gamma[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->khi[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->xi11[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetadx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetady[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetadz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
    }
  }
}

void initClack(Object *clack, Medium med, Pml *pml, int spx, int spy, int spz, int ranx, int rany, int ranz) {
  clack->med = med;
  spx = spx + pml->pl1.x - 1, spy = spy + pml->pl1.y - 1, spz = spz + pml->pl1.z - 1;//ok
  initCoord(&clack->sp, spx, spy, spz);
  initCoord(&clack->range, ranx, rany, ranz);
}