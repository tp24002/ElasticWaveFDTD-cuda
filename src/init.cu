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
  dif->dt = dif->dx / tmp / 100.;
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
  initCoord(&con->sp, spx, spy, spz);
  initCoord(&con->range, ranx, rany, ranz);
}

void initRange(Range *ran, int x, int y, int z, Pml pml) {
  // initCoord(&ran->sr.Txx, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  // initCoord(&ran->sr.Tyy, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  // initCoord(&ran->sr.Tzz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  // initCoord(&ran->tr.Txy, x + pml.pl1.x + pml.pl2.x - 1, y + pml.pl1.y + pml.pl2.y - 1, z + pml.pl1.z + pml.pl2.z    );
  // initCoord(&ran->tr.Tyz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y - 1, z + pml.pl1.z + pml.pl2.z - 1);
  // initCoord(&ran->tr.Tzx, x + pml.pl1.x + pml.pl2.x - 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z - 1);
  // initCoord(&ran->vr.Vx , x + pml.pl1.x + pml.pl2.x - 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  // initCoord(&ran->vr.Vy , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y - 1, z + pml.pl1.z + pml.pl2.z    );
  // initCoord(&ran->vr.Vz , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z - 1);
  initCoord(&ran->sr.Txx, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->sr.Tyy, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->sr.Tzz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->tr.Txy, x + pml.pl1.x + pml.pl2.x + 1, y + pml.pl1.y + pml.pl2.y + 1, z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->tr.Tyz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y + 1, z + pml.pl1.z + pml.pl2.z + 1);
  initCoord(&ran->tr.Tzx, x + pml.pl1.x + pml.pl2.x + 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z + 1);
  initCoord(&ran->vr.Vx , x + pml.pl1.x + pml.pl2.x + 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->vr.Vy , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y + 1, z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->vr.Vz , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z + 1);
} 

void initHostSigArr(SigArr *sa, SigRan sr) {
  int i, j;
  sa->Txx = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  sa->Txxx = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  sa->Txxy = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  sa->Txxz = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  for (i = 0; i < sr.Txx.x; i++) {
    sa->Txx[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    sa->Txxx[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    sa->Txxy[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    sa->Txxz[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    for (j = 0; j < sr.Txx.y; j++) {
      sa->Txx[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      sa->Txxx[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      sa->Txxy[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      sa->Txxz[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
    }
  }
  sa->Tyy = (double ***)malloc(sizeof(double **) * sr.Tyy.x);
  sa->Tyyx = (double ***)malloc(sizeof(double **) * sr.Tyy.x);
  sa->Tyyy = (double ***)malloc(sizeof(double **) * sr.Tyy.x);
  sa->Tyyz = (double ***)malloc(sizeof(double **) * sr.Tyy.x);
  for (i = 0; i < sr.Tyy.x; i++) {
    sa->Tyy[i] = (double **)malloc(sizeof(double *) * sr.Tyy.y);
    sa->Tyyx[i] = (double **)malloc(sizeof(double *) * sr.Tyy.y);
    sa->Tyyy[i] = (double **)malloc(sizeof(double *) * sr.Tyy.y);
    sa->Tyyz[i] = (double **)malloc(sizeof(double *) * sr.Tyy.y);
    for (j = 0; j < sr.Tyy.y; j++) {
      sa->Tyy[i][j] = (double *)malloc(sizeof(double) * sr.Tyy.x);
      sa->Tyyx[i][j] = (double *)malloc(sizeof(double) * sr.Tyy.x);
      sa->Tyyy[i][j] = (double *)malloc(sizeof(double) * sr.Tyy.x);
      sa->Tyyz[i][j] = (double *)malloc(sizeof(double) * sr.Tyy.x);
    }
  }
  sa->Tzz = (double ***)malloc(sizeof(double **) * sr.Tzz.x);
  sa->Tzzx = (double ***)malloc(sizeof(double **) * sr.Tzz.x);
  sa->Tzzy = (double ***)malloc(sizeof(double **) * sr.Tzz.x);
  sa->Tzzz = (double ***)malloc(sizeof(double **) * sr.Tzz.x);
  for (i = 0; i < sr.Tzz.x; i++) {
    sa->Tzz[i] = (double **)malloc(sizeof(double *) * sr.Tzz.y);
    sa->Tzzx[i] = (double **)malloc(sizeof(double *) * sr.Tzz.y);
    sa->Tzzy[i] = (double **)malloc(sizeof(double *) * sr.Tzz.y);
    sa->Tzzz[i] = (double **)malloc(sizeof(double *) * sr.Tzz.y);
    for (j = 0; j < sr.Tzz.y; j++) {
      sa->Tzz[i][j] = (double *)malloc(sizeof(double) * sr.Tzz.z);
      sa->Tzzx[i][j] = (double *)malloc(sizeof(double) * sr.Tzz.z);
      sa->Tzzy[i][j] = (double *)malloc(sizeof(double) * sr.Tzz.z);
      sa->Tzzz[i][j] = (double *)malloc(sizeof(double) * sr.Tzz.z);
    }
  }
}

void initHostTauArr(TauArr *ta, TauRan tr) {
  int i, j;
  ta->Txy = (double ***)malloc(sizeof(double **) * tr.Txy.x);
  ta->Txyx = (double ***)malloc(sizeof(double **) * tr.Txy.x);
  ta->Txyy = (double ***)malloc(sizeof(double **) * tr.Txy.x);
  for (i = 0; i < tr.Txy.x; i++) {
    ta->Txy[i] = (double **)malloc(sizeof(double *) * tr.Txy.y);
    ta->Txyx[i] = (double **)malloc(sizeof(double *) * tr.Txy.y);
    ta->Txyy[i] = (double **)malloc(sizeof(double *) * tr.Txy.y);
    for (j = 0; j < tr.Txy.y; j++) {
      ta->Txy[i][j] = (double *)malloc(sizeof(double) * tr.Txy.z);
      ta->Txyx[i][j] = (double *)malloc(sizeof(double) * tr.Txy.z);
      ta->Txyy[i][j] = (double *)malloc(sizeof(double) * tr.Txy.z);
    }
  }
  ta->Tyz = (double ***)malloc(sizeof(double **) * tr.Tyz.x);
  ta->Tyzy = (double ***)malloc(sizeof(double **) * tr.Tyz.x);
  ta->Tyzz = (double ***)malloc(sizeof(double **) * tr.Tyz.x);
  for (i = 0; i < tr.Tyz.x; i++) {
    ta->Tyz[i] = (double **)malloc(sizeof(double *) * tr.Tyz.y);
    ta->Tyzy[i] = (double **)malloc(sizeof(double *) * tr.Tyz.y);
    ta->Tyzz[i] = (double **)malloc(sizeof(double *) * tr.Tyz.y);
    for (j = 0; j < tr.Tyz.y; j++) {
      ta->Tyz[i][j] = (double *)malloc(sizeof(double) * tr.Tyz.z);
      ta->Tyzy[i][j] = (double *)malloc(sizeof(double) * tr.Tyz.z);
      ta->Tyzz[i][j] = (double *)malloc(sizeof(double) * tr.Tyz.z);
    }
  }
  ta->Tzx = (double ***)malloc(sizeof(double **) * tr.Tzx.x);
  ta->Tzxz = (double ***)malloc(sizeof(double **) * tr.Tzx.x);
  ta->Tzxx = (double ***)malloc(sizeof(double **) * tr.Tzx.x);
  for (i = 0; i < tr.Tzx.x; i++) {
    ta->Tzx[i] = (double **)malloc(sizeof(double *) * tr.Tzx.y);
    ta->Tzxz[i] = (double **)malloc(sizeof(double *) * tr.Tzx.y);
    ta->Tzxx[i] = (double **)malloc(sizeof(double *) * tr.Tzx.y);
    for (j = 0; j < tr.Tzx.y; j++) {
      ta->Tzx[i][j] = (double *)malloc(sizeof(double) * tr.Tzx.z);
      ta->Tzxz[i][j] = (double *)malloc(sizeof(double) * tr.Tzx.z);
      ta->Tzxx[i][j] = (double *)malloc(sizeof(double) * tr.Tzx.z);
    }
  }
}

void initHostVelArr(VelArr *va, VelRan vr) {
  int i, j;
  va->Vx = (double ***)malloc(sizeof(double **) * vr.Vx.x);
  va->Vxx = (double ***)malloc(sizeof(double **) * vr.Vx.x);
  va->Vxy = (double ***)malloc(sizeof(double **) * vr.Vx.x);
  va->Vxz = (double ***)malloc(sizeof(double **) * vr.Vx.x);
  for (i = 0; i < vr.Vx.x; i++) {
    va->Vx[i] = (double **)malloc(sizeof(double *) * vr.Vx.y);
    va->Vxx[i] = (double **)malloc(sizeof(double *) * vr.Vx.y);
    va->Vxy[i] = (double **)malloc(sizeof(double *) * vr.Vx.y);
    va->Vxz[i] = (double **)malloc(sizeof(double *) * vr.Vx.y);
    for (j = 0; j < vr.Vx.y; j++) {
      va->Vx[i][j] = (double *)malloc(sizeof(double) * vr.Vx.z);
      va->Vxx[i][j] = (double *)malloc(sizeof(double) * vr.Vx.z);
      va->Vxy[i][j] = (double *)malloc(sizeof(double) * vr.Vx.z);
      va->Vxz[i][j] = (double *)malloc(sizeof(double) * vr.Vx.z);
    }
  }
  va->Vy = (double ***)malloc(sizeof(double **) * vr.Vy.x);
  va->Vyx = (double ***)malloc(sizeof(double **) * vr.Vy.x);
  va->Vyy = (double ***)malloc(sizeof(double **) * vr.Vy.x);
  va->Vyz = (double ***)malloc(sizeof(double **) * vr.Vy.x);
  for (i = 0; i < vr.Vy.x; i++) {
    va->Vy[i] = (double **)malloc(sizeof(double *) * vr.Vy.y);
    va->Vyx[i] = (double **)malloc(sizeof(double *) * vr.Vy.y);
    va->Vyy[i] = (double **)malloc(sizeof(double *) * vr.Vy.y);
    va->Vyz[i] = (double **)malloc(sizeof(double *) * vr.Vy.y);
    for (j = 0; j < vr.Vy.y; j++) {
      va->Vy[i][j] = (double *)malloc(sizeof(double) * vr.Vy.z);
      va->Vyx[i][j] = (double *)malloc(sizeof(double) * vr.Vy.z);
      va->Vyy[i][j] = (double *)malloc(sizeof(double) * vr.Vy.z);
      va->Vyz[i][j] = (double *)malloc(sizeof(double) * vr.Vy.z);
    }
  }
  va->Vz = (double ***)malloc(sizeof(double **) * vr.Vz.x);
  va->Vzx = (double ***)malloc(sizeof(double **) * vr.Vz.x);
  va->Vzy = (double ***)malloc(sizeof(double **) * vr.Vz.x);
  va->Vzz = (double ***)malloc(sizeof(double **) * vr.Vz.x);
  for (i = 0; i < vr.Vz.x; i++) {
    va->Vz[i] = (double **)malloc(sizeof(double *) * vr.Vz.y);
    va->Vzx[i] = (double **)malloc(sizeof(double *) * vr.Vz.y);
    va->Vzy[i] = (double **)malloc(sizeof(double *) * vr.Vz.y);
    va->Vzz[i] = (double **)malloc(sizeof(double *) * vr.Vz.y);
    for (j = 0; j < vr.Vz.y; j++) {
      va->Vz[i][j] = (double *)malloc(sizeof(double) * vr.Vz.z);
      va->Vzx[i][j] = (double *)malloc(sizeof(double) * vr.Vz.z);
      va->Vzy[i][j] = (double *)malloc(sizeof(double) * vr.Vz.z);
      va->Vzz[i][j] = (double *)malloc(sizeof(double) * vr.Vz.z);
    }
  }
}


void initHostBefAft(BefAft *ba, Range ran) {
  initHostSigArr(&ba->sa, ran.sr);
  initHostTauArr(&ba->ta, ran.tr);
  initHostVelArr(&ba->va, ran.vr);
}

void initHostInpalse(Inpaluse *ip, SigRan sr, Pml pml, int mode, int x, int y, int z, double freq) {
  ip->freq = freq;
  ip->mode = mode;
  int i, j;
  ip->Txx = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ip->Tyy = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ip->Tzz = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  for (i = 0; i < sr.Txx.x; i++) {
    ip->Txx[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ip->Tyy[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ip->Tzz[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    for (j = 0; j < sr.Txx.y; j++) {
      ip->Txx[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ip->Tyy[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ip->Tzz[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
    }
  }
  // for (k = 0; k < sr.Txx.z; k++) {
  //   for (j = 0; j < sr.Txx.y; j++) {
  //     for (i = 0; i < sr.Txx.x; i++) {
  //       ip->Txx[i][j][k] = 0.;
  //       ip->Tyy[i][j][k] = 0.;
  //       ip->Tzz[i][j][k] = 0.;
  //     }
  //   }
  // }  
}

void initHostMedArr(MedArr *ma, SigRan sr) {
  int i, j;
  ma->ramda = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->mu = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->c11 = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->rho = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetaxx = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetaxy = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetaxz = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetayx = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetayy = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetayz = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetazx = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetazy = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetazz = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->gamma = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->khi = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->xi11 = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetadx = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetady = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  ma->zetadz = (double ***)malloc(sizeof(double **) * sr.Txx.x);
  for (i = 0; i < sr.Txx.x; i++) {
    ma->ramda[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->mu[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->c11[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->rho[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetaxx[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetaxy[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetaxz[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetayx[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetayy[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetayz[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetazx[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetazy[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetazz[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->gamma[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->khi[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->xi11[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetadx[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetady[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
    ma->zetadz[i] = (double **)malloc(sizeof(double *) * sr.Txx.y);
  }
  for (i = 0; i < sr.Txx.x; i++) {
    for (j = 0; j < sr.Txx.y; j++) {
      ma->ramda[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->mu[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->c11[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->rho[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetaxx[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetaxy[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetaxz[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetayx[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetayy[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetayz[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetazx[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetazy[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetazz[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->gamma[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->khi[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->xi11[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetadx[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetady[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
      ma->zetadz[i][j] = (double *)malloc(sizeof(double) * sr.Txx.z);
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
  initCoord(&clack->sp, spx, spy, spz);
  initCoord(&clack->range, ranx, rany, ranz);
}