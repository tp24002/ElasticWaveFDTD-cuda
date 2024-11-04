#define _USE_MATH_DEFINES
#include "../header/init.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/struct.h"

#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))


void initDimI3(DimI3 *co, int x, int y, int z) {
  co->x = x;
  co->y = y;
  co->z = z;
}

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

void initDiff(Diff *dif, Medium *med) {
  dif->dx = 0.001;
  dif->dy = 0.001;
  dif->dz = 0.001;
  double tmp = 0;
  
  for(int i = E_AIR; i < E_M_END; i++){
    tmp = MAX(sqrt((med[i].K + 4. / 3. * med[i].G) / med[i].rho),tmp);
  }
  dif->dt = dif->dx / (tmp * 1000);
}

void initPml(Pml *pml, Medium *med, Diff dif) {
  pml->ta = 4.;
  pml->fm = 3.574e4;
  double R = 1.e-20;
  double tmp = 0;
  double tmp_v = 0;//max
  initDimI3(&pml->pl1, 32, 32, 32);
  initDimI3(&pml->pl2, 32, 32, 32);
  //計算領域内最高速度
  for(int i = E_AIR; i < E_M_END; i++){
    tmp_v = MAX(sqrt((med[i].K + 4. / 3. * med[i].G) / med[i].rho), tmp_v);
  }
  //減衰係数最大値(PML層)
  tmp = ((tmp_v * (pml->ta + 1)) / (2. * (double)pml->pl1.x * dif.dx)) * log(1/R);
  pml->fm = MAX(tmp, pml->fm);
}

void initRange(Range *ran, DimI3 region, Pml pml) {
  // initDimI3(&ran->sr.Txx, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  // initDimI3(&ran->sr.Tyy, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  // initDimI3(&ran->sr.Tzz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  // initDimI3(&ran->tr.Txy, x + pml.pl1.x + pml.pl2.x - 1, y + pml.pl1.y + pml.pl2.y - 1, z + pml.pl1.z + pml.pl2.z    );
  // initDimI3(&ran->tr.Tyz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y - 1, z + pml.pl1.z + pml.pl2.z - 1);
  // initDimI3(&ran->tr.Tzx, x + pml.pl1.x + pml.pl2.x - 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z - 1);
  // initDimI3(&ran->vr.Vx , x + pml.pl1.x + pml.pl2.x - 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  // initDimI3(&ran->vr.Vy , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y - 1, z + pml.pl1.z + pml.pl2.z    );
  // initDimI3(&ran->vr.Vz , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z - 1);
  initDimI3(&ran->sr.Txx, region.x + pml.pl1.x + pml.pl2.x    , region.y + pml.pl1.y + pml.pl2.y    , region.z + pml.pl1.z + pml.pl2.z    );
  initDimI3(&ran->sr.Tyy, region.x + pml.pl1.x + pml.pl2.x    , region.y + pml.pl1.y + pml.pl2.y    , region.z + pml.pl1.z + pml.pl2.z    );
  initDimI3(&ran->sr.Tzz, region.x + pml.pl1.x + pml.pl2.x    , region.y + pml.pl1.y + pml.pl2.y    , region.z + pml.pl1.z + pml.pl2.z    );
  initDimI3(&ran->tr.Txy, region.x + pml.pl1.x + pml.pl2.x + 1, region.y + pml.pl1.y + pml.pl2.y + 1, region.z + pml.pl1.z + pml.pl2.z    );
  initDimI3(&ran->tr.Tyz, region.x + pml.pl1.x + pml.pl2.x    , region.y + pml.pl1.y + pml.pl2.y + 1, region.z + pml.pl1.z + pml.pl2.z + 1);
  initDimI3(&ran->tr.Tzx, region.x + pml.pl1.x + pml.pl2.x + 1, region.y + pml.pl1.y + pml.pl2.y    , region.z + pml.pl1.z + pml.pl2.z + 1);
  initDimI3(&ran->vr.Vx , region.x + pml.pl1.x + pml.pl2.x + 1, region.y + pml.pl1.y + pml.pl2.y    , region.z + pml.pl1.z + pml.pl2.z    );
  initDimI3(&ran->vr.Vy , region.x + pml.pl1.x + pml.pl2.x    , region.y + pml.pl1.y + pml.pl2.y + 1, region.z + pml.pl1.z + pml.pl2.z    );
  initDimI3(&ran->vr.Vz , region.x + pml.pl1.x + pml.pl2.x    , region.y + pml.pl1.y + pml.pl2.y    , region.z + pml.pl1.z + pml.pl2.z + 1);
} 

void initRandom(Object con, DimI3 *clack, int ratio) {
  if(con.range.x < 3 || con.range.y < 3 || con.range.z < 3) {
    printf("Cannot place defects.\n");
    return;
  }
  int count = 0;
  // コンクリートセル数
  int max_Patern = con.range.x * con.range.y * con.range.z;
  // 内部欠陥パターン数
  int max_ClackPatern = (con.range.x - 2) * (con.range.y - 2) * (con.range.z - 2);
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
    int rand1 = rand() % (con.range.x - 2) + 1;
    int rand2 = rand() % (con.range.y - 2) + 1;
    int rand3 = rand() % (con.range.z - 2) + 1;

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

__global__ void ZeroT(BefAft *ba, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;
  int id;
  if(i < imax && j < jmax && k < kmax) {
    id = k * imax * jmax + j * imax + i;
    ba[id].sa.Txx = 0;
    ba[id].sa.Txxx = 0;
    ba[id].sa.Txxy = 0;
    ba[id].sa.Txxz = 0;
    ba[id].sa.Tyy = 0;
    ba[id].sa.Tyyx = 0;
    ba[id].sa.Tyyy = 0;
    ba[id].sa.Tyyz = 0;
    ba[id].sa.Tzz = 0;
    ba[id].sa.Tzzx = 0;
    ba[id].sa.Tzzy = 0;
    ba[id].sa.Tzzz = 0;
  } 
}

__global__ void ZeroTxy(BefAft *ba, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->tr.Txy.x, jmax = ran->tr.Txy.y, kmax = ran->tr.Txy.z;
  int id;
  if(i < imax && j < jmax && k < kmax) {
    id = k * imax * jmax + j * imax + i;
    ba[id].ta.Txy = 0;
    ba[id].ta.Txyx = 0;
    ba[id].ta.Txyy = 0;
  } 
}

__global__ void ZeroTyz(BefAft *ba, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->tr.Tyz.x, jmax = ran->tr.Tyz.y, kmax = ran->tr.Tyz.z;
  int id;
  if(i < imax && j < jmax && k < kmax) {
    id = k * imax * jmax + j * imax + i;
    ba[id].ta.Tyz = 0;
    ba[id].ta.Tyzy = 0;
    ba[id].ta.Tyzz = 0;
  } 
}

__global__ void ZeroTzx(BefAft *ba, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->tr.Tzx.x, jmax = ran->tr.Tzx.y, kmax = ran->tr.Tzx.z;
  int id;
  if(i < imax && j < jmax && k < kmax) {
    id = k * imax * jmax + j * imax + i;
    ba[id].ta.Tzx = 0;
    ba[id].ta.Tzxz = 0;
    ba[id].ta.Tzxx = 0;
  } 
}

__global__ void ZeroVx(BefAft *ba, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->vr.Vx.x, jmax = ran->vr.Vx.y, kmax = ran->vr.Vx.z;
  int id;
  if(i < imax && j < jmax && k < kmax) {
    id = k * imax * jmax + j * imax + i;
    ba[id].va.Vx = 0;
    ba[id].va.Vxx = 0;
    ba[id].va.Vxy = 0;
    ba[id].va.Vxz = 0;
  } 
}

__global__ void ZeroVy(BefAft *ba, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->vr.Vy.x, jmax = ran->vr.Vy.y, kmax = ran->vr.Vy.z;
  int id;
  if(i < imax && j < jmax && k < kmax) {
    id = k * imax * jmax + j * imax + i;
    ba[id].va.Vy = 0;
    ba[id].va.Vyx = 0;
    ba[id].va.Vyy = 0;
    ba[id].va.Vyz = 0;
  } 
}

__global__ void ZeroVz(BefAft *ba, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->vr.Vz.x, jmax = ran->vr.Vz.y, kmax = ran->vr.Vz.z;
  int id;
  if(i < imax && j < jmax && k < kmax) {
    id = k * imax * jmax + j * imax + i;
    ba[id].va.Vz = 0;
    ba[id].va.Vzx = 0;
    ba[id].va.Vzy = 0;
    ba[id].va.Vzz = 0;
  } 
}
