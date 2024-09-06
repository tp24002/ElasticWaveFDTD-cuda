#define _USE_MATH_DEFINES
#include "../header/update.h"
#include "../header/struct.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/insert.h"
#include "../header/init.h"
#include "../header/extradition.h"

// 垂直応力

// 引数に(int imax, int jmax, int kmax)を用いている理由はない
// (SigRan sr)でも良い，ただ最初に書いたコードがこれだっただけ
// max:構造体にまとめていいかも
// 垂直応力更新並列関数
__global__ void TxxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Coord Txx) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = Txx.x, jmax = Txx.y, kmax = Txx.z;
  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax - 1) {
    aft->sa.Txxx[i][j][k] = (2. - ma->zetadx[i][j][k] * dif->dt) / (2. + ma->zetadx[i][j][k] * dif->dt) * bef->sa.Txxx[i][j][k]
      + 2. * (ma->c11[i][j][k] * dif->dt + ma->xi11[i][j][k]) / (2. + ma->zetadx[i][j][k] * dif->dt) * (aft->va.Vx[i][j][k] - aft->va.Vx[i - 1][j][k]) / dif->dx 
      - 2. * ma->xi11[i][j][k] / (2. + ma->zetadx[i][j][k] * dif->dt) * (bef->va.Vx[i][j][k] - bef->va.Vx[i - 1][j][k]) / dif->dx;

    aft->sa.Txxy[i][j][k] = (2. - ma->zetady[i][j][k] * dif->dt) / (2. + ma->zetady[i][j][k] * dif->dt) * bef->sa.Txxy[i][j][k]
      + 2. * (ma->ramda[i][j][k] * dif->dt + ma->khi[i][j][k]) / (2. + ma->zetady[i][j][k] * dif->dt) * (aft->va.Vy[i][j][k] - aft->va.Vy[i][j - 1][k]) / dif->dy
      - 2. * ma->khi[i][j][k] / (2. + ma->zetady[i][j][k] * dif->dt) * (bef->va.Vy[i][j][k] - bef->va.Vy[i][j - 1][k]) / dif->dy;

    aft->sa.Txxz[i][j][k] = (2. - ma->zetadz[i][j][k] * dif->dt) / (2. + ma->zetadz[i][j][k] * dif->dt) * bef->sa.Txxz[i][j][k]
      + 2. * (ma->ramda[i][j][k] * dif->dt + ma->khi[i][j][k]) / (2. + ma->zetadz[i][j][k] * dif->dt) * (aft->va.Vz[i][j][k] - aft->va.Vz[i][j][k - 1]) / dif->dz
      - 2. * ma->khi[i][j][k] / (2. + ma->zetadz[i][j][k] * dif->dt) * (bef->va.Vz[i][j][k] - bef->va.Vz[i][j][k - 1]) / dif->dz;
  }
}
// 垂直応力更新並列関数
__global__ void TyyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Coord Tyy) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
  
  int imax = Tyy.x, jmax = Tyy.y, kmax = Tyy.z;
  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax - 1) {
    aft->sa.Tyyx[i][j][k] = (2. - ma->zetadx[i][j][k] * dif->dt) / (2. + ma->zetadx[i][j][k] * dif->dt) * bef->sa.Tyyx[i][j][k]
      + 2. * (ma->ramda[i][j][k] * dif->dt + ma->khi[i][j][k]) / (2. + ma->zetadx[i][j][k] * dif->dt) * (aft->va.Vx[i][j][k] - aft->va.Vx[i - 1][j][k]) / dif->dx
      - 2. * ma->khi[i][j][k] / (2. + ma->zetadx[i][j][k] * dif->dt) * (bef->va.Vx[i][j][k] - bef->va.Vx[i - 1][j][k]) / dif->dx;

    aft->sa.Tyyy[i][j][k] = (2. - ma->zetady[i][j][k] * dif->dt) / (2. + ma->zetady[i][j][k] * dif->dt) * bef->sa.Tyyy[i][j][k]
      + 2. * (ma->c11[i][j][k] * dif->dt + ma->xi11[i][j][k]) / (2. + ma->zetady[i][j][k] * dif->dt) * (aft->va.Vy[i][j][k] - aft->va.Vy[i][j - 1][k]) / dif->dy
      - 2. * ma->xi11[i][j][k] / (2. + ma->zetady[i][j][k] * dif->dt) * (bef->va.Vy[i][j][k] - bef->va.Vy[i][j - 1][k]) / dif->dy;

    aft->sa.Tyyz[i][j][k] = (2. - ma->zetadz[i][j][k] * dif->dt) / (2. + ma->zetadz[i][j][k] * dif->dt) * bef->sa.Tyyz[i][j][k]
      + 2. * (ma->ramda[i][j][k] * dif->dt + ma->khi[i][j][k]) / (2. + ma->zetadz[i][j][k] * dif->dt) * (aft->va.Vz[i][j][k] - aft->va.Vz[i][j][k - 1]) / dif->dz
      - 2. * ma->khi[i][j][k] / (2. + ma->zetadz[i][j][k] * dif->dt)  * (bef->va.Vz[i][j][k] - bef->va.Vz[i][j][k - 1]) / dif->dz;
  }
}
// 垂直応力更新並列関数
__global__ void TzzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Coord Tzz){
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
  
  int imax = Tzz.x, jmax = Tzz.y, kmax = Tzz.z;
  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax - 1) {
    aft->sa.Tzzx[i][j][k] = (2. - ma->zetadx[i][j][k] * dif->dt) / (2. + ma->zetadx[i][j][k] * dif->dt) * bef->sa.Tzzx[i][j][k]
      + 2. * (ma->ramda[i][j][k] * dif->dt + ma->khi[i][j][k]) / (2. + ma->zetadx[i][j][k] * dif->dt) * (aft->va.Vx[i][j][k] - aft->va.Vx[i - 1][j][k]) / dif->dx
      - 2. * ma->khi[i][j][k] / (2. + ma->zetadx[i][j][k] * dif->dt) * (bef->va.Vx[i][j][k] - bef->va.Vx[i - 1][j][k]) / dif->dx;

    aft->sa.Tzzy[i][j][k] = (2. - ma->zetady[i][j][k] * dif->dt) / (2. + ma->zetady[i][j][k] * dif->dt) * bef->sa.Tzzy[i][j][k]
      + 2. * (ma->ramda[i][j][k] * dif->dt + ma->khi[i][j][k]) / (2. + ma->zetady[i][j][k] * dif->dt) * (aft->va.Vy[i][j][k] - aft->va.Vy[i][j - 1][k]) / dif->dy
      - 2. * ma->khi[i][j][k] / (2. + ma->zetady[i][j][k] * dif->dt)  * (bef->va.Vy[i][j][k] - bef->va.Vy[i][j - 1][k]) / dif->dy;

    aft->sa.Tzzz[i][j][k] = (2. - ma->zetadz[i][j][k] * dif->dt) / (2. + ma->zetadz[i][j][k] * dif->dt) * bef->sa.Tzzz[i][j][k]
      + 2. * (ma->c11[i][j][k] * dif->dt + ma->xi11[i][j][k]) / (2. + ma->zetadz[i][j][k] * dif->dt) * (aft->va.Vz[i][j][k] - aft->va.Vz[i][j][k - 1]) / dif->dz
      - 2. * ma->xi11[i][j][k] / (2. + ma->zetadz[i][j][k] * dif->dt) * (bef->va.Vz[i][j][k] - bef->va.Vz[i][j][k - 1]) / dif->dz;
  }
}
// T 0 padding
__global__ void ZeroT_XY(BefAft *aft, Coord ranmax, char check) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = ranmax.x, jmax = ranmax.y, kmax = ranmax.z;

  if(j > jmax || i > imax) {
    return;
  }
  // if (j <= jmax && i <= imax) {
  if(check == 'X') {
    aft->sa.Txxx[i][j][0] = 0.;
    aft->sa.Txxx[i][j][kmax] = 0.;
    aft->sa.Txxy[i][j][0] = 0.;
    aft->sa.Txxy[i][j][kmax] = 0.;
    aft->sa.Txxz[i][j][0] = 0.;
    aft->sa.Txxz[i][j][kmax] = 0.;
  } else if(check == 'Y') {
    aft->sa.Tyyx[i][j][0] = 0.;
    aft->sa.Tyyx[i][j][kmax] = 0.;
    aft->sa.Tyyy[i][j][0] = 0.;
    aft->sa.Tyyy[i][j][kmax] = 0.;
    aft->sa.Tyyz[i][j][0] = 0.;
    aft->sa.Tyyz[i][j][kmax] = 0.;
  } else if(check == 'Z') {
    aft->sa.Tzzx[i][j][0] = 0.;
    aft->sa.Tzzx[i][j][kmax] = 0.;
    aft->sa.Tzzy[i][j][0] = 0.;
    aft->sa.Tzzy[i][j][kmax] = 0.;
    aft->sa.Tzzz[i][j][0] = 0.;
    aft->sa.Tzzz[i][j][kmax] = 0.;
  } else {
    printf("ZeroT_XY");
  }
}
// T 0 padding
__global__ void ZeroT_YZ(BefAft *aft, Coord ranmax, char check) {
  int k = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = ranmax.x, jmax = ranmax.y, kmax = ranmax.z;

  if(k > kmax - 1 || j > jmax) {
    return;
  }
  // if (k <= kmax - 1 && j <= jmax) {
  if(check == 'X') {
    aft->sa.Txxx[0][j][k] = 0.;
    aft->sa.Txxx[imax][j][k] = 0.;
    aft->sa.Txxy[0][j][k] = 0.;
    aft->sa.Txxy[imax][j][k] = 0.;
    aft->sa.Txxz[0][j][k] = 0.;
    aft->sa.Txxz[imax][j][k] = 0.;
  } else if(check =='Y') {
    aft->sa.Tyyx[0][j][k] = 0.;
    aft->sa.Tyyx[imax][j][k] = 0.;
    aft->sa.Tyyy[0][j][k] = 0.;
    aft->sa.Tyyy[imax][j][k] = 0.;
    aft->sa.Tyyz[0][j][k] = 0.;
    aft->sa.Tyyz[imax][j][k] = 0.;
  } else if(check == 'Z') {
    aft->sa.Tzzx[0][j][k] = 0.;
    aft->sa.Tzzx[imax][j][k] = 0.;
    aft->sa.Tzzy[0][j][k] = 0.;
    aft->sa.Tzzy[imax][j][k] = 0.;
    aft->sa.Tzzz[0][j][k] = 0.;
    aft->sa.Tzzz[imax][j][k] = 0.;
  } else {
    printf("ZeroT_YZ");
  }
}
// T 0 padding
__global__ void ZeroT_ZX(BefAft *aft, Coord ranmax, char check) {
  int k = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int i = blockIdx.y * blockDim.y + threadIdx.y + 1;

  int imax = ranmax.x, jmax = ranmax.y, kmax = ranmax.z;

  if(k > kmax - 1 || i > imax - 1) {
    return;
  }
  // if (k <= Txxkmax - 1 && i <= Txximax - 1) {
  if(check == 'X') {
    aft->sa.Txxx[i][0][k] = 0.;
    aft->sa.Txxx[i][jmax][k] = 0.;
    aft->sa.Txxy[i][0][k] = 0.;
    aft->sa.Txxy[i][jmax][k] = 0.;
    aft->sa.Txxz[i][0][k] = 0.;
    aft->sa.Txxz[i][jmax][k] = 0.;
  } else if(check == 'Y') {
    aft->sa.Tyyx[i][0][k] = 0.;
    aft->sa.Tyyx[i][jmax][k] = 0.;
    aft->sa.Tyyy[i][0][k] = 0.;
    aft->sa.Tyyy[i][jmax][k] = 0.;
    aft->sa.Tyyz[i][0][k] = 0.;
    aft->sa.Tyyz[i][jmax][k] = 0.;
  } else if(check == 'Z'){
    aft->sa.Tzzx[i][0][k] = 0.;
    aft->sa.Tzzx[i][jmax][k] = 0.;
    aft->sa.Tzzy[i][0][k] = 0.;
    aft->sa.Tzzy[i][jmax][k] = 0.;
    aft->sa.Tzzz[i][0][k] = 0.;
    aft->sa.Tzzz[i][jmax][k] = 0.;
  } else {
    printf("ZeroT_ZX");
  }
}
// 全方向加算
__global__ void DirectionalAdd(BefAft *aft, Inpaluse *ip, Coord ranmax, char check) {
  // スレッドインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ranmax.x, jmax = ranmax.y, kmax = ranmax.z;

  if(i > imax || j > jmax || k > kmax) {
    return;
  }
  // if (i <= imax && j <= jmax && k <= kmax) {
  if(check == 'X') {
    aft->sa.Txx[i][j][k] = aft->sa.Txxx[i][j][k] + aft->sa.Txxy[i][j][k] + aft->sa.Txxz[i][j][k] + ip->Txx[i][j][k];
  } else if(check == 'Y') {
    aft->sa.Tyy[i][j][k] = aft->sa.Tyyx[i][j][k] + aft->sa.Tyyy[i][j][k] + aft->sa.Tyyz[i][j][k] + ip->Tyy[i][j][k];
  } else if(check == 'Z'){
    aft->sa.Tzz[i][j][k] = aft->sa.Tzzx[i][j][k] + aft->sa.Tzzy[i][j][k] + aft->sa.Tzzz[i][j][k] + ip->Tzz[i][j][k];
  } else {
    printf("error:DirectionalAdd");
  }
}
// Txxクラス的な(Blocks大丈夫かな？)
void Txx(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Inpaluse ip_h, Inpaluse *ip_d, int t, Coord threads) {
  char check = 'X';
  Coord ranmax;
  // host->deviceデータ転送
  copyBefAftToDevice(aft_d, aft_h, ran);
  copyBefAftToDevice(bef_d, bef_h, ran);
  copyMedArrToDevice(ma_d, &ma_h, ran);
  copyDiffToDevice(dif_d, &dif_h);
  copyInpaluseToDevice(ip_d, &ip_h, ran);

  int Txximax = ran.sr.Txx.x, Txxjmax = ran.sr.Txx.y, Txxkmax = ran.sr.Txx.z;
  initCoord(&ranmax, Txximax, Txxjmax, Txxkmax);

  dim3 threadsPerBlock(threads.x, threads.y, threads.z);  // ブロック内のスレッド数
  dim3 UpdateBlocks((Txximax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Txxjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Txxkmax - 1 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroXYBlocks(    (Txxjmax + threadsPerBlock.x - 1) / threadsPerBlock.x,     (Txximax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroYZBlocks((Txxkmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,     (Txxjmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroZXBlocks((Txxkmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Txximax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Txximax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Txxjmax + threadsPerBlock.y) / threadsPerBlock.y, 
                            (Txxkmax + threadsPerBlock.z) / threadsPerBlock.z);
  // Txx更新式
  TxxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran.sr.Txx);
  // 0 padding
  ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  //全方向加算
  DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ip_d, ranmax, check);
  copyBefAftToHost(aft_h, aft_d, ran);
  copyBefAftToHost(bef_h, bef_d, ran);
}
// Tyyクラス的な(Blocks大丈夫かな？)
void Tyy(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Inpaluse ip_h, Inpaluse *ip_d, int t, Coord threads) {
  char check = 'Y';
  Coord ranmax;
  // host->deviceデータ転送
  copyBefAftToDevice(aft_d, aft_h, ran);
  copyBefAftToDevice(bef_d, bef_h, ran);
  copyMedArrToDevice(ma_d, &ma_h, ran);
  copyDiffToDevice(dif_d, &dif_h);
  copyInpaluseToDevice(ip_d, &ip_h, ran);

  int Tyyimax = ran.sr.Tyy.x, Tyyjmax = ran.sr.Tyy.y, Tyykmax = ran.sr.Tyy.z;
  initCoord(&ranmax, Tyyimax, Tyyjmax, Tyykmax);

  dim3 threadsPerBlock(threads.x, threads.y, threads.z);  // ブロック内のスレッド数
  dim3 UpdateBlocks((Tyyimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tyyjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tyykmax - 1 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroXYBlocks(    (Tyyjmax + threadsPerBlock.x - 1) / threadsPerBlock.x,     (Tyyimax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroYZBlocks((Tyykmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,     (Tyyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroZXBlocks((Tyykmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tyyimax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Tyyimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Tyyjmax + threadsPerBlock.y) / threadsPerBlock.y,
                            (Tyykmax + threadsPerBlock.z) / threadsPerBlock.z);
  // Tyy更新式
  TyyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ranmax);
  // 0 padding
  ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  // 全方向加算
  DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ip_d, ranmax, check);
  // データ転送device to host
  copyBefAftToHost(aft_h, aft_d, ran);
  copyBefAftToHost(bef_h, bef_d, ran);
}
// Tzzクラス的な(Blocks大丈夫かな？)
void Tzz(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Inpaluse ip_h, Inpaluse *ip_d, int t, Coord threads) {
  char check = 'Z';
  Coord ranmax;
  // host->deviceデータ転送
  copyBefAftToDevice(aft_d, aft_h, ran);
  copyBefAftToDevice(bef_d, bef_h, ran);
  copyMedArrToDevice(ma_d, &ma_h, ran);
  copyDiffToDevice(dif_d, &dif_h);
  copyInpaluseToDevice(ip_d, &ip_h, ran);

  int Tzzimax = ran.sr.Tzz.x, Tzzjmax = ran.sr.Tzz.y, Tzzkmax = ran.sr.Tzz.z;
  initCoord(&ranmax, Tzzimax, Tzzjmax, Tzzkmax);
  dim3 threadsPerBlock(threads.x, threads.y, threads.z);  // ブロック内のスレッド数
  dim3 UpdateBlocks((Tzzimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tzzjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tzzkmax - 1 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroXYBlocks(    (Tzzjmax + threadsPerBlock.x - 1) / threadsPerBlock.x,     (Tzzimax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroYZBlocks((Tzzkmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,     (Tzzjmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroZXBlocks((Tzzkmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tzzimax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Tzzimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Tzzjmax + threadsPerBlock.y) / threadsPerBlock.y,
                            (Tzzkmax + threadsPerBlock.z) / threadsPerBlock.z);
  // Tzzの更新式
  TzzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ranmax);
  // 0 padding
  ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  // 全方向加算
  DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ip_d, ranmax, check);
  // データ転送device to host
  copyBefAftToHost(aft_h, aft_d, ran);
  copyBefAftToHost(bef_h, bef_d, ran);
}
// 垂直応力計算(main呼び出し関数)
void Sig(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Inpaluse ip_h, Inpaluse *ip_d, int t, Coord threads) {
  Txx(aft_h, bef_h, aft_d, bef_d, ma_h, ma_d, dif_h, dif_d, ran, ip_h, ip_d, t, threads);
  Tyy(aft_h, bef_h, aft_d, bef_d, ma_h, ma_d, dif_h, dif_d, ran, ip_h, ip_d, t, threads);
  Tzz(aft_h, bef_h, aft_d, bef_d, ma_h, ma_d, dif_h, dif_d, ran, ip_h, ip_d, t, threads);
}

// せん断応力

// せん断応力更新関数
__global__ void TxyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1; // 始点を+1

  int imax = tr.Txy.x, jmax = tr.Txy.y, kmax = tr.Txy.z;
  double Hzetadx, Hzetady, Hmu, Hgamma;

  if(i <= imax && j <= jmax && k <= kmax - 1) {
    //PML:減衰係数,計算領域:摩擦定数
    Hzetadx = 4. * pow((1. / ma->zetadx[i + 1][j + 1][k]) + (1. / ma->zetadx[i][j + 1][k]) + (1. / ma->zetadx[i + 1][j][k]) + (1. / ma->zetadx[i][j][k]), -1.);
    //PML:減衰係数,計算領域:摩擦定数
    Hzetady = 4. * pow((1. / ma->zetady[i + 1][j + 1][k]) + (1. / ma->zetady[i][j + 1][k]) + (1. / ma->zetady[i + 1][j][k]) + (1. / ma->zetady[i][j][k]), -1.);
    //第2ラメ，横弾性係数(剛性率)
    Hmu     = 4. * pow((1. /     ma->mu[i + 1][j + 1][k]) + (1. /     ma->mu[i][j + 1][k]) + (1. /     ma->mu[i + 1][j][k]) + (1. /     ma->mu[i][j][k]), -1.);
    //第１粘性定数
    Hgamma  = 4. * pow((1. /  ma->gamma[i + 1][j + 1][k]) + (1. /  ma->gamma[i][j + 1][k]) + (1. /  ma->gamma[i + 1][j][k]) + (1. /  ma->gamma[i][j][k]), -1.);
    aft->ta.Txyx[i][j][k] = (2. - Hzetadx * dif->dt) / (2. + Hzetadx * dif->dt) * bef->ta.Txyx[i][j][k]
      + 2. * (Hmu * dif->dt + Hgamma) / (2. + Hzetadx * dif->dt) * (aft->va.Vy[i + 1][j][k] - aft->va.Vy[i][j][k]) / dif->dx
      - 2. * Hgamma / (2. + Hzetadx * dif->dt) * (bef->va.Vy[i + 1][j][k] - bef->va.Vy[i][j][k]) / dif->dx;

    aft->ta.Txyy[i][j][k] = (2. - Hzetady * dif->dt) / (2. + Hzetady * dif->dt) * bef->ta.Txyy[i][j][k]
      + 2. * (Hmu * dif->dt + Hgamma) / (2. + Hzetady * dif->dt) * (aft->va.Vx[i][j + 1][k] - aft->va.Vx[i][j][k]) / dif->dy
      - 2. * Hgamma / (2. + Hzetady * dif->dt) * (bef->va.Vx[i][j + 1][k] - bef->va.Vx[i][j][k]) / dif->dy;
  }
}
// せん断応力更新関数
__global__ void TyzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  
  int imax = tr.Tyz.x, jmax = tr.Tyz.y, kmax = tr.Tyz.z;
  double Hzetady, Hzetadz, Hmu, Hgamma;

  if(i <= imax - 1 && j <= jmax && k <= kmax) {
    Hzetady = 4. * pow((1. / ma->zetady[i][j + 1][k + 1]) + (1. / ma->zetady[i][j + 1][k]) + (1. / ma->zetady[i][j][k + 1]) + (1. / ma->zetady[i][j][k]), -1.);
    Hzetadz = 4. * pow((1. / ma->zetadz[i][j + 1][k + 1]) + (1. / ma->zetadz[i][j + 1][k]) + (1. / ma->zetadz[i][j][k + 1]) + (1. / ma->zetadz[i][j][k]), -1.);
    Hmu     = 4. * pow((1. /     ma->mu[i][j + 1][k + 1]) + (1. /     ma->mu[i][j + 1][k]) + (1. /     ma->mu[i][j][k + 1]) + (1. /     ma->mu[i][j][k]), -1.);
    Hgamma  = 4. * pow((1. /  ma->gamma[i][j + 1][k + 1]) + (1. /  ma->gamma[i][j + 1][k]) + (1. /  ma->gamma[i][j][k + 1]) + (1. /  ma->gamma[i][j][k]), -1.);
    aft->ta.Tyzy[i][j][k] = (2. - Hzetady * dif->dt) / (2. + Hzetady * dif->dt) * bef->ta.Tyzy[i][j][k]
      + 2. * (Hmu * dif->dt + Hgamma) / (2. + Hzetady * dif->dt) * (aft->va.Vz[i][j + 1][k] - aft->va.Vz[i][j][k]) / dif->dy
      - 2. * Hgamma / (2. + Hzetady * dif->dt) * (bef->va.Vz[i][j + 1][k] - bef->va.Vz[i][j][k]) / dif->dy;

    aft->ta.Tyzz[i][j][k] = (2. - Hzetadz * dif->dt) / (2. + Hzetadz * dif->dt) * bef->ta.Tyzz[i][j][k]
      + 2. * (Hmu * dif->dt + Hgamma) / (2. + Hzetadz * dif->dt) * (aft->va.Vy[i][j][k + 1] - aft->va.Vy[i][j][k]) / dif->dz 
      - 2. * Hgamma / (2. + Hzetadz * dif->dt) * (bef->va.Vy[i][j][k + 1] - bef->va.Vy[i][j][k]) / dif->dz;
  }
}
// せん断応力更新関数
__global__ void TzxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  
  int imax = tr.Tzx.x, jmax = tr.Tzx.y, kmax = tr.Tzx.z;
  double Hzetadx, Hzetadz, Hmu, Hgamma;

  if(i <= imax && j <= jmax - 1 && k <= kmax) {
    Hzetadx = 4. * pow((1. / ma->zetadx[i + 1][j][k + 1]) + (1. / ma->zetadx[i + 1][j][k]) + (1. / ma->zetadx[i][j][k + 1]) + (1. / ma->zetadx[i][j][k]), -1.);
    Hzetadz = 4. * pow((1. / ma->zetadz[i + 1][j][k + 1]) + (1. / ma->zetadz[i + 1][j][k]) + (1. / ma->zetadz[i][j][k + 1]) + (1. / ma->zetadz[i][j][k]), -1.);
    Hmu     = 4. * pow((1. /     ma->mu[i + 1][j][k + 1]) + (1. /     ma->mu[i + 1][j][k]) + (1. /     ma->mu[i][j][k + 1]) + (1. /     ma->mu[i][j][k]), -1.);
    Hgamma  = 4. * pow((1. /  ma->gamma[i + 1][j][k + 1]) + (1. /  ma->gamma[i + 1][j][k]) + (1. /  ma->gamma[i][j][k + 1]) + (1. /  ma->gamma[i][j][k]), -1.);
    aft->ta.Tzxz[i][j][k] = (2. - Hzetadz * dif->dt) / (2. + Hzetadz * dif->dt) * bef->ta.Tzxz[i][j][k]
      + 2. * (Hmu * dif->dt + Hgamma) / (2. + Hzetadz * dif->dt) * (aft->va.Vx[i][j][k + 1] - aft->va.Vx[i][j][k]) / dif->dz
      - 2. * Hgamma / (2. + Hzetadz * dif->dt) * (bef->va.Vx[i][j][k + 1] - bef->va.Vx[i][j][k]) / dif->dz;

    aft->ta.Tzxx[i][j][k] = (2. - Hzetadx * dif->dt) / (2. + Hzetadx * dif->dt) * bef->ta.Tzxx[i][j][k]
      + 2. * (Hmu * dif->dt + Hgamma) / (2. + Hzetadx * dif->dt) * (aft->va.Vz[i + 1][j][k] - aft->va.Vz[i][j][k]) / dif->dx
      - 2. * Hgamma / (2. + Hzetadx * dif->dt) * (bef->va.Vz[i + 1][j][k] - bef->va.Vz[i][j][k]) / dif->dx;
  }
}

__global__ void ZeroTxy(BefAft *aft, Coord Tmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if(i <= Tmax.x && j <= Tmax.y) {
    aft->ta.Txyx[i][j][0] = 0.;
    aft->ta.Txyx[i][j][Tmax.z] = 0.;
    aft->ta.Txyy[i][j][0] = 0.;
    aft->ta.Txyy[i][j][Tmax.z] = 0.;
  }
}

__global__ void ZeroTyz(BefAft *aft, Coord Tmax) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;

  if(j <= Tmax.y && k <= Tmax.z) {
    aft->ta.Tyzy[0][j][k] = 0.;
    aft->ta.Tyzy[Tmax.x][j][k] = 0.;
    aft->ta.Tyzz[0][j][k] = 0.;
    aft->ta.Tyzz[Tmax.x][j][k] = 0.;
  }
}

__global__ void ZeroTzx(BefAft *aft, Coord Tmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;

  if(i <= Tmax.x && k <= Tmax.z) {
    aft->ta.Tzxx[i][0][k] = 0.;
    aft->ta.Tzxx[i][Tmax.z][k] = 0.;
    aft->ta.Tzxz[i][0][k] = 0.;
    aft->ta.Tzxz[i][Tmax.z][k] = 0.;
  }
}

__global__ void DirectionalAddT(BefAft *aft, Coord Tmax, char check) {
  // スレッドインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if(i > Tmax.x || j > Tmax.y || k > Tmax.z){
    return;
  }
  if(check == 'X') {
    aft->ta.Tyz[i][j][k] = aft->ta.Tyzy[i][j][k] + aft->ta.Tyzz[i][j][k];
  } else if(check == 'Y') {
    aft->ta.Tzx[i][j][k] = aft->ta.Tzxx[i][j][k] + aft->ta.Tzxz[i][j][k];
  } else if(check == 'Z') {
    aft->ta.Txy[i][j][k] = aft->ta.Txyx[i][j][k] + aft->ta.Txyy[i][j][k];
  } else {
    printf("error:DirectionalAddT");
  }
}
// Txyクラス的な
void Txy(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads) {

  // host->deviceデータ転送
  copyBefAftToDevice(aft_d, aft_h, ran);
  copyBefAftToDevice(bef_d, bef_h, ran);
  copyMedArrToDevice(ma_d, &ma_h, ran);
  copyDiffToDevice(dif_d, &dif_h);

  int Txyimax = ran.tr.Txy.x, Txyjmax = ran.tr.Txy.y, Txykmax = ran.tr.Txy.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Txyimax + threadsPerBlock.x - 1)     / threadsPerBlock.x,
                    (Txyjmax + threadsPerBlock.y - 1)     / threadsPerBlock.y,
                    (Txykmax - 1 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroXYBlocks((Txyimax + threadsPerBlock.x - 1) / threadsPerBlock.x,(Txyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Txyimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Txyjmax + threadsPerBlock.y) / threadsPerBlock.y,
                            (Txykmax + threadsPerBlock.z) / threadsPerBlock.z);                    
  TxyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran.tr);
  ZeroTxy<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran.tr.Txy);
  DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran.tr.Txy, 'Z');
  // データ転送device to host
  copyBefAftToHost(aft_h, aft_d, ran);
  copyBefAftToHost(bef_h, bef_d, ran);
}
// Tyzクラス的な
void Tyz(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads) {
  
  // host->deviceデータ転送
  copyBefAftToDevice(aft_d, aft_h, ran);
  copyBefAftToDevice(bef_d, bef_h, ran);
  copyMedArrToDevice(ma_d, &ma_h, ran);
  copyDiffToDevice(dif_d, &dif_h);

  int Tyzimax = ran.tr.Tyz.x, Tyzjmax = ran.tr.Tyz.y, Tyzkmax = ran.tr.Tyz.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Tyzimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tyzjmax + threadsPerBlock.y - 1)     / threadsPerBlock.y,
                    (Tyzkmax + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  dim3 ZeroYZBlocks((Tyzjmax + threadsPerBlock.x - 1) / threadsPerBlock.x,(Tyzkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Tyzimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Tyzjmax + threadsPerBlock.y) / threadsPerBlock.y,
                            (Tyzkmax + threadsPerBlock.z) / threadsPerBlock.z);                    
  TyzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran.tr);
  ZeroTyz<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran.tr.Tyz);
  DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran.tr.Tyz, 'X');
  // データ転送device to host
  copyBefAftToHost(aft_h, aft_d, ran);
  copyBefAftToHost(bef_h, bef_d, ran);
}
// Tzxクラス的な
void Tzx(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads) {
  
  // host->deviceデータ転送
  copyBefAftToDevice(aft_d, aft_h, ran);
  copyBefAftToDevice(bef_d, bef_h, ran);
  copyMedArrToDevice(ma_d, &ma_h, ran);
  copyDiffToDevice(dif_d, &dif_h);

  int Tzximax = ran.tr.Tzx.x, Tzxjmax = ran.tr.Tzx.y, Tzxkmax = ran.tr.Tzx.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Tzximax + threadsPerBlock.x - 1)     / threadsPerBlock.x,
                    (Tzxjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tzxkmax + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  dim3 ZeroZXBlocks((Tzximax + threadsPerBlock.x - 1) / threadsPerBlock.x,(Tzxkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);   
  dim3 DirectionalAddBlocks((Tzximax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Tzxjmax + threadsPerBlock.y) / threadsPerBlock.y, 
                            (Tzxkmax + threadsPerBlock.z) / threadsPerBlock.z);                  
  TzxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran.tr);
  ZeroTzx<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran.tr.Tzx);
  DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran.tr.Tzx , 'Y');
  // データ転送device to host
  copyBefAftToHost(aft_h, aft_d, ran);
  copyBefAftToHost(bef_h, bef_d, ran);
}
// せん断応力計算(main呼び出し関数)
void Tau(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads) {
  Txy(aft_h, bef_h, aft_d, bef_d, ma_h, ma_d, dif_h, dif_d, ran, threads);
  Tyz(aft_h, bef_h, aft_d, bef_d, ma_h, ma_d, dif_h, dif_d, ran, threads);
  Tzx(aft_h, bef_h, aft_d, bef_d, ma_h, ma_d, dif_h, dif_d, ran, threads);
}

// 粒子速度

// 粒子速度更新関数
__global__ void VxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, VelRan vr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = vr.Vx.x, jmax = vr.Vx.y, kmax = vr.Vx.z;
  double Azetaxx, Azetaxy, Azetaxz, Arho;

  if(i <= imax && j <= jmax - 1 && k <= kmax - 1) {
    Azetaxx = (ma->zetaxx[i + 1][j][k] + ma->zetaxx[i][j][k]) / 2.;
    Azetaxy = (ma->zetaxy[i + 1][j][k] + ma->zetaxy[i][j][k]) / 2.;
    Azetaxz = (ma->zetaxz[i + 1][j][k] + ma->zetaxz[i][j][k]) / 2.;
    Arho    = (   ma->rho[i + 1][j][k] +    ma->rho[i][j][k]) / 2.;
    aft->va.Vxx[i][j][k] = (2. * Arho - Azetaxx * dif->dt) / (2. * Arho + Azetaxx * dif->dt) * bef->va.Vxx[i][j][k]
      + 2. * dif->dt / (2. * Arho + Azetaxx * dif->dt) * (bef->sa.Txx[i + 1][j][k] - bef->sa.Txx[i][j][k]) / dif->dx;

    aft->va.Vxy[i][j][k] = (2. * Arho - Azetaxy * dif->dt) / (2. * Arho + Azetaxy * dif->dt) * bef->va.Vxy[i][j][k]
      + 2. * dif->dt / (2. * Arho + Azetaxy * dif->dt) * (bef->ta.Txy[i][j][k] - bef->ta.Txy[i][j - 1][k]) / dif->dy;

    aft->va.Vxz[i][j][k] = (2. * Arho - Azetaxz * dif->dt) / (2. * Arho + Azetaxz * dif->dt) * bef->va.Vxz[i][j][k]
      + 2. * dif->dt / (2. * Arho + Azetaxz * dif->dt) * (bef->ta.Tzx[i][j][k] - bef->ta.Tzx[i][j][k - 1]) / dif->dz;
  }
}
// 粒子速度更新関数
__global__ void VyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, VelRan vr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = vr.Vy.x, jmax = vr.Vy.y, kmax = vr.Vy.z;
  double Azetayx, Azetayy, Azetayz, Arho;

  if(i <= imax - 1 && j <= jmax && k <= kmax - 1) {
    Azetayx = (ma->zetayx[i][j + 1][k] + ma->zetayx[i][j][k]) / 2.;
    Azetayy = (ma->zetayy[i][j + 1][k] + ma->zetayy[i][j][k]) / 2.;
    Azetayz = (ma->zetayz[i][j + 1][k] + ma->zetayz[i][j][k]) / 2.;
    Arho    = (   ma->rho[i][j + 1][k] +    ma->rho[i][j][k]) / 2.;
    aft->va.Vyx[i][j][k] = (2. * Arho - Azetayx * dif->dt) / (2. * Arho + Azetayx * dif->dt) * bef->va.Vyx[i][j][k]
      + 2. * dif->dt / (2. * Arho + Azetayx * dif->dt) * (bef->ta.Txy[i][j][k] - bef->ta.Txy[i - 1][j][k]) / dif->dx;
    aft->va.Vyy[i][j][k] = (2. * Arho - Azetayy * dif->dt) / (2. * Arho + Azetayy * dif->dt) * bef->va.Vyy[i][j][k]
      + 2. * dif->dt / (2. * Arho + Azetayy * dif->dt) * (bef->sa.Tyy[i][j + 1][k] - bef->sa.Tyy[i][j][k]) / dif->dy;
    aft->va.Vyz[i][j][k] = (2. * Arho - Azetayz * dif->dt) / (2. * Arho + Azetayz * dif->dt) * bef->va.Vyz[i][j][k]
      + 2. * dif->dt / (2. * Arho + Azetayz * dif->dt) * (bef->ta.Tyz[i][j][k] - bef->ta.Tyz[i][j][k - 1]) / dif->dz;
  }
}
// 粒子速度更新関数
__global__ void VzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, VelRan vr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = vr.Vz.x, jmax = vr.Vz.y, kmax = vr.Vz.z;
  double Azetazx, Azetazy, Azetazz, Arho;

  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax) {
    Azetazx = (ma->zetazx[i][j][k + 1] + ma->zetazx[i][j][k]) / 2.;
    Azetazy = (ma->zetazy[i][j][k + 1] + ma->zetazy[i][j][k]) / 2.;
    Azetazz = (ma->zetazz[i][j][k + 1] + ma->zetazz[i][j][k]) / 2.;
    Arho    = (   ma->rho[i][j][k + 1] +    ma->rho[i][j][k]) / 2.;
    aft->va.Vzx[i][j][k] = (2. * Arho - Azetazx * dif->dt) / (2. * Arho + Azetazx * dif->dt) * bef->va.Vzx[i][j][k]
      + 2. * dif->dt / (2. * Arho + Azetazx * dif->dt) * (bef->ta.Tzx[i][j][k] - bef->ta.Tzx[i - 1][j][k]) / dif->dx;
    aft->va.Vzy[i][j][k] = (2. * Arho - Azetazy * dif->dt) / (2. * Arho + Azetazy * dif->dt) * bef->va.Vzy[i][j][k]
      + 2. * dif->dt / (2. * Arho + Azetazy * dif->dt) * (bef->ta.Tyz[i][j][k] - bef->ta.Tyz[i][j - 1][k]) / dif->dy;
    aft->va.Vzz[i][j][k] = (2. * Arho - Azetazz * dif->dt) / (2. * Arho + Azetazz * dif->dt) * bef->va.Vzz[i][j][k]
      + 2. * dif->dt / (2. * Arho + Azetazz * dif->dt) * (bef->sa.Tzz[i][j][k + 1] - bef->sa.Tzz[i][j][k]) / dif->dz;
  }
}

__global__ void ZeroVx_XY(BefAft *aft, Coord Vmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  if(i <= Vmax.x && j <= Vmax.y){
    aft->va.Vxx[i][j][0] = 0.;
    aft->va.Vxx[i][j][Vmax.z] = 0.;
    aft->va.Vxy[i][j][0] = 0.;
    aft->va.Vxy[i][j][Vmax.z] = 0.;
    aft->va.Vxz[i][j][0] = 0.;
    aft->va.Vxz[i][j][Vmax.z] = 0.;
  }
}

__global__ void ZeroVx_XZ(BefAft *aft, Coord Vmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;

  if(i <= Vmax.x && k <= Vmax.z) {
    aft->va.Vxx[i][0][k] = 0.;
    aft->va.Vxx[i][Vmax.y][k] = 0.;
    aft->va.Vxy[i][0][k] = 0.;
    aft->va.Vxy[i][Vmax.y][k] = 0.;
    aft->va.Vxz[i][0][k] = 0.;
    aft->va.Vxz[i][Vmax.y][k] = 0.;
  }
} 

__global__ void ZeroVy_YX(BefAft *aft, Coord Vmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  if(i <= Vmax.x && j <= Vmax.y){
    aft->va.Vyx[i][j][0] = 0.;
    aft->va.Vyx[i][j][Vmax.z] = 0.;
    aft->va.Vyy[i][j][0] = 0.;
    aft->va.Vyy[i][j][Vmax.z] = 0.;
    aft->va.Vyz[i][j][0] = 0.;
    aft->va.Vyz[i][j][Vmax.z] = 0.;
  }
}

__global__ void ZeroVy_YZ(BefAft *aft, Coord Vmax) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  if(j <= Vmax.y && k <= Vmax.z){
    aft->va.Vyx[0][j][k] = 0.;
    aft->va.Vyx[Vmax.x][j][k] = 0.;
    aft->va.Vyy[0][j][k] = 0.;
    aft->va.Vyy[Vmax.x][j][k] = 0.;
    aft->va.Vyz[0][j][k] = 0.;
    aft->va.Vyz[Vmax.x][j][k] = 0.;
  }
}

__global__ void ZeroVz_ZX(BefAft *aft, Coord Vmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  if(i <= Vmax.x && k <= Vmax.z){
    aft->va.Vzx[i][0][k] = 0.;
    aft->va.Vzx[i][Vmax.y][k] = 0.;
    aft->va.Vzy[i][0][k] = 0.;
    aft->va.Vzy[i][Vmax.y][k] = 0.;
    aft->va.Vzz[i][0][k] = 0.;
    aft->va.Vzz[i][Vmax.y][k] = 0.;
  }
}

__global__ void ZeroVz_ZY(BefAft *aft, Coord Vmax) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  if(j <= Vmax.y && k <= Vmax.z){
    aft->va.Vzx[0][j][k] = 0.;
    aft->va.Vzx[Vmax.x][j][k] = 0.;
    aft->va.Vzy[0][j][k] = 0.;
    aft->va.Vzy[Vmax.x][j][k] = 0.;
    aft->va.Vzz[0][j][k] = 0.;
    aft->va.Vzz[Vmax.x][j][k] = 0.;
  }
}

__global__ void DirectionalAddV(BefAft *aft, Coord Vmax, char check) {
  // スレッドインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if(i > Vmax.x || j > Vmax.y || k > Vmax.z){
    return;
  }
  if(check == 'X') {
    aft->va.Vx[i][j][k] = aft->va.Vxx[i][j][k] + aft->va.Vxy[i][j][k] + aft->va.Vxz[i][j][k];
  } else if(check == 'Y') {
    aft->va.Vy[i][j][k] = aft->va.Vyx[i][j][k] + aft->va.Vyy[i][j][k] + aft->va.Vyz[i][j][k];
  } else if(check == 'Z') {
    aft->va.Vx[i][j][k] = aft->va.Vzx[i][j][k] + aft->va.Vzy[i][j][k] + aft->va.Vzz[i][j][k];
  } else {
    printf("error:DirectionalAddV");
  }
}
// Vxクラス的な
void Vx(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads) {
  // BefAft->ran
  // Range->region,pml
  // pml->ma,dif
  
  // host->deviceデータ転送
  copyBefAftToDevice(aft_d, aft_h, ran);
  copyBefAftToDevice(bef_d, bef_h, ran);
  copyMedArrToDevice(ma_d, &ma_h, ran);
  copyDiffToDevice(dif_d, &dif_h);
  
  int Vximax = ran.vr.Vx.x, Vxjmax = ran.vr.Vx.y, Vxkmax = ran.vr.Vx.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Vximax + threadsPerBlock.x - 1)     / threadsPerBlock.x,
                    (Vxjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Vxkmax - 1 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroXYBlocks((Vximax + threadsPerBlock.x - 1) / threadsPerBlock.x, 
                    (Vxjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroXZBlocks((Vximax + threadsPerBlock.x - 1) / threadsPerBlock.x, 
                    (Vxkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);  
  dim3 DirectionalAddBlocks((Vximax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Vxjmax + threadsPerBlock.y) / threadsPerBlock.y, 
                            (Vxkmax + threadsPerBlock.z) / threadsPerBlock.z);
                
  VxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran.vr);
  ZeroVx_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran.vr.Vx);
  ZeroVx_XZ<<<ZeroXZBlocks, threadsPerBlock>>>(aft_d, ran.vr.Vx);
  DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran.vr.Vx , 'X');
  // device->hostデータ転送
  copyBefAftToHost(aft_h, aft_d, ran);
  copyBefAftToHost(bef_h, bef_d, ran);
}
// Vyクラス的な
void Vy(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads) {
  
  // host->deviceデータ転送
  copyBefAftToDevice(aft_d, aft_h, ran);
  copyBefAftToDevice(bef_d, bef_h, ran);
  copyMedArrToDevice(ma_d, &ma_h, ran);
  copyDiffToDevice(dif_d, &dif_h);

  int Vyimax = ran.vr.Vy.x, Vyjmax = ran.vr.Vy.y, Vykmax = ran.vr.Vy.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Vyimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Vyjmax + threadsPerBlock.y - 1)     / threadsPerBlock.y,
                    (Vykmax - 1 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroYXBlocks((Vyimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, 
                    (Vyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroYZBlocks((Vyjmax + threadsPerBlock.x - 1) / threadsPerBlock.x, 
                    (Vykmax + threadsPerBlock.y - 1) / threadsPerBlock.y);  
  dim3 DirectionalAddBlocks((Vyimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Vyjmax + threadsPerBlock.y) / threadsPerBlock.y, 
                            (Vykmax + threadsPerBlock.z) / threadsPerBlock.z);
  VyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran.vr);
  ZeroVy_YX<<<ZeroYXBlocks, threadsPerBlock>>>(aft_d, ran.vr.Vy);
  ZeroVy_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran.vr.Vy);
  //全方向加算
  DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran.vr.Vy , 'Y');
  // データ転送device to host
  copyBefAftToHost(aft_h, aft_d, ran);
  copyBefAftToHost(bef_h, bef_d, ran);
}
// Vzクラス的な
void Vz(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads) {

  // host->deviceデータ転送
  copyBefAftToDevice(aft_d, aft_h, ran);
  copyBefAftToDevice(bef_d, bef_h, ran);
  copyMedArrToDevice(ma_d, &ma_h, ran);
  copyDiffToDevice(dif_d, &dif_h);

  int Vzimax = ran.vr.Vz.x, Vzjmax = ran.vr.Vz.y, Vzkmax = ran.vr.Vz.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Vzimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Vzjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Vzkmax + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  dim3 ZeroZXBlocks((Vzimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, 
                    (Vzkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroZYBlocks((Vzjmax + threadsPerBlock.x - 1) / threadsPerBlock.x, 
                    (Vzkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);  
  dim3 DirectionalAddBlocks((Vzimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Vzjmax + threadsPerBlock.y) / threadsPerBlock.y, 
                            (Vzkmax + threadsPerBlock.z) / threadsPerBlock.z);                    
  VzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran.vr);
  ZeroVz_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran.vr.Vz);
  ZeroVz_ZY<<<ZeroZYBlocks, threadsPerBlock>>>(aft_d, ran.vr.Vz);
  //全方向加算
  DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran.vr.Vz , 'Z');
  // データ転送device to host
  copyBefAftToHost(aft_h, aft_d, ran);
  copyBefAftToHost(bef_h, bef_d, ran);
}
//粒子速度計算
void Vel(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads) {
  Vx(aft_h, bef_h, aft_d, bef_d, ma_h, ma_d, dif_h, dif_d, ran, threads);
  Vy(aft_h, bef_h, aft_d, bef_d, ma_h, ma_d, dif_h, dif_d, ran, threads);
  Vz(aft_h, bef_h, aft_d, bef_d, ma_h, ma_d, dif_h, dif_d, ran, threads);
}

void Acceleration(Coord_acc *Acc,BefAft *aft, BefAft *bef, Diff dif, Coord out){
  Acc->x = ((aft->va.Vx[out.x - 1][out.y][out.z] - bef->va.Vx[out.x - 1][out.y][out.z]) / dif.dt  + (aft->va.Vx[out.x][out.y][out.z] - bef->va.Vx[out.x][out.y][out.z]) / dif.dt) / 2;
  Acc->y = ((aft->va.Vy[out.x][out.y - 1][out.z] - bef->va.Vy[out.x][out.y - 1][out.z]) / dif.dt  + (aft->va.Vy[out.x][out.y][out.z] - bef->va.Vy[out.x][out.y][out.z]) / dif.dt) / 2;
  Acc->z = ((aft->va.Vz[out.x][out.y][out.z - 1] - bef->va.Vz[out.x][out.y][out.z - 1]) / dif.dt  + (aft->va.Vz[out.x][out.y][out.z] - bef->va.Vz[out.x][out.y][out.z]) / dif.dt) / 2;
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
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        bef->sa.Txx[i][j][k] = aft->sa.Txx[i][j][k];
        bef->sa.Txxx[i][j][k] = aft->sa.Txxx[i][j][k];
        bef->sa.Txxy[i][j][k] = aft->sa.Txxy[i][j][k];
        bef->sa.Txxz[i][j][k] = aft->sa.Txxz[i][j][k];
      }
    }
  }
  for (k = 0; k < Tyykmax; k++) {
    for (j = 0; j < Tyyjmax; j++) {
      for (i = 0; i < Tyyimax; i++) {
        bef->sa.Tyy[i][j][k] = aft->sa.Tyy[i][j][k];
        bef->sa.Tyyx[i][j][k] = aft->sa.Tyyx[i][j][k];
        bef->sa.Tyyy[i][j][k] = aft->sa.Tyyy[i][j][k];
        bef->sa.Tyyz[i][j][k] = aft->sa.Tyyz[i][j][k];
      }
    }
  }
  for (k = 0; k < Tzzkmax; k++) {
    for (j = 0; j < Tzzjmax; j++) {
      for (i = 0; i < Tzzimax; i++) {
        bef->sa.Tzz[i][j][k] = aft->sa.Tzz[i][j][k];
        bef->sa.Tzzx[i][j][k] = aft->sa.Tzzx[i][j][k];
        bef->sa.Tzzy[i][j][k] = aft->sa.Tzzy[i][j][k];
        bef->sa.Tzzz[i][j][k] = aft->sa.Tzzz[i][j][k];
      }
    }
  }
  for (i = 0; i < Txyimax; i++) {
    for (j = 0; j < Txyjmax; j++) {
      for (k = 0; k < Txykmax; k++) {
        bef->ta.Txy[i][j][k] = aft->ta.Txy[i][j][k];
        bef->ta.Txyx[i][j][k] = aft->ta.Txyx[i][j][k];
        bef->ta.Txyy[i][j][k] = aft->ta.Txyy[i][j][k];
      }
    }
  }
  for (i = 0; i < Tyzimax; i++) {
    for (j = 0; j < Tyzjmax; j++) {
      for (k = 0; k < Tyzkmax; k++) {
        bef->ta.Tyz[i][j][k] = aft->ta.Tyz[i][j][k];
        bef->ta.Tyzy[i][j][k] = aft->ta.Tyzy[i][j][k];
        bef->ta.Tyzz[i][j][k] = aft->ta.Tyzz[i][j][k];
      }
    }
  }
  for (i = 0; i < Tzximax; i++) {
    for (j = 0; j < Tzxjmax; j++) {
      for (k = 0; k < Tzxkmax; k++) {
        bef->ta.Tzx[i][j][k] = aft->ta.Tzx[i][j][k];
        bef->ta.Tzxz[i][j][k] = aft->ta.Tzxz[i][j][k];
        bef->ta.Tzxx[i][j][k] = aft->ta.Tzxx[i][j][k];
      }
    }
  }
  for (k = 0; k < Vxkmax; k++) {
    for (j = 0; j < Vxjmax; j++) {
      for (i = 0; i < Vximax; i++) {
        bef->va.Vx[i][j][k] = aft->va.Vx[i][j][k];
        bef->va.Vxx[i][j][k] = aft->va.Vxx[i][j][k];
        bef->va.Vxy[i][j][k] = aft->va.Vxy[i][j][k];
        bef->va.Vxz[i][j][k] = aft->va.Vxz[i][j][k];
      }
    }
  }
  for (k = 0; k < Vykmax; k++) {
    for (j = 0; j < Vyjmax; j++) {
      for (i = 0; i < Vyimax; i++) {
        bef->va.Vy[i][j][k] = aft->va.Vy[i][j][k];
        bef->va.Vyx[i][j][k] = aft->va.Vyx[i][j][k];
        bef->va.Vyy[i][j][k] = aft->va.Vyy[i][j][k];
        bef->va.Vyz[i][j][k] = aft->va.Vyz[i][j][k];
      }
    }
  }
  for (k = 0; k < Vzkmax; k++) {
    for (j = 0; j < Vzjmax; j++) {
      for (i = 0; i < Vzimax; i++) {
        bef->va.Vz[i][j][k] = aft->va.Vz[i][j][k];
        bef->va.Vzx[i][j][k] = aft->va.Vzx[i][j][k];
        bef->va.Vzy[i][j][k] = aft->va.Vzy[i][j][k];
        bef->va.Vzz[i][j][k] = aft->va.Vzz[i][j][k];
      }
    }
  }
}