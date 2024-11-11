#define _USE_MATH_DEFINES
#include "../header/update.h"
#include "../header/struct.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/init.h"
#include "../header/memory.h"

// 垂直応力

// 垂直応力更新並列関数
__global__ void TxxUpdate(SigArr aftsa, VelArr aftva, SigArr befsa, VelArr befva, MedArr *ma, Diff *dif, Range *ran, ImpulseArr *ipa) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;

  if(i < imax - 1 && j < jmax - 1 && k < kmax - 1) {
    // int idx    = Dgetid<<<1,1>>>(ran->sr.Txx, i, j, k);
    // int idx_i1 = Dgetid<<<1,1>>>(ran->sr.Txx, i + 1, j, k);
    // int idx_j1 = Dgetid<<<1,1>>>(ran->sr.Txx, i, j + 1, k);
    // int idx_k1 = Dgetid<<<1,1>>>(ran->sr.Txx, i, j, k + 1);
    int id    = k * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + i;
    int idvx  = k * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + i;
    int idvxi = k * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + (i + 1);
    int idvy  = k * ran->vr.Vy.x * ran->vr.Vy.y + j * ran->vr.Vy.x + i;
    int idvyj = k * ran->vr.Vy.x * ran->vr.Vy.y + (j + 1) * ran->vr.Vy.x + i;
    int idvz  = k * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i;
    int idvzk = (k + 1) * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i;

    aftsa.Txxx[id] = (2.0 - ma[id].zetadx * dif->dt) / (2.0 + ma[id].zetadx * dif->dt) * befsa.Txxx[id]
        + 2.0 * (ma[id].c11 * dif->dt + ma[id].xi11) / (2.0 + ma[id].zetadx * dif->dt) * (aftva.Vx[idvxi] - aftva.Vx[idvx]) / dif->dx
        - 2.0 * ma[id].xi11 / (2.0 + ma[id].zetadx * dif->dt) * (befva.Vx[idvxi] - befva.Vx[idvx]) / dif->dx;

    aftsa.Txxy[id] = (2.0 - ma[id].zetady * dif->dt) / (2.0 + ma[id].zetady * dif->dt) * befsa.Txxy[id]
        + 2.0 * (ma[id].ramda * dif->dt + ma[id].khi) / (2.0 + ma[id].zetady * dif->dt) * (aftva.Vy[idvyj] - aftva.Vy[idvy]) / dif->dy
        - 2.0 * ma[id].khi / (2.0 + ma[id].zetady * dif->dt) * (befva.Vy[idvyj] - befva.Vy[idvy]) / dif->dy;

    aftsa.Txxz[id] = (2.0 - ma[id].zetadz * dif->dt) / (2.0 + ma[id].zetadz * dif->dt) * befsa.Txxz[id]
        + 2.0 * (ma[id].ramda * dif->dt + ma[id].khi) / (2.0 + ma[id].zetadz * dif->dt) * (aftva.Vz[idvzk] - aftva.Vz[idvz]) / dif->dz
        - 2.0 * ma[id].khi / (2.0 + ma[id].zetadz * dif->dt) * (befva.Vz[idvzk] - befva.Vz[idvz]) / dif->dz;

    aftsa.Txx[id] = aftsa.Txxx[id] + aftsa.Txxy[id] + aftsa.Txxz[id] + ipa[id].Txx;
  }
}

// 垂直応力更新並列関数
__global__ void TyyUpdate(SigArr aftsa, VelArr aftva, SigArr befsa, VelArr befva, MedArr *ma, Diff *dif, Range *ran, ImpulseArr *ipa) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = ran->sr.Tyy.x, jmax = ran->sr.Tyy.y, kmax = ran->sr.Tyy.z;

  if(i < imax - 1 && j < jmax - 1 && k < kmax - 1) {
    // int idx    = Dgetid<<<1,1>>>(ran->sr.Tyy, i, j, k);
    // int idx_i1 = Dgetid<<<1,1>>>(ran->sr.Tyy, i + 1, j, k);
    // int idx_j1 = Dgetid<<<1,1>>>(ran->sr.Tyy, i, j + 1, k);
    // int idx_k1 = Dgetid<<<1,1>>>(ran->sr.Tyy, i, j, k + 1);
    int id    = k * ran->sr.Tyy.x * ran->sr.Tyy.y + j * ran->sr.Tyy.x + i;
    int idvx  = k * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + i;
    int idvxi = k * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + (i + 1);
    int idvy  = k * ran->vr.Vy.x * ran->vr.Vy.y + j * ran->vr.Vy.x + i;
    int idvyj = k * ran->vr.Vy.x * ran->vr.Vy.y + (j + 1) * ran->vr.Vy.x + i;
    int idvz  = k * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i;
    int idvzk = (k + 1) * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i;

    aftsa.Tyyx[id] = (2.0 - ma[id].zetadx * dif->dt) / (2.0 + ma[id].zetadx * dif->dt) * befsa.Tyyx[id]
        + 2.0 * (ma[id].ramda * dif->dt + ma[id].khi) / (2.0 + ma[id].zetadx * dif->dt) * (aftva.Vx[idvxi] - aftva.Vx[idvx]) / dif->dx
        - 2.0 * ma[id].khi / (2.0 + ma[id].zetadx * dif->dt) * (befva.Vx[idvxi] - befva.Vx[idvx]) / dif->dx;

    aftsa.Tyyy[id] = (2.0 - ma[id].zetady * dif->dt) / (2.0 + ma[id].zetady * dif->dt) * befsa.Tyyy[id]
        + 2.0 * (ma[id].c11 * dif->dt + ma[id].xi11) / (2.0 + ma[id].zetady * dif->dt) * (aftva.Vy[idvyj] - aftva.Vy[idvy]) / dif->dy
        - 2.0 * ma[id].xi11 / (2.0 + ma[id].zetady * dif->dt) * (befva.Vy[idvyj] - befva.Vy[idvy]) / dif->dy;

    aftsa.Tyyz[id] = (2.0 - ma[id].zetadz * dif->dt) / (2.0 + ma[id].zetadz * dif->dt) * befsa.Tyyz[id]
        + 2.0 * (ma[id].ramda * dif->dt + ma[id].khi) / (2.0 + ma[id].zetadz * dif->dt) * (aftva.Vz[idvzk] - aftva.Vz[idvz]) / dif->dz
        - 2.0 * ma[id].khi / (2.0 + ma[id].zetadz * dif->dt) * (befva.Vz[idvzk] - befva.Vz[idvz]) / dif->dz;

    aftsa.Tyy[id] = aftsa.Tyyx[id] + aftsa.Tyyy[id] + aftsa.Tyyz[id] + ipa[id].Tyy;
  }
}

// 垂直応力更新並列関数
__global__ void TzzUpdate(SigArr aftsa, VelArr aftva, SigArr befsa, VelArr befva, MedArr *ma, Diff *dif, Range *ran, ImpulseArr *ipa) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
  int imax = ran->sr.Tzz.x, jmax = ran->sr.Tzz.y, kmax = ran->sr.Tzz.z;
  
  if(i < imax - 1 && j < jmax - 1 && k < kmax - 1) {
    // int idx     = Dgetid<<<1,1>>>(ran->sr.Tzz, i, j, k);
    // int idx_im1 = Dgetid<<<1,1>>>(ran->sr.Tzz, i + 1, j, k);
    // int idx_jm1 = Dgetid<<<1,1>>>(ran->sr.Tzz, i, j + 1, k);
    // int idx_km1 = Dgetid<<<1,1>>>(ran->sr.Tzz, i, j, k + 1);
    int id    = k * ran->sr.Tzz.x * ran->sr.Tzz.y + j * ran->sr.Tzz.x + i;
    int idvx  = k * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + i;
    int idvxi = k * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + (i + 1);
    int idvy  = k * ran->vr.Vy.x * ran->vr.Vy.y + j * ran->vr.Vy.x + i;
    int idvyj = k * ran->vr.Vy.x * ran->vr.Vy.y + (j + 1) * ran->vr.Vy.x + i;
    int idvz  = k * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i;
    int idvzk = (k + 1) * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i;

    aftsa.Tzzx[id] = (2.0 - ma[id].zetadx * dif->dt) / (2.0 + ma[id].zetadx * dif->dt) * befsa.Tzzx[id]
        + 2.0 * (ma[id].ramda * dif->dt + ma[id].khi) / (2.0 + ma[id].zetadx * dif->dt) * (aftva.Vx[idvxi] - aftva.Vx[idvx]) / dif->dx
        - 2.0 * ma[id].khi / (2.0 + ma[id].zetadx * dif->dt) * (befva.Vx[idvxi] - befva.Vx[idvx]) / dif->dx;

    aftsa.Tzzy[id] = (2.0 - ma[id].zetady * dif->dt) / (2.0 + ma[id].zetady * dif->dt) * befsa.Tzzy[id]
        + 2.0 * (ma[id].ramda * dif->dt + ma[id].khi) / (2.0 + ma[id].zetady * dif->dt) * (aftva.Vy[idvyj] - aftva.Vy[idvy]) / dif->dy
        - 2.0 * ma[id].khi / (2.0 + ma[id].zetady * dif->dt) * (befva.Vy[idvyj] - befva.Vy[idvy]) / dif->dy;

    aftsa.Tzzz[id] = (2.0 - ma[id].zetadz * dif->dt) / (2.0 + ma[id].zetadz * dif->dt) * befsa.Tzzz[id]
        + 2.0 * (ma[id].c11 * dif->dt + ma[id].xi11) / (2.0 + ma[id].zetadz * dif->dt) * (aftva.Vz[idvzk] - aftva.Vz[idvz]) / dif->dz
        - 2.0 * ma[id].xi11 / (2.0 + ma[id].zetadz * dif->dt) * (befva.Vz[idvzk] - befva.Vz[idvz]) / dif->dz;

    aftsa.Tzz[id] = aftsa.Tzzx[id] + aftsa.Tzzy[id] + aftsa.Tzzz[id] + ipa[id].Tzz;
  }
}

// 全方向加算
__global__ void DirectionalAdd(SigArr aftsa, ImpulseArr *ipa, Range *ran, char check) {
  // 1Dインデックスの計算
  // 全範囲でもいい(端は計算されていないためずっと0)
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;

  // 1Dインデックス化
  // int id = Dgetid<<<1,1>>>(ran->sr.Txx, i, j, k);
  int id = k * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + i;
  if (i < imax - 1 && j < jmax - 1 && k < kmax - 1) {
    // 各方向に応じた計算を実行（ポインタ表記）
    if (check == 'X') {
      aftsa.Txx[id] = aftsa.Txxx[id] + aftsa.Txxy[id] + aftsa.Txxz[id] + ipa[id].Txx;
    } else if (check == 'Y') {
      aftsa.Tyy[id] = aftsa.Tyyx[id] + aftsa.Tyyy[id] + aftsa.Tyyz[id] + ipa[id].Tyy;
    } else if (check == 'Z') {
      aftsa.Tzz[id] = aftsa.Tzzx[id] + aftsa.Tzzy[id] + aftsa.Tzzz[id] + ipa[id].Tzz;
    } else {
      printf("error: DirectionalAdd\n");
    }
  }
}

__global__ void createImpulse(ImpulseArr *ipa, Impulse *ip, Diff *dif, Range *ran, int innum, int t) {
  int id;
  for(int num = 0; num < innum; num++) {
    // id = Dgetid<<<1,1>>>(ran->sr.Txx, ip[num].in.x, ip[num].in.y, ip[num].in.z);
    id = ip[num].in.z * ran->sr.Txx.x * ran->sr.Txx.y + ip[num].in.y * ran->sr.Txx.x + ip[num].in.x;
    if (ip->mode == E_SINE) {
      ipa[id].Txx = 0;
      ipa[id].Tyy = 0;
      ipa[id].Tzz = 8.e3 * 0.5 * sin(2. * M_PI * ip->freq * (double)t * dif->dt) / 2.;
    } else if (ip->mode == E_RCOS) {
      if (t < 1. / ip->freq / dif->dt) {
        ipa[id].Txx = 0;
        ipa[id].Tyy = 0;
        ipa[id].Tzz = 8.e3 * 0.5 * (1. - cos(2. * M_PI * ip->freq * (double)t * dif->dt)) / 2.;
      } else {
        ipa[id].Txx = 0;
        ipa[id].Tyy = 0;
        ipa[id].Tzz = 0;
      }
    } else {
      ipa[id].Txx = 0;
      ipa[id].Tyy = 0;
      ipa[id].Tzz = 0;
    }
  }
}

// Txxクラス的な(Blocks大丈夫かな？)
void Txx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, ImpulseArr *ipa_d, DimI3 threads) {
  // cudaError_t err;
  // char check = 'X';

  int Txximax = ran_h->sr.Txx.x, Txxjmax = ran_h->sr.Txx.y, Txxkmax = ran_h->sr.Txx.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z);  // ブロック内のスレッド数
  dim3 UpdateBlocks((Txximax - 2 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Txxjmax - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Txxkmax - 2 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // dim3 ZeroXYBlocks((Txximax     + threadsPerBlock.x - 1) / threadsPerBlock.x, (Txxjmax     + threadsPerBlock.y - 1) / threadsPerBlock.y);
  // dim3 ZeroYZBlocks((Txxjmax     + threadsPerBlock.x - 1) / threadsPerBlock.x, (Txxkmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  // dim3 ZeroZXBlocks((Txxkmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Txximax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Txximax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                            (Txxjmax + threadsPerBlock.y - 1) / threadsPerBlock.y, 
                            (Txxkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // Txx更新式
  TxxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d->sa, aft_d->va, bef_d->sa, bef_d->va, ma_d, dif_d, ran_d, ipa_d);
  // 0 padding
  // ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d, 'X');
  // ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d, 'X');
  // ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d, 'X');
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txx zero  : %s\n", cudaGetErrorString(err));
  //全方向加算
  // DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d->sa, ipa_d, ran_d, 'X');
  // cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txx add   : %s\n", cudaGetErrorString(err));

}
// Tyyクラス的な(Blocks大丈夫かな？)
void Tyy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, ImpulseArr *ipa_d, DimI3 threads) {
  // cudaError_t err;
  char check = 'Y';

  int Tyyimax = ran_h->sr.Tyy.x, Tyyjmax = ran_h->sr.Tyy.y, Tyykmax = ran_h->sr.Tyy.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z);  // ブロック内のスレッド数
  dim3 UpdateBlocks((Tyyimax - 2 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tyyjmax - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tyykmax - 2 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // dim3 ZeroXYBlocks((Tyyimax     + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tyyjmax     + threadsPerBlock.y - 1) / threadsPerBlock.y);
  // dim3 ZeroYZBlocks((Tyyjmax     + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tyykmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  // dim3 ZeroZXBlocks((Tyykmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tyyimax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Tyyimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                            (Tyyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                            (Tyykmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // Tyy更新式
  TyyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d->sa, aft_d->va, bef_d->sa, bef_d->va, ma_d, dif_d, ran_d, ipa_d);
  // 0 padding
  // ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  // ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  // ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyy zero  : %s\n", cudaGetErrorString(err));
  // 全方向加算
  // DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d->sa, ipa_d, ran_d, check);
  // cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyy add   : %s\n", cudaGetErrorString(err));
}
// Tzzクラス的な(Blocks大丈夫かな？)
void Tzz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, ImpulseArr *ipa_d, DimI3 threads) {
  // cudaError_t err;
  char check = 'Z';

  int Tzzimax = ran_h->sr.Tzz.x, Tzzjmax = ran_h->sr.Tzz.y, Tzzkmax = ran_h->sr.Tzz.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z);  // ブロック内のスレッド数
  dim3 UpdateBlocks((Tzzimax - 2 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tzzjmax - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tzzkmax - 2 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // dim3 ZeroXYBlocks((Tzzimax     + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tzzjmax     + threadsPerBlock.y - 1) / threadsPerBlock.y);
  // dim3 ZeroYZBlocks((Tzzjmax     + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tzzkmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  // dim3 ZeroZXBlocks((Tzzkmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tzzimax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Tzzimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                            (Tzzjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                            (Tzzkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // Tzzの更新式
  TzzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d->sa, aft_d->va, bef_d->sa, bef_d->va, ma_d, dif_d, ran_d, ipa_d);
  // 0 padding
  // ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  // ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  // ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzz zero  : %s\n", cudaGetErrorString(err));
  // 全方向加算
  // DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d->sa, ipa_d, ran_d, check);
  // cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzz add   : %s\n", cudaGetErrorString(err));
 
}
// 垂直応力計算(main呼び出し関数)
void Sig(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, ImpulseArr *ipa_d, DimI3 threads) {
  Txx(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, ipa_d, threads);
  Tyy(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, ipa_d, threads);
  Tzz(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, ipa_d, threads);
}

// せん断応力

// せん断応力更新関数
__global__ void TxyUpdate(TauArr aftta, VelArr aftva, TauArr befta, VelArr befva, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1; // 始点を+1

  // int imax = ran->tr.Txy.x, jmax = ran->tr.Txy.y, kmax = ran->tr.Txy.z;
  double Hzetadx, Hzetady, Hmu, Hgamma;

  if (i < ran->tr.Txy.x - 1 && j < ran->tr.Txy.y - 1 && k < ran->tr.Txy.z - 1) {
    // 各インデックスの計算
    // int idtxy = Dgetid<<<1,1>>>(ran->tr.Txy, i, j, k);

    // int id   = Dgetid<<<1,1>>>(ran->sr.Txx, i, j, k);
    // int idi  = Dgetid<<<1,1>>>(ran->sr.Txx, i - 1, j, k);
    // int idj  = Dgetid<<<1,1>>>(ran->sr.Txx, i, j - 1, k);
    // int idij = Dgetid<<<1,1>>>(ran->sr.Txx, i - 1, j - 1, k);

    // int idvx  = Dgetid<<<1,1>>>(ran->vr.Vx, i, j, k);
    // int idvxj = Dgetid<<<1,1>>>(ran->vr.Vx, i, j - 1, k);
    // int idvy  = Dgetid<<<1,1>>>(ran->vr.Vy, i, j, k);
    // int idvyi = Dgetid<<<1,1>>>(ran->vr.Vy, i - 1, j, k);
    int idtxy = k * ran->tr.Txy.x * ran->tr.Txy.y + j * ran->tr.Txy.x + i;

    int id = k * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + i;
    int idi = k * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + (i - 1);
    int idj = k * ran->sr.Txx.x * ran->sr.Txx.y + (j - 1) * ran->sr.Txx.x + i;
    int idij = k * ran->sr.Txx.x * ran->sr.Txx.y + (j - 1) * ran->sr.Txx.x + (i - 1);

    int idvx  = k * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + i;
    int idvxj = k * ran->vr.Vx.x * ran->vr.Vx.y + (j - 1) * ran->vr.Vx.x + i;
    int idvy  = k * ran->vr.Vy.x * ran->vr.Vy.y + j * ran->vr.Vy.x + i;
    int idvyi = k * ran->vr.Vy.x * ran->vr.Vy.y + j * ran->vr.Vy.x + (i - 1);


    // PML:減衰係数,計算領域:摩擦定数
    Hzetadx = 4. * pow((1. / ma[id].zetadx) + (1. / ma[idi].zetadx) + (1. / ma[idj].zetadx) + (1. / ma[idij].zetadx), -1.);
    Hzetady = 4. * pow((1. / ma[id].zetady) + (1. / ma[idi].zetady) + (1. / ma[idj].zetady) + (1. / ma[idij].zetady), -1.);
    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / ma[id].mu) + (1. / ma[idi].mu) + (1. / ma[idj].mu) + (1. / ma[idij].mu), -1.);
    // 第1粘性定数
    Hgamma = 4. * pow((1. / ma[id].gamma) + (1. / ma[idi].gamma) + (1. / ma[idj].gamma) + (1. / ma[idij].gamma), -1.);


    aftta.Txyx[idtxy] = (2.0 - Hzetadx * dif->dt) / (2.0 + Hzetadx * dif->dt) * befta.Txyx[idtxy]
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadx * dif->dt) * (aftva.Vy[idvy] - aftva.Vy[idvyi]) / dif->dx
        - 2.0 * Hgamma / (2.0 + Hzetadx * dif->dt) * (befva.Vy[idvy] - befva.Vy[idvyi]) / dif->dx;

    aftta.Txyy[idtxy] = (2.0 - Hzetady * dif->dt) / (2.0 + Hzetady * dif->dt) * befta.Txyy[idtxy]
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetady * dif->dt) * (aftva.Vx[idvx] - aftva.Vx[idvxj]) / dif->dy
        - 2.0 * Hgamma / (2.0 + Hzetady * dif->dt) * (befva.Vx[idvx] - befva.Vx[idvxj]) / dif->dy;
    aftta.Txy[idtxy] = aftta.Txyx[idtxy] + aftta.Txyy[idtxy];
  }
}

// せん断応力更新関数
__global__ void TyzUpdate(TauArr aftta, VelArr aftva, TauArr befta, VelArr befva, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  // int imax = ran->tr.Tyz.x, jmax = ran->tr.Tyz.y, kmax = ran->tr.Tyz.z;
  double Hzetady, Hzetadz, Hmu, Hgamma;

  if (i < ran->tr.Tyz.x - 1 && j < ran->tr.Tyz.y - 1 && k < ran->tr.Tyz.z - 1) {
    // 各インデックスの計算
    // int idtyz = Dgetid<<<1,1>>>(ran->tr.Tyz, i, j, k);

    // int id   = Dgetid<<<1,1>>>(ran->sr.Txx, i, j, k);
    // int idj  = Dgetid<<<1,1>>>(ran->sr.Txx, i, j - 1, k);
    // int idk  = Dgetid<<<1,1>>>(ran->sr.Txx, i, j, k - 1);
    // int idjk = Dgetid<<<1,1>>>(ran->sr.Txx, i, j - 1, k - 1);

    // int idvy  = Dgetid<<<1,1>>>(ran->vr.Vy, i, j, k);
    // int idvyk = Dgetid<<<1,1>>>(ran->vr.Vy, i, j, k - 1);
    // int idvz  = Dgetid<<<1,1>>>(ran->vr.Vz, i, j, k);
    // int idvzj = Dgetid<<<1,1>>>(ran->vr.Vz, i, j - 1, k);
    int idtyz = k * ran->tr.Tyz.x * ran->tr.Tyz.y + j * ran->tr.Tyz.x + i;

    int id = k * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + i;
    int idj = k * ran->sr.Txx.x * ran->sr.Txx.y + (j - 1) * ran->sr.Txx.x + i;
    int idk = (k - 1) * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + i;
    int idjk = (k - 1) * ran->sr.Txx.x * ran->sr.Txx.y + (j - 1) * ran->sr.Txx.x + i;

    int idvy = k * ran->vr.Vy.x * ran->vr.Vy.y + j * ran->vr.Vy.x + i;
    int idvyk = (k - 1) * ran->vr.Vy.x * ran->vr.Vy.y + j * ran->vr.Vy.x + i;
    int idvz = k * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i;
    int idvzj = k * ran->vr.Vz.x * ran->vr.Vz.y + (j - 1) * ran->vr.Vz.x + i;


    // PML:減衰係数,計算領域:摩擦定数
    Hzetady = 4. * pow((1. / ma[id].zetady) + (1. / ma[idj].zetady) + (1. / ma[idk].zetady) + (1. / ma[idjk].zetady), -1.);
    Hzetadz = 4. * pow((1. / ma[id].zetadz) + (1. / ma[idj].zetadz) + (1. / ma[idk].zetadz) + (1. / ma[idjk].zetadz), -1.);
    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / ma[id].mu) + (1. / ma[idj].mu) + (1. / ma[idk].mu) + (1. / ma[idjk].mu), -1.);
    // 第1粘性定数
    Hgamma = 4. * pow((1. / ma[id].gamma) + (1. / ma[idj].gamma) + (1. / ma[idk].gamma) + (1. / ma[idjk].gamma), -1.);

    aftta.Tyzy[idtyz] = (2.0 - Hzetady * dif->dt) / (2.0 + Hzetady * dif->dt) * befta.Tyzy[idtyz]
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetady * dif->dt) * (aftva.Vz[idvz] - aftva.Vz[idvzj]) / dif->dy
        - 2.0 * Hgamma / (2.0 + Hzetady * dif->dt) * (befva.Vz[idvz] - befva.Vz[idvzj]) / dif->dy;

    aftta.Tyzz[idtyz] = (2.0 - Hzetadz * dif->dt) / (2.0 + Hzetadz * dif->dt) * befta.Tyzz[idtyz]
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadz * dif->dt) * (aftva.Vy[idvy] - aftva.Vy[idvyk]) / dif->dz
        - 2.0 * Hgamma / (2.0 + Hzetadz * dif->dt) * (befva.Vy[idvy] - befva.Vy[idvyk]) / dif->dz;
    aftta.Tyz[idtyz] = aftta.Tyzy[idtyz] + aftta.Tyzz[idtyz];
    
  }
}

// せん断応力更新関数
__global__ void TzxUpdate(TauArr aftta, VelArr aftva, TauArr befta, VelArr befva, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  double Hzetadx, Hzetadz, Hmu, Hgamma;

  if (i < ran->tr.Tzx.x - 1 && j < ran->tr.Tzx.y - 1 && k < ran->tr.Tzx.z - 1) {
    // 各インデックスの計算
    // 求めるTzx
    int idtzx  = k * ran->tr.Tzx.x * ran->tr.Tzx.y + j * ran->tr.Tzx.x + i;
    // 格子点
    int id     = k * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + i;
    int idi    = k * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + (i - 1);
    int idk    = (k - 1) * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + i;
    int idki   = (k - 1) * ran->sr.Txx.x * ran->sr.Txx.y + j * ran->sr.Txx.x + (i - 1);
    // 速度点
    int idvx   = k * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + i;
    int idvxk  = (k - 1) * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + i;
    int idvz   = k * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i;
    int idvzi  = k * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + (i - 1);

    // PML:減衰係数,計算領域:摩擦定数
    Hzetadx = 4. * pow((1. / ma[id].zetadx) + (1. / ma[idi].zetadx) + (1. / ma[idk].zetadx) + (1. / ma[idki].zetadx), -1.);
    Hzetadz = 4. * pow((1. / ma[id].zetadz) + (1. / ma[idi].zetadz) + (1. / ma[idk].zetadz) + (1. / ma[idki].zetadz), -1.);
    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / ma[id].mu) + (1. / ma[idi].mu) + (1. / ma[idk].mu) + (1. / ma[idki].mu), -1.);
    // 第1粘性定数
    Hgamma = 4. * pow((1. / ma[id].gamma) + (1. / ma[idi].gamma) + (1. / ma[idk].gamma) + (1. / ma[idki].gamma), -1.);

    aftta.Tzxz[idtzx] = (2.0 - Hzetadz * dif->dt) / (2.0 + Hzetadz * dif->dt) * befta.Tzxz[idtzx]
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadz * dif->dt) * (aftva.Vx[idvx] - aftva.Vx[idvxk]) / dif->dz
        - 2.0 * Hgamma / (2.0 + Hzetadz * dif->dt) * (befva.Vx[idvx] - befva.Vx[idvxk]) / dif->dz;

    aftta.Tzxx[idtzx] = (2.0 - Hzetadx * dif->dt) / (2.0 + Hzetadx * dif->dt) * befta.Tzxx[idtzx]
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadx * dif->dt) * (aftva.Vz[idvz] - aftva.Vz[idvzi]) / dif->dx
        - 2.0 * Hgamma / (2.0 + Hzetadx * dif->dt) * (befva.Vz[idvz] - befva.Vz[idvzi]) / dif->dx;
    aftta.Tzx[idtzx] = aftta.Tzxx[idtzx] + aftta.Tzxz[idtzx];
  }
}

__global__ void DirectionalAddT(TauArr aftta, Range *ran, char check) {
  // スレッドインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int id;
  if (check == 'X' && i < ran->tr.Tyz.x - 1 && j < ran->tr.Tyz.y - 1 && k < ran->tr.Tyz.z - 1) {
    id = k * ran->tr.Tyz.x * ran->tr.Tyz.y + j * ran->tr.Tyz.x + i;
    aftta.Tyz[id] = aftta.Tyzy[id] + aftta.Tyzz[id];
  } else if (check == 'Y' && i < ran->tr.Tzx.x - 1 && j < ran->tr.Tzx.y - 1 && k < ran->tr.Tzx.z - 1) {
    id = k * ran->tr.Tzx.x * ran->tr.Tzx.y + j * ran->tr.Tzx.x + i;
    aftta.Tzx[id] = aftta.Tzxx[id] + aftta.Tzxz[id];
  } else if (check == 'Z' && i < ran->tr.Txy.x - 1 && j < ran->tr.Txy.y - 1 && k < ran->tr.Txy.z - 1) {
    id = k * ran->tr.Txy.x * ran->tr.Txy.y + j * ran->tr.Txy.x + i;
    aftta.Txy[id] = aftta.Txyx[id] + aftta.Txyy[id];
  }
}

// Txyクラス的な
void Txy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {
  // cudaError_t err;
  int Txyimax = ran_h->tr.Txy.x, Txyjmax = ran_h->tr.Txy.y, Txykmax = ran_h->tr.Txy.z;
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Txyimax - 2 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Txyjmax - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Txykmax - 2 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // dim3 ZeroXYBlocks((Txyimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Txyjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Txyimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                            (Txyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                            (Txykmax + threadsPerBlock.z - 1) / threadsPerBlock.z);                    
  TxyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d->ta, aft_d->va, bef_d->ta, bef_d->va, ma_d, dif_d, ran_d);
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txy update: %s\n", cudaGetErrorString(err));
  // ZeroTxy<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d);
  cudaDeviceSynchronize();

  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txy zero  : %s\n", cudaGetErrorString(err));
  // DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d->ta, ran_d, 'Z');
  // cudaDeviceSynchronize();

  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txy add   : %s\n", cudaGetErrorString(err));
}
// Tyzクラス的な
void Tyz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {
  // cudaError_t err;
  int Tyzimax = ran_h->tr.Tyz.x, Tyzjmax = ran_h->tr.Tyz.y, Tyzkmax = ran_h->tr.Tyz.z;
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Tyzimax - 2 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tyzjmax - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tyzkmax - 2 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // dim3 ZeroYZBlocks((Tyzjmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tyzkmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Tyzimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                            (Tyzjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                            (Tyzkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);                    
  TyzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d->ta, aft_d->va, bef_d->ta, bef_d->va, ma_d, dif_d, ran_d);
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyz update: %s\n", cudaGetErrorString(err));
  // ZeroTyz<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyz zero  : %s\n", cudaGetErrorString(err));
  // DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d->ta, ran_d, 'X');
  // cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyz add   : %s\n", cudaGetErrorString(err));
}
// Tzxクラス的な
void Tzx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {
  // cudaError_t err;
  int Tzximax = ran_h->tr.Tzx.x, Tzxjmax = ran_h->tr.Tzx.y, Tzxkmax = ran_h->tr.Tzx.z;
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Tzximax - 2 + threadsPerBlock.x - 1)     / threadsPerBlock.x,
                    (Tzxjmax - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tzxkmax - 2 + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  // dim3 ZeroZXBlocks((Tzxkmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Tzximax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);   
  dim3 DirectionalAddBlocks((Tzximax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                            (Tzxjmax + threadsPerBlock.y - 1) / threadsPerBlock.y, 
                            (Tzxkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);                  
  TzxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d->ta, aft_d->va, bef_d->ta, bef_d->va, ma_d, dif_d, ran_d);
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzx update: %s\n", cudaGetErrorString(err));
  // ZeroTzx<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzx zero  : %s\n", cudaGetErrorString(err));
  // DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d->ta, ran_d, 'Y');
  // cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzx add   : %s\n", cudaGetErrorString(err));

}
// せん断応力計算(main呼び出し関数)
void Tau(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {
  Txy(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
  Tyz(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
  Tzx(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
}

// 粒子速度

// 粒子速度更新関数
__global__ void VxUpdate(VelArr aftva, VelArr befva, SigArr befsa, TauArr befta, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  double Azetaxx, Azetaxy, Azetaxz, Arho;

  if(i < ran->vr.Vx.x - 1 && j < ran->vr.Vx.y - 1 && k < ran->vr.Vx.z - 1) {

    int idvx   =       k * ran->vr.Vx.x  * ran->vr.Vx.y  +       j * ran->vr.Vx.x  + i;
    int idtxx  =       k * ran->sr.Txx.x * ran->sr.Txx.y +       j * ran->sr.Txx.x + i;
    int idtxxi =       k * ran->sr.Txx.x * ran->sr.Txx.y +       j * ran->sr.Txx.x + (i - 1);
    int idtxy  =       k * ran->tr.Txy.x * ran->tr.Txy.y +       j * ran->tr.Txy.x + i;
    int idtxyj =       k * ran->tr.Txy.x * ran->tr.Txy.y + (j + 1) * ran->tr.Txy.x + i;
    int idtzx  =       k * ran->tr.Tzx.x * ran->tr.Tzx.y +       j * ran->tr.Tzx.x + i;
    int idtzxk = (k + 1) * ran->tr.Tzx.x * ran->tr.Tzx.y +       j * ran->tr.Tzx.x + i;

    Azetaxx = (ma[idtxx].zetaxx + ma[idtxxi].zetaxx) / 2.0;
    Azetaxy = (ma[idtxx].zetaxy + ma[idtxxi].zetaxy) / 2.0;
    Azetaxz = (ma[idtxx].zetaxz + ma[idtxxi].zetaxz) / 2.0;
    Arho    = (ma[idtxx].rho + ma[idtxxi].rho) / 2.0;

    aftva.Vxx[idvx] = (2.0 * Arho - Azetaxx * dif->dt) / (2.0 * Arho + Azetaxx * dif->dt) * befva.Vxx[idvx]
        + 2.0 * dif->dt / (2.0 * Arho + Azetaxx * dif->dt) * (befsa.Txx[idtxx] - befsa.Txx[idtxxi]) / dif->dx;

    aftva.Vxy[idvx] = (2.0 * Arho - Azetaxy * dif->dt) / (2.0 * Arho + Azetaxy * dif->dt) * befva.Vxy[idvx]
        + 2.0 * dif->dt / (2.0 * Arho + Azetaxy * dif->dt) * (befta.Txy[idtxyj] - befta.Txy[idtxy]) / dif->dy;

    aftva.Vxz[idvx] = (2.0 * Arho - Azetaxz * dif->dt) / (2.0 * Arho + Azetaxz * dif->dt) * befva.Vxz[idvx]
        + 2.0 * dif->dt / (2.0 * Arho + Azetaxz * dif->dt) * (befta.Tzx[idtzxk] - befta.Tzx[idtzx]) / dif->dz;

    aftva.Vx[idvx] = aftva.Vxx[idvx] + aftva.Vxy[idvx] + aftva.Vxz[idvx];
  }
}

// 粒子速度更新関数
__global__ void VyUpdate(VelArr aftva, VelArr befva, SigArr befsa, TauArr befta, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  // int imax = ran->vr.Vy.x, jmax = ran->vr.Vy.y, kmax = ran->vr.Vy.z;
  // int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;
  double Azetayx, Azetayy, Azetayz, Arho;

  if (i < ran->vr.Vy.x - 1 && j < ran->vr.Vy.y - 1 && k < ran->vr.Vy.z - 1) {
    // 各インデックスの計算
    int idvy   = k * ran->vr.Vy.x * ran->vr.Vy.y + j * ran->vr.Vy.x + i;
    int idtyy  = k * ran->sr.Tyy.x * ran->sr.Tyy.y + j * ran->sr.Tyy.x + i;
    int idtyyj = k * ran->sr.Tyy.x * ran->sr.Tyy.y + (j - 1) * ran->sr.Tyy.x + i;
    int idtxy  = k * ran->tr.Txy.x * ran->tr.Txy.y + j * ran->tr.Txy.x + i;
    int idtxyi = k * ran->tr.Txy.x * ran->tr.Txy.y + j * ran->tr.Txy.x + (i + 1);
    int idtyz  = k * ran->tr.Tyz.x * ran->tr.Tyz.y + j * ran->tr.Tyz.x + i;
    int idtyzk = (k + 1) * ran->tr.Tyz.x * ran->tr.Tyz.y + j * ran->tr.Tyz.x + i;


    // 各種パラメータの計算
    Azetayx = (ma[idtyy].zetayx + ma[idtyyj].zetayx) / 2.0;
    Azetayy = (ma[idtyy].zetayy + ma[idtyyj].zetayy) / 2.0;
    Azetayz = (ma[idtyy].zetayz + ma[idtyyj].zetayz) / 2.0;
    Arho    = (ma[idtyy].rho + ma[idtyyj].rho) / 2.0;

    // Vyxの更新
    aftva.Vyx[idvy] = (2.0 * Arho - Azetayx * dif->dt) / (2.0 * Arho + Azetayx * dif->dt) * befva.Vyx[idvy]
        + 2.0 * dif->dt / (2.0 * Arho + Azetayx * dif->dt) * (befta.Txy[idtxyi] - befta.Txy[idtxy]) / dif->dx;

    // Vyyの更新
    aftva.Vyy[idvy] = (2.0 * Arho - Azetayy * dif->dt) / (2.0 * Arho + Azetayy * dif->dt) * befva.Vyy[idvy]
        + 2.0 * dif->dt / (2.0 * Arho + Azetayy * dif->dt) * (befsa.Tyy[idtyy] - befsa.Tyy[idtyyj]) / dif->dy;

    // Vyzの更新
    aftva.Vyz[idvy] = (2.0 * Arho - Azetayz * dif->dt) / (2.0 * Arho + Azetayz * dif->dt) * befva.Vyz[idvy]
        + 2.0 * dif->dt / (2.0 * Arho + Azetayz * dif->dt) * (befta.Tyz[idtyzk] - befta.Tyz[idtyz]) / dif->dz;

    aftva.Vy[idvy] = aftva.Vyx[idvy] + aftva.Vyy[idvy] + aftva.Vyz[idvy];
  }
}

// 粒子速度更新関数
__global__ void VzUpdate(VelArr aftva, VelArr befva, SigArr befsa, TauArr befta, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  // int imax = ran->vr.Vz.x, jmax = ran->vr.Vz.y, kmax = ran->vr.Vz.z;
  // int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;

  double Azetazx, Azetazy, Azetazz, Arho;

  if(i < ran->vr.Vz.x - 1 && j < ran->vr.Vz.y - 1 && k < ran->vr.Vz.z - 1) {
    // 1D indexing for 3D arrays
    int idvz   = k * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i;
    int idtzz  = k * ran->sr.Tzz.x * ran->sr.Tzz.y + j * ran->sr.Tzz.x + i;
    int idtzzk = (k - 1) * ran->sr.Tzz.x * ran->sr.Tzz.y + j * ran->sr.Tzz.x + i;
    int idtyz  = k * ran->tr.Tyz.x * ran->tr.Tyz.y + j * ran->tr.Tyz.x + i;
    int idtyzj = k * ran->tr.Tyz.x * ran->tr.Tyz.y + (j + 1) * ran->tr.Tyz.x + i;
    int idtzx  = k * ran->tr.Tzx.x * ran->tr.Tzx.y + j * ran->tr.Tzx.x + i;
    int idtzxi = k * ran->tr.Tzx.x * ran->tr.Tzx.y + j * ran->tr.Tzx.x + (i + 1);

    Azetazx = (ma[idtzz].zetazx + ma[idtzzk].zetazx) / 2.0;
    Azetazy = (ma[idtzz].zetazy + ma[idtzzk].zetazy) / 2.0;
    Azetazz = (ma[idtzz].zetazz + ma[idtzzk].zetazz) / 2.0;
    Arho    = (ma[idtzz].rho + ma[idtzzk].rho) / 2.0;

    aftva.Vzx[idvz] = (2.0 * Arho - Azetazx * dif->dt) / (2.0 * Arho + Azetazx * dif->dt) * befva.Vzx[idvz]
        + 2.0 * dif->dt / (2.0 * Arho + Azetazx * dif->dt) * (befta.Tzx[idtzxi] - befta.Tzx[idtzx]) / dif->dx;

    aftva.Vzy[idvz] = (2.0 * Arho - Azetazy * dif->dt) / (2.0 * Arho + Azetazy * dif->dt) * befva.Vzy[idvz]
        + 2.0 * dif->dt / (2.0 * Arho + Azetazy * dif->dt) * (befta.Tyz[idtyzj] - befta.Tyz[idtyz]) / dif->dy;

    aftva.Vzz[idvz] = (2.0 * Arho - Azetazz * dif->dt) / (2.0 * Arho + Azetazz * dif->dt) * befva.Vzz[idvz]
        + 2.0 * dif->dt / (2.0 * Arho + Azetazz * dif->dt) * (befsa.Tzz[idtzz] - befsa.Tzz[idtzzk]) / dif->dz;

    aftva.Vz[idvz] = aftva.Vzx[idvz] + aftva.Vzy[idvz] + aftva.Vzz[idvz];
  }
}


__global__ void DirectionalAddV(VelArr aftva, Range *ran, char check) {
  // スレッドインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int id;

  if (check == 'X' && i < ran->vr.Vx.x - 1 && j < ran->vr.Vx.y - 1 && k < ran->vr.Vx.z - 1) {
    id = k * ran->vr.Vx.x * ran->vr.Vx.y + j * ran->vr.Vx.x + i; 
    aftva.Vx[id] = aftva.Vxx[id] + aftva.Vxy[id] + aftva.Vxz[id];
  } else if (check == 'Y' && i < ran->vr.Vy.x - 1 && j < ran->vr.Vy.y - 1 && k < ran->vr.Vy.z - 1) {
    id = k * ran->vr.Vy.x * ran->vr.Vy.y + j * ran->vr.Vy.x + i; 
    aftva.Vy[id] = aftva.Vyx[id] + aftva.Vyy[id] + aftva.Vyz[id];
  } else if (check == 'Z' && i < ran->vr.Vz.x && j < ran->vr.Vz.y && k < ran->vr.Vz.z) {
    id = k * ran->vr.Vz.x * ran->vr.Vz.y + j * ran->vr.Vz.x + i; 
    aftva.Vz[id] = aftva.Vzx[id] + aftva.Vzy[id] + aftva.Vzz[id];
  }
  
}

// Vxクラス的な
void Vx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {
  int Vximax = ran_h->vr.Vx.x, Vxjmax = ran_h->vr.Vx.y, Vxkmax = ran_h->vr.Vx.z;
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Vximax - 2 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Vxjmax - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Vxkmax - 2 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // dim3 ZeroXYBlocks((Vximax + threadsPerBlock.x - 1) / threadsPerBlock.x, 
  //                   (Vxjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  // dim3 ZeroXZBlocks((Vximax + threadsPerBlock.x - 1) / threadsPerBlock.x, 
  //                   (Vxkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);  
  dim3 DirectionalAddBlocks((Vximax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                            (Vxjmax + threadsPerBlock.y - 1) / threadsPerBlock.y, 
                            (Vxkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
  VxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d->va, bef_d->va, bef_d->sa, bef_d->ta, ma_d, dif_d, ran_d);
  // cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vx  update: %s\n", cudaGetErrorString(err));
  // ZeroVx_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d);
  // ZeroVx_XZ<<<ZeroXZBlocks, threadsPerBlock>>>(aft_d, ran_d);
  
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vx  zero  : %s\n", cudaGetErrorString(err));
  // DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d->va, ran_d, 'X');
  // cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vx  add   : %s\n", cudaGetErrorString(err));
  // cudaDeviceSynchronize();

}
// Vyクラス的な
void Vy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {

  int Vyimax = ran_h->vr.Vy.x, Vyjmax = ran_h->vr.Vy.y, Vykmax = ran_h->vr.Vy.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Vyimax - 2 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Vyjmax - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Vykmax - 2 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  // dim3 ZeroYXBlocks((Vyimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, 
  //                   (Vyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  // dim3 ZeroYZBlocks((Vyjmax + threadsPerBlock.x - 1) / threadsPerBlock.x, 
  //                   (Vykmax + threadsPerBlock.y - 1) / threadsPerBlock.y);  
  dim3 DirectionalAddBlocks((Vyimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                            (Vyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y, 
                            (Vykmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
  VyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d->va, bef_d->va, bef_d->sa, bef_d->ta, ma_d, dif_d, ran_d);
  // cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vy  update: %s\n", cudaGetErrorString(err));
  // ZeroVy_YX<<<ZeroYXBlocks, threadsPerBlock>>>(aft_d, ran_d);
  // ZeroVy_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d);
 
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vy  zero  : %s\n", cudaGetErrorString(err));

  //全方向加算
  // DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d->va, ran_d, 'Y');
  // cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vy  add   : %s\n", cudaGetErrorString(err));
  // cudaDeviceSynchronize();
}
// Vzクラス的な
void Vz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {

  int Vzimax = ran_h->vr.Vz.x, Vzjmax = ran_h->vr.Vz.y, Vzkmax = ran_h->vr.Vz.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Vzimax - 2 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Vzjmax - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Vzkmax - 2 + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  // dim3 ZeroZXBlocks((Vzimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, 
  //                   (Vzkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  // dim3 ZeroZYBlocks((Vzjmax + threadsPerBlock.x - 1) / threadsPerBlock.x, 
  //                   (Vzkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);  
  dim3 DirectionalAddBlocks((Vzimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                            (Vzjmax + threadsPerBlock.y - 1) / threadsPerBlock.y, 
                            (Vzkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);                    
  VzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d->va, bef_d->va, bef_d->sa, bef_d->ta, ma_d, dif_d, ran_d);

  // cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vz  update: %s\n", cudaGetErrorString(err));
  // ZeroVz_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d);
  // ZeroVz_ZY<<<ZeroZYBlocks, threadsPerBlock>>>(aft_d, ran_d);
 
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vz  zero  : %s\n", cudaGetErrorString(err));
  //全方向加算
  // DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d->va, ran_d, 'Z');
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vz  add   : %s\n", cudaGetErrorString(err));
  // cudaDeviceSynchronize();
}
//粒子速度計算
void Vel(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {
  Vx(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
  Vy(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
  Vz(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
}

__global__ void Acceleration(DimD3 *acc, VelArr aftva, VelArr befva, Diff *dif, DimI3 *out, Range *ran, int outnum) {
  for(int num = 0; num < outnum; num++) {
    // printf("%d,%d,%d\n", out[num].x, out[num].y, out[num].z);
    int idvx  = out[num].z * ran->vr.Vx.x * ran->vr.Vx.y + out[num].y * ran->vr.Vx.x + out[num].x;
    int idvxi = out[num].z * ran->vr.Vx.x * ran->vr.Vx.y + out[num].y * ran->vr.Vx.x + (out[num].x + 1);
    int idvy  = out[num].z * ran->vr.Vy.x * ran->vr.Vy.y + out[num].y * ran->vr.Vy.x + out[num].x;
    int idvyj = out[num].z * ran->vr.Vy.x * ran->vr.Vy.y + (out[num].y + 1) * ran->vr.Vy.x + out[num].x;
    int idvz  = out[num].z * ran->vr.Vz.x * ran->vr.Vz.y + out[num].y * ran->vr.Vz.x + out[num].x;
    int idvzk = (out[num].z + 1) * ran->vr.Vz.x * ran->vr.Vz.y + out[num].y * ran->vr.Vz.x + out[num].x;

    acc[num].x = ((aftva.Vx[idvx] - befva.Vx[idvx]) / dif->dt + (aftva.Vx[idvxi] - befva.Vx[idvxi]) / dif->dt) / 2;
    acc[num].y = ((aftva.Vy[idvy] - befva.Vy[idvy]) / dif->dt + (aftva.Vy[idvyj] - befva.Vy[idvyj]) / dif->dt) / 2;
    acc[num].z = ((aftva.Vz[idvz] - befva.Vz[idvz]) / dif->dt + (aftva.Vz[idvzk] - befva.Vz[idvzk]) / dif->dt) / 2;
  }
}

//更新
__global__ void swapTxx(SigArr aftsa, SigArr befsa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Txximax = ran->sr.Txx.x, Txxjmax = ran->sr.Txx.y, Txxkmax = ran->sr.Txx.z;
  if (i < Txximax && j < Txxjmax && k < Txxkmax) {
    int idx_Txx = k * Txximax * Txxjmax + j * Txximax + i;
    befsa.Txx[idx_Txx]  = aftsa.Txx[idx_Txx];
    befsa.Txxx[idx_Txx] = aftsa.Txxx[idx_Txx];
    befsa.Txxy[idx_Txx] = aftsa.Txxy[idx_Txx];
    befsa.Txxz[idx_Txx] = aftsa.Txxz[idx_Txx];
  }
}

__global__ void swapTyy(SigArr aftsa, SigArr befsa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tyyimax = ran->sr.Tyy.x, Tyyjmax = ran->sr.Tyy.y, Tyykmax = ran->sr.Tyy.z;
  
  if (i < Tyyimax && j < Tyyjmax && k < Tyykmax) {
    int idx_Tyy = k * Tyyimax * Tyyjmax + j * Tyyimax + i;
    befsa.Tyy[idx_Tyy]  = aftsa.Tyy[idx_Tyy];
    befsa.Tyyx[idx_Tyy] = aftsa.Tyyx[idx_Tyy];
    befsa.Tyyy[idx_Tyy] = aftsa.Tyyy[idx_Tyy];
    befsa.Tyyz[idx_Tyy] = aftsa.Tyyz[idx_Tyy];
  }
}

__global__ void swapTzz(SigArr aftsa, SigArr befsa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tzzimax = ran->sr.Tzz.x, Tzzjmax = ran->sr.Tzz.y, Tzzkmax = ran->sr.Tzz.z;

  if (i < Tzzimax && j < Tzzjmax && k < Tzzkmax) {
    int idx_Tzz = k * Tzzimax * Tzzjmax + j * Tzzimax + i;
    befsa.Tzz [idx_Tzz] = aftsa.Tzz[idx_Tzz];
    befsa.Tzzx[idx_Tzz] = aftsa.Tzzx[idx_Tzz];
    befsa.Tzzy[idx_Tzz] = aftsa.Tzzy[idx_Tzz];
    befsa.Tzzz[idx_Tzz] = aftsa.Tzzz[idx_Tzz];
  }
}

__global__ void swapTxy(TauArr aftta, TauArr befta, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Txyimax = ran->tr.Txy.x, Txyjmax = ran->tr.Txy.y, Txykmax = ran->tr.Txy.z;

  if (i < Txyimax && j < Txyjmax && k < Txykmax) {
    int idx_Txy = k * Txyimax * Txyjmax + j * Txyimax + i;
    befta.Txy [idx_Txy] = aftta.Txy[idx_Txy];
    befta.Txyx[idx_Txy] = aftta.Txyx[idx_Txy];
    befta.Txyy[idx_Txy] = aftta.Txyy[idx_Txy];
  }
}

__global__ void swapTyz(TauArr aftta, TauArr befta, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tyzimax = ran->tr.Tyz.x, Tyzjmax = ran->tr.Tyz.y, Tyzkmax = ran->tr.Tyz.z;

  if (i < Tyzimax && j < Tyzjmax && k < Tyzkmax) {
    int idx_Tyz = k * Tyzimax * Tyzjmax + j * Tyzimax + i;
    befta.Tyz [idx_Tyz] = aftta.Tyz[idx_Tyz];
    befta.Tyzy[idx_Tyz] = aftta.Tyzy[idx_Tyz];
    befta.Tyzz[idx_Tyz] = aftta.Tyzz[idx_Tyz];
  }
}

__global__ void swapTzx(TauArr aftta, TauArr befta, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tzximax = ran->tr.Tzx.x, Tzxjmax = ran->tr.Tzx.y, Tzxkmax = ran->tr.Tzx.z;
  
  if (i < Tzximax && j < Tzxjmax && k < Tzxkmax) {
    int idx_Tzx = k * Tzximax * Tzxjmax + j * Tzximax + i;
    befta.Tzx [idx_Tzx] = aftta.Tzx[idx_Tzx];
    befta.Tzxz[idx_Tzx] = aftta.Tzxz[idx_Tzx];
    befta.Tzxx[idx_Tzx] = aftta.Tzxx[idx_Tzx];
  }
}

__global__ void swapVx(VelArr aftva, VelArr befva, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Vximax = ran->vr.Vx.x, Vxjmax = ran->vr.Vx.y, Vxkmax = ran->vr.Vx.z;
  
  if (i < Vximax && j < Vxjmax && k < Vxkmax) {
    int idx_Vx = k * Vximax * Vxjmax + j * Vximax + i;
    befva.Vx [idx_Vx] = aftva.Vx[idx_Vx];
    befva.Vxx[idx_Vx] = aftva.Vxx[idx_Vx];
    befva.Vxy[idx_Vx] = aftva.Vxy[idx_Vx];
    befva.Vxz[idx_Vx] = aftva.Vxz[idx_Vx];
  }
}

__global__ void swapVy(VelArr aftva, VelArr befva, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Vyimax = ran->vr.Vy.x, Vyjmax = ran->vr.Vy.y, Vykmax = ran->vr.Vy.z;
  
  if (i < Vyimax && j < Vyjmax && k < Vykmax) {
    int idx_Vy = k * Vyimax * Vyjmax + j * Vyimax + i;
    befva.Vy [idx_Vy] = aftva.Vy[idx_Vy];
    befva.Vyx[idx_Vy] = aftva.Vyx[idx_Vy];
    befva.Vyy[idx_Vy] = aftva.Vyy[idx_Vy];
    befva.Vyz[idx_Vy] = aftva.Vyz[idx_Vy];
  }
}

__global__ void swapVz(VelArr aftva, VelArr befva, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Vzimax = ran->vr.Vz.x, Vzjmax = ran->vr.Vz.y, Vzkmax = ran->vr.Vz.z;
  
  if (i < Vzimax && j < Vzjmax && k < Vzkmax) {
    int idx_Vz = k * Vzimax * Vzjmax + j * Vzimax + i;
    befva.Vz [idx_Vz] = aftva.Vz[idx_Vz];
    befva.Vzx[idx_Vz] = aftva.Vzx[idx_Vz];
    befva.Vzy[idx_Vz] = aftva.Vzy[idx_Vz];
    befva.Vzz[idx_Vz] = aftva.Vzz[idx_Vz];
  }
}

void swapBefAft(BefAft *aft, BefAft *bef, Range *ran_h, Range *ran_d, DimI3 threads) {
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 SwapTxxBlocks((ran_h->sr.Txx.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h->sr.Txx.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h->sr.Txx.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTyyBlocks((ran_h->sr.Tyy.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h->sr.Tyy.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h->sr.Tyy.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTzzBlocks((ran_h->sr.Tzz.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h->sr.Tzz.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h->sr.Tzz.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTxyBlocks((ran_h->tr.Txy.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h->tr.Txy.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h->tr.Txy.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTyzBlocks((ran_h->tr.Tyz.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h->tr.Tyz.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h->tr.Tyz.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTzxBlocks((ran_h->tr.Tzx.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h->tr.Tzx.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h->tr.Tzx.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3  SwapVxBlocks((ran_h->vr.Vx.x  + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h->vr.Vx.y  + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h->vr.Vx.z  + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3  SwapVyBlocks((ran_h->vr.Vy.x  + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h->vr.Vy.y  + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h->vr.Vy.z  + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3  SwapVzBlocks((ran_h->vr.Vz.x  + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h->vr.Vz.y  + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h->vr.Vz.z  + threadsPerBlock.z - 1) / threadsPerBlock.z);
  swapTxx<<<SwapTxxBlocks, threadsPerBlock>>>(aft->sa, bef->sa, ran_d);
  swapTyy<<<SwapTyyBlocks, threadsPerBlock>>>(aft->sa, bef->sa, ran_d);
  swapTzz<<<SwapTzzBlocks, threadsPerBlock>>>(aft->sa, bef->sa, ran_d);
  swapTxy<<<SwapTxyBlocks, threadsPerBlock>>>(aft->ta, bef->ta, ran_d);
  swapTyz<<<SwapTyzBlocks, threadsPerBlock>>>(aft->ta, bef->ta, ran_d);
  swapTzx<<<SwapTzxBlocks, threadsPerBlock>>>(aft->ta, bef->ta, ran_d);
  swapVx<<<SwapVxBlocks, threadsPerBlock>>>(aft->va, bef->va, ran_d);
  swapVy<<<SwapVyBlocks, threadsPerBlock>>>(aft->va, bef->va, ran_d);
  swapVz<<<SwapVzBlocks, threadsPerBlock>>>(aft->va, bef->va, ran_d);
  cudaDeviceSynchronize();
}
