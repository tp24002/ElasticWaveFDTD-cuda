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

  int txxi = ran->sr.Txx.x, txxj = ran->sr.Txx.y, txxk = ran->sr.Txx.z;
  int vxi = ran->vr.Vx.x, vxj = ran->vr.Vx.y;
  int vyi = ran->vr.Vy.x, vyj = ran->vr.Vy.y;
  int vzi = ran->vr.Vz.x, vzj = ran->vr.Vz.y;
  if(i < txxi && j < txxj && k < txxk) {
    int idtxx = k * txxi * txxj + j * txxi + i;
    int idvx  = k * vxi * vxj + j * vxi + i;
    int idvxi = k * vxi * vxj + j * vxi + (i + 1);
    int idvy  = k * vyi * vyj + j * vyi + i;
    int idvyj = k * vyi * vyj + (j + 1) * vyi + i;
    int idvz  = k * vzi * vzj + j * vzi + i;
    int idvzk = (k + 1) * vzi * vzj + j * vzi + i;

    aftsa.Txxx[idtxx] = (2.0 - ma[idtxx].zetadx * dif->dt) / (2.0 + ma[idtxx].zetadx * dif->dt) * befsa.Txxx[idtxx]
        + 2.0 * (ma[idtxx].c11 * dif->dt + ma[idtxx].xi11) / (2.0 + ma[idtxx].zetadx * dif->dt) * (aftva.Vx[idvxi] - aftva.Vx[idvx]) / dif->dx
        - 2.0 * ma[idtxx].xi11 / (2.0 + ma[idtxx].zetadx * dif->dt) * (befva.Vx[idvxi] - befva.Vx[idvx]) / dif->dx;

    aftsa.Txxy[idtxx] = (2.0 - ma[idtxx].zetady * dif->dt) / (2.0 + ma[idtxx].zetady * dif->dt) * befsa.Txxy[idtxx]
        + 2.0 * (ma[idtxx].ramda * dif->dt + ma[idtxx].khi) / (2.0 + ma[idtxx].zetady * dif->dt) * (aftva.Vy[idvyj] - aftva.Vy[idvy]) / dif->dy
        - 2.0 * ma[idtxx].khi / (2.0 + ma[idtxx].zetady * dif->dt) * (befva.Vy[idvyj] - befva.Vy[idvy]) / dif->dy;

    aftsa.Txxz[idtxx] = (2.0 - ma[idtxx].zetadz * dif->dt) / (2.0 + ma[idtxx].zetadz * dif->dt) * befsa.Txxz[idtxx]
        + 2.0 * (ma[idtxx].ramda * dif->dt + ma[idtxx].khi) / (2.0 + ma[idtxx].zetadz * dif->dt) * (aftva.Vz[idvzk] - aftva.Vz[idvz]) / dif->dz
        - 2.0 * ma[idtxx].khi / (2.0 + ma[idtxx].zetadz * dif->dt) * (befva.Vz[idvzk] - befva.Vz[idvz]) / dif->dz;

    aftsa.Txx[idtxx] = aftsa.Txxx[idtxx] + aftsa.Txxy[idtxx] + aftsa.Txxz[idtxx] + ipa[idtxx].Txx;
  }
}

// 垂直応力更新並列関数
__global__ void TyyUpdate(SigArr aftsa, VelArr aftva, SigArr befsa, VelArr befva, MedArr *ma, Diff *dif, Range *ran, ImpulseArr *ipa) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int tyyi = ran->sr.Tyy.x, tyyj = ran->sr.Tyy.y, tyyk = ran->sr.Tyy.z;
  int vxi = ran->vr.Vx.x, vxj = ran->vr.Vx.y;
  int vyi = ran->vr.Vy.x, vyj = ran->vr.Vy.y;
  int vzi = ran->vr.Vz.x, vzj = ran->vr.Vz.y;
  if(i < tyyi && j < tyyj && k < tyyk) {
    int idtyy = k * tyyi * tyyj + j * tyyi + i;
    int idvx  = k * vxi * vxj + j * vxi + i;
    int idvxi = k * vxi * vxj + j * vxi + (i + 1);
    int idvy  = k * vyi * vyj + j * vyi + i;
    int idvyj = k * vyi * vyj + (j + 1) * vyi + i;
    int idvz  = k * vzi * vzj + j * vzi + i;
    int idvzk = (k + 1) * vzi * vzj + j * vzi + i;
    

    aftsa.Tyyx[idtyy] = (2.0 - ma[idtyy].zetadx * dif->dt) / (2.0 + ma[idtyy].zetadx * dif->dt) * befsa.Tyyx[idtyy]
        + 2.0 * (ma[idtyy].ramda * dif->dt + ma[idtyy].khi) / (2.0 + ma[idtyy].zetadx * dif->dt) * (aftva.Vx[idvxi] - aftva.Vx[idvx]) / dif->dx
        - 2.0 * ma[idtyy].khi / (2.0 + ma[idtyy].zetadx * dif->dt) * (befva.Vx[idvxi] - befva.Vx[idvx]) / dif->dx;

    aftsa.Tyyy[idtyy] = (2.0 - ma[idtyy].zetady * dif->dt) / (2.0 + ma[idtyy].zetady * dif->dt) * befsa.Tyyy[idtyy]
        + 2.0 * (ma[idtyy].c11 * dif->dt + ma[idtyy].xi11) / (2.0 + ma[idtyy].zetady * dif->dt) * (aftva.Vy[idvyj] - aftva.Vy[idvy]) / dif->dy
        - 2.0 * ma[idtyy].xi11 / (2.0 + ma[idtyy].zetady * dif->dt) * (befva.Vy[idvyj] - befva.Vy[idvy]) / dif->dy;

    aftsa.Tyyz[idtyy] = (2.0 - ma[idtyy].zetadz * dif->dt) / (2.0 + ma[idtyy].zetadz * dif->dt) * befsa.Tyyz[idtyy]
        + 2.0 * (ma[idtyy].ramda * dif->dt + ma[idtyy].khi) / (2.0 + ma[idtyy].zetadz * dif->dt) * (aftva.Vz[idvzk] - aftva.Vz[idvz]) / dif->dz
        - 2.0 * ma[idtyy].khi / (2.0 + ma[idtyy].zetadz * dif->dt) * (befva.Vz[idvzk] - befva.Vz[idvz]) / dif->dz;

    aftsa.Tyy[idtyy] = aftsa.Tyyx[idtyy] + aftsa.Tyyy[idtyy] + aftsa.Tyyz[idtyy] + ipa[idtyy].Tyy;
  }
}

// 垂直応力更新並列関数
__global__ void TzzUpdate(SigArr aftsa, VelArr aftva, SigArr befsa, VelArr befva, MedArr *ma, Diff *dif, Range *ran, ImpulseArr *ipa) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int tzzi = ran->sr.Tzz.x, tzzj = ran->sr.Tzz.y, tzzk = ran->sr.Tzz.z;
  int vxi = ran->vr.Vx.x, vxj = ran->vr.Vx.y;
  int vyi = ran->vr.Vy.x, vyj = ran->vr.Vy.y;
  int vzi = ran->vr.Vz.x, vzj = ran->vr.Vz.y;

  if(i < tzzi && j < tzzj && k < tzzk) {
    int idtzz = k * tzzi * tzzj + j * tzzi + i;
    int idvx  = k * vxi * vxj + j * vxi + i;
    int idvxi = k * vxi * vxj + j * vxi + (i + 1);
    int idvy  = k * vyi * vyj + j * vyi + i;
    int idvyj = k * vyi * vyj + (j + 1) * vyi + i;
    int idvz  = k * vzi * vzj + j * vzi + i;
    int idvzk = (k + 1) * vzi * vzj + j * vzi + i;

    aftsa.Tzzx[idtzz] = (2.0 - ma[idtzz].zetadx * dif->dt) / (2.0 + ma[idtzz].zetadx * dif->dt) * befsa.Tzzx[idtzz]
        + 2.0 * (ma[idtzz].ramda * dif->dt + ma[idtzz].khi) / (2.0 + ma[idtzz].zetadx * dif->dt) * (aftva.Vx[idvxi] - aftva.Vx[idvx]) / dif->dx
        - 2.0 * ma[idtzz].khi / (2.0 + ma[idtzz].zetadx * dif->dt) * (befva.Vx[idvxi] - befva.Vx[idvx]) / dif->dx;

    aftsa.Tzzy[idtzz] = (2.0 - ma[idtzz].zetady * dif->dt) / (2.0 + ma[idtzz].zetady * dif->dt) * befsa.Tzzy[idtzz]
        + 2.0 * (ma[idtzz].ramda * dif->dt + ma[idtzz].khi) / (2.0 + ma[idtzz].zetady * dif->dt) * (aftva.Vy[idvyj] - aftva.Vy[idvy]) / dif->dy
        - 2.0 * ma[idtzz].khi / (2.0 + ma[idtzz].zetady * dif->dt) * (befva.Vy[idvyj] - befva.Vy[idvy]) / dif->dy;

    aftsa.Tzzz[idtzz] = (2.0 - ma[idtzz].zetadz * dif->dt) / (2.0 + ma[idtzz].zetadz * dif->dt) * befsa.Tzzz[idtzz]
        + 2.0 * (ma[idtzz].c11 * dif->dt + ma[idtzz].xi11) / (2.0 + ma[idtzz].zetadz * dif->dt) * (aftva.Vz[idvzk] - aftva.Vz[idvz]) / dif->dz
        - 2.0 * ma[idtzz].xi11 / (2.0 + ma[idtzz].zetadz * dif->dt) * (befva.Vz[idvzk] - befva.Vz[idvz]) / dif->dz;

    aftsa.Tzz[idtzz] = aftsa.Tzzx[idtzz] + aftsa.Tzzy[idtzz] + aftsa.Tzzz[idtzz] + ipa[idtzz].Tzz;
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
  int Txximax = ran_h->sr.Txx.x, Txxjmax = ran_h->sr.Txx.y, Txxkmax = ran_h->sr.Txx.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z);  // ブロック内のスレッド数
  dim3 UpdateBlocks((Txximax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Txxjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Txxkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
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
  dim3 UpdateBlocks((Tyyimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tyyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tyykmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
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
  dim3 UpdateBlocks((Tzzimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tzzjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tzzkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
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
  int txyi = ran->tr.Txy.x, txyj = ran->tr.Txy.y, txyk = ran->tr.Txy.z;
  int txxi = ran->sr.Txx.x, txxj = ran->sr.Txx.y;
  int vxi = ran->vr.Vx.x, vxj = ran->vr.Vx.y;
  int vyi = ran->vr.Vy.x, vyj = ran->vr.Vy.y;
  if (i < txyi && j < txyj && k < txyk) {
    // 各インデックスの計算
    int idtxy = k * ran->tr.Txy.x * ran->tr.Txy.y + j * ran->tr.Txy.x + i;

    int idtxx   = k * txxi * txxj + j * txxi + i;
    int idtxxi  = k * txxi * txxj + j * txxi + (i - 1);
    int idtxxj  = k * txxi * txxj + (j - 1) * txxi + i;
    int idtxxij = k * txxi * txxj + (j - 1) * txxi + (i - 1);

    int idvx  = k * vxi * vxj + j * vxi + i;
    int idvxj = k * vxi * vxj + (j - 1) * vxi + i;
    int idvy  = k * vyi * vyj + j * vyi + i;
    int idvyi = k * vyi * vyj + j * vyi + (i - 1);


    // PML:減衰係数,計算領域:摩擦定数
    Hzetadx = 4. * pow((1. / ma[idtxx].zetadx) + (1. / ma[idtxxi].zetadx) + (1. / ma[idtxxj].zetadx) + (1. / ma[idtxxij].zetadx), -1.);
    Hzetady = 4. * pow((1. / ma[idtxx].zetady) + (1. / ma[idtxxi].zetady) + (1. / ma[idtxxj].zetady) + (1. / ma[idtxxij].zetady), -1.);
    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / ma[idtxx].mu) + (1. / ma[idtxxi].mu) + (1. / ma[idtxxj].mu) + (1. / ma[idtxxij].mu), -1.);
    // 第1粘性定数
    Hgamma = 4. * pow((1. / ma[idtxx].gamma) + (1. / ma[idtxxi].gamma) + (1. / ma[idtxxj].gamma) + (1. / ma[idtxxij].gamma), -1.);


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

  double Hzetady, Hzetadz, Hmu, Hgamma;

  int tyzi = ran->tr.Tyz.x, tyzj = ran->tr.Tyz.y, tyzk = ran->tr.Tyz.z;
  int txxi = ran->sr.Txx.x, txxj = ran->sr.Txx.y;
  int vyi = ran->vr.Vy.x, vyj = ran->vr.Vy.y;
  int vzi = ran->vr.Vz.x, vzj = ran->vr.Vz.y;
  if (i < tyzi && j < tyzj && k < tyzk) {
    // 各インデックスの計算
    int idtyz = k * tyzi * tyzj + j * tyzi + i;

    int idtxx   = k * txxi * txxj + j * txxi + i;
    int idtxxj  = k * txxi * txxj + (j - 1) * txxi + i;
    int idtxxk  = (k - 1) * txxi * txxj + j * txxi + i;
    int idtxxjk = (k - 1) * txxi * txxj + (j - 1) * txxi + i;

    int idvy  = k * vyi * vyj + j * vyi + i;
    int idvyk = (k - 1) * vyi * vyj + j * vyi + i;
    int idvz  = k * vzi * vzj + j * vzi + i;
    int idvzj = k * vzi * vzj + (j - 1) * vzi + i;


    // PML:減衰係数,計算領域:摩擦定数
    Hzetady = 4. * pow((1. / ma[idtxx].zetady) + (1. / ma[idtxxj].zetady) + (1. / ma[idtxxk].zetady) + (1. / ma[idtxxjk].zetady), -1.);
    Hzetadz = 4. * pow((1. / ma[idtxx].zetadz) + (1. / ma[idtxxj].zetadz) + (1. / ma[idtxxk].zetadz) + (1. / ma[idtxxjk].zetadz), -1.);
    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / ma[idtxx].mu) + (1. / ma[idtxxj].mu) + (1. / ma[idtxxk].mu) + (1. / ma[idtxxjk].mu), -1.);
    // 第1粘性定数
    Hgamma = 4. * pow((1. / ma[idtxx].gamma) + (1. / ma[idtxxj].gamma) + (1. / ma[idtxxk].gamma) + (1. / ma[idtxxjk].gamma), -1.);

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
  int tzxi = ran->tr.Tzx.x, tzxj = ran->tr.Tzx.y, tzxk = ran->tr.Tzx.z;
  int txxi = ran->sr.Txx.x, txxj = ran->sr.Txx.y;
  int vxi = ran->vr.Vx.x, vxj = ran->vr.Vx.y;
  int vzi = ran->vr.Vz.x, vzj = ran->vr.Vz.y;
  if (i < tzxi && j < tzxj && k < tzxk) {
    // 各インデックスの計算
    // 求めるTzx
    int idtzx  = k * tzxi * tzxj + j * tzxi + i;
    // 格子点
    int idtxx   = k * txxi * txxj + j * txxi + i;
    int idtxxi  = k * txxi * txxj + j * txxi + (i - 1);
    int idtxxk  = (k - 1) * txxi * txxj + j * txxi + i;
    int idtxxki = (k - 1) * txxi * txxj + j * txxi + (i - 1);
    // 速度点
    int idvx  = k * vxi * vxj + j * vxi + i;
    int idvxk = (k - 1) * vxi * vxj + j * vxi + i;
    int idvz  = k * vzi * vzj + j * vzi + i;
    int idvzi = k * vzi * vzj + j * vzi + (i - 1);

    // PML:減衰係数,計算領域:摩擦定数
    Hzetadx = 4. * pow((1. / ma[idtxx].zetadx) + (1. / ma[idtxxi].zetadx) + (1. / ma[idtxxk].zetadx) + (1. / ma[idtxxki].zetadx), -1.);
    Hzetadz = 4. * pow((1. / ma[idtxx].zetadz) + (1. / ma[idtxxi].zetadz) + (1. / ma[idtxxk].zetadz) + (1. / ma[idtxxki].zetadz), -1.);
    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / ma[idtxx].mu) + (1. / ma[idtxxi].mu) + (1. / ma[idtxxk].mu) + (1. / ma[idtxxki].mu), -1.);
    // 第1粘性定数
    Hgamma = 4. * pow((1. / ma[idtxx].gamma) + (1. / ma[idtxxi].gamma) + (1. / ma[idtxxk].gamma) + (1. / ma[idtxxki].gamma), -1.);

    aftta.Tzxz[idtzx] = (2.0 - Hzetadz * dif->dt) / (2.0 + Hzetadz * dif->dt) * befta.Tzxz[idtzx]
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadz * dif->dt) * (aftva.Vx[idvx] - aftva.Vx[idvxk]) / dif->dz
        - 2.0 * Hgamma / (2.0 + Hzetadz * dif->dt) * (befva.Vx[idvx] - befva.Vx[idvxk]) / dif->dz;

    aftta.Tzxx[idtzx] = (2.0 - Hzetadx * dif->dt) / (2.0 + Hzetadx * dif->dt) * befta.Tzxx[idtzx]
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadx * dif->dt) * (aftva.Vz[idvz] - aftva.Vz[idvzi]) / dif->dx
        - 2.0 * Hgamma / (2.0 + Hzetadx * dif->dt) * (befva.Vz[idvz] - befva.Vz[idvzi]) / dif->dx;
    aftta.Tzx[idtzx] = aftta.Tzxx[idtzx] + aftta.Tzxz[idtzx];
  }
}

// Txyクラス的な
void Txy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {
  // cudaError_t err;
  int Txyimax = ran_h->tr.Txy.x, Txyjmax = ran_h->tr.Txy.y, Txykmax = ran_h->tr.Txy.z;
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Txyimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Txyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Txykmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
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
  dim3 UpdateBlocks((Tyzimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tyzjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tyzkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
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
  dim3 UpdateBlocks((Tzximax + threadsPerBlock.x - 1)     / threadsPerBlock.x,
                    (Tzxjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tzxkmax + threadsPerBlock.z - 1)     / threadsPerBlock.z);
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
  int vxi = ran->vr.Vx.x, vxj = ran->vr.Vx.y, vxk = ran->vr.Vx.z;
  int txxi = ran->sr.Txx.x, txxj = ran->sr.Txx.y;
  int txyi = ran->tr.Txy.x, txyj = ran->tr.Txy.y;
  int tzxi = ran->tr.Tzx.x, tzxj = ran->tr.Tzx.y;
  if(i < vxi && j < vxj && k < vxk) {

    int idvx   = k * vxi * vxj + j * vxi + i;
    int idtxx  = k * txxi * txxj + j * txxi + i;
    int idtxxi = k * txxi * txxj + j * txxi + (i - 1);
    int idtxy  = k * txyi * txyj + j * txyi + i;
    int idtxyj = k * txyi * txyj + (j + 1) * txyi + i;
    int idtzx  = k * tzxi * tzxj + j * tzxi + i;
    int idtzxk = (k + 1) * tzxi * tzxj + j * tzxi + i;

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

  double Azetayx, Azetayy, Azetayz, Arho;
  int vyi = ran->vr.Vy.x, vyj = ran->vr.Vy.y, vyk = ran->vr.Vy.z;
  int tyyi = ran->sr.Tyy.x, tyyj = ran->sr.Tyy.y;
  int txyi = ran->tr.Txy.x, txyj = ran->tr.Txy.y;
  int tyzi = ran->tr.Tyz.x, tyzj = ran->tr.Tyz.y;
  if (i < vyi && j < vyj && k < vyk) {
    // 各インデックスの計算
    int idvy   = k * vyi * vyj + j * vyi + i;
    int idtyy  = k * tyyi * tyyj + j * tyyi + i;
    int idtyyj = k * tyyi * tyyj + (j - 1) * tyyi + i;
    int idtxy  = k * txyi * txyj + j * txyi + i;
    int idtxyi = k * txyi * txyj + j * txyi + (i + 1);
    int idtyz  = k * tyzi * tyzj + j * tyzi + i;
    int idtyzk = (k + 1) * tyzi * tyzj + j * tyzi + i;

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
  int vzi = ran->vr.Vz.x, vzj = ran->vr.Vz.y, vzk = ran->vr.Vz.z;
  int tzzi = ran->sr.Tzz.x, tzzj = ran->sr.Tzz.y;
  int tzxi = ran->tr.Tzx.x, tzxj = ran->tr.Tzx.y;
  int tyzi = ran->tr.Tyz.x, tyzj = ran->tr.Tyz.y;
  if(i < vzi && j < vzj && k < vzk) {
    // 1D indexing for 3D arrays
    int idvz   = k * vzi * vzj + j * vzi + i;
    int idtzz  = k * tzzi * tzzj + j * tzzi + i;
    int idtzzk = (k - 1) * tzzi * tzzj + j * tzzi + i;
    int idtyz  = k * tyzi * tyzj + j * tyzi + i;
    int idtyzj = k * tyzi * tyzj + (j + 1) * tyzi + i;
    int idtzx  = k * tzxi * tzxj + j * tzxi + i;
    int idtzxi = k * tzxi * tzxj + j * tzxi + (i + 1);

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

// Vxクラス的な
void Vx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads) {
  int Vximax = ran_h->vr.Vx.x, Vxjmax = ran_h->vr.Vx.y, Vxkmax = ran_h->vr.Vx.z;
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Vximax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Vxjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Vxkmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
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
  dim3 UpdateBlocks((Vyimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Vyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Vykmax + threadsPerBlock.z - 1) / threadsPerBlock.z);
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
  dim3 UpdateBlocks((Vzimax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Vzjmax + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Vzkmax + threadsPerBlock.z - 1)     / threadsPerBlock.z);
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
