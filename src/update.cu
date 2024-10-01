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
__global__ void TxxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;

  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax - 1) {
    int idx = k * imax * jmax + j * imax + i;
    int idx_im1 = k * imax * jmax + j * imax + (i - 1);
    int idx_jm1 = k * imax * jmax + (j - 1) * imax + i;
    int idx_km1 = (k - 1) * imax * jmax + j * imax + i;

    *(aft->sa.Txxx + idx) = (2.0 - *(ma->zetadx + idx) * dif->dt) / (2.0 + *(ma->zetadx + idx) * dif->dt) * *(bef->sa.Txxx + idx)
        + 2.0 * (*(ma->c11 + idx) * dif->dt + *(ma->xi11 + idx)) / (2.0 + *(ma->zetadx + idx) * dif->dt) * (*(aft->va.Vx + idx) - *(aft->va.Vx + idx_im1)) / dif->dx
        - 2.0 * *(ma->xi11 + idx) / (2.0 + *(ma->zetadx + idx) * dif->dt) * (*(bef->va.Vx + idx) - *(bef->va.Vx + idx_im1)) / dif->dx;

    *(aft->sa.Txxy + idx) = (2.0 - *(ma->zetady + idx) * dif->dt) / (2.0 + *(ma->zetady + idx) * dif->dt) * *(bef->sa.Txxy + idx)
        + 2.0 * (*(ma->ramda + idx) * dif->dt + *(ma->khi + idx)) / (2.0 + *(ma->zetady + idx) * dif->dt) * (*(aft->va.Vy + idx) - *(aft->va.Vy + idx_jm1)) / dif->dy
        - 2.0 * *(ma->khi + idx) / (2.0 + *(ma->zetady + idx) * dif->dt) * (*(bef->va.Vy + idx) - *(bef->va.Vy + idx_jm1)) / dif->dy;

    *(aft->sa.Txxz + idx) = (2.0 - *(ma->zetadz + idx) * dif->dt) / (2.0 + *(ma->zetadz + idx) * dif->dt) * *(bef->sa.Txxz + idx)
        + 2.0 * (*(ma->ramda + idx) * dif->dt + *(ma->khi + idx)) / (2.0 + *(ma->zetadz + idx) * dif->dt) * (*(aft->va.Vz + idx) - *(aft->va.Vz + idx_km1)) / dif->dz
        - 2.0 * *(ma->khi + idx) / (2.0 + *(ma->zetadz + idx) * dif->dt) * (*(bef->va.Vz + idx) - *(bef->va.Vz + idx_km1)) / dif->dz;
  }
}
// 垂直応力更新並列関数
__global__ void TyyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = ran->sr.Tyy.x, jmax = ran->sr.Tyy.y, kmax = ran->sr.Tyy.z;

  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax - 1) {
    int idx = k * imax * jmax + j * imax + i;
    int idx_im1 = k * imax * jmax + j * imax + (i - 1);
    int idx_jm1 = k * imax * jmax + (j - 1) * imax + i;
    int idx_km1 = (k - 1) * imax * jmax + j * imax + i;

    *(aft->sa.Tyyx + idx) = (2.0 - *(ma->zetadx + idx) * dif->dt) / (2.0 + *(ma->zetadx + idx) * dif->dt) * *(bef->sa.Tyyx + idx)
        + 2.0 * (*(ma->ramda + idx) * dif->dt + *(ma->khi + idx)) / (2.0 + *(ma->zetadx + idx) * dif->dt) * (*(aft->va.Vx + idx) - *(aft->va.Vx + idx_im1)) / dif->dx
        - 2.0 * *(ma->khi + idx) / (2.0 + *(ma->zetadx + idx) * dif->dt) * (*(bef->va.Vx + idx) - *(bef->va.Vx + idx_im1)) / dif->dx;

    *(aft->sa.Tyyy + idx) = (2.0 - *(ma->zetady + idx) * dif->dt) / (2.0 + *(ma->zetady + idx) * dif->dt) * *(bef->sa.Tyyy + idx)
        + 2.0 * (*(ma->c11 + idx) * dif->dt + *(ma->xi11 + idx)) / (2.0 + *(ma->zetady + idx) * dif->dt) * (*(aft->va.Vy + idx) - *(aft->va.Vy + idx_jm1)) / dif->dy
        - 2.0 * *(ma->xi11 + idx) / (2.0 + *(ma->zetady + idx) * dif->dt) * (*(bef->va.Vy + idx) - *(bef->va.Vy + idx_jm1)) / dif->dy;

    *(aft->sa.Tyyz + idx) = (2.0 - *(ma->zetadz + idx) * dif->dt) / (2.0 + *(ma->zetadz + idx) * dif->dt) * *(bef->sa.Tyyz + idx)
        + 2.0 * (*(ma->ramda + idx) * dif->dt + *(ma->khi + idx)) / (2.0 + *(ma->zetadz + idx) * dif->dt) * (*(aft->va.Vz + idx) - *(aft->va.Vz + idx_km1)) / dif->dz
        - 2.0 * *(ma->khi + idx) / (2.0 + *(ma->zetadz + idx) * dif->dt)  * (*(bef->va.Vz + idx) - *(bef->va.Vz + idx_km1)) / dif->dz;
  }
}
// 垂直応力更新並列関数
__global__ void TzzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = ran->sr.Tzz.x, jmax = ran->sr.Tzz.y, kmax = ran->sr.Tzz.z;
  
  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax - 1) {
    int idx = k * imax * jmax + j * imax + i;
    int idx_im1 = k * imax * jmax + j * imax + (i - 1);
    int idx_jm1 = k * imax * jmax + (j - 1) * imax + i;
    int idx_km1 = (k - 1) * imax * jmax + j * imax + i;

    *(aft->sa.Tzzx + idx) = (2.0 - *(ma->zetadx + idx) * dif->dt) / (2.0 + *(ma->zetadx + idx) * dif->dt) * *(bef->sa.Tzzx + idx)
        + 2.0 * (*(ma->ramda + idx) * dif->dt + *(ma->khi + idx)) / (2.0 + *(ma->zetadx + idx) * dif->dt) * (*(aft->va.Vx + idx) - *(aft->va.Vx + idx_im1)) / dif->dx
        - 2.0 * *(ma->khi + idx) / (2.0 + *(ma->zetadx + idx) * dif->dt) * (*(bef->va.Vx + idx) - *(bef->va.Vx + idx_im1)) / dif->dx;

    *(aft->sa.Tzzy + idx) = (2.0 - *(ma->zetady + idx) * dif->dt) / (2.0 + *(ma->zetady + idx) * dif->dt) * *(bef->sa.Tzzy + idx)
        + 2.0 * (*(ma->ramda + idx) * dif->dt + *(ma->khi + idx)) / (2.0 + *(ma->zetady + idx) * dif->dt) * (*(aft->va.Vy + idx) - *(aft->va.Vy + idx_jm1)) / dif->dy
        - 2.0 * *(ma->khi + idx) / (2.0 + *(ma->zetady + idx) * dif->dt)  * (*(bef->va.Vy + idx) - *(bef->va.Vy + idx_jm1)) / dif->dy;

    *(aft->sa.Tzzz + idx) = (2.0 - *(ma->zetadz + idx) * dif->dt) / (2.0 + *(ma->zetadz + idx) * dif->dt) * *(bef->sa.Tzzz + idx)
        + 2.0 * (*(ma->c11 + idx) * dif->dt + *(ma->xi11 + idx)) / (2.0 + *(ma->zetadz + idx) * dif->dt) * (*(aft->va.Vz + idx) - *(aft->va.Vz + idx_km1)) / dif->dz
        - 2.0 * *(ma->xi11 + idx) / (2.0 + *(ma->zetadz + idx) * dif->dt) * (*(bef->va.Vz + idx) - *(bef->va.Vz + idx_km1)) / dif->dz;
  }
}
// T 0 padding
__global__ void ZeroT_XY(BefAft *aft, Range *ran, char check) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  int k = 0;

  int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;

  if (j > jmax || i > imax) {
    return;
  }

  int idx_0 = k * imax * jmax + j * imax + i;
  int idx_kmax = (kmax - 1) * imax * jmax + j * imax + i;

  if (check == 'X') {
    *(aft->sa.Txxx + idx_0) = 0.0;
    *(aft->sa.Txxx + idx_kmax) = 0.0;
    *(aft->sa.Txxy + idx_0) = 0.0;
    *(aft->sa.Txxy + idx_kmax) = 0.0;
    *(aft->sa.Txxz + idx_0) = 0.0;
    *(aft->sa.Txxz + idx_kmax) = 0.0;
  } else if (check == 'Y') {
    *(aft->sa.Tyyx + idx_0) = 0.0;
    *(aft->sa.Tyyx + idx_kmax) = 0.0;
    *(aft->sa.Tyyy + idx_0) = 0.0;
    *(aft->sa.Tyyy + idx_kmax) = 0.0;
    *(aft->sa.Tyyz + idx_0) = 0.0;
    *(aft->sa.Tyyz + idx_kmax) = 0.0;
  } else if (check == 'Z') {
    *(aft->sa.Tzzx + idx_0) = 0.0;
    *(aft->sa.Tzzx + idx_kmax) = 0.0;
    *(aft->sa.Tzzy + idx_0) = 0.0;
    *(aft->sa.Tzzy + idx_kmax) = 0.0;
    *(aft->sa.Tzzz + idx_0) = 0.0;
    *(aft->sa.Tzzz + idx_kmax) = 0.0;
  } else {
    printf("ZeroT_XY");
  }
}
// T 0 padding
__global__ void ZeroT_YZ(BefAft *aft, Range *ran, char check) {
  int k = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int i = 0;

  int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;

  if (k > kmax - 1 || j > jmax) {
    return;
  }

  int idx_0 = k * imax * jmax + j * imax + i;
  int idx_imax = k * imax * jmax + j * imax + (imax - 1);

  if (check == 'X') {
    *(aft->sa.Txxx + idx_0) = 0.0;
    *(aft->sa.Txxx + idx_imax) = 0.0;
    *(aft->sa.Txxy + idx_0) = 0.0;
    *(aft->sa.Txxy + idx_imax) = 0.0;
    *(aft->sa.Txxz + idx_0) = 0.0;
    *(aft->sa.Txxz + idx_imax) = 0.0;
  } else if (check == 'Y') {
    *(aft->sa.Tyyx + idx_0) = 0.0;
    *(aft->sa.Tyyx + idx_imax) = 0.0;
    *(aft->sa.Tyyy + idx_0) = 0.0;
    *(aft->sa.Tyyy + idx_imax) = 0.0;
    *(aft->sa.Tyyz + idx_0) = 0.0;
    *(aft->sa.Tyyz + idx_imax) = 0.0;
  } else if (check == 'Z') {
    *(aft->sa.Tzzx + idx_0) = 0.0;
    *(aft->sa.Tzzx + idx_imax) = 0.0;
    *(aft->sa.Tzzy + idx_0) = 0.0;
    *(aft->sa.Tzzy + idx_imax) = 0.0;
    *(aft->sa.Tzzz + idx_0) = 0.0;
    *(aft->sa.Tzzz + idx_imax) = 0.0;
  } else {
    printf("ZeroT_YZ");
  }
}
// T 0 padding
__global__ void ZeroT_ZX(BefAft *aft, Range *ran, char check) {
  int k = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int j = 0;
  int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;

  if (k > kmax - 1 || i > imax - 1) {
    return;
  }

  int idx_0 = k * imax * jmax + j * imax + i;
  int idx_jmax = k * imax * jmax + (jmax - 1) * imax + i;

  if (check == 'X') {
    *(aft->sa.Txxx + idx_0) = 0.0;
    *(aft->sa.Txxx + idx_jmax) = 0.0;
    *(aft->sa.Txxy + idx_0) = 0.0;
    *(aft->sa.Txxy + idx_jmax) = 0.0;
    *(aft->sa.Txxz + idx_0) = 0.0;
    *(aft->sa.Txxz + idx_jmax) = 0.0;
  } else if (check == 'Y') {
    *(aft->sa.Tyyx + idx_0) = 0.0;
    *(aft->sa.Tyyx + idx_jmax) = 0.0;
    *(aft->sa.Tyyy + idx_0) = 0.0;
    *(aft->sa.Tyyy + idx_jmax) = 0.0;
    *(aft->sa.Tyyz + idx_0) = 0.0;
    *(aft->sa.Tyyz + idx_jmax) = 0.0;
  } else if (check == 'Z') {
    *(aft->sa.Tzzx + idx_0) = 0.0;
    *(aft->sa.Tzzx + idx_jmax) = 0.0;
    *(aft->sa.Tzzy + idx_0) = 0.0;
    *(aft->sa.Tzzy + idx_jmax) = 0.0;
    *(aft->sa.Tzzz + idx_0) = 0.0;
    *(aft->sa.Tzzz + idx_jmax) = 0.0;
  } else {
    printf("ZeroT_ZX");
  }
}
// 全方向加算
__global__ void DirectionalAdd(BefAft *aft, Impulse *ip, Range *ran, char check) {
  // 1Dインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->sr.Txx.x, jmax = ran->sr.Txx.y, kmax = ran->sr.Txx.z;

  // 1Dインデックス化
  int idx = k * imax * jmax + j * imax + i;
  // printf("%lf\n",*ip->Txx);
  if (i < imax && j < jmax && k < kmax) {
    // 各方向に応じた計算を実行（ポインタ表記）
    if (check == 'X') {
      *(aft->sa.Txx + idx) = *(aft->sa.Txxx + idx) + *(aft->sa.Txxy + idx) + *(aft->sa.Txxz + idx) + *(ip->Txx + idx);
    } else if (check == 'Y') {
      *(aft->sa.Tyy + idx) = *(aft->sa.Tyyx + idx) + *(aft->sa.Tyyy + idx) + *(aft->sa.Tyyz + idx) + *(ip->Tyy + idx);
    } else if (check == 'Z') {
      *(aft->sa.Tzz + idx) = *(aft->sa.Tzzx + idx) + *(aft->sa.Tzzy + idx) + *(aft->sa.Tzzz + idx) + *(ip->Tzz + idx);
    } else {
      printf("error: DirectionalAdd\n");
    }
  }
}


// Txxクラス的な(Blocks大丈夫かな？)
void Txx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Impulse *ip_d, Coord threads) {
  // cudaError_t err;
  // char check = 'X';

  int Txximax = ran_h->sr.Txx.x, Txxjmax = ran_h->sr.Txx.y, Txxkmax = ran_h->sr.Txx.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z);  // ブロック内のスレッド数
  dim3 UpdateBlocks((Txximax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Txxjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Txxkmax - 1 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroXYBlocks((Txxjmax     + threadsPerBlock.x - 1) / threadsPerBlock.x, (Txximax     + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroYZBlocks((Txxkmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Txxjmax     + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 ZeroZXBlocks((Txxkmax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x, (Txximax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Txximax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Txxjmax + threadsPerBlock.y) / threadsPerBlock.y, 
                            (Txxkmax + threadsPerBlock.z) / threadsPerBlock.z);
  // Txx更新式
  TxxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran_d);
  
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txx update: %s\n", cudaGetErrorString(err));
  // 0 padding
  ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d, 'X');
  ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d, 'X');
  ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d, 'X');
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txx zero  : %s\n", cudaGetErrorString(err));
  //全方向加算
  DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ip_d, ran_d, 'X');
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txx add   : %s\n", cudaGetErrorString(err));

}
// Tyyクラス的な(Blocks大丈夫かな？)
void Tyy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Impulse *ip_d, Coord threads) {
  // cudaError_t err;
  char check = 'Y';

  int Tyyimax = ran_h->sr.Tyy.x, Tyyjmax = ran_h->sr.Tyy.y, Tyykmax = ran_h->sr.Tyy.z;

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
  TyyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran_d);
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyy update: %s\n", cudaGetErrorString(err));
  // 0 padding
  ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyy zero  : %s\n", cudaGetErrorString(err));
  // 全方向加算
  DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ip_d, ran_d, check);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyy add   : %s\n", cudaGetErrorString(err));
}
// Tzzクラス的な(Blocks大丈夫かな？)
void Tzz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Impulse *ip_d, Coord threads) {
  // cudaError_t err;
  char check = 'Z';

  int Tzzimax = ran_h->sr.Tzz.x, Tzzjmax = ran_h->sr.Tzz.y, Tzzkmax = ran_h->sr.Tzz.z;

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
  TzzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran_d);
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzz update: %s\n", cudaGetErrorString(err));
  // 0 padding
  ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d, check);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzz zero  : %s\n", cudaGetErrorString(err));
  // 全方向加算
  DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ip_d, ran_d, check);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzz add   : %s\n", cudaGetErrorString(err));
 
}
// 垂直応力計算(main呼び出し関数)
void Sig(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Impulse *ip_d, Coord threads) {
  Txx(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, ip_d, threads);
  Tyy(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, ip_d, threads);
  Tzz(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, ip_d, threads);
}

// せん断応力

// せん断応力更新関数
__global__ void TxyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1; // 始点を+1

  int imax = ran->tr.Txy.x, jmax = ran->tr.Txy.y, kmax = ran->tr.Txy.z;
  double Hzetadx, Hzetady, Hmu, Hgamma;

  if (i <= imax && j <= jmax && k <= kmax - 1) {
    // 各インデックスの計算
    int idx      = k * imax * jmax + j * imax + i;
    int idx_ip1  = k * imax * jmax + j * imax + (i + 1);
    int idx_jp1  = k * imax * jmax + (j + 1) * imax + i;
    int idx_ijp1 = k * imax * jmax + (j + 1) * imax + (i + 1);

    // PML:減衰係数,計算領域:摩擦定数
    Hzetadx = 4. * pow((1. / *(ma->zetadx + idx_ip1)) + (1. / *(ma->zetadx + idx_jp1)) + (1. / *(ma->zetadx + idx_ijp1)) + (1. / *(ma->zetadx + idx)), -1.);
    Hzetady = 4. * pow((1. / *(ma->zetady + idx_ip1)) + (1. / *(ma->zetady + idx_jp1)) + (1. / *(ma->zetady + idx_ijp1)) + (1. / *(ma->zetady + idx)), -1.);

    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / *(ma->mu + idx_ip1)) + (1. / *(ma->mu + idx_jp1)) + (1. / *(ma->mu + idx_ijp1)) + (1. / *(ma->mu + idx)), -1.);

    // 第1粘性定数
    Hgamma = 4. * pow((1. / *(ma->gamma + idx_ip1)) + (1. / *(ma->gamma + idx_jp1)) + (1. / *(ma->gamma + idx_ijp1)) + (1. / *(ma->gamma + idx)), -1.);

    *(aft->ta.Txyx + idx) = (2.0 - Hzetadx * dif->dt) / (2.0 + Hzetadx * dif->dt) * *(bef->ta.Txyx + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadx * dif->dt) * (*(aft->va.Vy + idx_ip1) - *(aft->va.Vy + idx)) / dif->dx
        - 2.0 * Hgamma / (2.0 + Hzetadx * dif->dt) * (*(bef->va.Vy + idx_ip1) - *(bef->va.Vy + idx)) / dif->dx;

    *(aft->ta.Txyy + idx) = (2.0 - Hzetady * dif->dt) / (2.0 + Hzetady * dif->dt) * *(bef->ta.Txyy + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetady * dif->dt) * (*(aft->va.Vx + idx_jp1) - *(aft->va.Vx + idx)) / dif->dy
        - 2.0 * Hgamma / (2.0 + Hzetady * dif->dt) * (*(bef->va.Vx + idx_jp1) - *(bef->va.Vx + idx)) / dif->dy;
  }
}
// せん断応力更新関数
__global__ void TyzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->tr.Tyz.x, jmax = ran->tr.Tyz.y, kmax = ran->tr.Tyz.z;
  double Hzetady, Hzetadz, Hmu, Hgamma;

  if (i <= imax - 1 && j <= jmax && k <= kmax) {
    // 各インデックスの計算
    int idx      = k * imax * jmax + j * imax + i;
    int idx_jp1  = k * imax * jmax + (j + 1) * imax + i;
    int idx_kp1  = (k + 1) * imax * jmax + j * imax + i;
    int idx_jkp1 = (k + 1) * jmax * kmax + (j + 1) * kmax + i;

    // PML:減衰係数,計算領域:摩擦定数
    Hzetady = 4. * pow((1. / *(ma->zetady + idx_jkp1)) + (1. / *(ma->zetady + idx_jp1)) + (1. / *(ma->zetady + idx_kp1)) + (1. / *(ma->zetady + idx)), -1.);
    Hzetadz = 4. * pow((1. / *(ma->zetadz + idx_jkp1)) + (1. / *(ma->zetadz + idx_jp1)) + (1. / *(ma->zetadz + idx_kp1)) + (1. / *(ma->zetadz + idx)), -1.);
    
    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / *(ma->mu + idx_jkp1)) + (1. / *(ma->mu + idx_jp1)) + (1. / *(ma->mu + idx_kp1)) + (1. / *(ma->mu + idx)), -1.);

    // 第1粘性定数
    Hgamma = 4. * pow((1. / *(ma->gamma + idx_jkp1)) + (1. / *(ma->gamma + idx_jp1)) + (1. / *(ma->gamma + idx_kp1)) + (1. / *(ma->gamma + idx)), -1.);

    *(aft->ta.Tyzy + idx) = (2.0 - Hzetady * dif->dt) / (2.0 + Hzetady * dif->dt) * *(bef->ta.Tyzy + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetady * dif->dt) * (*(aft->va.Vz + idx_jp1) - *(aft->va.Vz + idx)) / dif->dy
        - 2.0 * Hgamma / (2.0 + Hzetady * dif->dt) * (*(bef->va.Vz + idx_jp1) - *(bef->va.Vz + idx)) / dif->dy;

    *(aft->ta.Tyzz + idx) = (2.0 - Hzetadz * dif->dt) / (2.0 + Hzetadz * dif->dt) * *(bef->ta.Tyzz + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadz * dif->dt) * (*(aft->va.Vy + idx_kp1) - *(aft->va.Vy + idx)) / dif->dz
        - 2.0 * Hgamma / (2.0 + Hzetadz * dif->dt) * (*(bef->va.Vy + idx_kp1) - *(bef->va.Vy + idx)) / dif->dz;
  }
}
// せん断応力更新関数
__global__ void TzxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->tr.Tzx.x, jmax = ran->tr.Tzx.y, kmax = ran->tr.Tzx.z;
  double Hzetadx, Hzetadz, Hmu, Hgamma;

  if (i <= imax && j <= jmax - 1 && k <= kmax) {
    // 各インデックスの計算
    int idx      = k * imax * jmax + j * imax + i;
    int idx_ip1  = k * imax * jmax + j * imax + (i + 1);
    int idx_kp1  = (k + 1) * imax * jmax + j * imax + i;
    int idx_kip1 = (k + 1) * imax * imax + j * imax + (i + 1);

    // PML:減衰係数,計算領域:摩擦定数
    Hzetadx = 4. * pow((1. / *(ma->zetadx + idx_kip1)) + (1. / *(ma->zetadx + idx_ip1)) + (1. / *(ma->zetadx + idx_kp1)) + (1. / *(ma->zetadx + idx)), -1.);
    Hzetadz = 4. * pow((1. / *(ma->zetadz + idx_kip1)) + (1. / *(ma->zetadz + idx_ip1)) + (1. / *(ma->zetadz + idx_kp1)) + (1. / *(ma->zetadz + idx)), -1.);

    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / *(ma->mu + idx_kip1)) + (1. / *(ma->mu + idx_ip1)) + (1. / *(ma->mu + idx_kp1)) + (1. / *(ma->mu + idx)), -1.);

    // 第1粘性定数
    Hgamma = 4. * pow((1. / *(ma->gamma + idx_kip1)) + (1. / *(ma->gamma + idx_ip1)) + (1. / *(ma->gamma + idx_kp1)) + (1. / *(ma->gamma + idx)), -1.);

    *(aft->ta.Tzxz + idx) = (2.0 - Hzetadz * dif->dt) / (2.0 + Hzetadz * dif->dt) * *(bef->ta.Tzxz + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadz * dif->dt) * (*(aft->va.Vx + idx_kp1) - *(aft->va.Vx + idx)) / dif->dz
        - 2.0 * Hgamma / (2.0 + Hzetadz * dif->dt) * (*(bef->va.Vx + idx_kp1) - *(bef->va.Vx + idx)) / dif->dz;

    *(aft->ta.Tzxx + idx) = (2.0 - Hzetadx * dif->dt) / (2.0 + Hzetadx * dif->dt) * *(bef->ta.Tzxx + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadx * dif->dt) * (*(aft->va.Vz + idx_ip1) - *(aft->va.Vz + idx)) / dif->dx
        - 2.0 * Hgamma / (2.0 + Hzetadx * dif->dt) * (*(bef->va.Vz + idx_ip1) - *(bef->va.Vz + idx)) / dif->dx;
  }
}

__global__ void ZeroTxy(BefAft *aft, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = 0;

  int imax = ran->tr.Txy.x, jmax = ran->tr.Txy.y, kmax = ran->tr.Txy.z;

  if (i <= imax && j <= jmax) {
    int idx_0 = k * imax * jmax + j * imax + i;
    int idx_kmax = (kmax - 1) * imax * jmax + j * imax + i;
    *(aft->ta.Txyx + idx_0) = 0.0;
    *(aft->ta.Txyx + idx_kmax) = 0.0;
    *(aft->ta.Txyy + idx_0) = 0.0;
    *(aft->ta.Txyy + idx_kmax) = 0.0;
  }
}

__global__ void ZeroTyz(BefAft *aft, Range *ran) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  int i = 0;

  int imax = ran->tr.Tyz.x;
  int jmax = ran->tr.Tyz.y;
  int kmax = ran->tr.Tyz.z;

  if (j <= jmax && k <= kmax) {
    int idx_0 = k * imax * jmax + j * imax + i;
    int idx_imax = k * imax * jmax + j * imax + (imax - 1);

    *(aft->ta.Tyzy + idx_0) = 0.0;
    *(aft->ta.Tyzy + idx_imax) = 0.0;
    *(aft->ta.Tyzz + idx_0) = 0.0;
    *(aft->ta.Tyzz + idx_imax) = 0.0;
  }
}

__global__ void ZeroTzx(BefAft *aft, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  int j = 0;

  int imax = ran->tr.Tyz.x, jmax = ran->tr.Tyz.y, kmax = ran->tr.Tyz.z;

  if (i <= imax && k <= kmax) {
    int idx_0 = k * imax * jmax + j * imax + i;
    int idx_jmax = k * imax * jmax + (jmax - 1) * imax + i;
    *(aft->ta.Tzxx + idx_0) = 0.0;
    *(aft->ta.Tzxx + idx_jmax) = 0.0;
    *(aft->ta.Tzxz + idx_0) = 0.0;
    *(aft->ta.Tzxz + idx_jmax) = 0.0;
  }
}

__global__ void DirectionalAddT(BefAft *aft, Range *ran, char check) {
  // スレッドインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->tr.Txy.x, jmax = ran->tr.Txy.y, kmax = ran->tr.Txy.z;

  if (i > imax || j > jmax || k > kmax) {
    return;
  }

  int idx = k * imax * jmax + j * imax + i;

  if (check == 'X') {
    *(aft->ta.Tyz + idx) = *(aft->ta.Tyzy + idx) + *(aft->ta.Tyzz + idx);
  } else if (check == 'Y') {
    *(aft->ta.Tzx + idx) = *(aft->ta.Tzxx + idx) + *(aft->ta.Tzxz + idx);
  } else if (check == 'Z') {
    *(aft->ta.Txy + idx) = *(aft->ta.Txyx + idx) + *(aft->ta.Txyy + idx);
  } else {
    printf("error: DirectionalAddT");
  }
}

// Txyクラス的な
void Txy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads) {
  // cudaError_t err;
  int Txyimax = ran_h->tr.Txy.x, Txyjmax = ran_h->tr.Txy.y, Txykmax = ran_h->tr.Txy.z;
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Txyimax + threadsPerBlock.x - 1)     / threadsPerBlock.x,
                    (Txyjmax + threadsPerBlock.y - 1)     / threadsPerBlock.y,
                    (Txykmax - 1 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroXYBlocks((Txyimax + threadsPerBlock.x - 1) / threadsPerBlock.x,(Txyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Txyimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Txyjmax + threadsPerBlock.y) / threadsPerBlock.y,
                            (Txykmax + threadsPerBlock.z) / threadsPerBlock.z);                    
  TxyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran_d);
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txy update: %s\n", cudaGetErrorString(err));
  ZeroTxy<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d);
  cudaDeviceSynchronize();

  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txy zero  : %s\n", cudaGetErrorString(err));
  DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran_d, 'Z');
  cudaDeviceSynchronize();

  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Txy add   : %s\n", cudaGetErrorString(err));
}

// Tyzクラス的な
void Tyz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads) {
  // cudaError_t err;
  int Tyzimax = ran_h->tr.Tyz.x, Tyzjmax = ran_h->tr.Tyz.y, Tyzkmax = ran_h->tr.Tyz.z;
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Tyzimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tyzjmax + threadsPerBlock.y - 1)     / threadsPerBlock.y,
                    (Tyzkmax + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  dim3 ZeroYZBlocks((Tyzjmax + threadsPerBlock.x - 1) / threadsPerBlock.x,(Tyzkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Tyzimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Tyzjmax + threadsPerBlock.y) / threadsPerBlock.y,
                            (Tyzkmax + threadsPerBlock.z) / threadsPerBlock.z);                    
  TyzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran_d);
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyz update: %s\n", cudaGetErrorString(err));
  ZeroTyz<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyz zero  : %s\n", cudaGetErrorString(err));
  DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran_d, 'X');
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tyz add   : %s\n", cudaGetErrorString(err));
}
// Tzxクラス的な
void Tzx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads) {
  // cudaError_t err;
  int Tzximax = ran_h->tr.Tzx.x, Tzxjmax = ran_h->tr.Tzx.y, Tzxkmax = ran_h->tr.Tzx.z;
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Tzximax + threadsPerBlock.x - 1)     / threadsPerBlock.x,
                    (Tzxjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tzxkmax + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  dim3 ZeroZXBlocks((Tzximax + threadsPerBlock.x - 1) / threadsPerBlock.x,(Tzxkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);   
  dim3 DirectionalAddBlocks((Tzximax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Tzxjmax + threadsPerBlock.y) / threadsPerBlock.y, 
                            (Tzxkmax + threadsPerBlock.z) / threadsPerBlock.z);                  
  TzxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran_d);
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzx update: %s\n", cudaGetErrorString(err));
  ZeroTzx<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d);
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzx zero  : %s\n", cudaGetErrorString(err));
  DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran_d, 'Y');
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Tzx add   : %s\n", cudaGetErrorString(err));

}
// せん断応力計算(main呼び出し関数)
void Tau(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads) {
  Txy(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
  Tyz(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
  Tzx(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
}

// 粒子速度

// 粒子速度更新関数
__global__ void VxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int imax = ran->vr.Vx.x, jmax = ran->vr.Vx.y, kmax = ran->vr.Vx.z;
  double Azetaxx, Azetaxy, Azetaxz, Arho;
  // printf("%d,%d,%d\n", imax ,jmax ,kmax);


  if(i <= imax && j <= jmax - 1 && k <= kmax - 1) {
    // 1D indexing for 3D arrays
    // printf("ok\n");
    int idx    = k * imax * jmax + j * imax + i;
    int idx_i1 = k * imax * jmax + j * imax + (i + 1);
    int idx_j1 = k * imax * jmax + (j - 1) * imax + i;
    int idx_k1 = (k - 1) * imax * jmax + j * imax + i;

    Azetaxx = (*(ma->zetaxx + idx_i1) + *(ma->zetaxx + idx)) / 2.;
    Azetaxy = (*(ma->zetaxy + idx_i1) + *(ma->zetaxy + idx)) / 2.;
    Azetaxz = (*(ma->zetaxz + idx_i1) + *(ma->zetaxz + idx)) / 2.;
    Arho    = (*(ma->rho + idx_i1) + *(ma->rho + idx)) / 2.;
    *(aft->va.Vxx + idx) = (2. * Arho - Azetaxx * dif->dt) / (2. * Arho + Azetaxx * dif->dt) * *(bef->va.Vxx + idx)
        + 2. * dif->dt / (2. * Arho + Azetaxx * dif->dt) * (*(bef->sa.Txx + idx_i1) - *(bef->sa.Txx + idx)) / dif->dx;

    *(aft->va.Vxy + idx) = (2. * Arho - Azetaxy * dif->dt) / (2. * Arho + Azetaxy * dif->dt) * *(bef->va.Vxy + idx)
        + 2. * dif->dt / (2. * Arho + Azetaxy * dif->dt) * (*(bef->ta.Txy + idx) - *(bef->ta.Txy + idx_j1)) / dif->dy;

    *(aft->va.Vxz + idx) = (2. * Arho - Azetaxz * dif->dt) / (2. * Arho + Azetaxz * dif->dt) * *(bef->va.Vxz + idx)
        + 2. * dif->dt / (2. * Arho + Azetaxz * dif->dt) * (*(bef->ta.Tzx + idx) - *(bef->ta.Tzx + idx_k1)) / dif->dz;
  }
}
// 粒子速度更新関数
__global__ void VyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = ran->vr.Vy.x, jmax = ran->vr.Vy.y, kmax = ran->vr.Vy.z;
  double Azetayx, Azetayy, Azetayz, Arho;

  if (i <= imax - 1 && j <= jmax && k <= kmax - 1) {
    // 各インデックスの計算
    int idx     = i * imax * jmax + j * imax + k;
    int idx_im1 = k * imax * jmax + j * imax + (i - 1);
    int idx_jp1 = k * imax * jmax + (j + 1) * imax + i;
    int idx_km1 = (k - 1) * imax * jmax + j * imax + i;

    // 各種パラメータの計算
    Azetayx = (*(ma->zetayx + idx_jp1) + *(ma->zetayx + idx)) / 2.0;
    Azetayy = (*(ma->zetayy + idx_jp1) + *(ma->zetayy + idx)) / 2.0;
    Azetayz = (*(ma->zetayz + idx_jp1) + *(ma->zetayz + idx)) / 2.0;
    Arho    = (*(ma->rho + idx_jp1) + *(ma->rho + idx)) / 2.0;

    // Vyxの更新
    *(aft->va.Vyx + idx) = (2.0 * Arho - Azetayx * dif->dt) / (2.0 * Arho + Azetayx * dif->dt) * *(bef->va.Vyx + idx)
        + 2.0 * dif->dt / (2.0 * Arho + Azetayx * dif->dt) * (*(bef->ta.Txy + idx) - *(bef->ta.Txy + idx_im1)) / dif->dx;

    // Vyyの更新
    *(aft->va.Vyy + idx) = (2.0 * Arho - Azetayy * dif->dt) / (2.0 * Arho + Azetayy * dif->dt) * *(bef->va.Vyy + idx)
        + 2.0 * dif->dt / (2.0 * Arho + Azetayy * dif->dt) * (*(bef->sa.Tyy + idx_jp1) - *(bef->sa.Tyy + idx)) / dif->dy;

    // Vyzの更新
    *(aft->va.Vyz + idx) = (2.0 * Arho - Azetayz * dif->dt) / (2.0 * Arho + Azetayz * dif->dt) * *(bef->va.Vyz + idx)
        + 2.0 * dif->dt / (2.0 * Arho + Azetayz * dif->dt) * (*(bef->ta.Tyz + idx) - *(bef->ta.Tyz + idx_km1)) / dif->dz;
  }
}
// 粒子速度更新関数
__global__ void VzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->vr.Vz.x, jmax = ran->vr.Vz.y, kmax = ran->vr.Vz.z;
  double Azetazx, Azetazy, Azetazz, Arho;



  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax) {
    // 1D indexing for 3D arrays
    int idx     = k * imax * jmax + j * imax + i;
    int idx_im1 = k * imax * jmax + j * imax + (i - 1);
    int idx_jm1 = k * imax * jmax + (j - 1) * imax + i;
    int idx_kp1 = (k + 1) * imax * jmax + j * imax + i;

    Azetazx = (ma->zetazx[idx_kp1] + ma->zetazx[idx]) / 2.;
    Azetazy = (ma->zetazy[idx_kp1] + ma->zetazy[idx]) / 2.;
    Azetazz = (ma->zetazz[idx_kp1] + ma->zetazz[idx]) / 2.;
    Arho    = (ma->rho[idx_kp1] + ma->rho[idx]) / 2.;
    aft->va.Vzx[idx] = (2. * Arho - Azetazx * dif->dt) / (2. * Arho + Azetazx * dif->dt) * bef->va.Vzx[idx]
        + 2. * dif->dt / (2. * Arho + Azetazx * dif->dt) * (bef->ta.Tzx[idx] - bef->ta.Tzx[idx_im1]) / dif->dx;

    aft->va.Vzy[idx] = (2. * Arho - Azetazy * dif->dt) / (2. * Arho + Azetazy * dif->dt) * bef->va.Vzy[idx]
        + 2. * dif->dt / (2. * Arho + Azetazy * dif->dt) * (bef->ta.Tyz[idx] - bef->ta.Tyz[idx_jm1]) / dif->dy;

    aft->va.Vzz[idx] = (2. * Arho - Azetazz * dif->dt) / (2. * Arho + Azetazz * dif->dt) * bef->va.Vzz[idx]
        + 2. * dif->dt / (2. * Arho + Azetazz * dif->dt) * (bef->sa.Tzz[idx_kp1] - bef->sa.Tzz[idx]) / dif->dz;
  }
}

__global__ void ZeroVx_XY(BefAft *aft, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = 0;
  int imax = ran->vr.Vx.x, jmax = ran->vr.Vx.y, kmax = ran->vr.Vx.z;


  if(i <= imax && j <= jmax) {
    // 1D indexing for 3D arrays
    int idx_0   = k * imax * jmax + j * imax + i;
    int idx_kmax = (kmax - 1) * imax * jmax + j * imax + i;

    aft->va.Vxx[idx_0] = 0.0;
    aft->va.Vxx[idx_kmax] = 0.0;

    aft->va.Vxy[idx_0] = 0.0;
    aft->va.Vxy[idx_kmax] = 0.0;

    aft->va.Vxz[idx_0] = 0.0;
    aft->va.Vxz[idx_kmax] = 0.0;
  }
}

__global__ void ZeroVx_XZ(BefAft *aft, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  int j = 0;
  int imax = ran->vr.Vx.x, jmax = ran->vr.Vx.y, kmax = ran->vr.Vx.z;

  if(i <= imax && k <= kmax) {
    // 1D indexing for 3D arrays
    int idx_0   = k * imax * jmax + j * imax + i;
    int idx_jmax = i * jmax * kmax + jmax * kmax + k;

    aft->va.Vxx[idx_0] = 0.0;
    aft->va.Vxx[idx_jmax] = 0.0;

    aft->va.Vxy[idx_0] = 0.0;
    aft->va.Vxy[idx_jmax] = 0.0;

    aft->va.Vxz[idx_0] = 0.0;
    aft->va.Vxz[idx_jmax] = 0.0;
  }
}

__global__ void ZeroVy_YX(BefAft *aft, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = 0;
  int imax = ran->vr.Vy.x, jmax = ran->vr.Vy.y, kmax = ran->vr.Vy.z;

  if (i <= imax && j <= jmax) {
    int idx_0 = k * imax * jmax + j * imax + i;
    int idx_zmax = (kmax - 1) * imax * jmax + j * imax + i;

    *(aft->va.Vyx + idx_0) = 0.0;
    *(aft->va.Vyx + idx_zmax) = 0.0;
    *(aft->va.Vyy + idx_0) = 0.0;
    *(aft->va.Vyy + idx_zmax) = 0.0;
    *(aft->va.Vyz + idx_0) = 0.0;
    *(aft->va.Vyz + idx_zmax) = 0.0;
  }
}

__global__ void ZeroVy_YZ(BefAft *aft, Range *ran) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  int i = 0;
  int imax = ran->vr.Vy.x, jmax = ran->vr.Vy.y, kmax = ran->vr.Vy.z;

  if(j <= jmax && k <= kmax) {
    // 1D indexing for 3D arrays
    int idx_0   = k * imax * jmax + j * imax + i;
    int idx_imax = k * imax * jmax + j * imax + (imax - 1);

    aft->va.Vyx[idx_0] = 0.0;
    aft->va.Vyx[idx_imax] = 0.0;

    aft->va.Vyy[idx_0] = 0.0;
    aft->va.Vyy[idx_imax] = 0.0;

    aft->va.Vyz[idx_0] = 0.0;
    aft->va.Vyz[idx_imax] = 0.0;
  }
}

__global__ void ZeroVz_ZX(BefAft *aft, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  int j = 0;
  int imax = ran->vr.Vz.x, jmax = ran->vr.Vz.y, kmax = ran->vr.Vz.z;

  if (i <= imax && k <= kmax) {
    int idx_0 = k * imax * jmax + j * imax + i;
    int idx_ymax = k * imax * jmax + (jmax - 1) * kmax + i;

    *(aft->va.Vzx + idx_0) = 0.0;
    *(aft->va.Vzx + idx_ymax) = 0.0;
    *(aft->va.Vzy + idx_0) = 0.0;
    *(aft->va.Vzy + idx_ymax) = 0.0;
    *(aft->va.Vzz + idx_0) = 0.0;
    *(aft->va.Vzz + idx_ymax) = 0.0;
  }
}

__global__ void ZeroVz_ZY(BefAft *aft, Range *ran) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  int i = 0;
  int imax = ran->vr.Vz.x, jmax = ran->vr.Vz.y, kmax = ran->vr.Vz.z;


  if(j <= jmax && k <= kmax) {
    // 1D indexing for 3D arrays
    int idx_0   = k * imax * jmax + j * imax + i;
    int idx_imax = k * imax * jmax + j * imax + (imax - 1);

    aft->va.Vzx[idx_0] = 0.0;
    aft->va.Vzx[idx_imax] = 0.0;

    aft->va.Vzy[idx_0] = 0.0;
    aft->va.Vzy[idx_imax] = 0.0;

    aft->va.Vzz[idx_0] = 0.0;
    aft->va.Vzz[idx_imax] = 0.0;
  }
}

__global__ void DirectionalAddV(BefAft *aft, Range *ran, char check) {
  // スレッドインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ran->vr.Vx.x, jmax = ran->vr.Vx.y, kmax = ran->vr.Vx.z;


  if (i > imax || j > jmax || k > kmax) {
    return;
  }

  int idx = k * imax * jmax + j * imax + i;

  if (check == 'X') {
    *(aft->va.Vx + idx) = *(aft->va.Vxx + idx) + *(aft->va.Vxy + idx) + *(aft->va.Vxz + idx);
  } else if (check == 'Y') {
    *(aft->va.Vy + idx) = *(aft->va.Vyx + idx) + *(aft->va.Vyy + idx) + *(aft->va.Vyz + idx);
  } else if (check == 'Z') {
    *(aft->va.Vz + idx) = *(aft->va.Vzx + idx) + *(aft->va.Vzy + idx) + *(aft->va.Vzz + idx);
  } else {
    printf("error: DirectionalAddV");
  }
}

// Vxクラス的な
void Vx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads) {
  int Vximax = ran_h->vr.Vx.x, Vxjmax = ran_h->vr.Vx.y, Vxkmax = ran_h->vr.Vx.z;
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
  VxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran_d);
  cudaDeviceSynchronize();
  // cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vx  update: %s\n", cudaGetErrorString(err));
  ZeroVx_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran_d);
  ZeroVx_XZ<<<ZeroXZBlocks, threadsPerBlock>>>(aft_d, ran_d);
  
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vx  zero  : %s\n", cudaGetErrorString(err));
  DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran_d, 'X');
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vx  add   : %s\n", cudaGetErrorString(err));

}
// Vyクラス的な
void Vy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads) {

  int Vyimax = ran_h->vr.Vy.x, Vyjmax = ran_h->vr.Vy.y, Vykmax = ran_h->vr.Vy.z;

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
  VyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran_d);
  // cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vy  update: %s\n", cudaGetErrorString(err));
  ZeroVy_YX<<<ZeroYXBlocks, threadsPerBlock>>>(aft_d, ran_d);
  ZeroVy_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran_d);
 
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vy  zero  : %s\n", cudaGetErrorString(err));

  //全方向加算
  DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran_d, 'Y');
  cudaDeviceSynchronize();
  cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vy  add   : %s\n", cudaGetErrorString(err));

}
// Vzクラス的な
void Vz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads) {

  int Vzimax = ran_h->vr.Vz.x, Vzjmax = ran_h->vr.Vz.y, Vzkmax = ran_h->vr.Vz.z;

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
  VzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran_d);
  cudaDeviceSynchronize();
  // cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vz  update: %s\n", cudaGetErrorString(err));
  ZeroVz_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran_d);
  ZeroVz_ZY<<<ZeroZYBlocks, threadsPerBlock>>>(aft_d, ran_d);
 
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vz  zero  : %s\n", cudaGetErrorString(err));
  //全方向加算
  DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran_d, 'Z');
  cudaDeviceSynchronize();
  // err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  // printf("CUDA kernel error Vz  add   : %s\n", cudaGetErrorString(err));

}
//粒子速度計算
void Vel(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads) {
  Vx(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
  Vy(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
  Vz(aft_d, bef_d, ma_d, dif_d, ran_d, ran_h, threads);
}

// __global__ void AccelerationCalculation(AccCoord *Acc, BefAft *aft, BefAft *bef, Diff *dif, Coord *out, Range *ran, int *outNum) {

//   int i = blockIdx.x * blockDim.x + threadIdx.x;

//   int ymax = ran->sr.Txx.y, zmax = ran->sr.Txx.z;
//   printf("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaammmmmmmmmmmm%d\n", *outNum);
//   if (i < *outNum) {

//     printf("ok\n");
//     int x = out[i].x;
//     int y = out[i].y;
//     int z = out[i].z;

//     // 1Dインデックスの計算
//     int idxX = (x - 1) * (ymax * zmax) +       y * zmax + z;
//     int idxY =       x * (ymax * zmax) + (y - 1) * zmax + z;
//     int idxZ =       x * (ymax * zmax) +       y * zmax + (z - 1);
//     int idx  =       x * (ymax * zmax) +       y * zmax + z;

//     Acc[i].x = ((*(aft->va.Vx + idxX) - *(bef->va.Vx + idxX)) / dif->dt + (*(aft->va.Vx + idx) - *(bef->va.Vx + idx)) / dif->dt) / 2;

//     Acc[i].y = ((*(aft->va.Vy + idxY) - *(bef->va.Vy + idxY)) / dif->dt + (*(aft->va.Vy + idx) - *(bef->va.Vy + idx)) / dif->dt) / 2;

//     Acc[i].z = ((*(aft->va.Vz + idxZ) - *(bef->va.Vz + idxZ)) / dif->dt + (*(aft->va.Vz + idx) - *(bef->va.Vz + idx)) / dif->dt) / 2;
//     printf("%d:%f,%f,%f\n", i, Acc[i].x, Acc[i].y, Acc[i].z);
//   }
// }

//更新
__global__ void swapTxx(BefAft *aft, BefAft *bef, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Txximax = ran->sr.Txx.x, Txxjmax = ran->sr.Txx.y, Txxkmax = ran->sr.Txx.z;
  int idx_Txx = k * Txximax * Txxjmax + j * Txximax + i;
  if (i < Txximax && j < Txxjmax && k < Txxkmax) {
    *(bef->sa.Txx  + idx_Txx) = *(aft->sa.Txx  + idx_Txx);
    *(bef->sa.Txxx + idx_Txx) = *(aft->sa.Txxx + idx_Txx);
    *(bef->sa.Txxy + idx_Txx) = *(aft->sa.Txxy + idx_Txx);
    *(bef->sa.Txxz + idx_Txx) = *(aft->sa.Txxz + idx_Txx);
  }
}

__global__ void swapTyy(BefAft *aft, BefAft *bef, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tyyimax = ran->sr.Tyy.x, Tyyjmax = ran->sr.Tyy.y, Tyykmax = ran->sr.Tyy.z;
  int idx_Tyy = k * Tyyimax * Tyyjmax + j * Tyyimax + i;


  if (i < Tyyimax && j < Tyyjmax && k < Tyykmax) {
    *(bef->sa.Tyy  + idx_Tyy) = *(aft->sa.Tyy  + idx_Tyy);
    *(bef->sa.Tyyx + idx_Tyy) = *(aft->sa.Tyyx + idx_Tyy);
    *(bef->sa.Tyyy + idx_Tyy) = *(aft->sa.Tyyy + idx_Tyy);
    *(bef->sa.Tyyz + idx_Tyy) = *(aft->sa.Tyyz + idx_Tyy);
  }
}

__global__ void swapTzz(BefAft *aft, BefAft *bef, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tzzimax = ran->sr.Tzz.x, Tzzjmax = ran->sr.Tzz.y, Tzzkmax = ran->sr.Tzz.z;
  int idx_Tzz = k * Tzzimax * Tzzjmax + j * Tzzimax + i;


  if (i < Tzzimax && j < Tzzjmax && k < Tzzkmax) {
    *(bef->sa.Tzz  + idx_Tzz) = *(aft->sa.Tzz  + idx_Tzz);
    *(bef->sa.Tzzx + idx_Tzz) = *(aft->sa.Tzzx + idx_Tzz);
    *(bef->sa.Tzzy + idx_Tzz) = *(aft->sa.Tzzy + idx_Tzz);
    *(bef->sa.Tzzz + idx_Tzz) = *(aft->sa.Tzzz + idx_Tzz);
  }
}

__global__ void swapTxy(BefAft *aft, BefAft *bef, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Txyimax = ran->tr.Txy.x, Txyjmax = ran->tr.Txy.y, Txykmax = ran->tr.Txy.z;
  int idx_Txy = k * Txyimax * Txyjmax + j * Txyimax + i;

  if (i < Txyimax && j < Txyjmax && k < Txykmax) {
    *(bef->ta.Txy  + idx_Txy) = *(aft->ta.Txy  + idx_Txy);
    *(bef->ta.Txyx + idx_Txy) = *(aft->ta.Txyx + idx_Txy);
    *(bef->ta.Txyy + idx_Txy) = *(aft->ta.Txyy + idx_Txy);
  }
}

__global__ void swapTyz(BefAft *aft, BefAft *bef, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tyzimax = ran->tr.Tyz.x, Tyzjmax = ran->tr.Tyz.y, Tyzkmax = ran->tr.Tyz.z;
  int idx_Tyz = k * Tyzimax * Tyzjmax + j * Tyzimax + i;

  if (i < Tyzimax && j < Tyzjmax && k < Tyzkmax) {
    *(bef->ta.Tyz  + idx_Tyz) = *(aft->ta.Tyz  + idx_Tyz);
    *(bef->ta.Tyzy + idx_Tyz) = *(aft->ta.Tyzy + idx_Tyz);
    *(bef->ta.Tyzz + idx_Tyz) = *(aft->ta.Tyzz + idx_Tyz);
  }
}

__global__ void swapTzx(BefAft *aft, BefAft *bef, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tzximax = ran->tr.Tzx.x, Tzxjmax = ran->tr.Tzx.y, Tzxkmax = ran->tr.Tzx.z;
  int idx_Tzx = k * Tzximax * Tzxjmax + j * Tzximax + i;

  if (i < Tzximax && j < Tzxjmax && k < Tzxkmax) {
    *(bef->ta.Tzx  + idx_Tzx) = *(aft->ta.Tzx  + idx_Tzx);
    *(bef->ta.Tzxz + idx_Tzx) = *(aft->ta.Tzxz + idx_Tzx);
    *(bef->ta.Tzxx + idx_Tzx) = *(aft->ta.Tzxx + idx_Tzx);
  }
}

__global__ void swapVx(BefAft *aft, BefAft *bef, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Vximax = ran->vr.Vx.x, Vxjmax = ran->vr.Vx.y, Vxkmax = ran->vr.Vx.z;
  int idx_Vx = k * Vximax * Vxjmax + j * Vximax + i;

  if (i < Vximax && j < Vxjmax && k < Vxkmax) {
    *(bef->va.Vx  + idx_Vx) = *(aft->va.Vx  + idx_Vx);
    *(bef->va.Vxx + idx_Vx) = *(aft->va.Vxx + idx_Vx);
    *(bef->va.Vxy + idx_Vx) = *(aft->va.Vxy + idx_Vx);
    *(bef->va.Vxz + idx_Vx) = *(aft->va.Vxz + idx_Vx);
  }
}

__global__ void swapVy(BefAft *aft, BefAft *bef, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Vyimax = ran->vr.Vy.x, Vyjmax = ran->vr.Vy.y, Vykmax = ran->vr.Vy.z;
  int idx_Vy = k * Vyimax * Vyjmax + j * Vyimax + i;

  if (i < Vyimax && j < Vyjmax && k < Vykmax) {
    *(bef->va.Vy  + idx_Vy) = *(aft->va.Vy  + idx_Vy);
    *(bef->va.Vyx + idx_Vy) = *(aft->va.Vyx + idx_Vy);
    *(bef->va.Vyy + idx_Vy) = *(aft->va.Vyy + idx_Vy);
    *(bef->va.Vyz + idx_Vy) = *(aft->va.Vyz + idx_Vy);
  }
}

__global__ void swapVz(BefAft *aft, BefAft *bef, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Vzimax = ran->vr.Vz.x, Vzjmax = ran->vr.Vz.y, Vzkmax = ran->vr.Vz.z;
  int idx_Vz = k * Vzimax * Vzjmax + j * Vzimax + i;

  if (i < Vzimax && j < Vzjmax && k < Vzkmax) {
    *(bef->va.Vz  + idx_Vz) = *(aft->va.Vz  + idx_Vz);
    *(bef->va.Vzx + idx_Vz) = *(aft->va.Vzx + idx_Vz);
    *(bef->va.Vzy + idx_Vz) = *(aft->va.Vzy + idx_Vz);
    *(bef->va.Vzz + idx_Vz) = *(aft->va.Vzz + idx_Vz);
  }
}

void swapBefAft(BefAft *aft, BefAft *bef, Range *ran_h, Range *ran_d, Coord threads) {
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
  swapTxx<<<SwapTxxBlocks, threadsPerBlock>>>(aft, bef, ran_d);
  swapTyy<<<SwapTyyBlocks, threadsPerBlock>>>(aft, bef, ran_d);
  swapTzz<<<SwapTzzBlocks, threadsPerBlock>>>(aft, bef, ran_d);
  swapTxy<<<SwapTxyBlocks, threadsPerBlock>>>(aft, bef, ran_d);
  swapTyz<<<SwapTyzBlocks, threadsPerBlock>>>(aft, bef, ran_d);
  swapTzx<<<SwapTzxBlocks, threadsPerBlock>>>(aft, bef, ran_d);
  swapVx<<<SwapVxBlocks, threadsPerBlock>>>(aft, bef, ran_d);
  swapVy<<<SwapVyBlocks, threadsPerBlock>>>(aft, bef, ran_d);
  swapVz<<<SwapVzBlocks, threadsPerBlock>>>(aft, bef, ran_d);
  cudaDeviceSynchronize();
}
