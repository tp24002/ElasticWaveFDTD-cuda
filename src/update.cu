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
__global__ void TxxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Coord Txx) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = Txx.x, jmax = Txx.y, kmax = Txx.z;

  // 


  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax - 1) {
    int idx = i * jmax * kmax + j * kmax + k;
    int idx_im1 = (i - 1) * jmax * kmax + j * kmax + k;
    int idx_jm1 = i * jmax * kmax + (j - 1) * kmax + k;
    int idx_km1 = i * jmax * kmax + j * kmax + (k - 1);

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
__global__ void TyyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Coord Tyy) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = Tyy.x, jmax = Tyy.y, kmax = Tyy.z;

  // 


  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax - 1) {
    int idx = i * jmax * kmax + j * kmax + k;
    int idx_im1 = (i - 1) * jmax * kmax + j * kmax + k;
    int idx_jm1 = i * jmax * kmax + (j - 1) * kmax + k;
    int idx_km1 = i * jmax * kmax + j * kmax + (k - 1);

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
__global__ void TzzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Coord Tzz) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
  
    int imax = Tzz.x, jmax = Tzz.y, kmax = Tzz.z;

    // 
    
    

    if(i <= imax - 1 && j <= jmax - 1 && k <= kmax - 1) {
        int idx = i * jmax * kmax + j * kmax + k;
        int idx_im1 = (i - 1) * jmax * kmax + j * kmax + k;
        int idx_jm1 = i * jmax * kmax + (j - 1) * kmax + k;
        int idx_km1 = i * jmax * kmax + j * kmax + (k - 1);

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
__global__ void ZeroT_XY(BefAft *aft, Coord ranmax, char check) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = ranmax.x, jmax = ranmax.y, kmax = ranmax.z;

  if (j > jmax || i > imax) {
    return;
  }

  // 


  int idx_0 = i * jmax * kmax + j * kmax;
  int idx_kmax = i * jmax * kmax + j * kmax + kmax;

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
__global__ void ZeroT_YZ(BefAft *aft, Coord ranmax, char check) {
  int k = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = ranmax.x, jmax = ranmax.y, kmax = ranmax.z;

  if (k > kmax - 1 || j > jmax) {
    return;
  }

  // 


  int idx_0 = j * kmax + k;
  int idx_imax = imax * jmax * kmax + j * kmax + k;

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
__global__ void ZeroT_ZX(BefAft *aft, Coord ranmax, char check) {
  int k = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int i = blockIdx.y * blockDim.y + threadIdx.y + 1;

  int imax = ranmax.x, jmax = ranmax.y, kmax = ranmax.z;

  if (k > kmax - 1 || i > imax - 1) {
    return;
  }

  // 


  int idx_0 = i * jmax * kmax + k;
  int idx_jmax = i * jmax * kmax + jmax * kmax + k;

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
__global__ void DirectionalAdd(BefAft *aft, Impulse *ip, Coord ranmax, char check) {
  // 1Dインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = ranmax.x, jmax = ranmax.y, kmax = ranmax.z;


  // 1Dインデックス化
  int idx = i * jmax * kmax + j * kmax + k;

  // 範囲外のスレッドは処理をスキップ
  if (i > imax || j > jmax || k > kmax) {
    return;
  }

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


// Txxクラス的な(Blocks大丈夫かな？)
__device__ void Txx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Impulse *ip_d, Coord threads) {
  char check = 'X';
  Coord ranmax;

  int Txximax = ran->sr.Txx.x, Txxjmax = ran->sr.Txx.y, Txxkmax = ran->sr.Txx.z;
  initDeviceCoord(&ranmax, Txximax, Txxjmax, Txxkmax);

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
  TxxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran->sr.Txx);
  // 0 padding
  ZeroT_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  ZeroT_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ranmax, check);
  ZeroT_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ranmax, check);

  //全方向加算
  DirectionalAdd<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ip_d, ranmax, check);
  cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  if (err != cudaSuccess) {
      printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
  }


}
// Tyyクラス的な(Blocks大丈夫かな？)
__device__ void Tyy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Impulse *ip_d, Coord threads) {
  char check = 'Y';
  Coord ranmax;

  int Tyyimax = ran->sr.Tyy.x, Tyyjmax = ran->sr.Tyy.y, Tyykmax = ran->sr.Tyy.z;
  initDeviceCoord(&ranmax, Tyyimax, Tyyjmax, Tyykmax);

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
  cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  if (err != cudaSuccess) {
    printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
  }

}
// Tzzクラス的な(Blocks大丈夫かな？)
__device__ void Tzz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Impulse *ip_d, Coord threads) {
  char check = 'Z';
  Coord ranmax;

  int Tzzimax = ran->sr.Tzz.x, Tzzjmax = ran->sr.Tzz.y, Tzzkmax = ran->sr.Tzz.z;
  initDeviceCoord(&ranmax, Tzzimax, Tzzjmax, Tzzkmax);
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
  cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  if (err != cudaSuccess) {
      printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
  }
 
}
// 垂直応力計算(main呼び出し関数)
__global__ void Sig(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Impulse *ip_d, Coord threads) {
  Txx(aft_d, bef_d, ma_d, dif_d, ran, ip_d, threads);
  Tyy(aft_d, bef_d, ma_d, dif_d, ran, ip_d, threads);
  Tzz(aft_d, bef_d, ma_d, dif_d, ran, ip_d, threads);
}

// せん断応力

// せん断応力更新関数
__global__ void TxyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1; // 始点を+1

  int imax = tr.Txy.x, jmax = tr.Txy.y, kmax = tr.Txy.z;
  double Hzetadx, Hzetady, Hmu, Hgamma;

  // 


  if (i <= imax && j <= jmax && k <= kmax - 1) {
    // 各インデックスの計算
    int idx = i * jmax * kmax + j * kmax + k;
    int idx_ip1 = (i + 1) * jmax * kmax + j * kmax + k;
    int idx_jp1 = i * jmax * kmax + (j + 1) * kmax + k;

    // PML:減衰係数,計算領域:摩擦定数
    Hzetadx = 4. * pow((1. / *(ma->zetadx + idx_ip1)) + (1. / *(ma->zetadx + idx_jp1)) + (1. / *(ma->zetadx + idx_ip1 - jmax * kmax)) + (1. / *(ma->zetadx + idx - jmax * kmax)), -1.);
    Hzetady = 4. * pow((1. / *(ma->zetady + idx_ip1)) + (1. / *(ma->zetady + idx_jp1)) + (1. / *(ma->zetady + idx_ip1 - jmax * kmax)) + (1. / *(ma->zetady + idx - jmax * kmax)), -1.);

    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / *(ma->mu + idx_ip1)) + (1. / *(ma->mu + idx_jp1)) + (1. / *(ma->mu + idx_ip1 - jmax * kmax)) + (1. / *(ma->mu + idx - jmax * kmax)), -1.);

    // 第1粘性定数
    Hgamma = 4. * pow((1. / *(ma->gamma + idx_ip1)) + (1. / *(ma->gamma + idx_jp1)) + (1. / *(ma->gamma + idx_ip1 - jmax * kmax)) + (1. / *(ma->gamma + idx - jmax * kmax)), -1.);

    *(aft->ta.Txyx + idx) = (2.0 - Hzetadx * dif->dt) / (2.0 + Hzetadx * dif->dt) * *(bef->ta.Txyx + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadx * dif->dt) * (*(aft->va.Vy + idx_ip1) - *(aft->va.Vy + idx)) / dif->dx
        - 2.0 * Hgamma / (2.0 + Hzetadx * dif->dt) * (*(bef->va.Vy + idx_ip1) - *(bef->va.Vy + idx)) / dif->dx;

    *(aft->ta.Txyy + idx) = (2.0 - Hzetady * dif->dt) / (2.0 + Hzetady * dif->dt) * *(bef->ta.Txyy + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetady * dif->dt) * (*(aft->va.Vx + idx_jp1) - *(aft->va.Vx + idx)) / dif->dy
        - 2.0 * Hgamma / (2.0 + Hzetady * dif->dt) * (*(bef->va.Vx + idx_jp1) - *(bef->va.Vx + idx)) / dif->dy;
  }
}

// せん断応力更新関数
__global__ void TyzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = tr.Tyz.x, jmax = tr.Tyz.y, kmax = tr.Tyz.z;
  double Hzetady, Hzetadz, Hmu, Hgamma;

  // 


  if (i <= imax - 1 && j <= jmax && k <= kmax) {
    // 各インデックスの計算
    int idx = i * jmax * kmax + j * kmax + k;
    int idx_jp1_kp1 = i * jmax * kmax + (j + 1) * kmax + (k + 1);
    int idx_jp1 = i * jmax * kmax + (j + 1) * kmax + k;
    int idx_kp1 = i * jmax * kmax + j * kmax + (k + 1);

    // PML:減衰係数,計算領域:摩擦定数
    Hzetady = 4. * pow((1. / *(ma->zetady + idx_jp1_kp1)) + (1. / *(ma->zetady + idx_jp1)) + (1. / *(ma->zetady + idx_kp1)) + (1. / *(ma->zetady + idx)), -1.);
    Hzetadz = 4. * pow((1. / *(ma->zetadz + idx_jp1_kp1)) + (1. / *(ma->zetadz + idx_jp1)) + (1. / *(ma->zetadz + idx_kp1)) + (1. / *(ma->zetadz + idx)), -1.);
    
    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / *(ma->mu + idx_jp1_kp1)) + (1. / *(ma->mu + idx_jp1)) + (1. / *(ma->mu + idx_kp1)) + (1. / *(ma->mu + idx)), -1.);

    // 第1粘性定数
    Hgamma = 4. * pow((1. / *(ma->gamma + idx_jp1_kp1)) + (1. / *(ma->gamma + idx_jp1)) + (1. / *(ma->gamma + idx_kp1)) + (1. / *(ma->gamma + idx)), -1.);

    *(aft->ta.Tyzy + idx) = (2.0 - Hzetady * dif->dt) / (2.0 + Hzetady * dif->dt) * *(bef->ta.Tyzy + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetady * dif->dt) * (*(aft->va.Vz + idx_jp1) - *(aft->va.Vz + idx)) / dif->dy
        - 2.0 * Hgamma / (2.0 + Hzetady * dif->dt) * (*(bef->va.Vz + idx_jp1) - *(bef->va.Vz + idx)) / dif->dy;

    *(aft->ta.Tyzz + idx) = (2.0 - Hzetadz * dif->dt) / (2.0 + Hzetadz * dif->dt) * *(bef->ta.Tyzz + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadz * dif->dt) * (*(aft->va.Vy + idx_kp1) - *(aft->va.Vy + idx)) / dif->dz
        - 2.0 * Hgamma / (2.0 + Hzetadz * dif->dt) * (*(bef->va.Vy + idx_kp1) - *(bef->va.Vy + idx)) / dif->dz;
  }
}

// せん断応力更新関数
__global__ void TzxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = tr.Tzx.x, jmax = tr.Tzx.y, kmax = tr.Tzx.z;
  double Hzetadx, Hzetadz, Hmu, Hgamma;

  


  if (i <= imax && j <= jmax - 1 && k <= kmax) {
    // 各インデックスの計算
    int idx = i * jmax * kmax + j * kmax + k;
    int idx_ip1 = (i + 1) * jmax * kmax + j * kmax + k;
    int idx_kp1 = i * jmax * kmax + j * kmax + (k + 1);
    int idx_ip1_kp1 = (i + 1) * jmax * kmax + j * kmax + (k + 1);

    // PML:減衰係数,計算領域:摩擦定数
    Hzetadx = 4. * pow((1. / *(ma->zetadx + idx_ip1_kp1)) + (1. / *(ma->zetadx + idx_ip1)) + (1. / *(ma->zetadx + idx_kp1)) + (1. / *(ma->zetadx + idx)), -1.);
    Hzetadz = 4. * pow((1. / *(ma->zetadz + idx_ip1_kp1)) + (1. / *(ma->zetadz + idx_ip1)) + (1. / *(ma->zetadz + idx_kp1)) + (1. / *(ma->zetadz + idx)), -1.);

    // 第2ラメ，横弾性係数(剛性率)
    Hmu = 4. * pow((1. / *(ma->mu + idx_ip1_kp1)) + (1. / *(ma->mu + idx_ip1)) + (1. / *(ma->mu + idx_kp1)) + (1. / *(ma->mu + idx)), -1.);

    // 第1粘性定数
    Hgamma = 4. * pow((1. / *(ma->gamma + idx_ip1_kp1)) + (1. / *(ma->gamma + idx_ip1)) + (1. / *(ma->gamma + idx_kp1)) + (1. / *(ma->gamma + idx)), -1.);

    *(aft->ta.Tzxz + idx) = (2.0 - Hzetadz * dif->dt) / (2.0 + Hzetadz * dif->dt) * *(bef->ta.Tzxz + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadz * dif->dt) * (*(aft->va.Vx + idx_kp1) - *(aft->va.Vx + idx)) / dif->dz
        - 2.0 * Hgamma / (2.0 + Hzetadz * dif->dt) * (*(bef->va.Vx + idx_kp1) - *(bef->va.Vx + idx)) / dif->dz;

    *(aft->ta.Tzxx + idx) = (2.0 - Hzetadx * dif->dt) / (2.0 + Hzetadx * dif->dt) * *(bef->ta.Tzxx + idx)
        + 2.0 * (Hmu * dif->dt + Hgamma) / (2.0 + Hzetadx * dif->dt) * (*(aft->va.Vz + idx_ip1) - *(aft->va.Vz + idx)) / dif->dx
        - 2.0 * Hgamma / (2.0 + Hzetadx * dif->dt) * (*(bef->va.Vz + idx_ip1) - *(bef->va.Vz + idx)) / dif->dx;
  }
}

__global__ void ZeroTxy(BefAft *aft, Coord Tmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = Tmax.x, jmax = Tmax.y, kmax = Tmax.z;

  


  if (i <= imax && j <= jmax) {
    int idx_0 = i * jmax * kmax + j * kmax;
    int idx_kmax = i * jmax * kmax + j * kmax + kmax;

    *(aft->ta.Txyx + idx_0) = 0.0;
    *(aft->ta.Txyx + idx_kmax) = 0.0;
    *(aft->ta.Txyy + idx_0) = 0.0;
    *(aft->ta.Txyy + idx_kmax) = 0.0;
  }
}

__global__ void ZeroTyz(BefAft *aft, Coord Tmax) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = Tmax.x;
  int jmax = Tmax.y;
  int kmax = Tmax.z;



  if (j <= jmax && k <= kmax) {
    int idx_0 = j * kmax + k;
    int idx_imax = imax * jmax * kmax + j * kmax + k;

    *(aft->ta.Tyzy + idx_0) = 0.0;
    *(aft->ta.Tyzy + idx_imax) = 0.0;
    *(aft->ta.Tyzz + idx_0) = 0.0;
    *(aft->ta.Tyzz + idx_imax) = 0.0;
  }
}

__global__ void ZeroTzx(BefAft *aft, Coord Tmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = Tmax.x;
  int kmax = Tmax.z;

  int jmax = Tmax.y;
  

  if (i <= imax && k <= kmax) {
    int idx_0 = i * jmax * kmax + k;
    int idx_jmax = i * jmax * kmax + Tmax.y * kmax + k;

    *(aft->ta.Tzxx + idx_0) = 0.0;
    *(aft->ta.Tzxx + idx_jmax) = 0.0;
    *(aft->ta.Tzxz + idx_0) = 0.0;
    *(aft->ta.Tzxz + idx_jmax) = 0.0;
  }
}

__global__ void DirectionalAddT(BefAft *aft, Coord Tmax, char check) {
  // スレッドインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = Tmax.x;
  int jmax = Tmax.y;
  int kmax = Tmax.z;

  if (i > imax || j > jmax || k > kmax) {
    return;
  }

  int idx = i * jmax * kmax + j * kmax + k;

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
__device__ void Txy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads) {

  int Txyimax = ran->tr.Txy.x, Txyjmax = ran->tr.Txy.y, Txykmax = ran->tr.Txy.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Txyimax + threadsPerBlock.x - 1)     / threadsPerBlock.x,
                    (Txyjmax + threadsPerBlock.y - 1)     / threadsPerBlock.y,
                    (Txykmax - 1 + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroXYBlocks((Txyimax + threadsPerBlock.x - 1) / threadsPerBlock.x,(Txyjmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Txyimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Txyjmax + threadsPerBlock.y) / threadsPerBlock.y,
                            (Txykmax + threadsPerBlock.z) / threadsPerBlock.z);                    
  TxyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran->tr);
  ZeroTxy<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran->tr.Txy);
 
  DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran->tr.Txy, 'Z');
  cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  if (err != cudaSuccess) {
      printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
  }
 
}
// Tyzクラス的な
__device__ void Tyz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads) {
  
  int Tyzimax = ran->tr.Tyz.x, Tyzjmax = ran->tr.Tyz.y, Tyzkmax = ran->tr.Tyz.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Tyzimax - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (Tyzjmax + threadsPerBlock.y - 1)     / threadsPerBlock.y,
                    (Tyzkmax + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  dim3 ZeroYZBlocks((Tyzjmax + threadsPerBlock.x - 1) / threadsPerBlock.x,(Tyzkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);
  dim3 DirectionalAddBlocks((Tyzimax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Tyzjmax + threadsPerBlock.y) / threadsPerBlock.y,
                            (Tyzkmax + threadsPerBlock.z) / threadsPerBlock.z);                    
  TyzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran->tr);
  ZeroTyz<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran->tr.Tyz);
 
  DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran->tr.Tyz, 'X');
  cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  if (err != cudaSuccess) {
      printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
  }


}
// Tzxクラス的な
__device__ void Tzx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads) {
  
  int Tzximax = ran->tr.Tzx.x, Tzxjmax = ran->tr.Tzx.y, Tzxkmax = ran->tr.Tzx.z;

  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 UpdateBlocks((Tzximax + threadsPerBlock.x - 1)     / threadsPerBlock.x,
                    (Tzxjmax - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (Tzxkmax + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  dim3 ZeroZXBlocks((Tzximax + threadsPerBlock.x - 1) / threadsPerBlock.x,(Tzxkmax + threadsPerBlock.y - 1) / threadsPerBlock.y);   
  dim3 DirectionalAddBlocks((Tzximax + threadsPerBlock.x) / threadsPerBlock.x,
                            (Tzxjmax + threadsPerBlock.y) / threadsPerBlock.y, 
                            (Tzxkmax + threadsPerBlock.z) / threadsPerBlock.z);                  
  TzxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran->tr);
  ZeroTzx<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran->tr.Tzx);
 
  DirectionalAddT<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran->tr.Tzx , 'Y');
  cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  if (err != cudaSuccess) {
      printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
  }

}
// せん断応力計算(main呼び出し関数)
__global__ void Tau(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads) {
  Txy(aft_d, bef_d, ma_d, dif_d, ran, threads);
  Tyz(aft_d, bef_d, ma_d, dif_d, ran, threads);
  Tzx(aft_d, bef_d, ma_d, dif_d, ran, threads);
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
    // 1D indexing for 3D arrays
    int idx_i1 = (i + 1) * jmax * kmax + j * kmax + k;
    int idx_i  = i * jmax * kmax + j * kmax + k;
    int idx_j1 = i * jmax * kmax + (j - 1) * kmax + k;
    int idx_k1 = i * jmax * kmax + j * kmax + (k - 1);

    Azetaxx = (ma->zetaxx[idx_i1] + ma->zetaxx[idx_i]) / 2.;
    Azetaxy = (ma->zetaxy[idx_i1] + ma->zetaxy[idx_i]) / 2.;
    Azetaxz = (ma->zetaxz[idx_i1] + ma->zetaxz[idx_i]) / 2.;
    Arho    = (ma->rho[idx_i1] + ma->rho[idx_i]) / 2.;

    aft->va.Vxx[idx_i] = (2. * Arho - Azetaxx * dif->dt) / (2. * Arho + Azetaxx * dif->dt) * bef->va.Vxx[idx_i]
        + 2. * dif->dt / (2. * Arho + Azetaxx * dif->dt) * (bef->sa.Txx[idx_i1] - bef->sa.Txx[idx_i]) / dif->dx;

    aft->va.Vxy[idx_i] = (2. * Arho - Azetaxy * dif->dt) / (2. * Arho + Azetaxy * dif->dt) * bef->va.Vxy[idx_i]
        + 2. * dif->dt / (2. * Arho + Azetaxy * dif->dt) * (bef->ta.Txy[idx_i] - bef->ta.Txy[idx_j1]) / dif->dy;

    aft->va.Vxz[idx_i] = (2. * Arho - Azetaxz * dif->dt) / (2. * Arho + Azetaxz * dif->dt) * bef->va.Vxz[idx_i]
        + 2. * dif->dt / (2. * Arho + Azetaxz * dif->dt) * (bef->ta.Tzx[idx_i] - bef->ta.Tzx[idx_k1]) / dif->dz;
  }
}
// 粒子速度更新関数
__global__ void VyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, VelRan vr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

  int imax = vr.Vy.x, jmax = vr.Vy.y, kmax = vr.Vy.z;
  double Azetayx, Azetayy, Azetayz, Arho;

  


  if (i <= imax - 1 && j <= jmax && k <= kmax - 1) {
    // 各インデックスの計算
    int idx = i * jmax * kmax + j * kmax + k;
    int idx_im1 = (i - 1) * jmax * kmax + j * kmax + k;
    int idx_jp1 = i * jmax * kmax + (j + 1) * kmax + k;
    int idx_km1 = i * jmax * kmax + j * kmax + (k - 1);

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
__global__ void VzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, VelRan vr) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = vr.Vz.x, jmax = vr.Vz.y, kmax = vr.Vz.z;
  double Azetazx, Azetazy, Azetazz, Arho;



  if(i <= imax - 1 && j <= jmax - 1 && k <= kmax) {
    // 1D indexing for 3D arrays
    int idx_i  = i * jmax * kmax + j * kmax + k;
    int idx_im1 = (i - 1) * jmax * kmax + j * kmax + k;
    int idx_jm1 = i * jmax * kmax + (j - 1) * kmax + k;
    int idx_kp1 = i * jmax * kmax + j * kmax + (k + 1);

    Azetazx = (ma->zetazx[idx_kp1] + ma->zetazx[idx_i]) / 2.;
    Azetazy = (ma->zetazy[idx_kp1] + ma->zetazy[idx_i]) / 2.;
    Azetazz = (ma->zetazz[idx_kp1] + ma->zetazz[idx_i]) / 2.;
    Arho    = (ma->rho[idx_kp1] + ma->rho[idx_i]) / 2.;

    aft->va.Vzx[idx_i] = (2. * Arho - Azetazx * dif->dt) / (2. * Arho + Azetazx * dif->dt) * bef->va.Vzx[idx_i]
        + 2. * dif->dt / (2. * Arho + Azetazx * dif->dt) * (bef->ta.Tzx[idx_i] - bef->ta.Tzx[idx_im1]) / dif->dx;

    aft->va.Vzy[idx_i] = (2. * Arho - Azetazy * dif->dt) / (2. * Arho + Azetazy * dif->dt) * bef->va.Vzy[idx_i]
        + 2. * dif->dt / (2. * Arho + Azetazy * dif->dt) * (bef->ta.Tyz[idx_i] - bef->ta.Tyz[idx_jm1]) / dif->dy;

    aft->va.Vzz[idx_i] = (2. * Arho - Azetazz * dif->dt) / (2. * Arho + Azetazz * dif->dt) * bef->va.Vzz[idx_i]
        + 2. * dif->dt / (2. * Arho + Azetazz * dif->dt) * (bef->sa.Tzz[idx_kp1] - bef->sa.Tzz[idx_i]) / dif->dz;
  }
}

__global__ void ZeroVx_XY(BefAft *aft, Coord Vmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

  int imax = Vmax.x, jmax = Vmax.y, kmax = Vmax.z;


  if(i <= imax && j <= jmax) {
    // 1D indexing for 3D arrays
    int idx_0   = i * jmax * kmax + j * kmax;
    int idx_kmax = i * jmax * kmax + j * kmax + kmax;

    aft->va.Vxx[idx_0] = 0.0;
    aft->va.Vxx[idx_kmax] = 0.0;

    aft->va.Vxy[idx_0] = 0.0;
    aft->va.Vxy[idx_kmax] = 0.0;

    aft->va.Vxz[idx_0] = 0.0;
    aft->va.Vxz[idx_kmax] = 0.0;
  }
}

__global__ void ZeroVx_XZ(BefAft *aft, Coord Vmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = Vmax.x, jmax = Vmax.y, kmax = Vmax.z;


  if(i <= imax && k <= kmax) {
    // 1D indexing for 3D arrays
    int idx_0   = i * jmax * kmax + k;
    int idx_jmax = i * jmax * kmax + jmax * kmax + k;

    aft->va.Vxx[idx_0] = 0.0;
    aft->va.Vxx[idx_jmax] = 0.0;

    aft->va.Vxy[idx_0] = 0.0;
    aft->va.Vxy[idx_jmax] = 0.0;

    aft->va.Vxz[idx_0] = 0.0;
    aft->va.Vxz[idx_jmax] = 0.0;
  }
}

__global__ void ZeroVy_YX(BefAft *aft, Coord Vmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = Vmax.x;
  int jmax = Vmax.y;
  int kmax = Vmax.z;

  if (i <= imax && j <= jmax) {
    int idx_0 = i * jmax * kmax + j * kmax;
    int idx_zmax = i * jmax * kmax + j * kmax + kmax;

    *(aft->va.Vyx + idx_0) = 0.0;
    *(aft->va.Vyx + idx_zmax) = 0.0;
    *(aft->va.Vyy + idx_0) = 0.0;
    *(aft->va.Vyy + idx_zmax) = 0.0;
    *(aft->va.Vyz + idx_0) = 0.0;
    *(aft->va.Vyz + idx_zmax) = 0.0;
  }
}

__global__ void ZeroVy_YZ(BefAft *aft, Coord Vmax) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = Vmax.x, jmax = Vmax.y, kmax = Vmax.z;


  if(j <= jmax && k <= kmax) {
    // 1D indexing for 3D arrays
    int idx_0   = j * kmax + k;
    int idx_imax = imax * jmax * kmax + j * kmax + k;

    aft->va.Vyx[idx_0] = 0.0;
    aft->va.Vyx[idx_imax] = 0.0;

    aft->va.Vyy[idx_0] = 0.0;
    aft->va.Vyy[idx_imax] = 0.0;

    aft->va.Vyz[idx_0] = 0.0;
    aft->va.Vyz[idx_imax] = 0.0;
  }
}

__global__ void ZeroVz_ZX(BefAft *aft, Coord Vmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int k = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = Vmax.x;
  int jmax = Vmax.y;
  int kmax = Vmax.z;

  if (i <= imax && k <= kmax) {
    int idx_0 = i * jmax * kmax + k;
    int idx_ymax = i * jmax * kmax + jmax * kmax + k;

    *(aft->va.Vzx + idx_0) = 0.0;
    *(aft->va.Vzx + idx_ymax) = 0.0;
    *(aft->va.Vzy + idx_0) = 0.0;
    *(aft->va.Vzy + idx_ymax) = 0.0;
    *(aft->va.Vzz + idx_0) = 0.0;
    *(aft->va.Vzz + idx_ymax) = 0.0;
  }
}

__global__ void ZeroVz_ZY(BefAft *aft, Coord Vmax) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;

  int imax = Vmax.x, jmax = Vmax.y, kmax = Vmax.z;


  if(j <= jmax && k <= kmax) {
    // 1D indexing for 3D arrays
    int idx_0   = j * kmax + k;
    int idx_imax = imax * jmax * kmax + j * kmax + k;

    aft->va.Vzx[idx_0] = 0.0;
    aft->va.Vzx[idx_imax] = 0.0;

    aft->va.Vzy[idx_0] = 0.0;
    aft->va.Vzy[idx_imax] = 0.0;

    aft->va.Vzz[idx_0] = 0.0;
    aft->va.Vzz[idx_imax] = 0.0;
  }
}

__global__ void DirectionalAddV(BefAft *aft, Coord Vmax, char check) {
  // スレッドインデックスの計算
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int imax = Vmax.x;
  int jmax = Vmax.y;
  int kmax = Vmax.z;

  if (i > imax || j > jmax || k > kmax) {
    return;
  }

  int idx = i * jmax * kmax + j * kmax + k;

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
__device__ void Vx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads) {
  
  int Vximax = ran->vr.Vx.x, Vxjmax = ran->vr.Vx.y, Vxkmax = ran->vr.Vx.z;

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
                
  VxUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran->vr);
  ZeroVx_XY<<<ZeroXYBlocks, threadsPerBlock>>>(aft_d, ran->vr.Vx);
  ZeroVx_XZ<<<ZeroXZBlocks, threadsPerBlock>>>(aft_d, ran->vr.Vx);
 
  DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran->vr.Vx , 'X');
  cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  if (err != cudaSuccess) {
      printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
  }

}
// Vyクラス的な
__device__ void Vy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads) {
  

  int Vyimax = ran->vr.Vy.x, Vyjmax = ran->vr.Vy.y, Vykmax = ran->vr.Vy.z;

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
  VyUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran->vr);
  ZeroVy_YX<<<ZeroYXBlocks, threadsPerBlock>>>(aft_d, ran->vr.Vy);
  ZeroVy_YZ<<<ZeroYZBlocks, threadsPerBlock>>>(aft_d, ran->vr.Vy);
 
  //全方向加算
  DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran->vr.Vy , 'Y');
  cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  if (err != cudaSuccess) {
      printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
  }

}
// Vzクラス的な
__device__ void Vz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads) {


  int Vzimax = ran->vr.Vz.x, Vzjmax = ran->vr.Vz.y, Vzkmax = ran->vr.Vz.z;

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
  VzUpdate<<<UpdateBlocks, threadsPerBlock>>>(aft_d, bef_d, ma_d, dif_d, ran->vr);
  ZeroVz_ZX<<<ZeroZXBlocks, threadsPerBlock>>>(aft_d, ran->vr.Vz);
  ZeroVz_ZY<<<ZeroZYBlocks, threadsPerBlock>>>(aft_d, ran->vr.Vz);
 
  //全方向加算
  DirectionalAddV<<<DirectionalAddBlocks, threadsPerBlock>>>(aft_d, ran->vr.Vz , 'Z');
  cudaError_t err = cudaGetLastError(); // カーネル呼び出し後にエラーチェック
  if (err != cudaSuccess) {
      printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
  }

}
//粒子速度計算
__global__ void Vel(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads) {
  Vx(aft_d, bef_d, ma_d, dif_d, ran, threads);
  Vy(aft_d, bef_d, ma_d, dif_d, ran, threads);
  Vz(aft_d, bef_d, ma_d, dif_d, ran, threads);
}

__global__ void AccelerationCalculation(AccCoord *Acc, BefAft *aft, BefAft *bef, Diff *dif, Coord *out, Range *ran, int *outNum, int *t, int *tmax) {
  int i = blockIdx.x * blockDim.x + threadIdx.x; // スレッドインデックスの計算

  int ymax = ran->sr.Txx.y, zmax = ran->sr.Txx.z;

  if (i < *outNum) {
    int x = out[i].x;
    int y = out[i].y;
    int z = out[i].z;

    // 1Dインデックスの計算
    int idxX = (x - 1) * (ymax * zmax) +       y * zmax + z;
    int idxY =       x * (ymax * zmax) + (y - 1) * zmax + z;
    int idxZ =       x * (ymax * zmax) +       y * zmax + (z - 1);
    int idx  =       x * (ymax * zmax) +       y * zmax + z;

    *(Acc->x + i * *tmax + *t) = ((*(aft->va.Vx + idxX) - *(bef->va.Vx + idxX)) / dif->dt +
                                (*(aft->va.Vx + idx) - *(bef->va.Vx + idx)) / dif->dt) / 2;

    *(Acc->y + i * *tmax + *t) = ((*(aft->va.Vy + idxY) - *(bef->va.Vy + idxY)) / dif->dt +
                                (*(aft->va.Vy + idx) - *(bef->va.Vy + idx)) / dif->dt) / 2;

    *(Acc->z + i * *tmax + *t) = ((*(aft->va.Vz + idxZ) - *(bef->va.Vz + idxZ)) / dif->dt +
                                (*(aft->va.Vz + idx) - *(bef->va.Vz + idx)) / dif->dt) / 2;
  }
}

//更新
__global__ void swapTxx(SigArr *aftSa, SigArr *befSa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Txximax = ran->sr.Txx.x, Txxjmax = ran->sr.Txx.y, Txxkmax = ran->sr.Txx.z;
  int idx_Txx = i * Txxjmax * Txxkmax + j * Txxkmax + k;

  if (i < Txximax && j < Txxjmax && k < Txxkmax) {
    *(befSa->Txx  + idx_Txx) = *(aftSa->Txx  + idx_Txx);
    *(befSa->Txxx + idx_Txx) = *(aftSa->Txxx + idx_Txx);
    *(befSa->Txxy + idx_Txx) = *(aftSa->Txxy + idx_Txx);
    *(befSa->Txxz + idx_Txx) = *(aftSa->Txxz + idx_Txx);
  }
}

__global__ void swapTyy(SigArr *aftSa, SigArr *befSa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tyyimax = ran->sr.Tyy.x, Tyyjmax = ran->sr.Tyy.y, Tyykmax = ran->sr.Tyy.z;
  int idx_Tyy = i * Tyyjmax * Tyykmax + j * Tyykmax + k;

  if (i < Tyyimax && j < Tyyjmax && k < Tyykmax) {
    *(befSa->Tyy  + idx_Tyy) = *(aftSa->Tyy  + idx_Tyy);
    *(befSa->Tyyx + idx_Tyy) = *(aftSa->Tyyx + idx_Tyy);
    *(befSa->Tyyy + idx_Tyy) = *(aftSa->Tyyy + idx_Tyy);
    *(befSa->Tyyz + idx_Tyy) = *(aftSa->Tyyz + idx_Tyy);
  }
}

__global__ void swapTzz(SigArr *aftSa, SigArr *befSa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tzzimax = ran->sr.Tzz.x, Tzzjmax = ran->sr.Tzz.y, Tzzkmax = ran->sr.Tzz.z;
  int idx_Tzz = i * Tzzjmax * Tzzkmax + j * Tzzkmax + k;

  if (i < Tzzimax && j < Tzzjmax && k < Tzzkmax) {
    *(befSa->Tzz  + idx_Tzz) = *(aftSa->Tzz  + idx_Tzz);
    *(befSa->Tzzx + idx_Tzz) = *(aftSa->Tzzx + idx_Tzz);
    *(befSa->Tzzy + idx_Tzz) = *(aftSa->Tzzy + idx_Tzz);
    *(befSa->Tzzz + idx_Tzz) = *(aftSa->Tzzz + idx_Tzz);
  }
}

__global__ void swapTxy(TauArr *aftTa, TauArr *befTa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Txyimax = ran->tr.Txy.x, Txyjmax = ran->tr.Txy.y, Txykmax = ran->tr.Txy.z;
  int idx_Txy = i * Txyjmax * Txykmax + j * Txykmax + k;

  if (i < Txyimax && j < Txyjmax && k < Txykmax) {
    *(befTa->Txy  + idx_Txy) = *(aftTa->Txy  + idx_Txy);
    *(befTa->Txyx + idx_Txy) = *(aftTa->Txyx + idx_Txy);
    *(befTa->Txyy + idx_Txy) = *(aftTa->Txyy + idx_Txy);
  }
}

__global__ void swapTyz(TauArr *aftTa, TauArr *befTa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tyzimax = ran->tr.Tyz.x, Tyzjmax = ran->tr.Tyz.y, Tyzkmax = ran->tr.Tyz.z;
  int idx_Tyz = i * Tyzjmax * Tyzkmax + j * Tyzkmax + k;

  if (i < Tyzimax && j < Tyzjmax && k < Tyzkmax) {
    *(befTa->Tyz  + idx_Tyz) = *(aftTa->Tyz  + idx_Tyz);
    *(befTa->Tyzy + idx_Tyz) = *(aftTa->Tyzy + idx_Tyz);
    *(befTa->Tyzz + idx_Tyz) = *(aftTa->Tyzz + idx_Tyz);
  }
}

__global__ void swapTzx(TauArr *aftTa, TauArr *befTa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Tzximax = ran->tr.Tzx.x, Tzxjmax = ran->tr.Tzx.y, Tzxkmax = ran->tr.Tzx.z;
  int idx_Tzx = i * Tzxjmax * Tzxkmax + j * Tzxkmax + k;

  if (i < Tzximax && j < Tzxjmax && k < Tzxkmax) {
    *(befTa->Tzx  + idx_Tzx) = *(aftTa->Tzx  + idx_Tzx);
    *(befTa->Tzxz + idx_Tzx) = *(aftTa->Tzxz + idx_Tzx);
    *(befTa->Tzxx + idx_Tzx) = *(aftTa->Tzxx + idx_Tzx);
  }
}

__global__ void swapVx(VelArr *aftVa, VelArr *befVa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Vximax = ran->vr.Vx.x, Vxjmax = ran->vr.Vx.y, Vxkmax = ran->vr.Vx.z;
  int idx_Vx = i * Vxjmax * Vxkmax + j * Vxkmax + k;

  if (i < Vximax && j < Vxjmax && k < Vxkmax) {
    *(befVa->Vx  + idx_Vx) = *(aftVa->Vx  + idx_Vx);
    *(befVa->Vxx + idx_Vx) = *(aftVa->Vxx + idx_Vx);
    *(befVa->Vxy + idx_Vx) = *(aftVa->Vxy + idx_Vx);
    *(befVa->Vxz + idx_Vx) = *(aftVa->Vxz + idx_Vx);
  }
}

__global__ void swapVy(VelArr *aftVa, VelArr *befVa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Vyimax = ran->vr.Vy.x, Vyjmax = ran->vr.Vy.y, Vykmax = ran->vr.Vy.z;
  int idx_Vy = i * Vyjmax * Vykmax + j * Vykmax + k;

  if (i < Vyimax && j < Vyjmax && k < Vykmax) {
    *(befVa->Vy  + idx_Vy) = *(aftVa->Vy  + idx_Vy);
    *(befVa->Vyx + idx_Vy) = *(aftVa->Vyx + idx_Vy);
    *(befVa->Vyy + idx_Vy) = *(aftVa->Vyy + idx_Vy);
    *(befVa->Vyz + idx_Vy) = *(aftVa->Vyz + idx_Vy);
  }
}

__global__ void swapVz(VelArr *aftVa, VelArr *befVa, Range *ran) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  int Vzimax = ran->vr.Vz.x, Vzjmax = ran->vr.Vz.y, Vzkmax = ran->vr.Vz.z;
  int idx_Vz = i * Vzjmax * Vzkmax + j * Vzkmax + k;

  if (i < Vzimax && j < Vzjmax && k < Vzkmax) {
    *(befVa->Vz  + idx_Vz) = *(aftVa->Vz  + idx_Vz);
    *(befVa->Vzx + idx_Vz) = *(aftVa->Vzx + idx_Vz);
    *(befVa->Vzy + idx_Vz) = *(aftVa->Vzy + idx_Vz);
    *(befVa->Vzz + idx_Vz) = *(aftVa->Vzz + idx_Vz);
  }
}

__global__ void swapBefAft(BefAft *aft, BefAft *bef, Range *ran, Coord threads) {
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 SwapTxxBlocks((ran->sr.Txx.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran->sr.Txx.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran->sr.Txx.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTyyBlocks((ran->sr.Tyy.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran->sr.Tyy.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran->sr.Tyy.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTzzBlocks((ran->sr.Tzz.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran->sr.Tzz.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran->sr.Tzz.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTxyBlocks((ran->tr.Txy.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran->tr.Txy.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran->tr.Txy.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTyzBlocks((ran->tr.Tyz.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran->tr.Tyz.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran->tr.Tyz.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 SwapTzxBlocks((ran->tr.Tzx.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran->tr.Tzx.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran->tr.Tzx.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3  SwapVxBlocks((ran->vr.Vx.x  + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran->vr.Vx.y  + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran->vr.Vx.z  + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3  SwapVyBlocks((ran->vr.Vy.x  + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran->vr.Vy.y  + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran->vr.Vy.z  + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3  SwapVzBlocks((ran->vr.Vz.x  + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran->vr.Vz.y  + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran->vr.Vz.z  + threadsPerBlock.z - 1) / threadsPerBlock.z);
  swapTxx<<<SwapTxxBlocks, threadsPerBlock>>>(&aft->sa, &bef->sa, ran);
  swapTyy<<<SwapTyyBlocks, threadsPerBlock>>>(&aft->sa, &bef->sa, ran);
  swapTzz<<<SwapTzzBlocks, threadsPerBlock>>>(&aft->sa, &bef->sa, ran);
  swapTxy<<<SwapTxyBlocks, threadsPerBlock>>>(&aft->ta, &bef->ta, ran);
  swapTyz<<<SwapTyzBlocks, threadsPerBlock>>>(&aft->ta, &bef->ta, ran);
  swapTzx<<<SwapTzxBlocks, threadsPerBlock>>>(&aft->ta, &bef->ta, ran);
  swapVx<<<SwapVxBlocks, threadsPerBlock>>>(&aft->va, &bef->va, ran);
  swapVy<<<SwapVyBlocks, threadsPerBlock>>>(&aft->va, &bef->va, ran);
  swapVz<<<SwapVzBlocks, threadsPerBlock>>>(&aft->va, &bef->va, ran);
}
