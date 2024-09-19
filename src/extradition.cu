#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>

#include "../header/struct.h"

////// メモリ確保

// 3d array メモリ確保
// void allocate3DArray(double ****array, int x, int y, int z) {
//     // 1. ホスト側のポインタをデバイス側に転送するためのメモリを確保
//     double ***temp_array;
//     double **temp_row;
//     double *temp_col;
//     cudaMalloc((void **)&temp_array, x * sizeof(double **));
//     for (int i = 0; i < x; i++) {
//         cudaMalloc((void **)&temp_row, y * sizeof(double *));
//         for (int j = 0; j < y; j++) {
//             cudaMalloc((void **)&temp_col, z * sizeof(double));
//             // デバイスメモリにtemp_row[j]を設定
//             cudaMemcpy(&temp_row[j], &temp_col, sizeof(double *), cudaMemcpyHostToDevice);
//         }
//         // デバイスメモリにtemp_array[i]を設定
//         cudaMemcpy(&temp_array[i], &temp_row, sizeof(double **), cudaMemcpyHostToDevice);
//     }
//     // ホスト側のポインタにデバイスのポインタを設定
//     *array = temp_array;
// }

void allocate3DArray(double ****array, int x, int y, int z) {
    double ***temp_array;
    double **temp_row;
    double *temp_col;

    // 1. 連続したデバイスメモリを一度に確保（3次元分）
    cudaError_t err = cudaMalloc((void **)&temp_col, x * y * z * sizeof(double));

    // 2. y方向のポインタ配列を連続して確保
    err = cudaMalloc((void **)&temp_row, x * y * sizeof(double*));

    // 3. x方向のポインタ配列を連続して確保
    err = cudaMalloc((void **)&temp_array, x * sizeof(double**));

    // 4. 各ポインタを設定
    for (int i = 0; i < x; i++) {
        // temp_rowのポインタをi番目の行に設定
        err = cudaMemcpy(&temp_array[i], &temp_row[i * y], sizeof(double*), cudaMemcpyHostToDevice);
        for (int j = 0; j < y; j++) {
            // temp_colのポインタをi番目、j番目に設定
            err = cudaMemcpy(&temp_row[i * y + j], &temp_col[(i * y + j) * z], sizeof(double), cudaMemcpyHostToDevice);
        }
    }
    if (err != cudaSuccess) {
        printf("cudaMalloc failed: %s\n", cudaGetErrorString(err));
        return;
    }
    // 5. ホスト側のポインタにデバイスのポインタを設定
    *array = temp_array;

}

// メモリ確保 BefAft
void allocateBefAft(BefAft *d_befAft, Range ran) {
    // メモリをGPUに割り当て
    cudaMalloc((void **)&(d_befAft->sa), sizeof(SigArr));
    cudaMalloc((void **)&(d_befAft->ta), sizeof(TauArr));
    cudaMalloc((void **)&(d_befAft->va), sizeof(VelArr));
    
    // SigArr のメンバのメモリ確保
    allocate3DArray(&d_befAft->sa.Txx , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_befAft->sa.Txxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_befAft->sa.Txxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_befAft->sa.Txxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_befAft->sa.Tyy , ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DArray(&d_befAft->sa.Tyyx, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DArray(&d_befAft->sa.Tyyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DArray(&d_befAft->sa.Tyyz, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DArray(&d_befAft->sa.Tzz , ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DArray(&d_befAft->sa.Tzzx, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DArray(&d_befAft->sa.Tzzy, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DArray(&d_befAft->sa.Tzzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    // TauArr のメンバのメモリ確保
    allocate3DArray(&d_befAft->ta.Txy , ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DArray(&d_befAft->ta.Txyx, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DArray(&d_befAft->ta.Txyy, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DArray(&d_befAft->ta.Tyz , ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DArray(&d_befAft->ta.Tyzy, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DArray(&d_befAft->ta.Tyzz, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DArray(&d_befAft->ta.Tzx , ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    allocate3DArray(&d_befAft->ta.Tzxz, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    allocate3DArray(&d_befAft->ta.Tzxx, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);

    // VelArr のメンバのメモリ確保
    allocate3DArray(&d_befAft->va.Vx , ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DArray(&d_befAft->va.Vxx, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DArray(&d_befAft->va.Vxy, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DArray(&d_befAft->va.Vxz, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DArray(&d_befAft->va.Vy , ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArray(&d_befAft->va.Vyx, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArray(&d_befAft->va.Vyy, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArray(&d_befAft->va.Vyz, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArray(&d_befAft->va.Vz , ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArray(&d_befAft->va.Vzx, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    allocate3DArray(&d_befAft->va.Vzy, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    allocate3DArray(&d_befAft->va.Vzz, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
}

// メモリ確保 MedArr
void allocateMedArr(MedArr *d_medArr, Range ran) {
    // 各メンバ変数のメモリをGPUに割り当て
    allocate3DArray(&d_medArr->ramda, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->mu, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->c11, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->rho, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    // ゼータ（zetaxx, zetaxy, など）のメモリをGPUに割り当て
    allocate3DArray(&d_medArr->zetaxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetaxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetaxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetayx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetayy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetayz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetazx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetazy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetazz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    // その他のメンバ変数のメモリをGPUに割り当て
    allocate3DArray(&d_medArr->gamma, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->khi, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->xi11, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetadx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetady, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_medArr->zetadz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
}

// メモリ確保 Impulse
void allocateImpulse(Impulse *d_impulse, Range ran) {
    cudaMalloc((void **)&(d_impulse->in), sizeof(Coord));
    allocate3DArray(&d_impulse->Txx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_impulse->Tyy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArray(&d_impulse->Tzz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
}


////// データ転送 ホスト->デバイス

// ホスト->デバイス 3d array
void copy3DArrayToDevice(double ***d_array, double ***h_array, int x, int y, int z) {
    cudaError_t err;

    // 1. 各次元ごとのメモリをデバイス側に割り当て
    for (int i = 0; i < x; i++) {
        double **d_row;
        double *d_col;

        // デバイス側で2次元目のポインタ用のメモリを確保
        err = cudaMalloc((void**)&d_row, y * sizeof(double *));
        if (err != cudaSuccess) {
            printf("Error allocating memory for d_row: %s\n", cudaGetErrorString(err));
            return;
        }

        for (int j = 0; j < y; j++) {
            // デバイス側で1次元目の配列（データ本体）用のメモリを確保
            err = cudaMalloc((void**)&d_col, z * sizeof(double));
            if (err != cudaSuccess) {
                printf("Error allocating memory for d_col: %s\n", cudaGetErrorString(err));
                return;
            }

            // ホストからデバイスへデータをコピー
            err = cudaMemcpy(d_col, h_array[i][j], z * sizeof(double), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) {
                printf("Error copying data to device: %s\n", cudaGetErrorString(err));
                return;
            }

            // 2次元目のポインタ配列に1次元目のポインタをセット
            err = cudaMemcpy(&d_row[j], &d_col, sizeof(double *), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) {
                printf("Error copying pointer to device1: %s\n", cudaGetErrorString(err));
                return;
            }
        }
        
        // デバイス側の3次元目のポインタ配列に2次元目のポインタをセット
        // err = cudaMemcpy(&d_array[i], &d_row, sizeof(double **), cudaMemcpyHostToDevice);
        err = cudaMemcpy(d_array + i, &d_row, sizeof(double **), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            printf("Error copying pointer to device2: %s\n", cudaGetErrorString(err));
            return;
        }
    }
}
// ホスト->デバイス BefAft
void copyBefAftToDevice(BefAft *d_befAft, BefAft *h_befAft, Range ran) {
    int Txxi = ran.sr.Txx.x, Txxj = ran.sr.Txx.y, Txxk = ran.sr.Txx.z;
    int Tyyi = ran.sr.Tyy.x, Tyyj = ran.sr.Tyy.y, Tyyk = ran.sr.Tyy.z;
    int Tzzi = ran.sr.Tzz.x, Tzzj = ran.sr.Tzz.y, Tzzk = ran.sr.Tzz.z;
    int Txyi = ran.tr.Txy.x, Txyj = ran.tr.Txy.y, Txyk = ran.tr.Txy.z;
    int Tyzi = ran.tr.Tyz.x, Tyzj = ran.tr.Tyz.y, Tyzk = ran.tr.Tyz.z;
    int Tzxi = ran.tr.Tzx.x, Tzxj = ran.tr.Tzx.y, Tzxk = ran.tr.Tzx.z;
    int Vxi  = ran.vr.Vx.x , Vxj  = ran.vr.Vx.y , Vxk  = ran.vr.Vx.z;
    int Vyi  = ran.vr.Vy.x , Vyj  = ran.vr.Vy.y , Vyk  = ran.vr.Vy.z;
    int Vzi  = ran.vr.Vz.x , Vzj  = ran.vr.Vz.y , Vzk  = ran.vr.Vz.z;
    // SigArr のメンバのデータコピー
    copy3DArrayToDevice(d_befAft->sa.Txx,  h_befAft->sa.Txx , Txxi, Txxj, Txxk);
    copy3DArrayToDevice(d_befAft->sa.Txxx, h_befAft->sa.Txxx, Txxi, Txxj, Txxk);
    copy3DArrayToDevice(d_befAft->sa.Txxy, h_befAft->sa.Txxy, Txxi, Txxj, Txxk);
    copy3DArrayToDevice(d_befAft->sa.Txxz, h_befAft->sa.Txxz, Txxi, Txxj, Txxk);
    copy3DArrayToDevice(d_befAft->sa.Tyy, h_befAft->sa.Tyy , Tyyi, Tyyj, Tyyk);
    copy3DArrayToDevice(d_befAft->sa.Tyyx, h_befAft->sa.Tyyx, Tyyi, Tyyj, Tyyk);
    copy3DArrayToDevice(d_befAft->sa.Tyyy, h_befAft->sa.Tyyy, Tyyi, Tyyj, Tyyk);
    copy3DArrayToDevice(d_befAft->sa.Tyyz, h_befAft->sa.Tyyz, Tyyi, Tyyj, Tyyk);
    copy3DArrayToDevice(d_befAft->sa.Tzz, h_befAft->sa.Tzz , Tzzi, Tzzj, Tzzk);
    copy3DArrayToDevice(d_befAft->sa.Tzzx, h_befAft->sa.Tzzx, Tzzi, Tzzj, Tzzk);
    copy3DArrayToDevice(d_befAft->sa.Tzzy, h_befAft->sa.Tzzy, Tzzi, Tzzj, Tzzk);
    copy3DArrayToDevice(d_befAft->sa.Tzzz, h_befAft->sa.Tzzz, Tzzi, Tzzj, Tzzk);
    // TauArr のメンバのデータコピー
    copy3DArrayToDevice(d_befAft->ta.Txy, h_befAft->ta.Txy , Txyi, Txyj, Txyk);
    copy3DArrayToDevice(d_befAft->ta.Txyx, h_befAft->ta.Txyx, Txyi, Txyj, Txyk);
    copy3DArrayToDevice(d_befAft->ta.Txyy, h_befAft->ta.Txyy, Txyi, Txyj, Txyk);
    copy3DArrayToDevice(d_befAft->ta.Tyz, h_befAft->ta.Tyz , Tyzi, Tyzj, Tyzk);
    copy3DArrayToDevice(d_befAft->ta.Tyzy, h_befAft->ta.Tyzy, Tyzi, Tyzj, Tyzk);
    copy3DArrayToDevice(d_befAft->ta.Tyzz, h_befAft->ta.Tyzz, Tyzi, Tyzj, Tyzk);
    copy3DArrayToDevice(d_befAft->ta.Tzx, h_befAft->ta.Tzx , Tzxi, Tzxj, Tzxk);
    copy3DArrayToDevice(d_befAft->ta.Tzxz, h_befAft->ta.Tzxz, Tzxi, Tzxj, Tzxk);
    copy3DArrayToDevice(d_befAft->ta.Tzxx, h_befAft->ta.Tzxx, Tzxi, Tzxj, Tzxk);

    // VelArr のメンバのデータコピー
    copy3DArrayToDevice(d_befAft->va.Vx, h_befAft->va.Vx , Vxi, Vxj, Vxk);
    copy3DArrayToDevice(d_befAft->va.Vxx, h_befAft->va.Vxx, Vxi, Vxj, Vxk);
    copy3DArrayToDevice(d_befAft->va.Vxy, h_befAft->va.Vxy, Vxi, Vxj, Vxk);
    copy3DArrayToDevice(d_befAft->va.Vxz, h_befAft->va.Vxz, Vxi, Vxj, Vxk);
    copy3DArrayToDevice(d_befAft->va.Vy, h_befAft->va.Vy , Vyi, Vyj, Vyk);
    copy3DArrayToDevice(d_befAft->va.Vyx, h_befAft->va.Vyx, Vyi, Vyj, Vyk);
    copy3DArrayToDevice(d_befAft->va.Vyy, h_befAft->va.Vyy, Vyi, Vyj, Vyk);
    copy3DArrayToDevice(d_befAft->va.Vyz, h_befAft->va.Vyz, Vyi, Vyj, Vyk);
    copy3DArrayToDevice(d_befAft->va.Vz, h_befAft->va.Vz , Vzi, Vzj, Vzk);
    copy3DArrayToDevice(d_befAft->va.Vzx, h_befAft->va.Vzx, Vzi, Vzj, Vzk);
    copy3DArrayToDevice(d_befAft->va.Vzy, h_befAft->va.Vzy, Vzi, Vzj, Vzk);
    copy3DArrayToDevice(d_befAft->va.Vzz, h_befAft->va.Vzz, Vzi, Vzj, Vzk);
}
// ホスト->デバイス Impulse
void copyImpulseToDevice(Impulse *d_impulse, Impulse *h_impulse, Range ran) {
    int Txxi = ran.sr.Txx.x, Txxj = ran.sr.Txx.y, Txxk = ran.sr.Txx.z;
    int Tyyi = ran.sr.Tyy.x, Tyyj = ran.sr.Tyy.y, Tyyk = ran.sr.Tyy.z;
    int Tzzi = ran.sr.Tzz.x, Tzzj = ran.sr.Tzz.y, Tzzk = ran.sr.Tzz.z;
    // 3D配列データの転送 (Txx, Tyy, Tzz)
    copy3DArrayToDevice(d_impulse->Txx, h_impulse->Txx, Txxi, Txxj, Txxk);
    copy3DArrayToDevice(d_impulse->Tyy, h_impulse->Tyy, Tyyi, Tyyj, Tyyk);
    copy3DArrayToDevice(d_impulse->Tzz, h_impulse->Tzz, Tzzi, Tzzj, Tzzk);

    cudaMemcpy(&(d_impulse->freq), &(h_impulse->freq), sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_impulse->mode), &(h_impulse->mode), sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_impulse->in), &(h_impulse->in), sizeof(Coord), cudaMemcpyHostToDevice);
    printf("host:%f\n",h_impulse->Tzz[h_impulse->in.x][h_impulse->in.y][h_impulse->in.z]);
}
// ホスト->デバイス MedArr
void copyMedArrToDevice(MedArr *d_medArr, MedArr *h_medArr, Range ran) {
    int xi = ran.sr.Txx.x, xj = ran.sr.Txx.y, xk = ran.sr.Txx.z;

    copy3DArrayToDevice(d_medArr->ramda, h_medArr->ramda, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->mu, h_medArr->mu, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->c11, h_medArr->c11, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->rho, h_medArr->rho, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetaxx, h_medArr->zetaxx, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetaxy, h_medArr->zetaxy, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetaxz, h_medArr->zetaxz, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetayx, h_medArr->zetayx, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetayy, h_medArr->zetayy, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetayz, h_medArr->zetayz, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetazx, h_medArr->zetazx, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetazy, h_medArr->zetazy, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetazz, h_medArr->zetazz, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->gamma, h_medArr->gamma, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->khi, h_medArr->khi, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->xi11, h_medArr->xi11, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetadx, h_medArr->zetadx, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetady, h_medArr->zetady, xi, xj, xk);
    copy3DArrayToDevice(d_medArr->zetadz, h_medArr->zetadz, xi, xj, xk);
}
// ホスト->デバイス Diff
void copyDiffToDevice(Diff *d_diff, Diff *h_diff) {
    cudaMemcpy(&(d_diff->dx), &(h_diff->dx), sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_diff->dy), &(h_diff->dy), sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_diff->dz), &(h_diff->dz), sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_diff->dt), &(h_diff->dt), sizeof(double), cudaMemcpyHostToDevice);
}

////// データ転送 デバイス->ホスト

// 3d array デバイス->ホスト
void copy3DArrayToHost(double ***h_array, double ***d_array, int x, int y, int z) {
    // 1. ホスト側の配列を受け取る準備
    double **d_row;
    double *d_col;

    // 2. 各次元ごとにデータをコピー
    for (int i = 0; i < x; i++) {
        // デバイス側の2次元ポインタを取得
        cudaError_t err = cudaMemcpy(&d_row, d_array + i, sizeof(double **), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            printf("Error copying pointer to host (d_row): %s\n", cudaGetErrorString(err));
            return;
        }

        for (int j = 0; j < y; j++) {
            // デバイス側の1次元ポインタを取得
            err = cudaMemcpy(&d_col, d_row + j, sizeof(double *), cudaMemcpyDeviceToHost);
            if (err != cudaSuccess) {
                printf("Error copying pointer to host (d_col): %s\n", cudaGetErrorString(err));
                return;
            }

            // デバイスからホストへデータコピー
            err = cudaMemcpy(h_array[i][j], d_col, z * sizeof(double), cudaMemcpyDeviceToHost);
            if (err != cudaSuccess) {
                printf("Error copying data to host: %s\n", cudaGetErrorString(err));
                return;
            }
        }
    }
}

// デバイス->ホスト BefAft
void copyBefAftToHost(BefAft *h_befAft, BefAft *d_befAft, Range ran) {
    int Txxi = ran.sr.Txx.x, Txxj = ran.sr.Txx.y, Txxk = ran.sr.Txx.z;
    int Tyyi = ran.sr.Tyy.x, Tyyj = ran.sr.Tyy.y, Tyyk = ran.sr.Tyy.z;
    int Tzzi = ran.sr.Tzz.x, Tzzj = ran.sr.Tzz.y, Tzzk = ran.sr.Tzz.z;
    int Txyi = ran.tr.Txy.x, Txyj = ran.tr.Txy.y, Txyk = ran.tr.Txy.z;
    int Tyzi = ran.tr.Tyz.x, Tyzj = ran.tr.Tyz.y, Tyzk = ran.tr.Tyz.z;
    int Tzxi = ran.tr.Tzx.x, Tzxj = ran.tr.Tzx.y, Tzxk = ran.tr.Tzx.z;
    int Vxi  = ran.vr.Vx.x , Vxj  = ran.vr.Vx.y , Vxk  = ran.vr.Vx.z;
    int Vyi  = ran.vr.Vy.x , Vyj  = ran.vr.Vy.y , Vyk  = ran.vr.Vy.z;
    int Vzi  = ran.vr.Vz.x , Vzj  = ran.vr.Vz.y , Vzk  = ran.vr.Vz.z;
    // SigArr のメンバのデータコピー
    copy3DArrayToHost(h_befAft->sa.Txx,  d_befAft->sa.Txx , Txxi, Txxj, Txxk);
    copy3DArrayToHost(h_befAft->sa.Txxx, d_befAft->sa.Txxx, Txxi, Txxj, Txxk);
    copy3DArrayToHost(h_befAft->sa.Txxy, d_befAft->sa.Txxy, Txxi, Txxj, Txxk);
    copy3DArrayToHost(h_befAft->sa.Txxz, d_befAft->sa.Txxz, Txxi, Txxj, Txxk);
    copy3DArrayToHost(h_befAft->sa.Tyy, d_befAft->sa.Tyy , Tyyi, Tyyj, Tyyk);
    copy3DArrayToHost(h_befAft->sa.Tyyx, d_befAft->sa.Tyyx, Tyyi, Tyyj, Tyyk);
    copy3DArrayToHost(h_befAft->sa.Tyyy, d_befAft->sa.Tyyy, Tyyi, Tyyj, Tyyk);
    copy3DArrayToHost(h_befAft->sa.Tyyz, d_befAft->sa.Tyyz, Tyyi, Tyyj, Tyyk);
    copy3DArrayToHost(h_befAft->sa.Tzz, d_befAft->sa.Tzz , Tzzi, Tzzj, Tzzk);
    copy3DArrayToHost(h_befAft->sa.Tzzx, d_befAft->sa.Tzzx, Tzzi, Tzzj, Tzzk);
    copy3DArrayToHost(h_befAft->sa.Tzzy, d_befAft->sa.Tzzy, Tzzi, Tzzj, Tzzk);
    copy3DArrayToHost(h_befAft->sa.Tzzz, d_befAft->sa.Tzzz, Tzzi, Tzzj, Tzzk);
    // TauArr のメンバのデータコピー
    copy3DArrayToHost(h_befAft->ta.Txy, d_befAft->ta.Txy , Txyi, Txyj, Txyk);
    copy3DArrayToHost(h_befAft->ta.Txyx, d_befAft->ta.Txyx, Txyi, Txyj, Txyk);
    copy3DArrayToHost(h_befAft->ta.Txyy, d_befAft->ta.Txyy, Txyi, Txyj, Txyk);
    copy3DArrayToHost(h_befAft->ta.Tyz, d_befAft->ta.Tyz , Tyzi, Tyzj, Tyzk);
    copy3DArrayToHost(h_befAft->ta.Tyzy, d_befAft->ta.Tyzy, Tyzi, Tyzj, Tyzk);
    copy3DArrayToHost(h_befAft->ta.Tyzz, d_befAft->ta.Tyzz, Tyzi, Tyzj, Tyzk);
    copy3DArrayToHost(h_befAft->ta.Tzx, d_befAft->ta.Tzx , Tzxi, Tzxj, Tzxk);
    copy3DArrayToHost(h_befAft->ta.Tzxz, d_befAft->ta.Tzxz, Tzxi, Tzxj, Tzxk);
    copy3DArrayToHost(h_befAft->ta.Tzxx, d_befAft->ta.Tzxx, Tzxi, Tzxj, Tzxk);

    // VelArr のメンバのデータコピー
    copy3DArrayToHost(h_befAft->va.Vx, d_befAft->va.Vx , Vxi, Vxj, Vxk);
    copy3DArrayToHost(h_befAft->va.Vxx, d_befAft->va.Vxx, Vxi, Vxj, Vxk);
    copy3DArrayToHost(h_befAft->va.Vxy, d_befAft->va.Vxy, Vxi, Vxj, Vxk);
    copy3DArrayToHost(h_befAft->va.Vxz, d_befAft->va.Vxz, Vxi, Vxj, Vxk);
    copy3DArrayToHost(h_befAft->va.Vy, d_befAft->va.Vy , Vyi, Vyj, Vyk);
    copy3DArrayToHost(h_befAft->va.Vyx, d_befAft->va.Vyx, Vyi, Vyj, Vyk);
    copy3DArrayToHost(h_befAft->va.Vyy, d_befAft->va.Vyy, Vyi, Vyj, Vyk);
    copy3DArrayToHost(h_befAft->va.Vyz, d_befAft->va.Vyz, Vyi, Vyj, Vyk);
    copy3DArrayToHost(h_befAft->va.Vz, d_befAft->va.Vz , Vzi, Vzj, Vzk);
    copy3DArrayToHost(h_befAft->va.Vzx, d_befAft->va.Vzx, Vzi, Vzj, Vzk);
    copy3DArrayToHost(h_befAft->va.Vzy, d_befAft->va.Vzy, Vzi, Vzj, Vzk);
    copy3DArrayToHost(h_befAft->va.Vzz, d_befAft->va.Vzz, Vzi, Vzj, Vzk);
}

// hostメモリ確保

void allocate3DArrayHost(double ****array, int x, int y, int z) {
    double ***temp_array;
    double **temp_row;
    double *temp_col;

    // 1. 連続したメモリ領域を確保（3次元分）
    temp_col = (double *)malloc(x * y * z * sizeof(double));
    if (temp_col == NULL) {
        printf("Memory allocation for data failed.\n");
        return;
    }

    // 2. y方向のポインタ配列を連続して確保
    temp_row = (double **)malloc(x * y * sizeof(double*));
    if (temp_row == NULL) {
        printf("Memory allocation for row pointers failed.\n");
        free(temp_col);  // 確保済みのメモリを解放
        return;
    }

    // 3. x方向のポインタ配列を連続して確保
    temp_array = (double ***)malloc(x * sizeof(double**));
    if (temp_array == NULL) {
        printf("Memory allocation for array pointers failed.\n");
        free(temp_row);  // 確保済みのメモリを解放
        free(temp_col);  // 確保済みのメモリを解放
        return;
    }

    // 4. 各ポインタを設定
    for (int i = 0; i < x; i++) {
        temp_array[i] = &temp_row[i * y];  // temp_rowのポインタを設定
        for (int j = 0; j < y; j++) {
            temp_row[i * y + j] = &temp_col[(i * y + j) * z];  // temp_colのポインタを設定
        }
    }

    // 5. ホスト側のポインタに確保したメモリを設定
    *array = temp_array;
}

// メモリ確保 BefAft
void allocateBefAftHost(BefAft *h_befAft, Range ran) {
    // メモリをGPUに割り当て
    cudaMalloc((void **)&(h_befAft->sa), sizeof(SigArr));
    cudaMalloc((void **)&(h_befAft->sa), sizeof(SigArr));
    cudaMalloc((void **)&(h_befAft->ta), sizeof(TauArr));
    cudaMalloc((void **)&(h_befAft->va), sizeof(VelArr));
    
    // SigArr のメンバのメモリ確保
    allocate3DArrayHost(&h_befAft->sa.Txx , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&h_befAft->sa.Txxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&h_befAft->sa.Txxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&h_befAft->sa.Txxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&h_befAft->sa.Tyy , ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DArrayHost(&h_befAft->sa.Tyyx, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DArrayHost(&h_befAft->sa.Tyyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DArrayHost(&h_befAft->sa.Tyyz, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DArrayHost(&h_befAft->sa.Tzz , ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DArrayHost(&h_befAft->sa.Tzzx, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DArrayHost(&h_befAft->sa.Tzzy, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DArrayHost(&h_befAft->sa.Tzzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    // TauArr のメンバのメモリ確保
    allocate3DArrayHost(&h_befAft->ta.Txy , ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DArrayHost(&h_befAft->ta.Txyx, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DArrayHost(&h_befAft->ta.Txyy, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DArrayHost(&h_befAft->ta.Tyz , ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DArrayHost(&h_befAft->ta.Tyzy, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DArrayHost(&h_befAft->ta.Tyzz, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DArrayHost(&h_befAft->ta.Tzx , ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    allocate3DArrayHost(&h_befAft->ta.Tzxz, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    allocate3DArrayHost(&h_befAft->ta.Tzxx, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);

    // VelArr のメンバのメモリ確保
    allocate3DArrayHost(&h_befAft->va.Vx , ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DArrayHost(&h_befAft->va.Vxx, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DArrayHost(&h_befAft->va.Vxy, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DArrayHost(&h_befAft->va.Vxz, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DArrayHost(&h_befAft->va.Vy , ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArrayHost(&h_befAft->va.Vyx, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArrayHost(&h_befAft->va.Vyy, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArrayHost(&h_befAft->va.Vyz, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArrayHost(&h_befAft->va.Vz , ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DArrayHost(&h_befAft->va.Vzx, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    allocate3DArrayHost(&h_befAft->va.Vzy, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    allocate3DArrayHost(&h_befAft->va.Vzz, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
}

// メモリ確保 MedArr
void allocateMedArrHost(MedArr *d_medArr, Range ran) {
    // 各メンバ変数のメモリをGPUに割り当て
    allocate3DArrayHost(&d_medArr->ramda, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->mu, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->c11, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->rho, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    // ゼータ（zetaxx, zetaxy, など）のメモリをGPUに割り当て
    allocate3DArrayHost(&d_medArr->zetaxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetaxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetaxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetayx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetayy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetayz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetazx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetazy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetazz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    // その他のメンバ変数のメモリをGPUに割り当て
    allocate3DArrayHost(&d_medArr->gamma, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->khi, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->xi11, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetadx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetady, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_medArr->zetadz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
}

// メモリ確保 Impulse
void allocateImpulseHost(Impulse *d_impulse, Range ran) {
    cudaMalloc((void **)&(d_impulse->in), sizeof(Coord));
    allocate3DArrayHost(&d_impulse->Txx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_impulse->Tyy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DArrayHost(&d_impulse->Tzz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
}
