#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>

#include "../header/struct.h"

void allocate3DArray(double ****array, int x, int y, int z) {
    // 1. ホスト側のポインタをデバイス側に転送するためのメモリを確保
    double ***temp_array;
    cudaMalloc((void **)&temp_array, x * sizeof(double **));
    for (int i = 0; i < x; i++) {
        double **temp_row;
        cudaMalloc((void **)&temp_row, y * sizeof(double *));
        for (int j = 0; j < y; j++) {
            double *temp_col;
            cudaMalloc((void **)&temp_col, z * sizeof(double));
            // デバイスメモリにtemp_row[j]を設定
            cudaMemcpy(&temp_row[j], &temp_col, sizeof(double *), cudaMemcpyHostToDevice);
        }
        // デバイスメモリにtemp_array[i]を設定
        cudaMemcpy(&temp_array[i], &temp_row, sizeof(double **), cudaMemcpyHostToDevice);
    }
    // ホスト側のポインタにデバイスのポインタを設定
    *array = temp_array;
}

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

void copy3DArrayToDevice(double ***d_array, double ***h_array, int x, int y, int z) {
    // 1. ホスト側の配列をデバイス側にコピーするための準備
    double **d_row;
    double *d_col;

    // 2. 各次元ごとにデータをコピー
    for (int i = 0; i < x; i++) {
        // デバイス側の2次元ポインタを取得
        cudaMemcpy(&d_row, &d_array[i], sizeof(double **), cudaMemcpyDeviceToHost);

        for (int j = 0; j < y; j++) {
            // デバイス側の1次元ポインタを取得
            cudaMemcpy(&d_col, &d_row[j], sizeof(double *), cudaMemcpyDeviceToHost);

            // ホストからデバイスへデータコピー
            cudaMemcpy(d_col, h_array[i][j], z * sizeof(double), cudaMemcpyHostToDevice);
        }
    }
}

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

void copy3DArrayToHost(double ***h_array, double ***d_array, int x, int y, int z) {
    // 1. ホスト側の配列を受け取る準備
    double **d_row;
    double *d_col;

    // 2. 各次元ごとにデータをコピー
    for (int i = 0; i < x; i++) {
        // デバイス側の2次元ポインタを取得
        cudaMemcpy(&d_row, &d_array[i], sizeof(double **), cudaMemcpyDeviceToHost);

        for (int j = 0; j < y; j++) {
            // デバイス側の1次元ポインタを取得
            cudaMemcpy(&d_col, &d_row[j], sizeof(double *), cudaMemcpyDeviceToHost);

            // デバイスからホストへデータコピー
            cudaMemcpy(h_array[i][j], d_col, z * sizeof(double), cudaMemcpyDeviceToHost);
        }
    }
}

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