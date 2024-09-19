#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>

#include "../header/struct.h"
#include "../header/memory.h"

// hostメモリ確保

void allocate3DarrayHost(double ****array, int x, int y ,int z) {
    *array = (double ***)malloc(x * sizeof(double **));
    for (int i = 0; i < x; i++) {
        (*array)[i] = (double **)malloc(y * sizeof(double *));
        for (int j = 0; j < y; j++) {
            (*array)[i][j] = (double *)malloc(z * sizeof(double));
        }
    }
}

void allocateHostBefAft(BefAft **befaft_host_ptr, Range ran) {
    BefAft *befaft_host = (BefAft *)malloc(sizeof(BefAft));
    // SigArr
    printf("no\n");
    allocate3DarrayHost(&befaft_host->sa.Txx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&befaft_host->sa.Tyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DarrayHost(&befaft_host->sa.Tzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    allocate3DarrayHost(&befaft_host->sa.Txxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&befaft_host->sa.Txxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&befaft_host->sa.Txxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    allocate3DarrayHost(&befaft_host->sa.Tyyx, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DarrayHost(&befaft_host->sa.Tyyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DarrayHost(&befaft_host->sa.Tyyz, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);

    allocate3DarrayHost(&befaft_host->sa.Tzzx, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DarrayHost(&befaft_host->sa.Tzzy, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DarrayHost(&befaft_host->sa.Tzzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    // TauArr
    allocate3DarrayHost(&befaft_host->ta.Txy, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DarrayHost(&befaft_host->ta.Tyz, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DarrayHost(&befaft_host->ta.Tzx, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);

    allocate3DarrayHost(&befaft_host->ta.Txyx, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DarrayHost(&befaft_host->ta.Txyy, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);

    allocate3DarrayHost(&befaft_host->ta.Tyzy, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DarrayHost(&befaft_host->ta.Tyzz, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);

    allocate3DarrayHost(&befaft_host->ta.Tzxx, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    allocate3DarrayHost(&befaft_host->ta.Tzxz, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);

    // VelArr
    allocate3DarrayHost(&befaft_host->va.Vx, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DarrayHost(&befaft_host->va.Vy, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DarrayHost(&befaft_host->va.Vz, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);

    allocate3DarrayHost(&befaft_host->va.Vxx, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DarrayHost(&befaft_host->va.Vxy, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DarrayHost(&befaft_host->va.Vxz, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);

    allocate3DarrayHost(&befaft_host->va.Vyx, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DarrayHost(&befaft_host->va.Vyy, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DarrayHost(&befaft_host->va.Vyz, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);

    allocate3DarrayHost(&befaft_host->va.Vzx, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    allocate3DarrayHost(&befaft_host->va.Vzy, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    allocate3DarrayHost(&befaft_host->va.Vzz, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);

    *befaft_host_ptr = befaft_host;
}

void allocateHostMedArr(MedArr **medarr_host_ptr, Range ran) {
    MedArr *medarr_host = (MedArr *)malloc(sizeof(MedArr));
    allocate3DarrayHost(&medarr_host->ramda, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->mu, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->c11, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->rho, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    
    allocate3DarrayHost(&medarr_host->zetaxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetaxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetaxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetayx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetayy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetayz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetazx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetazy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetazz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    allocate3DarrayHost(&medarr_host->gamma, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->khi, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->xi11, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetadx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetady, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&medarr_host->zetadz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    *medarr_host_ptr = medarr_host;
}
// ポインタ以外の部分確保できているかわからない
void allocateHostImpulse(Impulse **impulse_host_ptr, Range ran) {
    Impulse *impulse_host = (Impulse *)malloc(sizeof(Impulse));
    allocate3DarrayHost(&impulse_host->Txx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayHost(&impulse_host->Tyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DarrayHost(&impulse_host->Tzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    *impulse_host_ptr = impulse_host;
}

// deviceメモリ確保

// void allocate3DarrayDevice(&double ***array, int x, int y, int z) {
//     double ***d_level1;
//     double  **d_level2;
//     double   *d_level3;

//     cudaError_t err;

//     err = cudaMalloc((void ***)&d_level1, x * sizeof(double **));
//     if(err != cudaSuccess) {
//         printf("Malloc error(Device Level1): %s\n",cudaGetErrorString(err));
//         return;
//     }
//     for(int i = 0; i < x; i++) {
//         err = cudaMalloc((void ***)&d_level2, y * sizeof(double *));
//         if(err != cudaSuccess) {
//             printf("Malloc error(Device Level2): %s\n",cudaGetErrorString(err));
//             return;
//         }
//         for(int j = 0; j < y; j++) {
//             err = cudaMalloc((void **)&d_level3, z * sizeof(double));
//             if(err != cudaSuccess) {
//                 printf("Malloc error(Device Level3): %s\n",cudaGetErrorString(err));
//                 return;
//             }
//             err = cudaMemcpy(&d_level2[j], &d_level3, sizeof(double *), cudaMemcpyHostToDevice);
//             if(err != cudaSuccess) {
//                 printf("Malloc error(Device Level3 to Level2): %s\n",cudaGetErrorString(err));
//                 return;
//             }
//         }
//         err = cudaMemcpy(&d_level1[i], &d_level2, sizeof(double **), cudaMemcpyHostToDevice);
//         if(err != cudaSuccess) {
//             printf("Memcpy error(Device Level2 to Level1): %s\n",cudaGetErrorString(err));
//             return;
//         }
//     }
//     err = cudaMemcpy(&array, &d_level1, sizeof(double ***), cudaMemcpyHostToDevice);
//     if(err != cudaSuccess) {
//         printf("Memcpy error(Device Level1 to Array): %s\n",cudaGetErrorString(err));
//         return;
//     }
// }

void allocate3DarrayDevice(double ****array, int x, int y, int z) {
    double *data;
    cudaError_t err;

    // データ部分のメモリをデバイス上に確保
    err = cudaMalloc((void **)&data, x * y * z * sizeof(double));
    if(err != cudaSuccess) {
        printf("Malloc error(Device Data): %s\n", cudaGetErrorString(err));
        return;
    }

    // ポインタ配列をデバイス上に確保
    double **d_level2;
    err = cudaMalloc((void **)&d_level2, x * y * sizeof(double *));
    if(err != cudaSuccess) {
        printf("Malloc error(Device Level2): %s\n", cudaGetErrorString(err));
        return;
    }

    double ***d_level1;
    err = cudaMalloc((void **)&d_level1, x * sizeof(double **));
    if(err != cudaSuccess) {
        printf("Malloc error(Device Level1): %s\n", cudaGetErrorString(err));
        return;
    }

    // ホスト側でポインタ配列を作成
    double **h_level2 = (double **)malloc(x * y * sizeof(double *));
    double ***h_level1 = (double ***)malloc(x * sizeof(double **));

    for (int i = 0; i < x; i++) {
        h_level1[i] = d_level2 + i * y; // デバイス上のアドレス
        for (int j = 0; j < y; j++) {
            h_level2[i * y + j] = data + (i * y * z) + (j * z); // デバイス上のアドレス
        }
    }

    // ポインタ配列をデバイスにコピー
    err = cudaMemcpy(d_level2, h_level2, x * y * sizeof(double *), cudaMemcpyHostToDevice);
    if(err != cudaSuccess) {
        printf("Memcpy error(Level2): %s\n", cudaGetErrorString(err));
        return;
    }

    err = cudaMemcpy(d_level1, h_level1, x * sizeof(double **), cudaMemcpyHostToDevice);
    if(err != cudaSuccess) {
        printf("Memcpy error(Level1): %s\n", cudaGetErrorString(err));
        return;
    }

    // 呼び出し元にデバイス上のトップレベルのポインタを返す
    *array = d_level1;

    // ホスト側のメモリを解放
    free(h_level2);
    free(h_level1);
}



void allocateDeviceBefAft(BefAft **befaft_device_ptr, Range ran) {
    BefAft *befaft_device = (BefAft*)malloc(sizeof(BefAft));
    // cudaMalloc((void **)befaft_device, sizeof(BefAft));

    // SigArr
    allocate3DarrayDevice(&befaft_device->sa.Txx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&befaft_device->sa.Tyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DarrayDevice(&befaft_device->sa.Tzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    allocate3DarrayDevice(&befaft_device->sa.Txxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&befaft_device->sa.Txxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&befaft_device->sa.Txxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    allocate3DarrayDevice(&befaft_device->sa.Tyyx, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DarrayDevice(&befaft_device->sa.Tyyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DarrayDevice(&befaft_device->sa.Tyyz, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);

    allocate3DarrayDevice(&befaft_device->sa.Tzzx, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DarrayDevice(&befaft_device->sa.Tzzy, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    allocate3DarrayDevice(&befaft_device->sa.Tzzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    // TauArr
    allocate3DarrayDevice(&befaft_device->ta.Txy, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DarrayDevice(&befaft_device->ta.Tyz, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DarrayDevice(&befaft_device->ta.Tzx, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);

    allocate3DarrayDevice(&befaft_device->ta.Txyx, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    allocate3DarrayDevice(&befaft_device->ta.Txyy, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);

    allocate3DarrayDevice(&befaft_device->ta.Tyzy, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    allocate3DarrayDevice(&befaft_device->ta.Tyzz, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);

    allocate3DarrayDevice(&befaft_device->ta.Tzxx, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    allocate3DarrayDevice(&befaft_device->ta.Tzxz, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);

    // VelArr
    allocate3DarrayDevice(&befaft_device->va.Vx, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DarrayDevice(&befaft_device->va.Vy, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DarrayDevice(&befaft_device->va.Vz, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);

    allocate3DarrayDevice(&befaft_device->va.Vxx, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DarrayDevice(&befaft_device->va.Vxy, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    allocate3DarrayDevice(&befaft_device->va.Vxz, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);

    allocate3DarrayDevice(&befaft_device->va.Vyx, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DarrayDevice(&befaft_device->va.Vyy, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    allocate3DarrayDevice(&befaft_device->va.Vyz, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);

    allocate3DarrayDevice(&befaft_device->va.Vzx, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    allocate3DarrayDevice(&befaft_device->va.Vzy, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    allocate3DarrayDevice(&befaft_device->va.Vzz, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);

    *befaft_device_ptr = befaft_device;
}

void allocateDeviceMedArr(MedArr **medarr_device_ptr, Range ran) {
    MedArr *medarr_device = (MedArr*)malloc(sizeof(MedArr));
    // cudaMalloc((void **)medarr_device, sizeof(MedArr));

    allocate3DarrayDevice(&medarr_device->ramda, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->mu, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->c11, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->rho, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    
    allocate3DarrayDevice(&medarr_device->zetaxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetaxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetaxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetayx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetayy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetayz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetazx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetazy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetazz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    allocate3DarrayDevice(&medarr_device->gamma, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->khi, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->xi11, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetadx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetady, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&medarr_device->zetadz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    *medarr_device_ptr = medarr_device;
}

// ポインタ以外の部分確保できているかわからない
void allocateDeviceImpulse(Impulse **impulse_device_ptr, Range ran) {
    Impulse *impulse_device = (Impulse *)malloc(sizeof(Impulse));
    allocate3DarrayDevice(&impulse_device->Txx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    allocate3DarrayDevice(&impulse_device->Tyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    allocate3DarrayDevice(&impulse_device->Tzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    *impulse_device_ptr = impulse_device;
}
// ポインタ以外の部分確保できているかわからない
void allocateDeviceDiff(Diff **diff_device_ptr, Range ran) {
    Diff *diff_device = (Diff *)malloc(sizeof(Diff));

    *diff_device_ptr = diff_device;
}

// データ転送 host to device

void copy3DarrayHostToDevice(double ***array_host, double ***array_device, int x, int y, int z) {
    cudaError_t err;
    double **h_level2_device;
    double  *h_level3_device;
    for (int i = 0; i < x; i++) {
        err = cudaMemcpy(&h_level2_device, &(array_device[i]), sizeof(double **), cudaMemcpyDeviceToHost);
        if(err != cudaSuccess) {
            printf("Memcpy error host to device(level2): %s\n", cudaGetErrorString(err));
            return;
        }
        for (int j = 0; j < y; j++) {
            err = cudaMemcpy(&h_level3_device, &(h_level2_device[j]), sizeof(double *), cudaMemcpyDeviceToHost);
            if(err != cudaSuccess) {
                printf("Memcpy error host to device(level3): %s\n", cudaGetErrorString(err));
                return;
            }
            // ホストからデバイスへのデータ転送
            err = cudaMemcpy(h_level3_device, array_host[i][j], z * sizeof(double), cudaMemcpyHostToDevice);
            if(err != cudaSuccess) {
                printf("Memcpy error host to device(Txx.Tx data): %s\n", cudaGetErrorString(err));
                return;
            }
        }
    }
}

// void copy3DarrayHostToDevice(double ***array_host, double ***array_device, int x, int y, int z) {
//     cudaError_t err;
//     for (int i = 0; i < x; i++) {
//         double **h_level2_device = (double **)malloc(y * sizeof(double *));
//         // デバイス上のポインタ配列をホスト側にコピー
//         err = cudaMemcpy(h_level2_device, array_device[i], y * sizeof(double *), cudaMemcpyDeviceToHost);
//         if (err != cudaSuccess) {
//             printf("Memcpy error host to device(level2): %s\n", cudaGetErrorString(err));
//             free(h_level2_device);
//             return;
//         }
//         for (int j = 0; j < y; j++) {
//             double *d_level3 = h_level2_device[j]; // デバイス上のポインタ
//             // ホストからデバイスへのデータ転送
//             err = cudaMemcpy(d_level3, array_host[i][j], z * sizeof(double), cudaMemcpyHostToDevice);
//             if (err != cudaSuccess) {
//                 printf("Memcpy error host to device(Txx.Tx data): %s\n", cudaGetErrorString(err));
//                 free(h_level2_device);
//                 return;
//             }
//         }
//         free(h_level2_device);
//     }
// }


// void copyBefAftHostToDevice(BefAft *befaft_host, BefAft *befaft_device, Range ran) {
//     cudaError_t err;
//     BefAft h_befaft_device;
//     err = cudaMemcpy(&h_befaft_device, befaft_device, sizeof(BefAft), cudaMemcpyDeviceToHost);
//     if(err != cudaSuccess) {
//         printf("Memcpy error host to device(level1): %s\n", cudaGetErrorString(err));
//         return;
//     }
    
//     // Copy function for SigArr
//     copy3DarrayHostToDevice(befaft_host->sa.Txx , h_befaft_device.sa.Txx , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Txxx, h_befaft_device.sa.Txxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Txxy, h_befaft_device.sa.Txxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Txxz, h_befaft_device.sa.Txxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Tyy , h_befaft_device.sa.Tyy , ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Tyyx, h_befaft_device.sa.Tyyx, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Tyyy, h_befaft_device.sa.Tyyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Tyyz, h_befaft_device.sa.Tyyz, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Tzz , h_befaft_device.sa.Tzz , ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Tzzx, h_befaft_device.sa.Tzzx, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Tzzy, h_befaft_device.sa.Tzzy, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
//     copy3DarrayHostToDevice(befaft_host->sa.Tzzz, h_befaft_device.sa.Tzzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

//     // Copy function for TauArr
//     copy3DarrayHostToDevice(befaft_host->ta.Txy , h_befaft_device.ta.Txy , ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
//     copy3DarrayHostToDevice(befaft_host->ta.Txyx, h_befaft_device.ta.Txyx, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
//     copy3DarrayHostToDevice(befaft_host->ta.Txyy, h_befaft_device.ta.Txyy, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
//     copy3DarrayHostToDevice(befaft_host->ta.Tyz , h_befaft_device.ta.Tyz , ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
//     copy3DarrayHostToDevice(befaft_host->ta.Tyzy, h_befaft_device.ta.Tyzy, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
//     copy3DarrayHostToDevice(befaft_host->ta.Tyzz, h_befaft_device.ta.Tyzz, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
//     copy3DarrayHostToDevice(befaft_host->ta.Tzx , h_befaft_device.ta.Tzx , ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
//     copy3DarrayHostToDevice(befaft_host->ta.Tzxz, h_befaft_device.ta.Tzxz, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
//     copy3DarrayHostToDevice(befaft_host->ta.Tzxx, h_befaft_device.ta.Tzxx, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);

//     // Copy function for VelArr
//     copy3DarrayHostToDevice(befaft_host->va.Vx , h_befaft_device.va.Vx , ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vxx, h_befaft_device.va.Vxx, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vxy, h_befaft_device.va.Vxy, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vxz, h_befaft_device.va.Vxz, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vy , h_befaft_device.va.Vy , ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vyx, h_befaft_device.va.Vyx, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vyy, h_befaft_device.va.Vyy, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vyz, h_befaft_device.va.Vyz, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vz , h_befaft_device.va.Vz , ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vzx, h_befaft_device.va.Vzx, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vzy, h_befaft_device.va.Vzy, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
//     copy3DarrayHostToDevice(befaft_host->va.Vzz, h_befaft_device.va.Vzz, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
// }

void copyBefAftHostToDevice(BefAft *befaft_host, BefAft *befaft_device, Range ran) {
    cudaError_t err;
    // 構造体のサイズ分だけデバイスメモリをホスト側にコピー
    BefAft h_befaft_device;
    err = cudaMemcpy(&h_befaft_device, befaft_device, sizeof(BefAft), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        printf("Memcpy error host to devicehhhhhhhhhh(level1): %s\n", cudaGetErrorString(err));
        return;
    }
    
    // Copy function for SigArr
    copy3DarrayHostToDevice(befaft_host->sa.Txx , h_befaft_device.sa.Txx , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(befaft_host->sa.Txxx, h_befaft_device.sa.Txxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(befaft_host->sa.Txxy, h_befaft_device.sa.Txxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(befaft_host->sa.Txxz, h_befaft_device.sa.Txxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tyy , h_befaft_device.sa.Tyy , ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tyyx, h_befaft_device.sa.Tyyx, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tyyy, h_befaft_device.sa.Tyyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tyyz, h_befaft_device.sa.Tyyz, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tzz , h_befaft_device.sa.Tzz , ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tzzx, h_befaft_device.sa.Tzzx, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tzzy, h_befaft_device.sa.Tzzy, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tzzz, h_befaft_device.sa.Tzzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    // Copy function for TauArr
    copy3DarrayHostToDevice(befaft_host->ta.Txy , h_befaft_device.ta.Txy , ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    copy3DarrayHostToDevice(befaft_host->ta.Txyx, h_befaft_device.ta.Txyx, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    copy3DarrayHostToDevice(befaft_host->ta.Txyy, h_befaft_device.ta.Txyy, ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tyz , h_befaft_device.ta.Tyz , ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tyzy, h_befaft_device.ta.Tyzy, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tyzz, h_befaft_device.ta.Tyzz, ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tzx , h_befaft_device.ta.Tzx , ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tzxz, h_befaft_device.ta.Tzxz, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tzxx, h_befaft_device.ta.Tzxx, ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);

    // Copy function for VelArr
    copy3DarrayHostToDevice(befaft_host->va.Vx , h_befaft_device.va.Vx , ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    copy3DarrayHostToDevice(befaft_host->va.Vxx, h_befaft_device.va.Vxx, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    copy3DarrayHostToDevice(befaft_host->va.Vxy, h_befaft_device.va.Vxy, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    copy3DarrayHostToDevice(befaft_host->va.Vxz, h_befaft_device.va.Vxz, ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    copy3DarrayHostToDevice(befaft_host->va.Vy , h_befaft_device.va.Vy , ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    copy3DarrayHostToDevice(befaft_host->va.Vyx, h_befaft_device.va.Vyx, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    copy3DarrayHostToDevice(befaft_host->va.Vyy, h_befaft_device.va.Vyy, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    copy3DarrayHostToDevice(befaft_host->va.Vyz, h_befaft_device.va.Vyz, ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    copy3DarrayHostToDevice(befaft_host->va.Vz , h_befaft_device.va.Vz , ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    copy3DarrayHostToDevice(befaft_host->va.Vzx, h_befaft_device.va.Vzx, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    copy3DarrayHostToDevice(befaft_host->va.Vzy, h_befaft_device.va.Vzy, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    copy3DarrayHostToDevice(befaft_host->va.Vzz, h_befaft_device.va.Vzz, ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
}

//
void copyMedArrHostToDevice(MedArr *medarr_host, MedArr *medarr_device, Range ran) {
    cudaError_t err;
    MedArr h_medarr_device;
    err = cudaMemcpy(&h_medarr_device, medarr_device, sizeof(MedArr), cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) {
        printf("Memcpy error host to device(level1): %s\n", cudaGetErrorString(err));
        return;
    }
    copy3DarrayHostToDevice(medarr_device->ramda , h_medarr_device.ramda , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->mu    , h_medarr_device.mu    , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->c11   , h_medarr_device.c11   , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->rho   , h_medarr_device.rho   , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    
    copy3DarrayHostToDevice(medarr_device->zetaxx, h_medarr_device.zetaxx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetaxy, h_medarr_device.zetaxy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetaxz, h_medarr_device.zetaxz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetayx, h_medarr_device.zetayx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetayy, h_medarr_device.zetayy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetayz, h_medarr_device.zetayz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetazx, h_medarr_device.zetazx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetazy, h_medarr_device.zetazy, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetazz, h_medarr_device.zetazz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);

    copy3DarrayHostToDevice(medarr_device->gamma , h_medarr_device.gamma , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->khi   , h_medarr_device.khi   , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->xi11  , h_medarr_device.xi11  , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetadx, h_medarr_device.zetadx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetady, h_medarr_device.zetady, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(medarr_device->zetadz, h_medarr_device.zetadz, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
}

void copyImpulseHostToDevice(Impulse *impulse_host, Impulse *impulse_device, Range ran) {
    cudaError_t err;
    Impulse h_impulse_device;
    err = cudaMemcpy(&h_impulse_device, impulse_device, sizeof(Impulse), cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) {
        printf("Memcpy error host to device(level1): %s\n", cudaGetErrorString(err));
        return;
    }

    copy3DarrayHostToDevice(impulse_host->Txx, h_impulse_device.Txx, ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(impulse_host->Tyy, h_impulse_device.Tyy, ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    copy3DarrayHostToDevice(impulse_host->Tzz, h_impulse_device.Tzz, ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
}

void copyDiffHostToDevice(Diff *diff_host, Diff *diff_device) {
    cudaError_t err;
    Diff h_diff_device;
    err = cudaMemcpy(&h_diff_device, diff_device, sizeof(Diff), cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) {
        printf("Memcpy error host to device(level1): %s\n", cudaGetErrorString(err));
        return;
    }
}


// データ転送 device to host
void copy3DarrayDeviceToHost(double ***array_device, double ***array_host, int x, int y, int z) {
    // デバイス側のポインタをホスト側にコピー
    cudaError_t err;
    // double ***h_level1_device;
    double  **h_level2_device;
    double   *h_level3_device;
    
    // cudaError_t err = cudaMemcpy(&h_level1_device, array_device, sizeof(double ***), cudaMemcpyDeviceToHost);
    // if(err != cudaSuccess) {
    //     printf("Memcpy error device to host(tau_device): %s\n", cudaGetErrorString(err));
    //     return;
    // }

    // Txx.Tx のデータ転送
    for (int i = 0; i < x; i++) {
        err = cudaMemcpy(&h_level2_device, &(array_device[i]), sizeof(double **), cudaMemcpyDeviceToHost);
        if(err != cudaSuccess) {
            printf("Memcpy error device to host(Txx.Tx level2): %s\n", cudaGetErrorString(err));
            return;
        }
        for (int j = 0; j < y; j++) {
            err = cudaMemcpy(&h_level3_device, &(h_level2_device[j]), sizeof(double *), cudaMemcpyDeviceToHost);
            if(err != cudaSuccess) {
                printf("Memcpy error device to host(Txx.Tx level3): %s\n", cudaGetErrorString(err));
                return;
            }
            // デバイスからホストへのデータ転送
            err = cudaMemcpy(array_host[i][j], h_level3_device, z * sizeof(double), cudaMemcpyDeviceToHost);
            if(err != cudaSuccess) {
                printf("Memcpy error device to host(Txx.Tx data): %s\n", cudaGetErrorString(err));
                return;
            }
        }
    }
}

