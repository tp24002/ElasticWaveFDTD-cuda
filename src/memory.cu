#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>

#include "../header/struct.h"
#include "../header/memory.h"


AccCoord* allocateHostAccCoord(int outNum) {
    AccCoord *acccoordptr;
    acccoordptr = (AccCoord*)malloc(outNum * sizeof(AccCoord));

    // acccoordptr->x = (double*)malloc(outNum * sizeof(double));
    // acccoordptr->y = (double*)malloc(outNum * sizeof(double));
    // acccoordptr->z = (double*)malloc(outNum * sizeof(double));

    // printf("%p\n",acccoordptr);
    // printf("%p\n",&acccoordptr->y);
    return acccoordptr;

}
////////////////////////////////////////////////////////////////////////////////////////////
// deviceメモリ確保

BefAft* allocateDeviceBefAft(Range *ran) {
    BefAft *d_befaftptr;// = (BefAft*)malloc(sizeof(BefAft));
    BefAft *h_befaftptr;

    cudaMalloc(&d_befaftptr, sizeof(BefAft));
    h_befaftptr = (BefAft*)malloc(sizeof(BefAft));

    // SigArr
    cudaMalloc((void**)&(h_befaftptr->sa.Txx ), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->sa.Txxx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->sa.Txxy), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->sa.Txxz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    cudaMalloc((void**)&(h_befaftptr->sa.Tyy ), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->sa.Tyyx), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->sa.Tyyy), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->sa.Tyyz), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));

    cudaMalloc((void**)&(h_befaftptr->sa.Tzz ), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->sa.Tzzx), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->sa.Tzzy), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->sa.Tzzz), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));

    // TauArr
    cudaMalloc((void**)&(h_befaftptr->ta.Txy ), ran->tr.Txy.x * ran->tr.Txy.y * ran->tr.Txy.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->ta.Txyx), ran->tr.Txy.x * ran->tr.Txy.y * ran->tr.Txy.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->ta.Txyy), ran->tr.Txy.x * ran->tr.Txy.y * ran->tr.Txy.z * sizeof(double));

    cudaMalloc((void**)&(h_befaftptr->ta.Tyz ), ran->tr.Tyz.x * ran->tr.Tyz.y * ran->tr.Tyz.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->ta.Tyzy), ran->tr.Tyz.x * ran->tr.Tyz.y * ran->tr.Tyz.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->ta.Tyzz), ran->tr.Tyz.x * ran->tr.Tyz.y * ran->tr.Tyz.z * sizeof(double));

    cudaMalloc((void**)&(h_befaftptr->ta.Tzx ), ran->tr.Tzx.x * ran->tr.Tzx.y * ran->tr.Tzx.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->ta.Tzxz), ran->tr.Tzx.x * ran->tr.Tzx.y * ran->tr.Tzx.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->ta.Tzxx), ran->tr.Tzx.x * ran->tr.Tzx.y * ran->tr.Tzx.z * sizeof(double));

        // VelArr
    cudaMalloc((void**)&(h_befaftptr->va.Vx  ), ran->vr.Vx.x * ran->vr.Vx.y * ran->vr.Vx.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->va.Vxx ), ran->vr.Vx.x * ran->vr.Vx.y * ran->vr.Vx.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->va.Vxy ), ran->vr.Vx.x * ran->vr.Vx.y * ran->vr.Vx.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->va.Vxz ), ran->vr.Vx.x * ran->vr.Vx.y * ran->vr.Vx.z * sizeof(double));

    cudaMalloc((void**)&(h_befaftptr->va.Vy  ), ran->vr.Vy.x * ran->vr.Vy.y * ran->vr.Vy.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->va.Vyx ), ran->vr.Vy.x * ran->vr.Vy.y * ran->vr.Vy.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->va.Vyy ), ran->vr.Vy.x * ran->vr.Vy.y * ran->vr.Vy.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->va.Vyz ), ran->vr.Vy.x * ran->vr.Vy.y * ran->vr.Vy.z * sizeof(double));

    cudaMalloc((void**)&(h_befaftptr->va.Vz  ), ran->vr.Vz.x * ran->vr.Vz.y * ran->vr.Vz.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->va.Vzx ), ran->vr.Vz.x * ran->vr.Vz.y * ran->vr.Vz.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->va.Vzy ), ran->vr.Vz.x * ran->vr.Vz.y * ran->vr.Vz.z * sizeof(double));
    cudaMalloc((void**)&(h_befaftptr->va.Vzz ), ran->vr.Vz.x * ran->vr.Vz.y * ran->vr.Vz.z * sizeof(double));
    // Copy SigArr (Stress tensors)
    cudaMemcpy(&(d_befaftptr->sa.Txx), &h_befaftptr->sa.Txx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->sa.Txxx), &h_befaftptr->sa.Txxx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->sa.Txxy), &h_befaftptr->sa.Txxy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->sa.Txxz), &h_befaftptr->sa.Txxz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&(d_befaftptr->sa.Tyy), &h_befaftptr->sa.Tyy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->sa.Tyyx), &h_befaftptr->sa.Tyyx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->sa.Tyyy), &h_befaftptr->sa.Tyyy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->sa.Tyyz), &h_befaftptr->sa.Tyyz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&(d_befaftptr->sa.Tzz), &h_befaftptr->sa.Tzz, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->sa.Tzzx), &h_befaftptr->sa.Tzzx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->sa.Tzzy), &h_befaftptr->sa.Tzzy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->sa.Tzzz), &h_befaftptr->sa.Tzzz, sizeof(double*), cudaMemcpyHostToDevice);

    // Copy TauArr (Shear stress tensors)
    cudaMemcpy(&(d_befaftptr->ta.Txy), &h_befaftptr->ta.Txy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->ta.Txyx), &h_befaftptr->ta.Txyx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->ta.Txyy), &h_befaftptr->ta.Txyy, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&(d_befaftptr->ta.Tyz), &h_befaftptr->ta.Tyz, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->ta.Tyzy), &h_befaftptr->ta.Tyzy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->ta.Tyzz), &h_befaftptr->ta.Tyzz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&(d_befaftptr->ta.Tzx), &h_befaftptr->ta.Tzx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->ta.Tzxz), &h_befaftptr->ta.Tzxz, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->ta.Tzxx), &h_befaftptr->ta.Tzxx, sizeof(double*), cudaMemcpyHostToDevice);

    // Copy VelArr (Velocities)
    cudaMemcpy(&(d_befaftptr->va.Vx), &h_befaftptr->va.Vx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->va.Vxx), &h_befaftptr->va.Vxx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->va.Vxy), &h_befaftptr->va.Vxy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->va.Vxz), &h_befaftptr->va.Vxz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&(d_befaftptr->va.Vy), &h_befaftptr->va.Vy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->va.Vyx), &h_befaftptr->va.Vyx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->va.Vyy), &h_befaftptr->va.Vyy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->va.Vyz), &h_befaftptr->va.Vyz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&(d_befaftptr->va.Vz), &h_befaftptr->va.Vz, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->va.Vzx), &h_befaftptr->va.Vzx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->va.Vzy), &h_befaftptr->va.Vzy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_befaftptr->va.Vzz), &h_befaftptr->va.Vzz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaError_t err = cudaGetLastError();
    printf("allocate device befaft: %s\n", cudaGetErrorString(err));

    return d_befaftptr;
}

MedArr* allocateDeviceMedArr(Range *ran) {
    MedArr *d_medarrptr;
    MedArr *h_medarrptr;

    cudaMalloc(&d_medarrptr, sizeof(MedArr));
    h_medarrptr = (MedArr*)malloc(sizeof(MedArr));

    cudaMalloc((void**)&(h_medarrptr->ramda), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->mu), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->c11), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->rho), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    cudaMalloc((void**)&(h_medarrptr->zetaxx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->zetaxy), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->zetaxz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    
    cudaMalloc((void**)&(h_medarrptr->zetayx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->zetayy), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->zetayz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    
    cudaMalloc((void**)&(h_medarrptr->zetazx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->zetazy), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->zetazz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    cudaMalloc((void**)&(h_medarrptr->gamma), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->khi), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->xi11), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    cudaMalloc((void**)&(h_medarrptr->zetadx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->zetady), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_medarrptr->zetadz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    cudaMemcpy(&d_medarrptr->ramda, &h_medarrptr->ramda, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->mu, &h_medarrptr->mu, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->c11, &h_medarrptr->c11, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->rho, &h_medarrptr->rho, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->zetaxx, &h_medarrptr->zetaxx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetaxy, &h_medarrptr->zetaxy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetaxz, &h_medarrptr->zetaxz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->zetayx, &h_medarrptr->zetayx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetayy, &h_medarrptr->zetayy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetayz, &h_medarrptr->zetayz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->zetazx, &h_medarrptr->zetazx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetazy, &h_medarrptr->zetazy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetazz, &h_medarrptr->zetazz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->gamma, &h_medarrptr->gamma, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->khi, &h_medarrptr->khi, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->xi11, &h_medarrptr->xi11, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->zetadx, &h_medarrptr->zetadx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetady, &h_medarrptr->zetady, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetadz, &h_medarrptr->zetadz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaError_t err = cudaGetLastError();
    printf("allocate device medarr: %s\n", cudaGetErrorString(err));
    // デバイスメモリのポインタを返す
    return d_medarrptr;

}

Impulse* allocateDeviceImpulse(Range *ran) {
    Impulse *d_impulseptr;
    Impulse *h_impulseptr; // ホスト側の一時構造体

    cudaMalloc((void**)&d_impulseptr, sizeof(Impulse)); // Impulse構造体のためのメモリをデバイスに確保
    h_impulseptr = (Impulse*)malloc(sizeof(Impulse));

    cudaMalloc((void**)&(h_impulseptr->Txx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(h_impulseptr->Tyy), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));
    cudaMalloc((void**)&(h_impulseptr->Tzz), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));

    // メモリのポインタをデバイスに戻す
    cudaMemcpy(&(d_impulseptr->Txx), &h_impulseptr->Txx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_impulseptr->Tyy), &h_impulseptr->Tyy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_impulseptr->Tzz), &h_impulseptr->Tzz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaError_t err = cudaGetLastError();
    printf("allocate device impulse: %s\n", cudaGetErrorString(err));
    return d_impulseptr;
}


AccCoord* allocateDeviceAccCoord(int outNum) {
    // AccCoord *h_acccoordptr;
    AccCoord *d_acccoordptr;
    
    cudaMalloc(&d_acccoordptr, outNum * sizeof(AccCoord));  // デバイスメモリにAccCoord配列を確保
    // h_acccoordptr = (AccCoord*)malloc(outNum * sizeof(AccCoord));  // ホストメモリにAccCoord配列を確保

    // // 各AccCoord構造体のx, y, zに対してメモリ確保
    // cudaMalloc((void**)&(h_acccoordptr->x), outNum * sizeof(double));
    // cudaMalloc((void**)&(h_acccoordptr->y), outNum * sizeof(double));
    // cudaMalloc((void**)&(h_acccoordptr->z), outNum * sizeof(double));

    // // ホストからデバイスにx, y, zのポインタをコピー
    // cudaMemcpy(&(d_acccoordptr->x), &(h_acccoordptr->x), sizeof(double*), cudaMemcpyHostToDevice);
    // cudaMemcpy(&(d_acccoordptr->y), &(h_acccoordptr->y), sizeof(double*), cudaMemcpyHostToDevice);
    // cudaMemcpy(&(d_acccoordptr->z), &(h_acccoordptr->z), sizeof(double*), cudaMemcpyHostToDevice);

    cudaError_t err = cudaGetLastError();
    printf("allocate device acccoord: %s\n", cudaGetErrorString(err));
    return d_acccoordptr;
}


//////////////
// データ転送

void AccCoordDeviceToHost(AccCoord *acc_d, AccCoord *acc_h, int outNum) {
    // AccCoord acc_ptr;
    cudaError_t err = cudaMemcpy(acc_h, acc_d, outNum * sizeof(double), cudaMemcpyDeviceToHost);

    // cudaMalloc(&acc_ptr.x, outNum * sizeof(double));
    // cudaMalloc(&acc_ptr.y, outNum * sizeof(double));
    // cudaMalloc(&acc_ptr.z, outNum * sizeof(double));

    // cudaMemcpy(&acc_ptr.x, &acc_d->x, sizeof(double*), cudaMemcpyDeviceToHost);
    // cudaMemcpy(&acc_ptr.y, &acc_d->y, sizeof(double*), cudaMemcpyDeviceToHost);
    // cudaMemcpy(&acc_ptr.z, &acc_d->z, sizeof(double*), cudaMemcpyDeviceToHost);

    // cudaMemcpy(acc_h->x, acc_ptr.x, outNum * sizeof(double), cudaMemcpyDeviceToHost);
    // cudaMemcpy(acc_h->y, acc_ptr.y, outNum * sizeof(double), cudaMemcpyDeviceToHost);
    // err = cudaMemcpy(acc_h->z, acc_ptr.z, outNum * sizeof(double), cudaMemcpyDeviceToHost);
    // cudaError_t err = cudaGetLastError();

    printf("acc device to host:%s\n", cudaGetErrorString(err));
}

