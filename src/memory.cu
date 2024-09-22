#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>

#include "../header/struct.h"
#include "../header/memory.h"
////////////////////////////////////////////////////////////////////////////////////////////
// hostメモリ確保
void allocateHostAccCoord(AccCoord **acccoordptr, int tmax, int outNum) {
    *acccoordptr = (AccCoord*)malloc(outNum * sizeof(AccCoord));
    for(int i = 0; i < outNum; i++) {
        (*acccoordptr)[i].x = (double*)malloc(tmax * sizeof(double));
        (*acccoordptr)[i].y = (double*)malloc(tmax * sizeof(double));
        (*acccoordptr)[i].z = (double*)malloc(tmax * sizeof(double));
    }
    if((*acccoordptr)[0].z == NULL) {
        printf("err\n");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////
// deviceメモリ確保

BefAft* allocateDeviceBefAft(Range *ran) {
    BefAft *befaftptr = (BefAft*)malloc(sizeof(BefAft));

    // SigArr
    cudaMalloc((void**)&(befaftptr->sa.Txx ), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->sa.Txxx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->sa.Txxy), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->sa.Txxz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    cudaMalloc((void**)&(befaftptr->sa.Tyy ), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->sa.Tyyx), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->sa.Tyyy), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->sa.Tyyz), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));

    cudaMalloc((void**)&(befaftptr->sa.Tzz ), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->sa.Tzzx), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->sa.Tzzy), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->sa.Tzzz), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));

    // TauArr
    cudaMalloc((void**)&(befaftptr->ta.Txy ), ran->tr.Txy.x * ran->tr.Txy.y * ran->tr.Txy.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->ta.Txyx), ran->tr.Txy.x * ran->tr.Txy.y * ran->tr.Txy.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->ta.Txyy), ran->tr.Txy.x * ran->tr.Txy.y * ran->tr.Txy.z * sizeof(double));

    cudaMalloc((void**)&(befaftptr->ta.Tyz ), ran->tr.Tyz.x * ran->tr.Tyz.y * ran->tr.Tyz.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->ta.Tyzy), ran->tr.Tyz.x * ran->tr.Tyz.y * ran->tr.Tyz.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->ta.Tyzz), ran->tr.Tyz.x * ran->tr.Tyz.y * ran->tr.Tyz.z * sizeof(double));

    cudaMalloc((void**)&(befaftptr->ta.Tzx ), ran->tr.Tzx.x * ran->tr.Tzx.y * ran->tr.Tzx.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->ta.Tzxz), ran->tr.Tzx.x * ran->tr.Tzx.y * ran->tr.Tzx.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->ta.Tzxx), ran->tr.Tzx.x * ran->tr.Tzx.y * ran->tr.Tzx.z * sizeof(double));

        // VelArr
    cudaMalloc((void**)&(befaftptr->va.Vx  ), ran->vr.Vx.x * ran->vr.Vx.y * ran->vr.Vx.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->va.Vxx ), ran->vr.Vx.x * ran->vr.Vx.y * ran->vr.Vx.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->va.Vxy ), ran->vr.Vx.x * ran->vr.Vx.y * ran->vr.Vx.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->va.Vxz ), ran->vr.Vx.x * ran->vr.Vx.y * ran->vr.Vx.z * sizeof(double));

    cudaMalloc((void**)&(befaftptr->va.Vy  ), ran->vr.Vy.x * ran->vr.Vy.y * ran->vr.Vy.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->va.Vyx ), ran->vr.Vy.x * ran->vr.Vy.y * ran->vr.Vy.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->va.Vyy ), ran->vr.Vy.x * ran->vr.Vy.y * ran->vr.Vy.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->va.Vyz ), ran->vr.Vy.x * ran->vr.Vy.y * ran->vr.Vy.z * sizeof(double));

    cudaMalloc((void**)&(befaftptr->va.Vz  ), ran->vr.Vz.x * ran->vr.Vz.y * ran->vr.Vz.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->va.Vzx ), ran->vr.Vz.x * ran->vr.Vz.y * ran->vr.Vz.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->va.Vzy ), ran->vr.Vz.x * ran->vr.Vz.y * ran->vr.Vz.z * sizeof(double));
    cudaMalloc((void**)&(befaftptr->va.Vzz ), ran->vr.Vz.x * ran->vr.Vz.y * ran->vr.Vz.z * sizeof(double));

    return befaftptr;
}

MedArr* allocateDeviceMedArr(Range *ran) {
    MedArr* medarrptr = (MedArr*)malloc(sizeof(MedArr));

    cudaMalloc((void**)&(medarrptr->ramda), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->mu), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->c11), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->rho), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    cudaMalloc((void**)&(medarrptr->zetaxx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->zetaxy), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->zetaxz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    
    cudaMalloc((void**)&(medarrptr->zetayx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->zetayy), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->zetayz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    
    cudaMalloc((void**)&(medarrptr->zetazx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->zetazy), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->zetazz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    cudaMalloc((void**)&(medarrptr->gamma), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->khi), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->xi11), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    cudaMalloc((void**)&(medarrptr->zetadx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->zetady), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(medarrptr->zetadz), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));

    return medarrptr;
}

Impulse* allocateDeviceImpulse(Range *ran) {
    Impulse* impulseptr = (Impulse*)malloc(sizeof(Impulse));

    cudaMalloc((void**)&(impulseptr->Txx), ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    cudaMalloc((void**)&(impulseptr->Tyy), ran->sr.Tyy.x * ran->sr.Tyy.y * ran->sr.Tyy.z * sizeof(double));
    cudaMalloc((void**)&(impulseptr->Tzz), ran->sr.Tzz.x * ran->sr.Tzz.y * ran->sr.Tzz.z * sizeof(double));

    return impulseptr;
}

void allocateDeviceAccCoord(AccCoord **acccoordptr, int tmax, int outNum) {
    cudaMalloc((void**)acccoordptr, outNum * sizeof(AccCoord));
    cudaError_t err;
    for(int i = 0; i < outNum; i++) {
        cudaMalloc((void**)&((*acccoordptr)[i].x), tmax * sizeof(double));
        cudaMalloc((void**)&((*acccoordptr)[i].y), tmax * sizeof(double));
        err = cudaMalloc((void**)&((*acccoordptr)[i].z), tmax * sizeof(double));
    }
    printf("allocate device acc: %s\n", cudaGetErrorString(err));
}

// データ転送 device to host
void AccCoordDeviceToHost(AccCoord *acc_d, AccCoord *acc_h, int outNum, int tmax) {
    cudaError_t err;
    for(int i = 0; i < outNum; i++) {
        err = cudaMemcpy(acc_h[i].x, acc_d[i].x, tmax * sizeof(double), cudaMemcpyDeviceToHost);
        err = cudaMemcpy(acc_h[i].y, acc_d[i].y, tmax * sizeof(double), cudaMemcpyDeviceToHost);
        err = cudaMemcpy(acc_h[i].z, acc_d[i].z, tmax * sizeof(double), cudaMemcpyDeviceToHost);
    }
    
    
    printf("device to host acc: %s", cudaGetErrorString(err));
}
