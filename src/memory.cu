#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>

#include "../header/struct.h"
#include "../header/memory.h"

///////////////////////////////
// hostメモリ確保

MedArr* allocateHostMedArr(Range *ran) {
    MedArr *medarrptr;
    medarrptr = (MedArr*)malloc(sizeof(MedArr));

    medarrptr->c11    = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->gamma  = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->khi    = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->mu     = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->ramda  = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->rho    = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->xi11   = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetadx = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetady = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetadz = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetaxx = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetaxy = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetaxz = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetayx = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetayy = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetayz = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetazx = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetazy = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    medarrptr->zetazz = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    return medarrptr;
}

Impulse* allocateHostImpulse(Range *ran) {
    Impulse *impulseptr;
    impulseptr = (Impulse*)malloc(sizeof(Impulse));

    impulseptr->Txx = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    impulseptr->Tyy = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    impulseptr->Tzz = (double*)malloc(ran->sr.Txx.x * ran->sr.Txx.y * ran->sr.Txx.z * sizeof(double));
    return impulseptr;
}

AccCoord* allocateHostAccCoord(int outNum) {
    AccCoord *acccoordptr;
    acccoordptr = (AccCoord*)malloc(outNum * sizeof(AccCoord));
    return acccoordptr;
}

Coord* allocateHostCoord(int outNum) {
    Coord *coordptr;
    coordptr = (Coord*)malloc(outNum * sizeof(Coord));
    return coordptr;
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
    cudaError_t err = cudaMemcpy(&(d_befaftptr->va.Vzz), &h_befaftptr->va.Vzz, sizeof(double*), cudaMemcpyHostToDevice);

    // printf("allocate device BefAft: %s\n", cudaGetErrorString(err));

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
    cudaError_t err = cudaMemcpy(&d_medarrptr->zetadz, &h_medarrptr->zetadz, sizeof(double*), cudaMemcpyHostToDevice);

    // printf("allocate device MedArr: %s\n", cudaGetErrorString(err));
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
    cudaError_t err = cudaMemcpy(&(d_impulseptr->Tzz), &h_impulseptr->Tzz, sizeof(double*), cudaMemcpyHostToDevice);

    // printf("allocate device impulse: %s\n", cudaGetErrorString(err));
    return d_impulseptr;
}

AccCoord* allocateDeviceAccCoord(int outNum) {
    AccCoord *d_acccoordptr;
    
    cudaError_t err = cudaMalloc(&d_acccoordptr, outNum * sizeof(AccCoord));

    return d_acccoordptr;
}

Coord* allocateDeviceCoord(int outNum) {
    Coord* d_coordptr;
    cudaError_t err = cudaMalloc(&d_coordptr, outNum * sizeof(Coord));

    return d_coordptr;
}

//////////////
// データ転送
// host to device
void MedArrHostToDevice(MedArr *ma_h, MedArr *ma_d, Range ran) {
    MedArr ma_tmp;
    int size = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    double *ramda;
    double *mu;
    double *c11;
    double *rho;
    double *zetaxx;
    double *zetaxy;
    double *zetaxz;
    double *zetayx;
    double *zetayy;
    double *zetayz;
    double *zetazx;
    double *zetazy;
    double *zetazz;
    double *gamma;
    double *khi;
    double *xi11;
    double *zetadx;
    double *zetady;
    double *zetadz;
    cudaMalloc(&ramda, size * sizeof(double));
    cudaMalloc(&mu, size * sizeof(double));
    cudaMalloc(&c11, size * sizeof(double));
    cudaMalloc(&rho, size * sizeof(double));
    cudaMalloc(&zetaxx, size * sizeof(double));
    cudaMalloc(&zetaxy, size * sizeof(double));
    cudaMalloc(&zetaxz, size * sizeof(double));
    cudaMalloc(&zetayx, size * sizeof(double));
    cudaMalloc(&zetayy, size * sizeof(double));
    cudaMalloc(&zetayz, size * sizeof(double));
    cudaMalloc(&zetazx, size * sizeof(double));
    cudaMalloc(&zetazy, size * sizeof(double));
    cudaMalloc(&zetazz, size * sizeof(double));
    cudaMalloc(&gamma, size * sizeof(double));
    cudaMalloc(&khi, size * sizeof(double));
    cudaMalloc(&xi11, size * sizeof(double));
    cudaMalloc(&zetadx, size * sizeof(double));
    cudaMalloc(&zetady, size * sizeof(double));
    cudaMalloc(&zetadz, size * sizeof(double));

    cudaMemcpy(ramda, ma_h->ramda, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(mu, ma_h->mu, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(c11, ma_h->c11, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(rho, ma_h->rho, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetaxx, ma_h->zetaxx, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetaxy, ma_h->zetaxy, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetaxz, ma_h->zetaxz, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetayx, ma_h->zetayx, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetayy, ma_h->zetayy, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetayz, ma_h->zetayz, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetazx, ma_h->zetazx, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetazy, ma_h->zetazy, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetazz, ma_h->zetazz, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gamma, ma_h->gamma, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(khi, ma_h->khi, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(xi11, ma_h->xi11, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetadx, ma_h->zetadx, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetady, ma_h->zetady, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(zetadz, ma_h->zetadz, size * sizeof(double), cudaMemcpyHostToDevice);
    ma_tmp.ramda = ramda;
    ma_tmp.mu = mu;
    ma_tmp.c11 = c11;
    ma_tmp.rho = rho;
    ma_tmp.zetaxx = zetaxx;
    ma_tmp.zetaxy = zetaxy;
    ma_tmp.zetaxz = zetaxz;
    ma_tmp.zetayx = zetayx;
    ma_tmp.zetayy = zetayy;
    ma_tmp.zetayz = zetayz;
    ma_tmp.zetazx = zetazx;
    ma_tmp.zetazy = zetazy;
    ma_tmp.zetazz = zetazz;
    ma_tmp.gamma = gamma;
    ma_tmp.khi = khi;
    ma_tmp.xi11 = xi11;
    ma_tmp.zetadx = zetadx;
    ma_tmp.zetady = zetady;
    ma_tmp.zetadz = zetadz;

    cudaError_t err = cudaMemcpy(ma_d, &ma_tmp, sizeof(MedArr), cudaMemcpyHostToDevice); 
    // printf("host to device MedArr: %s\n", cudaGetErrorString(err));
}

void ImpulseHostToDevice(Impulse *ip_h, Impulse *ip_d, Range ran) {
    Impulse ip_tmp;
    int size = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    double *Txx;
    double *Tyy;
    double *Tzz;
    cudaMalloc(&Txx, size * sizeof(double));
    cudaMalloc(&Tyy, size * sizeof(double));
    cudaMalloc(&Tzz, size * sizeof(double));


    cudaMemcpy(Txx, ip_h->Txx, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Tyy, ip_h->Tyy, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Tzz, ip_h->Tzz, size * sizeof(double), cudaMemcpyHostToDevice);


    ip_tmp.Txx = Txx;
    ip_tmp.Tyy = Tyy;
    ip_tmp.Tzz = Tzz;
    ip_tmp.freq = ip_h->freq;
    ip_tmp.mode = ip_h->mode;
    ip_tmp.in = ip_h->in;
    // ip_tmp.in.x = ip_h->in.x;
    // ip_tmp.in.y = ip_h->in.y;
    // ip_tmp.in.z = ip_h->in.z;

    cudaError_t err = cudaMemcpy(ip_d, &ip_tmp, sizeof(Impulse), cudaMemcpyHostToDevice);
    // printf("host to device Impulse: %s\n", cudaGetErrorString(err));
    // cudaFree(Txx);
    // cudaFree(Tyy);
    // cudaFree(Tzz);
}

void RangeHostToDevice(Range *ran_h, Range *ran_d) {
    cudaError_t err = cudaMemcpy(ran_d, ran_h, sizeof(Range), cudaMemcpyHostToDevice);

    // printf("host to device Range: %s\n", cudaGetErrorString(err));
}

void DiffHostToDevice(Diff *dif_h, Diff *dif_d) {
    cudaError_t err = cudaMemcpy(dif_d, dif_h, sizeof(Diff), cudaMemcpyHostToDevice);

    // printf("host to device Diff: %s\n", cudaGetErrorString(err));
}

void CoordHostToDevice(Coord *out_h, Coord *out_d, int outNum) {
    cudaError_t err = cudaMemcpy(out_d, out_h, outNum * sizeof(Coord), cudaMemcpyHostToDevice);
    // printf("host to device Coord: %s\n", cudaGetErrorString(err));
}
// device to host

void AccCoordDeviceToHost(AccCoord *acc_d, AccCoord *acc_h, int outNum) {
    // AccCoord acc_ptr;
    cudaError_t err = cudaMemcpy(acc_h, acc_d, outNum * sizeof(AccCoord), cudaMemcpyDeviceToHost);
    // printf("device to host AccCoord: %s\n", cudaGetErrorString(err));
}

void BefAftDeviceToHost(BefAft *ba_d, BefAft *ba_h, Range ran) {
    int N = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    // Step 1: デバイスからホストへ構造体をコピー
    cudaMemcpy(ba_h, ba_d, sizeof(BefAft), cudaMemcpyDeviceToHost);

    // Step 2: デバイス側のポインタを保存
    double *device_Txx = ba_h->sa.Txx;
    double *device_Tyy = ba_h->sa.Tyy;
    double *device_Tzz = ba_h->sa.Tzz;

    double *device_Txy = ba_h->ta.Txy;
    double *device_Tyz = ba_h->ta.Tyz;
    double *device_Tzx = ba_h->ta.Tzx;

    double *device_Vx = ba_h->va.Vx;
    double *device_Vy = ba_h->va.Vy;
    double *device_Vz = ba_h->va.Vz;

    // Step 3: ホスト側でメモリを確保
    ba_h->sa.Txx = (double *)malloc(N * sizeof(double));
    ba_h->sa.Tyy = (double *)malloc(N * sizeof(double));
    ba_h->sa.Tzz = (double *)malloc(N * sizeof(double));

    ba_h->ta.Txy = (double *)malloc(N * sizeof(double));
    ba_h->ta.Tyz = (double *)malloc(N * sizeof(double));
    ba_h->ta.Tzx = (double *)malloc(N * sizeof(double));

    ba_h->va.Vx = (double *)malloc(N * sizeof(double));
    ba_h->va.Vy = (double *)malloc(N * sizeof(double));
    ba_h->va.Vz = (double *)malloc(N * sizeof(double));

    // Step 4: デバイスからホストへデータをコピー
    cudaMemcpy(ba_h->sa.Txx, device_Txx, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->sa.Tyy, device_Tyy, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->sa.Tzz, device_Tzz, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaMemcpy(ba_h->ta.Txy, device_Txy, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->ta.Tyz, device_Tyz, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->ta.Tzx, device_Tzx, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaMemcpy(ba_h->va.Vx, device_Vx, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->va.Vy, device_Vy, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->va.Vz, device_Vz, N * sizeof(double), cudaMemcpyDeviceToHost);
}
