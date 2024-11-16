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

MedArr* allocateHostMedArr(Range ran) {
    MedArr *medarrptr;
    int cell = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    medarrptr = (MedArr*)malloc(cell * sizeof(MedArr));
    // memset(medarrptr, 0, cell * sizeof(MedArr));
    return medarrptr;
}

BefAft* allocateHostBefAft(Range ran) {
    BefAft *ba;
    ba = (BefAft*)malloc(sizeof(BefAft));
    ba->sa.Txx = (double*)malloc(ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
    ba->sa.Tyy = (double*)malloc(ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));
    ba->sa.Tzz = (double*)malloc(ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));
    ba->ta.Txy = (double*)malloc(ran.tr.Txy.x * ran.tr.Txy.y * ran.tr.Txy.z * sizeof(double));
    ba->ta.Tyz = (double*)malloc(ran.tr.Tyz.x * ran.tr.Tyz.y * ran.tr.Tyz.z * sizeof(double));
    ba->ta.Tzx = (double*)malloc(ran.tr.Tzx.x * ran.tr.Tzx.y * ran.tr.Tzx.z * sizeof(double));
    ba->va.Vx  = (double*)malloc(ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double));
    ba->va.Vy  = (double*)malloc(ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double));
    ba->va.Vz  = (double*)malloc(ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double));
    return ba;
}

Impulse* allocateHostImpulse(int innum) {
    Impulse *impulseptr;
    impulseptr = (Impulse*)malloc(innum * sizeof(Impulse));
    return impulseptr;
}

DimI3* allocateHostDimI3(int outnum) {
    DimI3 *DI;
    DI = (DimI3*)malloc(outnum * sizeof(DimI3));
    return DI;
}

DimD3* allocateHostDimD3(int outnum) {
    DimD3 *DD;
    DD = (DimD3*)malloc(outnum * sizeof(DimD3));
    return DD;
}

////////////////////////////////////////////////////////////////////////////////////////////
// deviceメモリ確保
MedArr* allocateDeviceMedArr(Range ran) {
    MedArr *med;
    int cell = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    cudaError_t err = cudaMalloc(&med, cell * sizeof(MedArr));
    printf("allocateDeviceMedArr:%s\n", cudaGetErrorString(err));
    return med;
}

BefAft* allocateDeviceBefAft(Range ran) {
    BefAft *ba;
    
    ba = (BefAft*)malloc(sizeof(BefAft));

    int celltxx = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    int celltyy = ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z;
    int celltzz = ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z;
    int celltxy = ran.tr.Txy.x * ran.tr.Txy.y * ran.tr.Txy.z;
    int celltyz = ran.tr.Tyz.x * ran.tr.Tyz.y * ran.tr.Tyz.z;
    int celltzx = ran.tr.Tzx.x * ran.tr.Tzx.y * ran.tr.Tzx.z;
    int cellvx  = ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z;
    int cellvy  = ran.vr.Vx.x * ran.vr.Vy.y * ran.vr.Vy.z;
    int cellvz  = ran.vr.Vx.x * ran.vr.Vz.y * ran.vr.Vz.z;
    // SigArr
    cudaError_t err = cudaMalloc(&ba->sa.Txx , celltxx * sizeof(double));
    printf("allocateDeviceBefAft:%s\n", cudaGetErrorString(err));
    cudaMalloc(&ba->sa.Txxx, celltxx * sizeof(double));
    cudaMalloc(&ba->sa.Txxy, celltxx * sizeof(double));
    cudaMalloc(&ba->sa.Txxz, celltxx * sizeof(double));

    cudaMalloc(&ba->sa.Tyy , celltyy * sizeof(double));
    cudaMalloc(&ba->sa.Tyyx, celltyy * sizeof(double));
    cudaMalloc(&ba->sa.Tyyy, celltyy * sizeof(double));
    cudaMalloc(&ba->sa.Tyyz, celltyy * sizeof(double));

    cudaMalloc(&ba->sa.Tzz , celltzz * sizeof(double));
    cudaMalloc(&ba->sa.Tzzx, celltzz * sizeof(double));
    cudaMalloc(&ba->sa.Tzzy, celltzz * sizeof(double));
    cudaMalloc(&ba->sa.Tzzz, celltzz * sizeof(double));

    // TauArr
    cudaMalloc(&ba->ta.Txy , celltxy * sizeof(double));
    cudaMalloc(&ba->ta.Txyx, celltxy * sizeof(double));
    cudaMalloc(&ba->ta.Txyy, celltxy * sizeof(double));

    cudaMalloc(&ba->ta.Tyz , celltyz * sizeof(double));
    cudaMalloc(&ba->ta.Tyzy, celltyz * sizeof(double));
    cudaMalloc(&ba->ta.Tyzz, celltyz * sizeof(double));

    cudaMalloc(&ba->ta.Tzx , celltzx * sizeof(double));
    cudaMalloc(&ba->ta.Tzxz, celltzx * sizeof(double));
    cudaMalloc(&ba->ta.Tzxx, celltzx * sizeof(double));

     // VelArr
    cudaMalloc(&ba->va.Vx  , cellvx * sizeof(double));
    cudaMalloc(&ba->va.Vxx , cellvx * sizeof(double));
    cudaMalloc(&ba->va.Vxy , cellvx * sizeof(double));
    cudaMalloc(&ba->va.Vxz , cellvx * sizeof(double));

    cudaMalloc(&ba->va.Vy  , cellvy * sizeof(double));
    cudaMalloc(&ba->va.Vyx , cellvy * sizeof(double));
    cudaMalloc(&ba->va.Vyy , cellvy * sizeof(double));
    cudaMalloc(&ba->va.Vyz , cellvy * sizeof(double));

    cudaMalloc(&ba->va.Vz  , cellvz * sizeof(double));
    cudaMalloc(&ba->va.Vzx , cellvz * sizeof(double));
    cudaMalloc(&ba->va.Vzy , cellvz * sizeof(double));
    cudaMalloc(&ba->va.Vzz , cellvz * sizeof(double));

    return ba;
}

// ImpulseArr* allocateDeviceImpulseArr(Range ran) {
//     ImpulseArr *ip;
//     ip = (ImpulseArr*)malloc(sizeof(ImpulseArr));
//     cudaError_t err = cudaMalloc(&ip->Txx, ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
//     printf("allocateDeviceImpulseArr:%s\n", cudaGetErrorString(err));
//     cudaMalloc(&ip->Tyy, ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));
//     cudaMalloc(&ip->Tzz, ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));
//     return ip;
// }

ImpulseArr* allocateDeviceImpulseArr(Range ran) {
    ImpulseArr *ipa;
    int cell = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    ipa = (ImpulseArr*)malloc(sizeof(ImpulseArr));
    cudaError_t err = cudaMalloc(&ipa, cell * sizeof(ImpulseArr));
    printf("allocateDeviceImpulseArr:%s\n", cudaGetErrorString(err));
    return ipa;
}

Impulse* allocateDeviceImpulse(int innum) {
    Impulse *ip; // ホスト側の一時構造体
    cudaError_t err = cudaMalloc(&ip, innum * sizeof(Impulse));
    printf("allocateDeviceImpulse:%s\n", cudaGetErrorString(err));
    return ip;
}

DimI3* allocateDeviceDimI3(int outNum) {
    DimI3 *DI;
    cudaError_t err = cudaMalloc(&DI, outNum * sizeof(DimD3));
    printf("allocateDeviceDimI3:%s\n", cudaGetErrorString(err));
    return DI;
}

DimD3* allocateDeviceDimD3(int outNum) {
    DimD3 *DD;
    cudaError_t err = cudaMalloc(&DD, outNum * sizeof(DimD3));
    printf("allocateDeviceDimD3:%s\n", cudaGetErrorString(err));
    return DD;
}


//////////////
// データ転送
// host to device
void RangeHostToDevice(Range *ran_d, Range *ran_h) {
    cudaMemcpy(ran_d, ran_h, sizeof(Range), cudaMemcpyHostToDevice);
}

void DiffHostToDevice(Diff *dif_d, Diff *dif_h) {
    cudaMemcpy(dif_d, dif_h, sizeof(Diff), cudaMemcpyHostToDevice);
}

void MedArrHostToDevice(MedArr *ma_d, MedArr *ma_h, Range ran) {
    int cell = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    cudaMemcpy(ma_d, ma_h, cell * sizeof(MedArr), cudaMemcpyHostToDevice);
}

void ImpulseHostToDevice(Impulse *ip_d, Impulse *ip_h, int innum) {
    cudaMemcpy(ip_d, ip_h, innum * sizeof(Impulse), cudaMemcpyHostToDevice);
}

void DimI3HostToDevice(DimI3 *di_d, DimI3 *di_h, int outnum) {
    cudaMemcpy(di_d, di_h, outnum * sizeof(DimI3), cudaMemcpyHostToDevice);
}

// device to host

void BefAftDeviceToHost(BefAft *ba_h, BefAft *ba_d, Range ran) {
    cudaMemcpy(ba_h->sa.Txx, ba_d->sa.Txx, ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->sa.Tyy, ba_d->sa.Tyy, ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->sa.Tzz, ba_d->sa.Tzz, ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->ta.Txy, ba_d->ta.Txy, ran.tr.Txy.x * ran.tr.Txy.y * ran.tr.Txy.z * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->ta.Tyz, ba_d->ta.Tyz, ran.tr.Tyz.x * ran.tr.Tyz.y * ran.tr.Tyz.z * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->ta.Tzx, ba_d->ta.Tzx, ran.tr.Tzx.x * ran.tr.Tzx.y * ran.tr.Tzx.z * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->va.Vx , ba_d->va.Vx , ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->va.Vy , ba_d->va.Vy , ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->va.Vz , ba_d->va.Vz , ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double), cudaMemcpyDeviceToHost);
}

void DimD3DeviceToHost(DimD3 *acc_h, DimD3 *acc_d, int outNum) {
    cudaMemcpy(acc_h, acc_d, outNum * sizeof(DimD3), cudaMemcpyDeviceToHost);
}

