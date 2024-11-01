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
    medarrptr = (MedArr*)malloc(ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(MedArr));
    return medarrptr;
}

// BefAft* allocateHostBefAft(Range *ran) {
//     BefAft *baptr;
//     baptr = (BefAft*)malloc(sizeof(BefAft));
//     // SigArr
//     baptr->sa.Txx = (double*)malloc(ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
//     baptr->sa.Txxx = (double*)malloc(ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
//     baptr->sa.Txxy = (double*)malloc(ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
//     baptr->sa.Txxz = (double*)malloc(ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
//     baptr->sa.Tyy = (double*)malloc(ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));
//     baptr->sa.Tyyx = (double*)malloc(ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));
//     baptr->sa.Tyyy = (double*)malloc(ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));
//     baptr->sa.Tyyz = (double*)malloc(ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));
//     baptr->sa.Tzz = (double*)malloc(ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));
//     baptr->sa.Tzzx = (double*)malloc(ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));
//     baptr->sa.Tzzy = (double*)malloc(ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));
//     baptr->sa.Tzzz = (double*)malloc(ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));
//     // TauArr
//     baptr->ta.Txy = (double*)malloc(ran.tr.Txy.x * ran.tr.Txy.y * ran.tr.Txy.z * sizeof(double));
//     baptr->ta.Txyx = (double*)malloc(ran.tr.Txy.x * ran.tr.Txy.y * ran.tr.Txy.z * sizeof(double));
//     baptr->ta.Txyy = (double*)malloc(ran.tr.Txy.x * ran.tr.Txy.y * ran.tr.Txy.z * sizeof(double));
//     baptr->ta.Tyz = (double*)malloc(ran.tr.Tyz.x * ran.tr.Tyz.y * ran.tr.Tyz.z * sizeof(double));
//     baptr->ta.Tyzy = (double*)malloc(ran.tr.Tyz.x * ran.tr.Tyz.y * ran.tr.Tyz.z * sizeof(double));
//     baptr->ta.Tyzz = (double*)malloc(ran.tr.Tyz.x * ran.tr.Tyz.y * ran.tr.Tyz.z * sizeof(double));
//     baptr->ta.Tzx = (double*)malloc(ran.tr.Tzx.x * ran.tr.Tzx.y * ran.tr.Tzx.z * sizeof(double));
//     baptr->ta.Tzxz = (double*)malloc(ran.tr.Tzx.x * ran.tr.Tzx.y * ran.tr.Tzx.z * sizeof(double));
//     baptr->ta.Tzxx = (double*)malloc(ran.tr.Tzx.x * ran.tr.Tzx.y * ran.tr.Tzx.z * sizeof(double));
//     // VelArr
//     baptr->va.Vx = (double*)malloc(ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double));
//     baptr->va.Vxx = (double*)malloc(ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double));
//     baptr->va.Vxy = (double*)malloc(ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double));
//     baptr->va.Vxz = (double*)malloc(ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double));
//     baptr->va.Vy = (double*)malloc(ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double));
//     baptr->va.Vyx = (double*)malloc(ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double));
//     baptr->va.Vyy = (double*)malloc(ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double));
//     baptr->va.Vyz = (double*)malloc(ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double));
//     baptr->va.Vz = (double*)malloc(ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double));
//     baptr->va.Vzx = (double*)malloc(ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double));
//     baptr->va.Vzy = (double*)malloc(ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double));
//     baptr->va.Vzz = (double*)malloc(ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double));
//     return baptr;
// }

Impulse* allocateHostImpulse(int innum) {
    Impulse *impulseptr;
    impulseptr = (Impulse*)malloc(innum * sizeof(Impulse));
    // impulseptr->Txx = (double*)malloc(ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
    // impulseptr->Tyy = (double*)malloc(ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
    // impulseptr->Tzz = (double*)malloc(ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
    // impulseptr->in  = (DimI3*)malloc(innum * sizeof(DimI3));
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
    cudaError_t err = cudaMalloc(&med, ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(MedArr));
    printf("allocateDeviceMedArr:%s\n", cudaGetErrorString(err));
    return med;
}

BefAft* allocateDeviceBefAft(Range ran) {
    BefAft *ba;
    
    ba = (BefAft*)malloc(sizeof(BefAft));

    // SigArr
    cudaError_t err = cudaMalloc(&ba->sa.Txx , ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
    printf("allocateDeviceBefAft:%s\n", cudaGetErrorString(err));
    cudaMalloc(&ba->sa.Txxx, ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
    cudaMalloc(&ba->sa.Txxy, ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));
    cudaMalloc(&ba->sa.Txxz, ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(double));

    cudaMalloc(&ba->sa.Tyy , ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));
    cudaMalloc(&ba->sa.Tyyx, ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));
    cudaMalloc(&ba->sa.Tyyy, ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));
    cudaMalloc(&ba->sa.Tyyz, ran.sr.Tyy.x * ran.sr.Tyy.y * ran.sr.Tyy.z * sizeof(double));

    cudaMalloc(&ba->sa.Tzz , ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));
    cudaMalloc(&ba->sa.Tzzx, ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));
    cudaMalloc(&ba->sa.Tzzy, ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));
    cudaMalloc(&ba->sa.Tzzz, ran.sr.Tzz.x * ran.sr.Tzz.y * ran.sr.Tzz.z * sizeof(double));

    // TauArr
    cudaMalloc(&ba->ta.Txy , ran.tr.Txy.x * ran.tr.Txy.y * ran.tr.Txy.z * sizeof(double));
    cudaMalloc(&ba->ta.Txyx, ran.tr.Txy.x * ran.tr.Txy.y * ran.tr.Txy.z * sizeof(double));
    cudaMalloc(&ba->ta.Txyy, ran.tr.Txy.x * ran.tr.Txy.y * ran.tr.Txy.z * sizeof(double));

    cudaMalloc(&ba->ta.Tyz , ran.tr.Tyz.x * ran.tr.Tyz.y * ran.tr.Tyz.z * sizeof(double));
    cudaMalloc(&ba->ta.Tyzy, ran.tr.Tyz.x * ran.tr.Tyz.y * ran.tr.Tyz.z * sizeof(double));
    cudaMalloc(&ba->ta.Tyzz, ran.tr.Tyz.x * ran.tr.Tyz.y * ran.tr.Tyz.z * sizeof(double));

    cudaMalloc(&ba->ta.Tzx , ran.tr.Tzx.x * ran.tr.Tzx.y * ran.tr.Tzx.z * sizeof(double));
    cudaMalloc(&ba->ta.Tzxz, ran.tr.Tzx.x * ran.tr.Tzx.y * ran.tr.Tzx.z * sizeof(double));
    cudaMalloc(&ba->ta.Tzxx, ran.tr.Tzx.x * ran.tr.Tzx.y * ran.tr.Tzx.z * sizeof(double));

     // VelArr
    cudaMalloc(&ba->va.Vx  , ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double));
    cudaMalloc(&ba->va.Vxx , ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double));
    cudaMalloc(&ba->va.Vxy , ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double));
    cudaMalloc(&ba->va.Vxz , ran.vr.Vx.x * ran.vr.Vx.y * ran.vr.Vx.z * sizeof(double));

    cudaMalloc(&ba->va.Vy  , ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double));
    cudaMalloc(&ba->va.Vyx , ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double));
    cudaMalloc(&ba->va.Vyy , ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double));
    cudaMalloc(&ba->va.Vyz , ran.vr.Vy.x * ran.vr.Vy.y * ran.vr.Vy.z * sizeof(double));

    cudaMalloc(&ba->va.Vz  , ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double));
    cudaMalloc(&ba->va.Vzx , ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double));
    cudaMalloc(&ba->va.Vzy , ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double));
    cudaMalloc(&ba->va.Vzz , ran.vr.Vz.x * ran.vr.Vz.y * ran.vr.Vz.z * sizeof(double));

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
    ipa = (ImpulseArr*)malloc(sizeof(ImpulseArr));
    cudaError_t err = cudaMalloc(&ipa, ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(ImpulseArr));
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
    cudaMemcpy(ma_d, ma_h, ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z * sizeof(MedArr), cudaMemcpyHostToDevice);
}

void ImpulseHostToDevice(Impulse *ip_d, Impulse *ip_h, int innum) {
    cudaMemcpy(ip_d, ip_h, innum * sizeof(Impulse), cudaMemcpyHostToDevice);
}

void DimI3HostToDevice(DimI3 *di_d, DimI3 *di_h, int outnum) {
    cudaMemcpy(di_d, di_h, outnum * sizeof(DimI3), cudaMemcpyHostToDevice);
}

// device to host

void DimD3DeviceToHost(DimD3 *acc_h, DimD3 *acc_d, int outNum) {
    cudaMemcpy(acc_h, acc_d, outNum * sizeof(DimD3), cudaMemcpyDeviceToHost);
}

