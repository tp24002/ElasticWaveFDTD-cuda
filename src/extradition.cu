#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>

#include "../header/struct.h"


void MemoryBefAftToDevice(BefAft *h, BefAft *d, Range ran) {
    // cudaMemcpy(&d, &h, sizeof(BefAft), cudaMemcpyHostToDevice);
    cudaMemcpy(&d->sa, &h->sa, sizeof(SigArr), cudaMemcpyHostToDevice);
    cudaMemcpy(&d->ta, &h->ta, sizeof(TauArr), cudaMemcpyHostToDevice);
    cudaMemcpy(&d->va, &h->va, sizeof(VelArr), cudaMemcpyHostToDevice);
    for(int i = 0; i <= ran.sr.Txx.x; i++) {
        for(int j = 0; j <= ran.sr.Txx.y; j++){
            cudaMemcpy(&d->sa.Txx, &h->sa.Txx, sizeof(double) *  (ran.sr.Txx.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->sa.Txxx, &h->sa.Txxx, sizeof(double) *  (ran.sr.Txx.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->sa.Txxy, &h->sa.Txxy, sizeof(double) *  (ran.sr.Txx.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->sa.Txxz, &h->sa.Txxz, sizeof(double) *  (ran.sr.Txx.z + 1), cudaMemcpyHostToDevice);
        }
    }
    for(int i = 0; i <= ran.sr.Tyy.x; i++) {
        for(int j = 0; j <= ran.sr.Tyy.y; j++){
            cudaMemcpy(&d->sa.Tyy[i][j], &h->sa.Tyy[i][j], sizeof(double) *  (ran.sr.Tyy.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->sa.Tyyx[i][j], &h->sa.Tyyx[i][j], sizeof(double) *  (ran.sr.Tyy.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->sa.Tyyy[i][j], &h->sa.Tyyy[i][j], sizeof(double) *  (ran.sr.Tyy.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->sa.Tyyz[i][j], &h->sa.Tyyz[i][j], sizeof(double) *  (ran.sr.Tyy.z + 1), cudaMemcpyHostToDevice);
        }
    }
    for(int i = 0; i <= ran.sr.Tzz.x; i++) {
        for(int j = 0; j <= ran.sr.Tzz.y; j++){
            cudaMemcpy(&d->sa.Tzz[i][j], &h->sa.Tzz[i][j], sizeof(double) *  (ran.sr.Tzz.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->sa.Tzzx[i][j], &h->sa.Tzzx[i][j], sizeof(double) *  (ran.sr.Tzz.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->sa.Tzzy[i][j], &h->sa.Tzzy[i][j], sizeof(double) *  (ran.sr.Tzz.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->sa.Tzzz[i][j], &h->sa.Tzzz[i][j], sizeof(double) *  (ran.sr.Tzz.z + 1), cudaMemcpyHostToDevice);
        }
    }
    //////
    for(int i = 0; i <= ran.tr.Txy.x; i++) {
        for(int j = 0; j <= ran.tr.Txy.y; j++){
            cudaMemcpy(&d->ta.Txy[i][j], &h->ta.Txy[i][j], sizeof(double) *  (ran.tr.Txy.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->ta.Txyx[i][j], &h->ta.Txyx[i][j], sizeof(double) *  (ran.tr.Txy.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->ta.Txyy[i][j], &h->ta.Txyy[i][j], sizeof(double) *  (ran.tr.Txy.z + 1), cudaMemcpyHostToDevice);
        }
    }
    for(int i = 0; i <= ran.tr.Tyz.x; i++) {
        for(int j = 0; j <= ran.tr.Tyz.y; j++){
            cudaMemcpy(&d->ta.Tyz[i][j], &h->ta.Tyz[i][j], sizeof(double) *  (ran.tr.Tyz.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->ta.Tyzy[i][j], &h->ta.Tyzy[i][j], sizeof(double) *  (ran.tr.Tyz.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->ta.Tyzz[i][j], &h->ta.Tyzz[i][j], sizeof(double) *  (ran.tr.Tyz.z + 1), cudaMemcpyHostToDevice);
        }
    }
    for(int i = 0; i <= ran.tr.Tzx.x; i++) {
        for(int j = 0; j <= ran.tr.Tzx.y; j++){
            cudaMemcpy(&d->ta.Tzx[i][j], &h->ta.Tzx[i][j], sizeof(double) *  (ran.tr.Tzx.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->ta.Tzxz[i][j], &h->ta.Tzxz[i][j], sizeof(double) *  (ran.tr.Tzx.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->ta.Tzxx[i][j], &h->ta.Tzxx[i][j], sizeof(double) *  (ran.tr.Tzx.z + 1), cudaMemcpyHostToDevice);
        }
    }
    ///////
    for(int i = 0; i <= ran.vr.Vx.x; i++) {
        for(int j = 0; j <= ran.vr.Vx.y; j++){
            cudaMemcpy(&d->va.Vx, &h->va.Vx, sizeof(double) * ( ran.vr.Vx.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->va.Vxx, &h->va.Vxx, sizeof(double) * ( ran.vr.Vx.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->va.Vxy, &h->va.Vxy, sizeof(double) * ( ran.vr.Vx.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->va.Vxz, &h->va.Vxz, sizeof(double) * ( ran.vr.Vx.z + 1), cudaMemcpyHostToDevice);
        }
    }
    for(int i = 0; i <= ran.vr.Vy.x; i++) {
        for(int j = 0; j <= ran.vr.Vy.y; j++){
            cudaMemcpy(&d->va.Vy, &h->va.Vy, sizeof(double) * ( ran.vr.Vy.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->va.Vyx, &h->va.Vyx, sizeof(double) * ( ran.vr.Vy.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->va.Vyy, &h->va.Vyy, sizeof(double) * ( ran.vr.Vy.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->va.Vyz, &h->va.Vyz, sizeof(double) * ( ran.vr.Vy.z + 1), cudaMemcpyHostToDevice);
        }
    }
    for(int i = 0; i <= ran.vr.Vz.x; i++) {
        for(int j = 0; j <= ran.vr.Vz.y; j++) {
            cudaMemcpy(&d->va.Vz, &h->va.Vz, sizeof(double) * ( ran.vr.Vz.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->va.Vzx, &h->va.Vzx, sizeof(double) * ( ran.vr.Vz.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->va.Vzy, &h->va.Vzy, sizeof(double) * ( ran.vr.Vz.z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&d->va.Vzz, &h->va.Vzz, sizeof(double) * ( ran.vr.Vz.z + 1), cudaMemcpyHostToDevice);
        }
    }
}

void MemoryBefAftToHost(BefAft *h, BefAft *d, Range ran) {
    cudaMemcpy(&h, &d, sizeof(BefAft), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h->sa, &d->sa, sizeof(SigArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h->ta, &d->ta, sizeof(TauArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h->va, &d->va, sizeof(VelArr), cudaMemcpyDeviceToHost);
    for(int i = 0; i <= ran.sr.Txx.x; i++) {
        for(int j = 0; j <= ran.sr.Txx.y; j++){
            cudaMemcpy(&h->sa.Txx, &d->sa.Txx, sizeof(double) *  (ran.sr.Txx.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Txxx, &d->sa.Txxx, sizeof(double) *  (ran.sr.Txx.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Txxy, &d->sa.Txxy, sizeof(double) *  (ran.sr.Txx.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Txxz, &d->sa.Txxz, sizeof(double) *  (ran.sr.Txx.z + 1), cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i <= ran.sr.Tyy.x; i++) {
        for(int j = 0; j <= ran.sr.Tyy.y; j++){
            cudaMemcpy(&h->sa.Tyy[i][j], &d->sa.Tyy[i][j], sizeof(double) *  (ran.sr.Tyy.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tyyx[i][j], &d->sa.Tyyx[i][j], sizeof(double) *  (ran.sr.Tyy.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tyyy[i][j], &d->sa.Tyyy[i][j], sizeof(double) *  (ran.sr.Tyy.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tyyz[i][j], &d->sa.Tyyz[i][j], sizeof(double) *  (ran.sr.Tyy.z + 1), cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i <= ran.sr.Tzz.x; i++) {
        for(int j = 0; j <= ran.sr.Tzz.y; j++){
            cudaMemcpy(&h->sa.Tzz[i][j], &d->sa.Tzz[i][j], sizeof(double) *  (ran.sr.Tzz.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tzzx[i][j], &d->sa.Tzzx[i][j], sizeof(double) *  (ran.sr.Tzz.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tzzy[i][j], &d->sa.Tzzy[i][j], sizeof(double) *  (ran.sr.Tzz.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tzzz[i][j], &d->sa.Tzzz[i][j], sizeof(double) *  (ran.sr.Tzz.z + 1), cudaMemcpyDeviceToHost);
        }
    }
    //////
    for(int i = 0; i <= ran.tr.Txy.x; i++) {
        for(int j = 0; j <= ran.tr.Txy.y; j++){
            cudaMemcpy(&h->ta.Txy[i][j], &d->ta.Txy[i][j], sizeof(double) *  (ran.tr.Txy.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Txyx[i][j], &d->ta.Txyx[i][j], sizeof(double) *  (ran.tr.Txy.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Txyy[i][j], &d->ta.Txyy[i][j], sizeof(double) *  (ran.tr.Txy.z + 1), cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i <= ran.tr.Tyz.x; i++) {
        for(int j = 0; j <= ran.tr.Tyz.y; j++){
            cudaMemcpy(&h->ta.Tyz[i][j], &d->ta.Tyz[i][j], sizeof(double) *  (ran.tr.Tyz.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Tyzy[i][j], &d->ta.Tyzy[i][j], sizeof(double) *  (ran.tr.Tyz.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Tyzz[i][j], &d->ta.Tyzz[i][j], sizeof(double) *  (ran.tr.Tyz.z + 1), cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i <= ran.tr.Tzx.x; i++) {
        for(int j = 0; j <= ran.tr.Tzx.y; j++){
            cudaMemcpy(&h->ta.Tzx[i][j], &d->ta.Tzx[i][j], sizeof(double) *  (ran.tr.Tzx.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Tzxz[i][j], &d->ta.Tzxz[i][j], sizeof(double) *  (ran.tr.Tzx.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Tzxx[i][j], &d->ta.Tzxx[i][j], sizeof(double) *  (ran.tr.Tzx.z + 1), cudaMemcpyDeviceToHost);
        }
    }
    ///////
    for(int i = 0; i <= ran.vr.Vx.x; i++) {
        for(int j = 0; j <= ran.vr.Vx.y; j++){
            cudaMemcpy(&h->va.Vx, &d->va.Vx, sizeof(double) * ( ran.vr.Vx.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vxx, &d->va.Vxx, sizeof(double) * ( ran.vr.Vx.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vxy, &d->va.Vxy, sizeof(double) * ( ran.vr.Vx.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vxz, &d->va.Vxz, sizeof(double) * ( ran.vr.Vx.z + 1), cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i <= ran.vr.Vy.x; i++) {
        for(int j = 0; j <= ran.vr.Vy.y; j++){
            cudaMemcpy(&h->va.Vy, &d->va.Vy, sizeof(double) * ( ran.vr.Vy.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vyx, &d->va.Vyx, sizeof(double) * ( ran.vr.Vy.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vyy, &d->va.Vyy, sizeof(double) * ( ran.vr.Vy.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vyz, &d->va.Vyz, sizeof(double) * ( ran.vr.Vy.z + 1), cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i <= ran.vr.Vz.x; i++) {
        for(int j = 0; j <= ran.vr.Vz.y; j++) {
            cudaMemcpy(&h->va.Vz, &d->va.Vz, sizeof(double) * ( ran.vr.Vz.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vzx, &d->va.Vzx, sizeof(double) * ( ran.vr.Vz.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vzy, &d->va.Vzy, sizeof(double) * ( ran.vr.Vz.z + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vzz, &d->va.Vzz, sizeof(double) * ( ran.vr.Vz.z + 1), cudaMemcpyDeviceToHost);
        }
    }
}

void onceHtoD(MedArr *ma_h, MedArr *ma_d, Diff *dif_h, Diff *dif_d, Range *ran_h, Range *ran_d) {
    int x = ran_h->sr.Txx.x, y = ran_h->sr.Txx.y, z = ran_h->sr.Txx.z;

    cudaMemcpy(&ma_d, &ma_h, sizeof(MedArr), cudaMemcpyHostToDevice);
    for (int i = 0; i <= x; i++) {
        for (int j = 0; j <= y; j++) {
            cudaMemcpy(&ma_d->ramda[i][j], &ma_h->ramda[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->mu[i][j], &ma_h->mu[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->c11[i][j], &ma_h->c11[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->rho[i][j], &ma_h->rho[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetaxx[i][j], &ma_h->zetaxx[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetaxy[i][j], &ma_h->zetaxy[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetaxz[i][j], &ma_h->zetaxz[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetayx[i][j], &ma_h->zetayx[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetayy[i][j], &ma_h->zetayy[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetayz[i][j], &ma_h->zetayz[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetazx[i][j], &ma_h->zetazx[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetazy[i][j], &ma_h->zetazy[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetazz[i][j], &ma_h->zetazz[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->gamma[i][j], &ma_h->gamma[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->khi[i][j], &ma_h->khi[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->xi11[i][j], &ma_h->xi11[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetadx[i][j], &ma_h->zetadx[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetady[i][j], &ma_h->zetady[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(&ma_d->zetadz[i][j], &ma_h->zetadz[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
        }
    }
    
    cudaMemcpy(&dif_d, &dif_h, sizeof(Diff), cudaMemcpyHostToDevice);
    cudaMemcpy(&dif_d->dt, &dif_h->dt, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&dif_d->dx, &dif_h->dx, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&dif_d->dy, &dif_h->dy, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&dif_d->dz, &dif_h->dz, sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(&ran_d, &ran_h, sizeof(Range), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr, &ran_h->sr, sizeof(SigRan), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr, &ran_h->tr, sizeof(TauRan), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->vr, &ran_h->vr, sizeof(VelRan), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Txx, &ran_h->sr.Txx, sizeof(Coord), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tyy, &ran_h->sr.Tyy, sizeof(Coord), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tzz, &ran_h->sr.Tzz, sizeof(Coord), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Txy, &ran_h->tr.Txy, sizeof(Coord), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tyz, &ran_h->tr.Tyz, sizeof(Coord), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tzx, &ran_h->tr.Tzx, sizeof(Coord), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Txx.x, &ran_h->sr.Txx.x, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Txx.y, &ran_h->sr.Txx.y, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Txx.z, &ran_h->sr.Txx.z, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tyy.x, &ran_h->sr.Tyy.x, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tyy.y, &ran_h->sr.Tyy.y, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tyy.z, &ran_h->sr.Tyy.z, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tzz.x, &ran_h->sr.Tzz.x, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tzz.y, &ran_h->sr.Tzz.y, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tzz.z, &ran_h->sr.Tzz.z, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Txy.x, &ran_h->tr.Txy.x, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Txy.y, &ran_h->tr.Txy.y, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Txy.z, &ran_h->tr.Txy.z, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tyz.x, &ran_h->tr.Tyz.x, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tyz.y, &ran_h->tr.Tyz.y, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tyz.z, &ran_h->tr.Tyz.z, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tzx.x, &ran_h->tr.Tzx.x, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tzx.y, &ran_h->tr.Tzx.y, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tzx.z, &ran_h->tr.Tzx.z, sizeof(int), cudaMemcpyHostToDevice);
}

void MemoryInpaluseToDevice(Inpaluse *ip_h, Inpaluse *ip_d, Range ran) {
    int x = ran.sr.Txx.x, y = ran.sr.Txx.y, z = ran.sr.Txx.z;
    // cudaMemcpy(ip_d, &ip_h, sizeof(Inpaluse), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->freq, &ip_h->freq, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->mode, &ip_h->mode, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->in, &ip_h->in, sizeof(Coord), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->in.x, &ip_h->in.x, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->in.y, &ip_h->in.y, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->in.z, &ip_h->in.z, sizeof(int), cudaMemcpyHostToDevice);
    for(int i = 0; i <= x; i++) {
        for(int j = 0; j <= y; j++){
            cudaMemcpy(ip_d->Txx[i][j], &ip_h->Txx[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(ip_d->Tyy[i][j], &ip_h->Tyy[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
            cudaMemcpy(ip_d->Tzz[i][j], &ip_h->Tzz[i][j], sizeof(double) * (z + 1), cudaMemcpyHostToDevice);
        }
    }
}

void loopDtoH(BefAft *aft_h, BefAft *aft_d, BefAft *bef_h, BefAft *bef_d, Range ran) {
    MemoryBefAftToHost(bef_h, bef_d, ran);
    MemoryBefAftToHost(aft_h, aft_d, ran);
}




void allocate3DArray(double ****array, int x, int y, int z) {
    cudaMalloc((void ***)array, x * sizeof(double **));
    for (int i = 0; i < x; i++) {
        cudaMalloc((void **)&(array[i]), y * sizeof(double *));

        for (int j = 0; j < y; j++) {
            cudaMalloc((void **)&(array[i][j]), z * sizeof(double));

        }
    }
}

void allocateBefAft(BefAft *d_befAft, int x, int y, int z) {
    // メモリをGPUに割り当て
    cudaError_t err = cudaMalloc((void **)&(d_befAft->sa), sizeof(SigArr));
    cudaMalloc((void **)&(d_befAft->ta), sizeof(TauArr));
    cudaMalloc((void **)&(d_befAft->va), sizeof(VelArr));
    if (err != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
    }
    // SigArr のメンバのメモリ確保
    allocate3DArray(&d_befAft->sa.Txx, x, y, z);
    allocate3DArray(&d_befAft->sa.Txxx, x, y, z);
    allocate3DArray(&d_befAft->sa.Txxy, x, y, z);
    allocate3DArray(&d_befAft->sa.Txxz, x, y, z);
    allocate3DArray(&d_befAft->sa.Tyy, x, y, z);
    allocate3DArray(&d_befAft->sa.Tyyx, x, y, z);
    allocate3DArray(&d_befAft->sa.Tyyy, x, y, z);
    allocate3DArray(&d_befAft->sa.Tyyz, x, y, z);
    allocate3DArray(&d_befAft->sa.Tzz, x, y, z);
    allocate3DArray(&d_befAft->sa.Tzzx, x, y, z);
    allocate3DArray(&d_befAft->sa.Tzzy, x, y, z);
    allocate3DArray(&d_befAft->sa.Tzzz, x, y, z);

    // TauArr のメンバのメモリ確保
    allocate3DArray(&d_befAft->ta.Txy, x, y, z);
    allocate3DArray(&d_befAft->ta.Txyx, x, y, z);
    allocate3DArray(&d_befAft->ta.Txyy, x, y, z);
    allocate3DArray(&d_befAft->ta.Tyz, x, y, z);
    allocate3DArray(&d_befAft->ta.Tyzy, x, y, z);
    allocate3DArray(&d_befAft->ta.Tyzz, x, y, z);
    allocate3DArray(&d_befAft->ta.Tzx, x, y, z);
    allocate3DArray(&d_befAft->ta.Tzxz, x, y, z);
    allocate3DArray(&d_befAft->ta.Tzxx, x, y, z);

    // VelArr のメンバのメモリ確保
    allocate3DArray(&d_befAft->va.Vx, x, y, z);
    allocate3DArray(&d_befAft->va.Vxx, x, y, z);
    allocate3DArray(&d_befAft->va.Vxy, x, y, z);
    allocate3DArray(&d_befAft->va.Vxz, x, y, z);
    allocate3DArray(&d_befAft->va.Vy, x, y, z);
    allocate3DArray(&d_befAft->va.Vyx, x, y, z);
    allocate3DArray(&d_befAft->va.Vyy, x, y, z);
    allocate3DArray(&d_befAft->va.Vyz, x, y, z);
    allocate3DArray(&d_befAft->va.Vz, x, y, z);
    allocate3DArray(&d_befAft->va.Vzx, x, y, z);
    allocate3DArray(&d_befAft->va.Vzy, x, y, z);
    allocate3DArray(&d_befAft->va.Vzz, x, y, z);
}

void copy3DArrayToDevice(double ***d_array, double ***h_array, int x, int y, int z) {
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            cudaMemcpy(d_array[i][j], h_array[i][j], z * sizeof(double), cudaMemcpyHostToDevice);
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
    printf("no\n");
    // SigArr のメンバのデータコピー
    copy3DArrayToDevice(d_befAft->sa.Txx,  h_befAft->sa.Txx , Txxi, Txxj, Txxk);
    printf("oooooo\n");
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
    printf("okkk\n");
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