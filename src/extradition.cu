#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../header/init.h"
#include "../header/insert.h"
#include "../header/struct.h"
#include "../header/update.h"
#include "../header/para_in.h"

void onceHtoD(MedArr *ma_h, MedArr *ma_d, Diff *dif_h, Diff *dif_d, Range *ran_h, Range *ran_d) {
    int num_elements = ran_h->sr.Txx.x * ran_h->sr.Txx.y * ran_h->sr.Txx.z;
    cudaMemcpy(&ma_d, &ma_h, sizeof(MedArr), cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->ramda, &ma_h->ramda, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->mu, &ma_h->mu, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->c11, &ma_h->c11, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->rho, &ma_h->rho, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetaxx, &ma_h->zetaxx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetaxy, &ma_h->zetaxy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetaxz, &ma_h->zetaxz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetayx, &ma_h->zetayx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetayy, &ma_h->zetayy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetayz, &ma_h->zetayz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetazx, &ma_h->zetazx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetazy, &ma_h->zetazy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetazz, &ma_h->zetazz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->gamma, &ma_h->gamma, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->khi, &ma_h->khi, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->xi11, &ma_h->xi11, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetadx, &ma_h->zetadx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetady, &ma_h->zetady, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ma_d->zetadz, &ma_h->zetadz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);

    cudaMemcpy(&dif_d, &dif_h, sizeof(Diff), cudaMemcpyHostToDevice);
    cudaMemcpy(&dif_d->dt, &dif_h->dt, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&dif_d->dx, &dif_h->dx, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&dif_d->dy, &dif_h->dy, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&dif_d->dz, &dif_h->dz, sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(&ran_d, &ran_h, sizeof(Range), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr, &ran_h->sr, sizeof(SigRan), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr, &ran_h->tr, sizeof(TauRan), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->vr, &ran_h->vr, sizeof(VelRan), cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Txx, &ran_h->sr.Txx, sizeof(Coord) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tyy, &ran_h->sr.Tyy, sizeof(Coord) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tzz, &ran_h->sr.Tzz, sizeof(Coord) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Txy, &ran_h->tr.Txy, sizeof(Coord) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tyz, &ran_h->tr.Tyz, sizeof(Coord) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tzx, &ran_h->tr.Tzx, sizeof(Coord) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Txx.x, &ran_h->sr.Txx.x, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Txx.y, &ran_h->sr.Txx.y, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Txx.z, &ran_h->sr.Txx.z, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tyy.x, &ran_h->sr.Tyy.x, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tyy.y, &ran_h->sr.Tyy.y, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tyy.z, &ran_h->sr.Tyy.z, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tzz.x, &ran_h->sr.Tzz.x, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tzz.y, &ran_h->sr.Tzz.y, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->sr.Tzz.z, &ran_h->sr.Tzz.z, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Txy.x, &ran_h->tr.Txy.x, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Txy.y, &ran_h->tr.Txy.y, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Txy.z, &ran_h->tr.Txy.z, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tyz.x, &ran_h->tr.Tyz.x, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tyz.y, &ran_h->tr.Tyz.y, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tyz.z, &ran_h->tr.Tyz.z, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tzx.x, &ran_h->tr.Tzx.x, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tzx.y, &ran_h->tr.Tzx.y, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ran_d->tr.Tzx.z, &ran_h->tr.Tzx.z, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
}

void loopHtoD(Inpaluse *ip_h, Inpaluse *ip_d, BefAft *aft_h, BefAft *aft_d, BefAft *bef_h, BefAft *bef_d, Range ran) {
    int num_elements = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    cudaMemcpy(ip_d, &ip_h, sizeof(Inpaluse), cudaMemcpyHostToDevice);
    cudaMemcpy(ip_d->Txx, &ip_h->Txx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(ip_d->Tyy, &ip_h->Tyy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(ip_d->Tzz, &ip_h->Tzz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->freq, &ip_h->freq, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->mode, &ip_h->mode, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->in, &ip_h->in, sizeof(Coord), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->in.x, &ip_h->in.x, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->in.y, &ip_h->in.y, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&ip_d->in.z, &ip_h->in.z, sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(&aft_d, &aft_h, sizeof(BefAft), cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->sa, &aft_h->sa, sizeof(SigArr), cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->ta, &aft_h->ta, sizeof(TauArr), cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va, &aft_h->va, sizeof(VelArr), cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Txx, aft_h->sa.Txx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Txxx, aft_h->sa.Txxx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Txxy, aft_h->sa.Txxy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Txxz, aft_h->sa.Txxz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Tyy, aft_h->sa.Tyy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Tyyx, aft_h->sa.Tyyx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Tyyy, aft_h->sa.Tyyy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Tyyz, aft_h->sa.Tyyz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Tzz, aft_h->sa.Tzz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Tzzx, aft_h->sa.Tzzx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Tzzy, aft_h->sa.Tzzy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(aft_d->sa.Tzzz, aft_h->sa.Tzzz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    
    cudaMemcpy(&aft_d->ta.Txy, &aft_h->ta.Txy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->ta.Txyx, &aft_h->ta.Txyx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->ta.Txyy, &aft_h->ta.Txyy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->ta.Tyz, &aft_h->ta.Tyz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->ta.Tyzy, &aft_h->ta.Tyzy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->ta.Tyzz, &aft_h->ta.Tyzz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->ta.Tzx, &aft_h->ta.Tzx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->ta.Tzxz, &aft_h->ta.Tzxz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->ta.Tzxx, &aft_h->ta.Tzxx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);

    cudaMemcpy(&aft_d->va.Vx, &aft_h->va.Vx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vxx, &aft_h->va.Vxx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vxy, &aft_h->va.Vxy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vxz, &aft_h->va.Vxz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vy, &aft_h->va.Vy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vyx, &aft_h->va.Vyx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vyy, &aft_h->va.Vyy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vyz, &aft_h->va.Vyz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vz, &aft_h->va.Vz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vzx, &aft_h->va.Vzx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vzy, &aft_h->va.Vzy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&aft_d->va.Vzz, &aft_h->va.Vzz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);

    cudaMemcpy(&bef_d, &bef_h, sizeof(BefAft), cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->sa, &bef_h->sa, sizeof(SigArr), cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->ta, &bef_h->ta, sizeof(TauArr), cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va, &bef_h->va, sizeof(VelArr), cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Txx, bef_h->sa.Txx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Txxx, bef_h->sa.Txxx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Txxy, bef_h->sa.Txxy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Txxz, bef_h->sa.Txxz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Tyy, bef_h->sa.Tyy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Tyyx, bef_h->sa.Tyyx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Tyyy, bef_h->sa.Tyyy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Tyyz, bef_h->sa.Tyyz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Tzz, bef_h->sa.Tzz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Tzzx, bef_h->sa.Tzzx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Tzzy, bef_h->sa.Tzzy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(bef_d->sa.Tzzz, bef_h->sa.Tzzz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    
    cudaMemcpy(&bef_d->ta.Txy, &bef_h->ta.Txy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->ta.Txyx, &bef_h->ta.Txyx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->ta.Txyy, &bef_h->ta.Txyy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->ta.Tyz, &bef_h->ta.Tyz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->ta.Tyzy, &bef_h->ta.Tyzy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->ta.Tyzz, &bef_h->ta.Tyzz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->ta.Tzx, &bef_h->ta.Tzx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->ta.Tzxz, &bef_h->ta.Tzxz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->ta.Tzxx, &bef_h->ta.Tzxx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);

    cudaMemcpy(&bef_d->va.Vx, &bef_h->va.Vx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vxx, &bef_h->va.Vxx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vxy, &bef_h->va.Vxy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vxz, &bef_h->va.Vxz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vy, &bef_h->va.Vy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vyx, &bef_h->va.Vyx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vyy, &bef_h->va.Vyy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vyz, &bef_h->va.Vyz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vz, &bef_h->va.Vz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vzx, &bef_h->va.Vzx, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vzy, &bef_h->va.Vzy, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
    cudaMemcpy(&bef_d->va.Vzz, &bef_h->va.Vzz, sizeof(double ***) * num_elements, cudaMemcpyHostToDevice);
}

void loopDtoH(BefAft *aft_h, BefAft *aft_d, BefAft *bef_h, BefAft *bef_d, Range ran) {
    int num_elements = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    cudaMemcpy(&aft_h, &aft_d, sizeof(BefAft), cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->sa, &aft_d->sa, sizeof(SigArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->ta, &aft_d->ta, sizeof(TauArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va, &aft_d->va, sizeof(VelArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Txx, aft_d->sa.Txx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Txxx, aft_d->sa.Txxx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Txxy, aft_d->sa.Txxy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Txxz, aft_d->sa.Txxz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Tyy, aft_d->sa.Tyy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Tyyx, aft_d->sa.Tyyx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Tyyy, aft_d->sa.Tyyy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Tyyz, aft_d->sa.Tyyz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Tzz, aft_d->sa.Tzz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Tzzx, aft_d->sa.Tzzx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Tzzy, aft_d->sa.Tzzy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(aft_h->sa.Tzzz, aft_d->sa.Tzzz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    
    cudaMemcpy(&aft_h->ta.Txy, &aft_d->ta.Txy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->ta.Txyx, &aft_d->ta.Txyx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->ta.Txyy, &aft_d->ta.Txyy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->ta.Tyz, &aft_d->ta.Tyz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->ta.Tyzy, &aft_d->ta.Tyzy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->ta.Tyzz, &aft_d->ta.Tyzz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->ta.Tzx, &aft_d->ta.Tzx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->ta.Tzxz, &aft_d->ta.Tzxz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->ta.Tzxx, &aft_d->ta.Tzxx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);

    cudaMemcpy(&aft_h->va.Vx, &aft_d->va.Vx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vxx, &aft_d->va.Vxx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vxy, &aft_d->va.Vxy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vxz, &aft_d->va.Vxz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vy, &aft_d->va.Vy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vyx, &aft_d->va.Vyx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vyy, &aft_d->va.Vyy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vyz, &aft_d->va.Vyz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vz, &aft_d->va.Vz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vzx, &aft_d->va.Vzx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vzy, &aft_d->va.Vzy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&aft_h->va.Vzz, &aft_d->va.Vzz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);

    cudaMemcpy(&bef_h, &bef_d, sizeof(BefAft), cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->sa, &bef_d->sa, sizeof(SigArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->ta, &bef_d->ta, sizeof(TauArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va, &bef_d->va, sizeof(VelArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Txx, bef_d->sa.Txx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Txxx, bef_d->sa.Txxx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Txxy, bef_d->sa.Txxy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Txxz, bef_d->sa.Txxz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Tyy, bef_d->sa.Tyy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Tyyx, bef_d->sa.Tyyx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Tyyy, bef_d->sa.Tyyy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Tyyz, bef_d->sa.Tyyz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Tzz, bef_d->sa.Tzz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Tzzx, bef_d->sa.Tzzx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Tzzy, bef_d->sa.Tzzy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(bef_h->sa.Tzzz, bef_d->sa.Tzzz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    
    cudaMemcpy(&bef_h->ta.Txy, &bef_d->ta.Txy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->ta.Txyx, &bef_d->ta.Txyx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->ta.Txyy, &bef_d->ta.Txyy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->ta.Tyz, &bef_d->ta.Tyz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->ta.Tyzy, &bef_d->ta.Tyzy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->ta.Tyzz, &bef_d->ta.Tyzz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->ta.Tzx, &bef_d->ta.Tzx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->ta.Tzxz, &bef_d->ta.Tzxz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->ta.Tzxx, &bef_d->ta.Tzxx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);

    cudaMemcpy(&bef_h->va.Vx, &bef_d->va.Vx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vxx, &bef_d->va.Vxx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vxy, &bef_d->va.Vxy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vxz, &bef_d->va.Vxz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vy, &bef_d->va.Vy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vyx, &bef_d->va.Vyx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vyy, &bef_d->va.Vyy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vyz, &bef_d->va.Vyz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vz, &bef_d->va.Vz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vzx, &bef_d->va.Vzx, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vzy, &bef_d->va.Vzy, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
    cudaMemcpy(&bef_h->va.Vzz, &bef_d->va.Vzz, sizeof(double ***) * num_elements, cudaMemcpyDeviceToHost);
}