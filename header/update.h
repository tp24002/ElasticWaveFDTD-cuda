#pragma once
#include "./struct.h"

// 垂直応力

__global__ void TxxUpdate(SigArr aftsa, VelArr aftva, SigArr befsa, VelArr befva, MedArr *ma, Diff *dif, Range *ran);
__global__ void TyyUpdate(SigArr aftsa, VelArr aftva, SigArr befsa, VelArr befva, MedArr *ma, Diff *dif, Range *ran);
__global__ void TzzUpdate(SigArr aftsa, VelArr aftva, SigArr befsa, VelArr befva, MedArr *ma, Diff *dif, Range *ran);
__global__ void DirectionalAdd(SigArr aftsa, ImpulseArr *ipa, Range *ran, char check);
__global__ void createImpulse(ImpulseArr *ipa, Impulse *ip, Diff *dif, Range *ran, int innum, int t);
void Txx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, ImpulseArr *ipa_d, DimI3 threads);
void Tyy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, ImpulseArr *ipa_d, DimI3 threads);
void Tzz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, ImpulseArr *ipa_d, DimI3 threads);

void Sig(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, ImpulseArr *ipa_d, DimI3 threads);

// せん断応力

__global__ void TxyUpdate(TauArr aftta, VelArr aftva, TauArr befta, VelArr befva, MedArr *ma, Diff *dif, Range *ran);
__global__ void TyzUpdate(TauArr aftta, VelArr aftva, TauArr befta, VelArr befva, MedArr *ma, Diff *dif, Range *ran);
__global__ void TzxUpdate(TauArr aftta, VelArr aftva, TauArr befta, VelArr befva, MedArr *ma, Diff *dif, Range *ran);
__global__ void DirectionalAddT(TauArr aftta, Range *ran, char check);
void Txy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads);
void Tyz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads);
void Tzx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads);

void Tau(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads);

// 粒子速度

__global__ void VxUpdate(VelArr aftva, VelArr befva, SigArr befsa, TauArr befta, MedArr *ma, Diff *dif, Range *ran);
__global__ void VyUpdate(VelArr aftva, VelArr befva, SigArr befsa, TauArr befta, MedArr *ma, Diff *dif, Range *ran);
__global__ void VzUpdate(VelArr aftva, VelArr befva, SigArr befsa, TauArr befta, MedArr *ma, Diff *dif, Range *ran);
__global__ void DirectionalAddV(VelArr aftva, Range *ran, char check);
void Vx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads);
void Vy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads);
void Vz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads);

void Vel(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, DimI3 threads);


__global__ void Acceleration(DimD3 *acc, VelArr aftva, VelArr befva, Diff *dif, DimI3 *out, Range *ran, int outnum);


__global__ void swapTxx(SigArr aftsa, SigArr befsa, Range *ran);
__global__ void swapTyy(SigArr aftsa, SigArr befsa, Range *ran);
__global__ void swapTzz(SigArr aftsa, SigArr befsa, Range *ran);
__global__ void swapTxy(TauArr aftta, TauArr befta, Range *ran);
__global__ void swapTyz(TauArr aftta, TauArr befta, Range *ran);
__global__ void swapTzx(TauArr aftta, TauArr befta, Range *ran);
__global__ void swapVx(VelArr aftva, VelArr befva, Range *ran);
__global__ void swapVy(VelArr aftva, VelArr befva, Range *ran);
__global__ void swapVz(VelArr aftva, VelArr befva, Range *ran);
void swapBefAft(BefAft *aft, BefAft *bef, Range *ran_h, Range *ran_d, DimI3 threads);