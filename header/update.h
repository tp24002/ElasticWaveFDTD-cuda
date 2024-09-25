#pragma once
#include "./struct.h"

// 垂直応力

__global__ void TxxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Coord Txx);
__global__ void TyyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Coord Tyy);
__global__ void TzzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Coord Tzz);
__global__ void ZeroT_XY(BefAft *aft, Coord ranmax, char check);
__global__ void ZeroT_YZ(BefAft *aft, Coord ranmax, char check);
__global__ void ZeroT_ZX(BefAft *aft, Coord ranmax, char check);
__global__ void DirectionalAdd(BefAft *aft, Impulse *ip, Coord ranmax, char check);
__device__ void Txx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Impulse *ip_d, Coord threads);
__device__ void Tyy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Impulse *ip_d, Coord threads);
__device__ void Tzz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Impulse *ip_d, Coord threads);

__global__ void Sig(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Impulse *ip_d, Coord threads);

// せん断応力

__global__ void TxyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr);
__global__ void TyzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr);
__global__ void TzxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr);
__global__ void ZeroTxy(BefAft *aft, Coord Tmax);
__global__ void ZeroTyz(BefAft *aft, Coord Tmax);
__global__ void ZeroTzx(BefAft *aft, Coord Tmax);
__global__ void DirectionalAddT(BefAft *aft, Coord Tmax, char check);
__device__ void Txy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads);
__device__ void Tyz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads);
__device__ void Tzx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads);

__global__ void Tau(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads);

// 粒子速度

__global__ void VxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, VelRan vr);
__global__ void VyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, VelRan vr);
__global__ void VzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, VelRan vr);
__global__ void ZeroVx_XY(BefAft *aft, Coord Vmax);
__global__ void ZeroVx_XZ(BefAft *aft, Coord Vmax);
__global__ void ZeroVy_YX(BefAft *aft, Coord Vmax);
__global__ void ZeroVy_YZ(BefAft *aft, Coord Vmax);
__global__ void ZeroVz_ZX(BefAft *aft, Coord Vmax);
__global__ void ZeroVz_ZY(BefAft *aft, Coord Vmax);
__global__ void DirectionalAddV(BefAft *aft, Coord Vmax, char check);
__device__ void Vx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads);
__device__ void Vy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads);
__device__ void Vz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads);

__global__ void Vel(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran, Coord threads);


__global__ void AccelerationCalculation(AccCoord *Acc, BefAft *aft, BefAft *bef, Diff *dif, Coord *out, Range *ran, int *outNum);


__global__ void swapTxx(SigArr *aftSa, SigArr *befSa, Range *ran);
__global__ void swapTyy(SigArr *aftSa, SigArr *befSa, Range *ran);
__global__ void swapTzz(SigArr *aftSa, SigArr *befSa, Range *ran);
__global__ void swapTxy(TauArr *aftTa, TauArr *befTa, Range *ran);
__global__ void swapTyz(TauArr *aftTa, TauArr *befTa, Range *ran);
__global__ void swapTzx(TauArr *aftTa, TauArr *befTa, Range *ran);
__global__ void swapVx(VelArr *aftVa, VelArr *befVa, Range *ran);
__global__ void swapVy(VelArr *aftVa, VelArr *befVa, Range *ran);
__global__ void swapVz(VelArr *aftVa, VelArr *befVa, Range *ran);
__global__ void swapBefAft(BefAft *aft, BefAft *bef, Range *ran, Coord threads);