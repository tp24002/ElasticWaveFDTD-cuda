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
void Txx(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Impulse ip_h, Impulse *ip_d, int t, Coord threads);
void Tyy(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Impulse ip_h, Impulse *ip_d, int t, Coord threads);
void Tzz(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Impulse ip_h, Impulse *ip_d, int t, Coord threads);

void Sig(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Impulse ip_h, Impulse *ip_d, int t, Coord threads);

// せん断応力

__global__ void TxyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr);
__global__ void TyzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr);
__global__ void TzxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, TauRan tr);
__global__ void ZeroTxy(BefAft *aft, Coord Tmax);
__global__ void ZeroTyz(BefAft *aft, Coord Tmax);
__global__ void ZeroTzx(BefAft *aft, Coord Tmax);
__global__ void DirectionalAddT(BefAft *aft, Coord Tmax, char check);
void Txy(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads);
void Tyz(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads);
void Tzx(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads);

void Tau(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads);

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
void Vx(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads);
void Vy(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads);
void Vz(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads);

void Vel(BefAft *aft_h, BefAft *bef_h, BefAft *aft_d, BefAft *bef_d, MedArr ma_h, MedArr *ma_d, Diff dif_h, Diff *dif_d, Range ran, Coord threads);


__global__ void Acceleration(Coord_acc **Acc,BefAft *aft, BefAft *bef, Diff dif, Coord *out,int outnum, int t);
__global__ void swapBefAft(BefAft *aft, BefAft *bef, Range ran);