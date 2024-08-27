#pragma once
#include "./struct.h"

// 垂直応力

__device__ void TxxUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, int imax, int jmax, int kmax);
__device__ void TyyUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, int imax, int jmax, int kmax);
__device__ void TzzUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, int imax, int jmax, int kmax);
__device__ void ZeroT_XY(BefAft *aft, int imax, int jmax, int kmax, char check);
__device__ void ZeroT_YZ(BefAft *aft, int imax, int jmax, int kmax, char check);
__device__ void ZeroT_ZX(BefAft *aft, int imax, int jmax, int kmax, char check);
__device__ void DirectionalAdd(BefAft *aft, Inpaluse ip, int imax, int jmax, int kmax, char check);
__device__ void Txx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);
__device__ void Tyy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);
__device__ void Tzz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);

__global__ void Sig(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);

// せん断応力

__device__ void TxyUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__device__ void TyzUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__device__ void TzxUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__device__ void Txy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__device__ void Tyz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__device__ void Tzx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);

__global__ void Tau(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);

// 粒子速度

__device__ void VxUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__device__ void VyUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__device__ void VzUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__device__ void Vx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__device__ void Vy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__device__ void Vz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);

__global__ void Vel(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);


void Acc(Coord_acc *A,BefAft *aft, BefAft *bef, Diff dif, Coord out);
void swapBefAft(BefAft *aft, BefAft *bef, Range ran);