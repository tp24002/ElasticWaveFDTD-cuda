#pragma once
#include "./struct.h"

// 垂直応力

__global__ void TxxUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, int imax, int jmax, int kmax);
__global__ void TyyUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, int imax, int jmax, int kmax);
__global__ void TzzUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, int imax, int jmax, int kmax);
__global__ void ZeroT_XY(BefAft *aft, int imax, int jmax, int kmax, char check);
__global__ void ZeroT_YZ(BefAft *aft, int imax, int jmax, int kmax, char check);
__global__ void ZeroT_ZX(BefAft *aft, int imax, int jmax, int kmax, char check);
__global__ void DirectionalAdd(BefAft *aft, Inpaluse ip, int imax, int jmax, int kmax, char check);
void Txx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t, Coord threads);
void Tyy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t, Coord threads);
void Tzz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t, Coord threads);

void Sig(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t, Coord threads);

// せん断応力

__global__ void TxyUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__global__ void TyzUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__global__ void TzxUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__global__ void ZeroTxy(BefAft *aft, Coord Tmax);
__global__ void ZeroTyz(BefAft *aft, Coord Tmax);
__global__ void ZeroTzx(BefAft *aft, Coord Tmax);
__global__ void DirectionalAddT(BefAft *aft, Coord Tmax, char check);
void Txy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr, Coord threads);
void Tyz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr, Coord threads);
void Tzx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr, Coord threads);

void Tau(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr, Coord threads);

// 粒子速度

__global__ void VxUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__global__ void VyUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__global__ void VzUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__global__ void ZeroVx_XY(BefAft *aft, Coord Vmax);
__global__ void ZeroVx_XZ(BefAft *aft, Coord Vmax);
__global__ void ZeroVx_YX(BefAft *aft, Coord Vmax);
__global__ void ZeroVx_YZ(BefAft *aft, Coord Vmax);
__global__ void ZeroVx_ZX(BefAft *aft, Coord Vmax);
__global__ void ZeroVx_ZY(BefAft *aft, Coord Vmax);
__global__ void DirectionalAddV(BefAft *aft, Coord Vmax, char check);
void Vx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr, Coord threads);
void Vy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr, Coord threads);
void Vz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr, Coord threads);

void Vel(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr, Coord threads);


void Acc(Coord_acc *A,BefAft *aft, BefAft *bef, Diff dif, Coord out);
void swapBefAft(BefAft *aft, BefAft *bef, Range ran);