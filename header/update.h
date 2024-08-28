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
void Txx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);
void Tyy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);
void Tzz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);

void Sig(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);

// せん断応力

__global__ void TxyUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__global__ void TyzUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
__global__ void TzxUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
void Txy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
void Tyz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
void Tzx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);

void Tau(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);

// 粒子速度

__global__ void VxUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__global__ void VyUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
__global__ void VzUpdate(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
void Vx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
void Vy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
void Vz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);

void Vel(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);


void Acc(Coord_acc *A,BefAft *aft, BefAft *bef, Diff dif, Coord out);
void swapBefAft(BefAft *aft, BefAft *bef, Range ran);