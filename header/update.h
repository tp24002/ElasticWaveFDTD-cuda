#pragma once
#include "./struct.h"

// 垂直応力

__global__ void TxxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran);
__global__ void TyyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran);
__global__ void TzzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran);
__global__ void DirectionalAdd(BefAft *aft, Impulse *ip, Range *ran, char check);
void Txx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Impulse *ip_d, Coord threads);
void Tyy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Impulse *ip_d, Coord threads);
void Tzz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Impulse *ip_d, Coord threads);

void Sig(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Impulse *ip_d, Coord threads);

// せん断応力

__global__ void TxyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran);
__global__ void TyzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran);
__global__ void TzxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran);
__global__ void DirectionalAddT(BefAft *aft, Range *ran, char check);
void Txy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads);
void Tyz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads);
void Tzx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads);

void Tau(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads);

// 粒子速度

__global__ void VxUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran);
__global__ void VyUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran);
__global__ void VzUpdate(BefAft *aft, BefAft *bef, MedArr *ma, Diff *dif, Range *ran);
__global__ void DirectionalAddV(BefAft *aft, Range *ran, char check);
void Vx(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads);
void Vy(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads);
void Vz(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads);

void Vel(BefAft *aft_d, BefAft *bef_d, MedArr *ma_d, Diff *dif_d, Range *ran_d, Range *ran_h, Coord threads);


// __global__ void AccelerationCalculation(AccCoord *Acc, BefAft *aft, BefAft *bef, Diff *dif, Coord *out, Range *ran, int *outNum);


__global__ void swapTxx(BefAft *aft, BefAft *bef, Range *ran);
__global__ void swapTyy(BefAft *aft, BefAft *bef, Range *ran);
__global__ void swapTzz(BefAft *aft, BefAft *bef, Range *ran);
__global__ void swapTxy(BefAft *aft, BefAft *bef, Range *ran);
__global__ void swapTyz(BefAft *aft, BefAft *bef, Range *ran);
__global__ void swapTzx(BefAft *aft, BefAft *bef, Range *ran);
__global__ void swapVx(BefAft *aft, BefAft *bef, Range *ran);
__global__ void swapVy(BefAft *aft, BefAft *bef, Range *ran);
__global__ void swapVz(BefAft *aft, BefAft *bef, Range *ran);
void swapBefAft(BefAft *aft, BefAft *bef, Range *ran_h, Range *ran_d, Coord threads);