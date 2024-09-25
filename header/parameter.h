#pragma once
#include "./struct.h"

__global__ void StaticVariable(Medium *med, Pml *pml, Range *ran, Diff *dif, Object *air, Object *con, Object *clack, MedArr *ma, int *tmax, int *outNum);
__global__ void DynamicVariable(BefAft *bef, BefAft *aft, AccCoord *acc, MedArr *ma, Impulse *ip, Range *ran, Medium *med, Object *air, Object *con, Object *clack, Pml *pml, Diff *dif, Coord *out, int *outNum);
__device__ void insertBefAft(BefAft *ba, Range *ran);
__device__ void insertAccCoord(AccCoord *acc, int *outNum);
__device__ void insertAir(MedArr *ma, Object *air, Range *ran);
__device__ void insertConcrete(MedArr *ma, Object *con, Range *ran);
__device__ void insertClack(MedArr *ma, Object *clack, Range *ran);
__device__ void insertPml(MedArr *ma, Pml *pml, Range *ran);
__global__ void insertImpulse(Impulse *ip, Diff *dif, int *t, Range *ran);