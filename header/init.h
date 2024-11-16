#pragma once
#include "./struct.h"


void initDimI3(DimI3 *co, int x, int y, int z);
void initMedium(Medium *med);
void initDiff(Diff *dif, Medium *med);
void initPml(Pml *pml, Medium *med, Diff dif);
void initRange(Range *ran, DimI3 region, Pml pml);
void initRandom(Object con, DimI3 *clack, int ratio);
__global__ void ZeroT(SigArr sa, Range *ran);
__global__ void ZeroTxy(TauArr ta, Range *ran);
__global__ void ZeroTyz(TauArr ta, Range *ran);
__global__ void ZeroTzx(TauArr ta, Range *ran);
__global__ void ZeroVx(VelArr va, Range *ran);
__global__ void ZeroVy(VelArr va, Range *ran);
__global__ void ZeroVz(VelArr va, Range *ran);
