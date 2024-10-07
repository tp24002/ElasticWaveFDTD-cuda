#pragma once
#include "./struct.h"

void initCoord(Coord *co, int x, int y, int z);
void initMedium(Medium *med);
void initDiff(Diff *dif, Medium *med);
void initPml(Pml *pml, Medium *med, Diff dif);
void initRange(Range *ran, Coord region, Pml pml);
void initRandom(Object con, Coord *clack, int ratio);
__global__ void ZeroT(BefAft *ba, Range *ran);
__global__ void ZeroTxy(BefAft *ba, Range *ran);
__global__ void ZeroTyz(BefAft *ba, Range *ran);
__global__ void ZeroTzx(BefAft *ba, Range *ran);
__global__ void ZeroVx(BefAft *ba, Range *ran);
__global__ void ZeroVy(BefAft *ba, Range *ran);
__global__ void ZeroVz(BefAft *ba, Range *ran);
