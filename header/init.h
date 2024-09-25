#pragma once
#include "./struct.h"

__host__ void initHostCoord(Coord *co, int x, int y, int z);
__device__ void initDeviceCoord(Coord *co, int x, int y, int z);
__global__ void initMedium(Medium *med);
__global__ void initDiff(Diff *dif, Medium *med);
__global__ void initPml(Pml *pml, Medium *med, Diff dif);
__global__ void initRange(Range *ran, Coord region, Pml pml);
void initRandom(Coord con_size, Coord *clack, int ratio);