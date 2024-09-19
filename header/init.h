#pragma once
#include "./struct.h"

void initMedium(Medium *med);
void initCoord(Coord *co, int x, int y, int z);
void initDiff(Diff *dif, Medium *med);
void initPml(Pml *pml, Medium *med, Diff dif);
void initRange(Range *ran, Coord region, Pml pml);
void initrandom(Coord con_size, Coord *clack, int ratio);