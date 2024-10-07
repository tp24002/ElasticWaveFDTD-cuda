#pragma once
#include "./struct.h"

void StaticVariable(Medium *med, Pml *pml, Range *ran, Diff *dif, Object *air, Object *con, Object *clack, int *tmax, int *outNum);
void DynamicVariable(AccCoord *acc, MedArr *ma, Impulse *ip, Range ran, Object air, Object con, Object clack, Pml pml, Coord *out, int outNum);
void insertAccCoord(AccCoord *acc, int outNum);
void insertAir(MedArr *ma, Object air, Range ran);
void insertConcrete(MedArr *ma, Object con, Range ran);
void insertClack(MedArr *ma, Object clack, Range ran);
void insertPml(MedArr *ma, Pml pml, Range ran);
void insertImpulse(Impulse *ip, Diff dif, int t, Range ran);