#pragma once
#include "./struct.h"

void StaticVariable(Medium *med, Pml *pml, Range *ran, Diff *dif, Object *air, Object *con, Object *clack, int *tmax, int *outnum, int *innum);
void DynamicVariable(DimD3 *acc, MedArr *ma, Impulse *ip, Range ran, Object air, Object con, Object clack, Pml pml, DimI3 *out, int outnum);
void insertDimD3(DimD3 *acc, int outnum);
void insertAir(MedArr *ma, Object air, Range ran);
void insertConcrete(MedArr *ma, Object con, Range ran);
void insertClack(MedArr *ma, Object clack, Range ran);
void insertPml(MedArr *ma, Pml pml, Range ran);
void insertImpulse(Impulse *ip, Diff dif, int t, Range ran, double freq, int mode, int innum);