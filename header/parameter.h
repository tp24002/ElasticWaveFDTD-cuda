#pragma once
#include "./struct.h"


void StaticVariable(Range *ran, Pml *pml, Diff *dif, Object *con, Object *clack, Medium *med, MedArr *ma, int *tmax);
void DynamicVariable(MedArr *ma, Impulse *ip, Range ran, Medium *med, Object con, Object clack, Pml pml, Diff dif, int tmax);
void insertAir(MedArr *ma, SigRan sr, Medium air);
void insertConcrete(MedArr *ma, Object con);
void insertClack(MedArr *ma, Object clack);
void insertPml(MedArr *ma, SigRan sr, Pml pml);
void insertImpulse(Impulse *ip, Diff dif, int tmax);
