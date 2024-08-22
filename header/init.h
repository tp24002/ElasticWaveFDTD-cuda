#pragma once
#include "./struct.h"

void initMedium(Medium *med);
void initCoord(Coord *co, int x, int y, int z);
void initDiff(Diff *dif, Medium *med);
void initPml(Pml *pml, Medium *med, Diff dif);
void initConcrete(Object *con, Medium med, Pml pml, int spx, int spy, int spz, int ranx, int rany, int ranz);
void initSteel(Object *steel, Medium med, Pml *pml, int spx, int spy, int spz, int ranx, int rany, int ranz);
void initRange(Range *ran, int x, int y, int z, Pml pml);
void initSigArr(SigArr *sa, SigRan sr);
void initTauArr(TauArr *ta, TauRan tr);
void initVelArr(VelArr *va, VelRan vr);
void initBefAft(BefAft *ba, Range ran);
void initInpalse(Inpaluse *ip, SigRan sr, Pml pml, int mode, int x, int y, int z, double freq);
void initMedArr(MedArr *ma, SigRan sr);
void initClack(Object *clack, Medium med, Pml *pml, int spx, int spy, int spz, int ranx, int rany, int ranz);