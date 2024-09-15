#pragma once
#include "./struct.h"

void initMedium(Medium *med);
void initCoord(Coord *co, int x, int y, int z);
void initDiff(Diff *dif, Medium *med);
void initPml(Pml *pml, Medium *med, Diff dif);
void initConcrete(Object *con, Medium med, Pml pml, int spx, int spy, int spz, int ranx, int rany, int ranz);
void initRange(Range *ran, int x, int y, int z, Pml pml);
void initHostSigArr(SigArr *sa, SigRan sr);
void initHostTauArr(TauArr *ta, TauRan tr);
void initHostVelArr(VelArr *va, VelRan vr);
void initHostBefAft(BefAft *ba, Range ran);
void initHostInpalse(Inpaluse *ip, SigRan sr, Pml pml, int mode, int x, int y, int z, double freq);
void initHostMedArr(MedArr *ma, SigRan sr);
void initrandom(Coord con_size, Coord *clack, int ratio);
void initClack(Object *clack, Medium med, Pml *pml, int spx, int spy, int spz, int ranx, int rany, int ranz);