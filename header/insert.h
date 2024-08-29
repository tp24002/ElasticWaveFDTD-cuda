#pragma once
#include "./struct.h"

void insertAir(MedArr *ma, SigRan sr, Medium air);
void insertConcrete(MedArr *ma, Object con);
void insertClack(MedArr *ma, Object *clack, int clackNum);
void insertPml(MedArr *ma, SigRan sr, Pml pml);
void insertInpulse(Inpaluse *ip, Diff dif, int t);
void zeroPadSig(SigArr *sa, SigRan sr);
void zeroPadTau(TauArr *ta, TauRan tr);
void zeroPadVel(VelArr *va, VelRan vr);
void zeroPadding(BefAft *ba, Range ran);
