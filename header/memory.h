#pragma once
#include "./struct.h"
// hostメモリ確保
void allocateHostAccCoord(AccCoord **acccoordptr, int tmax, int outNum);
// deviceメモリ確保
BefAft* allocateDeviceBefAft(Range *ran);
MedArr* allocateDeviceMedArr(Range *ran);
Impulse* allocateDeviceImpulse(Range *ran);
void allocateDeviceAccCoord(AccCoord **acccoordptr, int tmax, int outNum);
// device to host
void AccCoordDeviceToHost(AccCoord *acc_d, AccCoord *acc_h, int outNum, int tmax);
