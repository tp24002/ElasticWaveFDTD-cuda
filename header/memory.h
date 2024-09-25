#pragma once
#include "./struct.h"
// hostメモリ確保
AccCoord* allocateHostAccCoord(int outNum);
// deviceメモリ確保
BefAft* allocateDeviceBefAft(Range *ran);
MedArr* allocateDeviceMedArr(Range *ran);
Impulse* allocateDeviceImpulse(Range *ran);
AccCoord* allocateDeviceAccCoord(int outNum);
// device to host
// void AccCoordDeviceToHost(AccCoord *acc_d, AccCoord *acc_h, int outNum, int tmax);
void AccCoordDeviceToHost(AccCoord *acc_d, AccCoord *acc_h, int outNum);