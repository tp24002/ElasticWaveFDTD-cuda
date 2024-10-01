#pragma once
#include "./struct.h"
// hostメモリ確保
MedArr* allocateHostMedArr(Range *ran);
Impulse* allocateHostImpulse(Range *ran);
AccCoord* allocateHostAccCoord(int outNum);
Coord* allocateHostCoord(int outNum);
// deviceメモリ確保
BefAft* allocateDeviceBefAft(Range *ran);
MedArr* allocateDeviceMedArr(Range *ran);
Impulse* allocateDeviceImpulse(Range *ran);
AccCoord* allocateDeviceAccCoord(int outNum);
Coord* allocateDeviceCoord(int outNum);
// データ転送 ホスト->デバイス
void MedArrHostToDevice(MedArr *ma_h, MedArr *ma_d, Range ran);
void ImpulseHostToDevice(Impulse *ip_h, Impulse *ip_d, Range ran);
void RangeHostToDevice(Range *ran_h, Range *ran_d);
void DiffHostToDevice(Diff *dif_h, Diff *dif_d);
void CoordHostToDevice(Coord *out_h, Coord *out_d, int outNum);
// データ転送　デバイス->ホスト
void AccCoordDeviceToHost(AccCoord *acc_d, AccCoord *acc_h, int outNum);
void BefAftDeviceToHost(BefAft *ba_d, BefAft *ba_h, Range ran);
