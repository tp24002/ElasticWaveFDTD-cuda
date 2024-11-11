#pragma once
#include "./struct.h"
// hostメモリ確保
MedArr* allocateHostMedArr(Range ran);
BefAft* allocateHostBefAft(Range ran);
Impulse* allocateHostImpulse(int innum);
DimI3* allocateHostDimI3(int outnum);
DimD3* allocateHostDimD3(int outNum);
// deviceメモリ確保
MedArr* allocateDeviceMedArr(Range ran);
BefAft* allocateDeviceBefAft(Range ran);
ImpulseArr* allocateDeviceImpulseArr(Range ran);
Impulse* allocateDeviceImpulse(int innum);
DimI3* allocateDeviceDimI3(int outNum);
DimD3* allocateDeviceDimD3(int outNum);
// データ転送 ホスト->デバイス
void RangeHostToDevice(Range *ran_d, Range *ran_h);
void DiffHostToDevice(Diff *dif_d, Diff *dif_h);
void MedArrHostToDevice(MedArr *ma_d, MedArr *ma_h, Range ran);
void ImpulseHostToDevice(Impulse *ip_d, Impulse *ip_h, int innum);
void DimI3HostToDevice(DimI3 *di_d, DimI3 *di_h, int outnum);
// データ転送　デバイス->ホスト
void BefAftDeviceToHost(BefAft *ba_h, BefAft *ba_d, Range ran);
void DimD3DeviceToHost(DimD3 *acc_h, DimD3 *acc_d, int outNum);
