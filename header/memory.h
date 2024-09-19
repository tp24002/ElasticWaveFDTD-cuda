#pragma once
#include "./struct.h"
// hostメモリ確保
void allocate3DarrayHost(double ****array, int x, int y ,int z);
void allocateHostBefAft(BefAft **befaft_host_ptr, Range ran);
void allocateHostMedArr(MedArr **medarr_host_ptr, Range ran);
void allocateHostImpulse(Impulse **impulse_host_ptr, Range ran);
// deviceメモリ確保
// void allocate3DarrayDevice(double ***array, int x, int y, int z);
void allocate3DarrayDevice(double ****array, int x, int y, int z);
void allocateDeviceBefAft(BefAft **befaft_device_ptr, Range ran);
void allocateDeviceMedArr(MedArr **medarr_device_ptr, Range ran);
void allocateDeviceImpulse(Impulse **impulse_device_ptr, Range ran);
void allocateDeviceDiff(Diff **diff_device_ptr, Range ran);
// host to device
void copy3DarrayHostToDevice(double ***array_host, double ***array_device, int x, int y, int z);
void copyBefAftHostToDevice(BefAft *befaft_host, BefAft *befaft_device, Range ran);
void copyMedArrHostToDevice(MedArr *medarr_host, MedArr *medarr_device, Range ran);
void copyImpulseHostToDevice(Impulse *impulse_host, Impulse *impulse_device, Range ran);
void copyDiffHostToDevice(Diff *diff_host, Diff *diff_device);
// device to host
void copy3DarrayDeviceToHost(double ***array_device, double ***array_host, int x, int y, int z);
