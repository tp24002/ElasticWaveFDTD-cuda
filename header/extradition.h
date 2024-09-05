#pragma once

#include "./struct.h"

// メモリ確保
void allocate3DArray(double ****array, int x, int y, int z);
void allocateBefAft(BefAft *d_befAft, Range ran);
// データ転送 to Device
void copy3DArrayToDevice(double ***d_array, double ***h_array, int x, int y, int z);
void copyBefAftToDevice(BefAft *d_befAft, BefAft *h_befAft, Range ran);
// データ転送 to Host
void copy3DArrayToHost(double ***h_array, double ***d_array, int x, int y, int z);
void copyBefAftToHost(BefAft *h_befAft, BefAft *d_befAft, Range ran);