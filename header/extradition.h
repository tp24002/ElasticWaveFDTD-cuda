#pragma once

#include "./struct.h"

void MemoryBefAftToDevice(BefAft *h, BefAft *d, Range ran);
void MemoryBefAftToHost(BefAft *h, BefAft *d, Range ran);
void onceHtoD(MedArr *ma_h, MedArr *ma_d, Diff *dif_h, Diff *dif_d, Range *ran_h, Range *ran_d);
void MemoryInpaluseToDevice(Inpaluse *ip_h, Inpaluse *ip_d, Range ran);
void allocate3DArray(double ****array, int x, int y, int z);
void allocateBefAft(BefAft *d_befAft, int x, int y, int z);
void copy3DArrayToDevice(double ***d_array, double ***h_array, int x, int y, int z);
void copyBefAftToDevice(BefAft *d_befAft, BefAft *h_befAft, Range ran);