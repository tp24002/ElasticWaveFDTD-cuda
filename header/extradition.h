#include "./struct.h"

void onceHtoD(MedArr *ma_h, MedArr *ma_d, Diff *dif_h, Diff *dif_d, Range *ran_h, Range *ran_d);
void loopHtoD(Inpaluse *ip_h, Inpaluse *ip_d, BefAft *aft_h, BefAft *aft_d, BefAft *bef_h, BefAft *bef_d, Range ran);
void loopDtoH(BefAft *aft_h, BefAft *aft_d, BefAft *bef_h, BefAft *bef_d, Range ran);