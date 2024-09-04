#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../header/init.h"
#include "../header/insert.h"
#include "../header/struct.h"
#include "../header/update.h"
#include "../header/para_in.h"


// device構造体中身(メンバ)メモリ確保関数
void MemoryDeviceBefAft(BefAft *d, Range ran) {
    cudaMalloc((void **)&(d->sa), sizeof(SigArr));
    cudaMalloc((void **)&(d->ta), sizeof(TauArr));
    cudaMalloc((void **)&(d->va), sizeof(VelArr));


}

void MemoryDeviceSigArr(SigArr *sa, SigRan sr) {
    int x = sr.Txx.x, y = sr.Txx.y, z = sr.Txx.z;
    cudaMalloc((void **)&sa, x * sizeof(double **));
    for(int i = 0; i < sr.Txx.x; i++) {

    }
}