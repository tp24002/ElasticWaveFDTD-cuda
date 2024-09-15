#define _USE_MATH_DEFINES
#include "../header/print.h"
#include "../header/struct.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/insert.h"
#include "../header/init.h"
#include "../header/extradition.h"


__global__ void CheckPrint(double out) {
    printf("%f\n",out);
}
