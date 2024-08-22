#include <iostream>
#include <stdio.h>

__global__ void helloWorldKernel() {
    printf("Hello, World from GPU!\n");
}

int main() {

    helloWorldKernel<<<2, 2>>>();

    return 0;
}
