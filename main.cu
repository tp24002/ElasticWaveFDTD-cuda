#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "./header/init.h"
#include "./header/struct.h"
#include "./header/update.h"
#include "./header/parameter.h"
#include "./header/memory.h"

void progressBar(int now, int max);

int main(void) {
  // 静的変数
  Medium *med_d;
  Object *air_d;
  Object *con_d;
  Object *clack_d;
  Range ran_h, *ran_d;
  Pml *pml_d;
  Diff *dif_d;

  // 動的変数
  MedArr *ma_d;
  BefAft *bef_d;
  BefAft *aft_d;
  Impulse *ip_d;
  AccCoord *acc_h, *acc_d;
  int tmax_h, *tmax_d;
  int t_h, *t_d;
  // int RegionArea;

  FILE *fp1;
  // FILE *fp1,*fp2,*fp3,*fp4;
  char fn1[256];
  // char fn1[256],fn2[256],fn3[256],fn4[256];
  // int tmp = 0;

  Coord *out;
  // Coord center;
  int outNum_h, *outNum_d;
  // // int make_models; // 作成するモデルの数
  // int model_count = 0; // いくつ目のモデルを作成中か
  // int ratio;
  // int max_Patern; // コンクリートのセル数
  // // int max_ClackPatern; // 欠陥を配置できる最大のパターン数
  // int clack_count; // 割合による欠陥数
  Coord threads;

  // スレッド数
  initHostCoord(&threads, 4, 4, 8);

  // デバイス変数
  // 静的変数メモリ確保
  cudaMalloc((void**)&med_d  , E_M_END * sizeof(Medium));
  cudaMalloc((void**)&pml_d  , sizeof(Pml));
  cudaMalloc((void**)&ran_d  , sizeof(Range));
  cudaMalloc((void**)&dif_d  , sizeof(Diff));
  cudaMalloc((void**)&air_d  , sizeof(Object));
  cudaMalloc((void**)&con_d  , sizeof(Object));
  cudaMalloc((void**)&clack_d, sizeof(Object));
  cudaMalloc((void**)&ma_d   , sizeof(MedArr));
  cudaMalloc((void**)&tmax_d, sizeof(int));
  cudaMalloc((void**)&outNum_d, sizeof(int));



  
  // データ格納
  StaticVariable<<<1,1>>>(med_d, pml_d, ran_d, dif_d, air_d, con_d, clack_d, ma_d, tmax_d, outNum_d);
  // Range デバイス変数->ホスト変数
  cudaMemcpy(&ran_h, ran_d, sizeof(Range), cudaMemcpyDeviceToHost);
  cudaMemcpy(&tmax_h, tmax_d, sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&outNum_h, outNum_d, sizeof(int), cudaMemcpyDeviceToHost);
  
  // 動的変数デバイスメモリ確保
  // cudaMalloc(&bef_d, sizeof(BefAft));
  // cudaMalloc(&aft_d, sizeof(BefAft));
  // cudaMalloc(&ma_d, sizeof(MedArr));
  // cudaMalloc(&ip_d, sizeof(Impulse));
  // cudaMalloc(&acc_d, outNum_h * sizeof(AccCoord));
  cudaMalloc((void**)&out, outNum_h * sizeof(Coord));


  // 動的変数デバイスメモリ確保
  bef_d = allocateDeviceBefAft(&ran_h);
  aft_d = allocateDeviceBefAft(&ran_h);
  ma_d = allocateDeviceMedArr(&ran_h);
  ip_d = allocateDeviceImpulse(&ran_h);
  acc_d = allocateDeviceAccCoord(outNum_h);
  
  

  DynamicVariable<<<1,1>>>(bef_d, aft_d, acc_d, ma_d, ip_d, ran_d, med_d, air_d, con_d, clack_d, pml_d, dif_d, out, outNum_d);

  // 動的変数ホストメモリ確保
  acc_h = allocateHostAccCoord(outNum_h);

  // 出力
  printf("time:%d\n", tmax_h);
  printf("range:%d,%d,%d(in pml)\n", ran_h.sr.Txx.x, ran_h.sr.Txx.y, ran_h.sr.Txx.z);

  ///////////clack
  // ratio = 10;
  // max_Patern = con_size.x * con_size.y * con_size.z;
  // // max_ClackPatern = (con_size.x - 2) * (con_size.y - 2) * (con_size.z - 2);
  // clack_count = max_Patern * ratio / 100;
  // if(ratio != 0){
  //   clack = (Object *)malloc(sizeof(Object) * clack_count);
  //   printf("clackhalfok\n");
  //   initClack(clack,med[E_AIR], &pml, clack_st.x, clack_st.y, clack_st.z, clack_size.x, clack_size.y, clack_size.z);
  //   printf("ratio:%d\n", ratio);
  //   insertClack(&ma_h, clack, ratio);
  // }
  // if(ratio != 0){
  //   // model_count++;
  //   sprintf(fn1, "./clack/ratio%d/clack_%d.csv", ratio, (model_count + 1));
  //   fp1 = fopen(fn1, "wb");
  //   fprintf(fp1, "sp.x,sp.y,sp.z,ln.x,ln.y,ln,z\n");
  //   // for(int i = 0; i < clack_count; i++){
  //   //   fprintf(fp1, "%d,%d,%d,", clack[i].sp.x, clack[i].sp.y, clack[i].sp.z, clack[i].range.x,clack[i].range.y, clack[i].range.z);
  //   // }
  // }


  //ファイル名出力
  // printf("%.*s\n", (int) sizeof fn1, fn1);
  // fp1 = fopen(fn1, "wb");

  // double test;
  int ratio = 10;
  int model_count = 0;
  sprintf(fn1, "./clack/ratio%d/clack_%d.csv", ratio, (model_count + 1));
  fp1 = fopen(fn1, "w");
  dim3 threadsPerBlock(threads.x, threads.y, threads.z); // 1ブロックあたりのスレッド数
  dim3 AccBlocks((outNum_h - 1 + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (outNum_h - 1 + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (outNum_h + threadsPerBlock.z - 1)     / threadsPerBlock.z);
  for (t_h = 0; t_h < tmax_h; t_h++) {
    cudaMemcpy(&t_d, &t_h, sizeof(int), cudaMemcpyHostToDevice);
    // printf("host to device t_h,t_d\n");
    insertImpulse<<<1,1>>>(ip_d, dif_d, t_d, ran_d);
    // printf("insert Impluse\n");
    // printf("%f\n",ip_h.Tzz[ip_h.in.x][ip_h.in.y][ip_h.in.z]);

    // printf("%f\n",ip_d.freq);
    Vel<<<1,1>>>(aft_d, bef_d, ma_d, dif_d, ran_d, threads);
    // printf("Vel OK\n");
    Sig<<<1,1>>>(aft_d, bef_d, ma_d, dif_d, ran_d, ip_d, threads);
    // printf("Sig OK\n");
    Tau<<<1,1>>>(aft_d, bef_d, ma_d, dif_d, ran_d, threads);
    // printf("Tau OK\n");

    // 加速度算出＆書き込み
    AccelerationCalculation<<<AccBlocks, threadsPerBlock>>>(acc_d, aft_d, bef_d, dif_d, out, ran_d, outNum_d);
    printf("Tau OK\n");

    AccCoordDeviceToHost(acc_d, acc_h, outNum_h);

    for(int j = 0; j < outNum_h; j++){
      fprintf(fp1,"%le,%le,%le,", acc_h[j].x, acc_h[j].y, acc_h[j].z);
    }
    fprintf(fp1, "\n");

    // printf("acc calcu\n");
    swapBefAft<<<1, 1>>>(aft_d, bef_d, ran_d, threads);
    // printf("swap befaft\n");
    // progressBar(t_h, tmax_h);
  }
  fclose(fp1);
  printf("loop end\n");
  return 0;
}

void progressBar(int now, int max) {
  int bar_width = 50;
  double progress = (double)(now + 1) / (double)max;
  int bar_length = (int)(progress * bar_width);
  printf("Progress: [");
  for (int j = 0; j < bar_length; j++) {
    printf("=");
  }
  for (int j = bar_length; j < bar_width; j++) {
    printf(" ");
  }
  printf("] %.2f%%\r", progress * 100);
  fflush(stdout);
}
