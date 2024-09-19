#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "./header/init.h"
#include "./header/insert.h"
#include "./header/struct.h"
#include "./header/update.h"
#include "./header/parameter.h"
#include "./header/memory.h"

void progressBar(int now, int max);

int main(void) {
  // 宣言
  Medium med[E_M_END];
  Object *con;
  Object *clack;
  Range ran;
  Pml *pml;
  //host data(cpu)
  Diff *dif_h, *dif_d;
  MedArr *ma_h, *ma_d;
  BefAft *bef_h, *bef_d;
  BefAft *aft_h, *aft_d;
  Impulse *ip_h, *ip_d;
  int tmax;


  // Coord out1,out2,out3,out4;
  FILE *fp1;
  // FILE *fp1,*fp2,*fp3,*fp4;
  // char fn1[256];
  // char fn1[256],fn2[256],fn3[256],fn4[256];
  // int tmp = 0;

  Coord out[16];
  // Coord center;
  int outNum = 10;
  // // int make_models; // 作成するモデルの数
  // int model_count = 0; // いくつ目のモデルを作成中か
  // int ratio;
  // int max_Patern; // コンクリートのセル数
  // // int max_ClackPatern; // 欠陥を配置できる最大のパターン数
  // int clack_count; // 割合による欠陥数
  Coord threads;
  Coord_acc **Acc_h, **Acc_d;

  // スレッド数
  initCoord(&threads, 4, 4, 8);
  // ran = (Range *)malloc(sizeof(Range));
  pml = (Pml *)malloc(sizeof(Pml));
  dif_h = (Diff *)malloc(sizeof(Diff));
  con = (Object *)malloc(sizeof(Object));
  clack = (Object *)malloc(sizeof(Object));
  ma_h = (MedArr *)malloc(sizeof(MedArr));

  // データ格納
  StaticVariable(&ran, pml, dif_h, con, clack, med, ma_h, &tmax);

  allocateHostBefAft(&bef_h, ran);
  allocateHostBefAft(&aft_h, ran);
  zeroPadding(bef_h, ran);
  zeroPadding(aft_h, ran);
  printf("aaa\n");
  printf("a:%f\n", bef_h->sa.Txx[73][73][73]);
  printf("nono\n");

  allocateHostMedArr(&ma_h, ran);
  printf("nono\n");
  allocateHostImpulse(&ip_h, ran);
  DynamicVariable(ma_h, ip_h, ran, med, *con, *clack, *pml, *dif_h, tmax);
  
  // 出力
  printf("time:%d\n", tmax);
  printf("range:%d,%d,%d(in pml)\n", ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
  printf("pml:%d,%d,%d\n", pml->pl1.x, pml->pl1.y, pml->pl1.z);
  printf("dif_xyz:%f,%f,%f\n", dif_h->dx, dif_h->dy, dif_h->dz);
  printf("dif_time:%f\n", dif_h->dt);
  printf("in:%d,%d,%d\n", ip_h->in.x, ip_h->in.y, ip_h->in.z);
  if(ip_h->mode == E_SINE){
    printf("sin:%f\n", ip_h->freq);
  } else if(ip_h->mode == E_RCOS){
    printf("cos:%f\n", ip_h->freq);
  }

  // 関数化推奨
  // 加速度メモリ確保
  Acc_h = (Coord_acc **)malloc(sizeof(Coord_acc *) * outNum);
  for (int i = 0; i < outNum; i++) {
    Acc_h[i] = (Coord_acc *)malloc(tmax * sizeof(Coord_acc));
  }
  cudaMalloc((void **)&Acc_d, outNum * sizeof(Coord_acc *));
  for (int i = 0; i < outNum; i++) {
    Coord_acc *temp_d;
    cudaMalloc((void **)&temp_d, tmax * sizeof(Coord_acc));
    cudaMemcpy(Acc_h[i], temp_d, tmax * sizeof(Coord_acc), cudaMemcpyHostToDevice);
    cudaMemcpy(&Acc_d[i], &temp_d, sizeof(Coord_acc *), cudaMemcpyHostToDevice);
  }
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

  // device構造体本体のメモリ確保
  cudaMalloc((void **)&aft_d, sizeof(BefAft));
  cudaMalloc((void **)&bef_d, sizeof(BefAft));
  cudaMalloc((void **)&ma_d, sizeof(MedArr));
  cudaMalloc((void **)&dif_d, sizeof(Diff));
  cudaMalloc((void **)&ip_d, sizeof(Impulse));
  // device構造体中身(メンバ)のメモリ確保関数
  allocateDeviceBefAft(&aft_d, ran);
  allocateDeviceBefAft(&bef_d, ran);
  printf("aloocate BefAft ok\n");
  allocateDeviceMedArr(&ma_d, ran);
  printf("aloocate MedArr ok\n");
  allocateDeviceImpulse(&ip_d, ran);
  printf("aloocate Impulse ok\n");
  // double test;
  for (int t = 0; t < tmax; t++) {

    
    // printf("%f\n",ip_h.Tzz[ip_h.in.x][ip_h.in.y][ip_h.in.z]);

    // printf("%f\n",ip_d.freq);
    Vel(aft_h, bef_h, aft_d, bef_d, *ma_h, ma_d, *dif_h, dif_d, ran, threads);
    Sig(aft_h, bef_h, aft_d, bef_d, *ma_h, ma_d, *dif_h, dif_d, ran, *ip_h, ip_d, t, threads);
    Tau(aft_h, bef_h, aft_d, bef_d, *ma_h, ma_d, *dif_h, dif_d, ran, threads);


    // 加速度算出＆書き込み
    Acceleration<<<1, 1>>>(Acc_d, aft_d, bef_d, *dif_h, out, outNum, t);

    // swapBefAft<<<1, 1>>>(&aft_d, &bef_d, ran);
    // progressBar(t, tmax);
  }

  cudaError_t err = cudaMemcpy(Acc_d, Acc_h, outNum * sizeof(Coord_acc *), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    printf("acc d to h Error: %s\n", cudaGetErrorString(err));
    return;
  }
  for (int j = 0; j < tmax; j++) {
    for (int i = 0; i < outNum; i++) {
      fprintf(fp1,"%le,%le,%le,", Acc_h[i][j].x, Acc_h[i][j].y, Acc_h[i][j].z);
    }
    fprintf(fp1, "\n");
    progressBar(j, tmax);
  }
  fclose(fp1);
  printf("loop end.\n");
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
