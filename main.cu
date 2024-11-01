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
void AccelerationCalculation(DimD3 *Acc, BefAft *aft, BefAft *bef, Diff dif, DimI3 out, Range ran);

int main(void) {
  int tmax_h;
  int outnum_h, innum_h;
  // 静的変数
  Medium med_h[E_M_END];
  Object air_h;
  Object con_h;
  Object clack_h;
  Pml pml_h;
  Range ran_h, *ran_d;
  Diff dif_h, *dif_d;

  // 動的変数(計算領域の大きさで大きさ決定)
  MedArr *ma_h, *ma_d; // 転送後即解放
  BefAft *bef_d;
  BefAft *aft_d;
  ImpulseArr *ipa_d;
  // 動的変数(入出力地点数で大きさ決定)
  Impulse *ip_h, *ip_d;
  DimD3 *acc_h, *acc_d;

  DimI3 *out_h, *out_d;
  // DimI3 center;
  // // int make_models; // 作成するモデルの数
  // int model_count = 0; // いくつ目のモデルを作成中か
  // int ratio;
  // int max_Patern; // コンクリートのセル数
  // // int max_ClackPatern; // 欠陥を配置できる最大のパターン数
  // int clack_count; // 割合による欠陥数
  DimI3 threads;

  // スレッド数
  initDimI3(&threads, 4, 4, 8);
  
  // データ格納
  StaticVariable(med_h, &pml_h, &ran_h, &dif_h, &air_h, &con_h, &clack_h, &tmax_h, &outnum_h, &innum_h);
  // ホスト動的変数
  ma_h  = allocateHostMedArr(ran_h);
  ip_h  = allocateHostImpulse(innum_h);
  acc_h = allocateHostDimD3(outnum_h);
  out_h = allocateHostDimI3(outnum_h);
  DynamicVariable(acc_h, ma_h, ip_h, ran_h, air_h, con_h, clack_h, pml_h, out_h, outnum_h);
  
  // デバイス動的変数
  ma_d  = allocateDeviceMedArr(ran_h);
  bef_d = allocateDeviceBefAft(ran_h);
  aft_d = allocateDeviceBefAft(ran_h);
  ipa_d = allocateDeviceImpulseArr(ran_h);
  ip_d  = allocateDeviceImpulse(innum_h);
  acc_d = allocateDeviceDimD3(outnum_h);
  out_d = allocateDeviceDimI3(outnum_h);
  
  // デバイス静的変数
  cudaMalloc(&ran_d, sizeof(Range));
  cudaMalloc(&dif_d, sizeof(Diff));

  // ホスト->デバイス　データ転送
  RangeHostToDevice(ran_d, &ran_h);
  DiffHostToDevice(dif_d, &dif_h);
  MedArrHostToDevice(ma_d, ma_h, ran_h);
  ImpulseHostToDevice(ip_d, ip_h, innum_h);
  DimI3HostToDevice(out_d, out_h, outnum_h);

  free(ma_h);

  dim3 threadsPerBlock(threads.x, threads.y, threads.z);  // ブロック内のスレッド数
  dim3 ZeroTBlocks((ran_h.sr.Txx.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (ran_h.sr.Txx.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                   (ran_h.sr.Txx.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroTxyBlocks((ran_h.tr.Txy.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h.tr.Txy.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h.tr.Txy.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroTyzBlocks((ran_h.tr.Tyz.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h.tr.Tyz.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h.tr.Tyz.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroTzxBlocks((ran_h.tr.Tzx.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (ran_h.tr.Tzx.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (ran_h.tr.Tzx.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroVxBlocks((ran_h.vr.Vx.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (ran_h.vr.Vx.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (ran_h.vr.Vx.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroVyBlocks((ran_h.vr.Vy.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (ran_h.vr.Vy.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (ran_h.vr.Vy.z + threadsPerBlock.z - 1) / threadsPerBlock.z);
  dim3 ZeroVzBlocks((ran_h.vr.Vz.x + threadsPerBlock.x - 1) / threadsPerBlock.x,
                    (ran_h.vr.Vz.y + threadsPerBlock.y - 1) / threadsPerBlock.y,
                    (ran_h.vr.Vz.z + threadsPerBlock.z - 1) / threadsPerBlock.z);

  // 出力
  // Medium
  for(int i = 0; i < E_M_END; i++) {
    printf("Med:%f\n", med_h[i].rho);
  }
  // Pml
  printf("pml.fm:%f\n", pml_h.fm);
  printf("pml.ta:%f\n", pml_h.ta);
  // Range
  printf("Range:%d,%d,%d\n", ran_h.sr.Txx.x, ran_h.sr.Txx.y, ran_h.sr.Txx.z);
  // Diff
  printf("sp   diff:%f,%f,%f\n", dif_h.dx, dif_h.dy, dif_h.dz);
  printf("time diff:%e\n", dif_h.dt);
  // Impulse
  for(int i = 0; i < innum_h; i++) {
    if(ip_h[i].mode == E_SINE) {
      printf("ip:%lf(sin)\n", ip_h[i].freq);
    } else {
      printf("ip:%lf(cos)\n", ip_h[i].freq);
    }
    printf("in[%d]:%d,%d,%d\n", i, ip_h[i].in.x, ip_h[i].in.y, ip_h[i].in.z);
  }
  
  for(int i = 0; i < outnum_h; i++) {
    printf("out[%d]:%d,%d,%d\n", i, out_h[i].x, out_h[i].y, out_h[i].z);
  }
  printf("time:%d\n", tmax_h);

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
  // fp1 = fopen(fn1, "wb");

  // double test;
  // int ratio = 10;
  // int model_count = 0;
  // sprintf(fn1, "./clack/ratio%d/clack_%d.csv", ratio, (model_count + 1));
  // printf("%.*s\n", (int) sizeof fn1, fn1);
  // fp1 = fopen(fn1, "w");

  // 0 padding
  ZeroT<<<ZeroTBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroTxy<<<ZeroTxyBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroTyz<<<ZeroTyzBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroTzx<<<ZeroTzxBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroVx<<<ZeroVxBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroVy<<<ZeroVyBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroVz<<<ZeroVzBlocks,threadsPerBlock>>>(aft_d, ran_d);
  cudaDeviceSynchronize();
  FILE *fp;
  char fn1[256];
  sprintf(fn1, "./clack.csv");
  fp = fopen(fn1, "w");
  for (int t_h = 0; t_h < tmax_h; t_h++) {
    // printf("ok\n");
    createImpulse<<<1,1>>>(ipa_d, ip_d, dif_d, ran_d, innum_h, t_h);////////
    cudaDeviceSynchronize();
    Vel(aft_d, bef_d, ma_d, dif_d, ran_d, &ran_h, threads);//////////////////////////////
    Sig(aft_d, bef_d, ma_d, dif_d, ran_d, &ran_h, ipa_d, threads);
    Tau(aft_d, bef_d, ma_d, dif_d, ran_d, &ran_h, threads);
    Acceleration<<<1,1>>>(acc_d, aft_d->va, bef_d->va, dif_d, out_d, ran_d, outnum_h);
    DimD3DeviceToHost(acc_h, acc_d, outnum_h);
    swapBefAft(aft_d, bef_d, &ran_h, ran_d, threads);
    for(int i = 0; i < outnum_h; i++) {
      fprintf(fp, "%.6lf,%.6lf,%.6lf,", acc_h[i].x, acc_h[i].y, acc_h[i].z);
    }
    fprintf(fp, "\n");
    progressBar(t_h, tmax_h);
  }
  fclose(fp);
  printf("\nloop end.\n");

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

