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
void AccelerationCalculation(AccCoord *Acc, BefAft aft, BefAft bef, Diff dif, Coord out, Range ran);

int main(void) {
  // 静的変数
  Medium med_h[E_M_END];
  Object air_h;
  Object con_h;
  Object clack_h;
  Range ran_h, *ran_d;
  Pml pml_h;
  Diff dif_h, *dif_d;
  int tmax_h;
  int t_h;
  int outNum_h;


  // 動的変数(計算領域の大きさで大きさ決定)
  MedArr *ma_h, *ma_d;
  BefAft bef_h, *bef_d;
  BefAft aft_h, *aft_d;
  Impulse *ip_h, *ip_d;
  
  AccCoord *acc_h;

  // int RegionArea;

  FILE *fp1;
  // FILE *fp1,*fp2,*fp3,*fp4;
  char fn1[256];
  // char fn1[256],fn2[256],fn3[256],fn4[256];
  // int tmp = 0;

  Coord *out_h;
  // Coord center;
  // // int make_models; // 作成するモデルの数
  // int model_count = 0; // いくつ目のモデルを作成中か
  // int ratio;
  // int max_Patern; // コンクリートのセル数
  // // int max_ClackPatern; // 欠陥を配置できる最大のパターン数
  // int clack_count; // 割合による欠陥数
  Coord threads;

  // スレッド数
  initCoord(&threads, 4, 4, 8);
  
  // データ格納
  StaticVariable(med_h, &pml_h, &ran_h, &dif_h, &air_h, &con_h, &clack_h, &tmax_h, &outNum_h);
  // ホスト動的変数
  ma_h  = allocateHostMedArr(&ran_h);
  ip_h  = allocateHostImpulse(&ran_h);
  acc_h = allocateHostAccCoord(outNum_h);
  out_h = allocateHostCoord(outNum_h);
  
  // デバイス動的変数
  ma_d  = allocateDeviceMedArr(&ran_h);
  bef_d = allocateDeviceBefAft(&ran_h);
  aft_d = allocateDeviceBefAft(&ran_h);
  ip_d  = allocateDeviceImpulse(&ran_h);

  // デバイス静的変数
  cudaMalloc(&ran_d, sizeof(Range));
  cudaMalloc(&dif_d, sizeof(Diff));

  DynamicVariable(acc_h, ma_h, ip_h, ran_h, air_h, con_h, clack_h, pml_h, out_h, outNum_h);

  // ホスト->デバイス　データ転送
  RangeHostToDevice(&ran_h, ran_d);
  DiffHostToDevice(&dif_h, dif_d);

  MedArrHostToDevice(ma_h, ma_d, ran_h);
  ImpulseHostToDevice(ip_h, ip_d, ran_h);

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
  // printf("")
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

  
  // int idx, idy, idz;
  // idx = ip_h->in.z * ran_h.sr.Txx.x * ran_h.sr.Txx.y + ip_h->in.y * ran_h.sr.Txx.x + ip_h->in.x;
  // idy = out_h[0].z * ran_h.sr.Txx.x * ran_h.sr.Txx.y + out_h[0].y * ran_h.sr.Txx.x + out_h[0].x;
  // idz = out_h[1].z * ran_h.sr.Txx.x * ran_h.sr.Txx.y + out_h[1].y * ran_h.sr.Txx.x + out_h[1].x;
  // for(int i = 0; i < ran_h.sr.Txx.x; i++) {
  //   for(int j = 0; j < ran_h.sr.Txx.y; j++) {
  //     id = j * ran_h.sr.Txx.x * ran_h.sr.Txx.y + ip_h->in.y * ran_h.sr.Txx.x + i;
  //     if(id == idx) {
  //       printf("*");
  //     } else if(id == idy) {
  //       printf("x");
  //     } else if(id == idz) {
  //       printf("x");
  //     } else {
  //       printf("o");
  //     }
  //   }
  //   printf("\n");
  // }

  //ファイル名出力
  // printf("%.*s\n", (int) sizeof fn1, fn1);
  // fp1 = fopen(fn1, "wb");

  // double test;
  int ratio = 10;
  int model_count = 0;
  sprintf(fn1, "./clack/ratio%d/clack_%d.csv", ratio, (model_count + 1));
  fp1 = fopen(fn1, "w");

  // 0 padding
  ZeroT<<<ZeroTBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroTxy<<<ZeroTxyBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroTyz<<<ZeroTyzBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroTzx<<<ZeroTzxBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroVx<<<ZeroVxBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroVy<<<ZeroVyBlocks,threadsPerBlock>>>(aft_d, ran_d);
  ZeroVz<<<ZeroVzBlocks,threadsPerBlock>>>(aft_d, ran_d);
  for (t_h = 0; t_h < tmax_h; t_h++) {
    // 入力情報作成
    insertImpulse(ip_h, dif_h, t_h, ran_h);
    ImpulseHostToDevice(ip_h, ip_d, ran_h);

    Vel(aft_d, bef_d, ma_d, dif_d, ran_d, &ran_h, threads);
    Sig(aft_d, bef_d, ma_d, dif_d, ran_d, &ran_h, ip_d, threads);
    Tau(aft_d, bef_d, ma_d, dif_d, ran_d, &ran_h, threads);

    BefAftDeviceToHost(aft_d, &aft_h, ran_h);
    BefAftDeviceToHost(bef_d, &bef_h, ran_h);

    for(int j = 0; j < outNum_h; j++) {
      AccelerationCalculation(&acc_h[j], aft_h, bef_h, dif_h, out_h[j], ran_h);
      fprintf(fp1,"%le,%le,%le,", acc_h[j].x, acc_h[j].y, acc_h[j].z);
    }
    fprintf(fp1, "\n");

    swapBefAft(aft_d, bef_d, &ran_h, ran_d, threads);
    progressBar(t_h, tmax_h);
  }
  printf("%lf\n", aft_h.sa.Tzz[0]);
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


void AccelerationCalculation(AccCoord *Acc, BefAft aft, BefAft bef, Diff dif, Coord out, Range ran) {

  int xmax = ran.sr.Txx.x, ymax = ran.sr.Txx.y;

  int x = out.x;
  int y = out.y;
  int z = out.z;

  // 1Dインデックスの計算
  int idx   = z * (xmax * ymax) + y * xmax + x;
  int idxX  = z * (xmax * ymax) + y * xmax + (x + 1);
  int idxY  = z * (xmax * ymax) + (y + 1) * xmax + x;
  int idxZ  = (z + 1) * (xmax * ymax) + y * xmax + x;

  Acc->x = ((*(aft.va.Vx + idxX) - *(bef.va.Vx + idxX)) / dif.dt + (*(aft.va.Vx + idx) - *(bef.va.Vx + idx)) / dif.dt) / 2;

  Acc->y = ((*(aft.va.Vy + idxY) - *(bef.va.Vy + idxY)) / dif.dt + (*(aft.va.Vy + idx) - *(bef.va.Vy + idx)) / dif.dt) / 2;

  Acc->z = ((*(aft.va.Vz + idxZ) - *(bef.va.Vz + idxZ)) / dif.dt + (*(aft.va.Vz + idx) - *(bef.va.Vz + idx)) / dif.dt) / 2;
  // printf("acc:%f,%f,%f\n", Acc.x, Acc.y, Acc.z);
}