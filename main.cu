#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "./header/init.h"
#include "./header/insert.h"
#include "./header/struct.h"
#include "./header/update.h"
#include "./header/para_in.h"

void progressBar(int now, int max);

int main(void) {
  // 宣言
  Medium med[E_M_END];
  Object con;
  Object *clack;
  Range ran_host;
  Range ran_device;
  Pml pml;
  Diff dif_host;
  Diff dif_device;
  MedArr ma_host;
  MedArr ma_device;
  //hosu data(cpu)
  BefAft bef_host;
  BefAft aft_host;
  //device data(gpu)
  BefAft bef_device;
  BefAft aft_device;
  Inpaluse ip_host;
  Inpaluse ip_device;
  // Coord out1,out2,out3,out4;
  FILE *fp1;
  // FILE *fp1,*fp2,*fp3,*fp4;
  char fn1[256];
  // char fn1[256],fn2[256],fn3[256],fn4[256];
  // int tmp = 0;
  Coord region;
  Coord con_size;
  Coord con_st;
  Coord clack_size;
  Coord clack_st;
  Coord out[256];
  Coord center;
  int outNum;
  int tmax;
  // int make_models; // 作成するモデルの数
  int model_count = 0; // いくつ目のモデルを作成中か
  int ratio;
  int max_Patern; // コンクリートのセル数
  // int max_ClackPatern; // 欠陥を配置できる最大のパターン数
  int clack_count; // 割合による欠陥数

  para_in(&region,&center,&con_st,&con_size,&clack_st,&clack_size,out,&ip_host,&outNum,&tmax);
  //入力地点出力
  printf("in:%d,%d,%d\n", ip_host.in.x, ip_host.in.y, ip_host.in.z);
  printf("center:%d,%d,%d\n",center.x,center.y,center.z);
  initMedium(med);
  printf("med[E_CON].gamma = %le\n", med[E_CON].gamma);
  printf("med[E_CON].khi = %le\n", med[E_CON].khi);
  initDiff(&dif_host, med);
  printf("dif_host.dt = %le\n", dif_host.dt);
  initPml(&pml, med, dif_host);//pml 32*2
  printf("pml.fm = %le\n", pml.fm);
  initConcrete(&con, med[E_CON], pml, con_st.x, con_st.y, con_st.z, con_size.x, con_size.y, con_size.z);//x:+5,-5 y:+3,-3 z:+2,-0
  //コンクリ部分出力
  printf("start:%d,%d,%d\n",con.sp.x,con.sp.y,con.sp.z);
  printf("size:%d,%d,%d\n",con.range.x,con.range.y,con.range.z);
  initRange(&ran_host, region.x, region.y, region.z, pml);//PMLOK
  initRange(&ran_device, region.x, region.y, region.z, pml);//PMLOK
  // hostメモリ確保
  initHostBefAft(&bef_host, ran_host);
  initHostBefAft(&aft_host, ran_host);

  // deviceメモリ確保
  initDeviceBefAft(&bef_device, ran_device);
  initDeviceBefAft(&aft_device, ran_device);

  initHostMedArr(&ma_host, ran_host.sr);
  initDeviceMedArr(&ma_device, ran_host.sr);

  initHostInpalse(&ip_host, ran_host.sr, pml, ip_host.mode, ip_host.in.x, ip_host.in.y, ip_host.in.z, ip_host.freq);//
  printf("ok\n");
  initDeviceInpalse(&ip_device, ran_device.sr, pml, ip_device.mode, ip_device.in.x, ip_device.in.y, ip_device.in.z, ip_device.freq);//
  printf("in:%d,%d,%d\n", ip_host.in.x, ip_host.in.y, ip_host.in.z);
  for(int outnum = 0; outnum < outNum; outnum++){
    printf("out:%d,%d,%d\n",out[outnum].x,out[outnum].y,out[outnum].z);
  }
  if(ip_host.mode == E_SINE){
    printf("sin:%f\n", ip_host.freq);
  } else if(ip_host.mode == E_RCOS){
    printf("cos:%f\n", ip_host.freq);
  }

  insertAir(&ma_host, ran_host.sr, med[E_AIR]);
  insertConcrete(&ma_host, con);//maに格納
  ratio = 10;
  max_Patern = con_size.x * con_size.y * con_size.z;
  // max_ClackPatern = (con_size.x - 2) * (con_size.y - 2) * (con_size.z - 2);
  clack_count = max_Patern * ratio / 100;
  if(ratio != 0){
    clack = (Object *)malloc(sizeof(Object) * clack_count);
    initClack(clack,med[E_AIR], &pml, clack_st.x, clack_st.y, clack_st.z, clack_size.x, clack_size.y, clack_size.z);
    printf("ratio:%d\n", ratio);
    insertClack(&ma_host, clack, ratio);
  }

  insertPml(&ma_host, ran_host.sr, pml);
  zeroPadding(&bef_host, ran_host);
  zeroPadding(&aft_host, ran_host);

  if(ratio != 0){
    // model_count++;
    sprintf(fn1, "./clack/ratio%d/clack_%d.csv", ratio, (model_count + 1));
    fp1 = fopen(fn1, "wb");
    fprintf(fp1, "sp.x,sp.y,sp.z,ln.x,ln.y,ln,z\n");
    // for(int i = 0; i < clack_count; i++){
    //   fprintf(fp1, "%d,%d,%d,", clack[i].sp.x, clack[i].sp.y, clack[i].sp.z, clack[i].range.x,clack[i].range.y, clack[i].range.z);
    // }
  } else {
    sprintf(fn1, "./%d_%d_%d_cos/first_con.csv",tmax, region.x, con_size.x);
  }
  fclose(fp1);


  //ファイル名出力
  printf("%.*s\n", (int) sizeof fn1, fn1);
  // fp1 = fopen(fn1, "wb");

  Coord_acc A[256];
  // int count,counter;
  

  // 変数引き渡しcpu->gpu
  cudaMemcpy(&aft_host, &aft_device, sizeof(BefAft), cudaMemcpyDeviceToHost);

  cudaMemcpy(&bef_host, &bef_device, sizeof(BefAft), cudaMemcpyDeviceToHost);
  cudaMemcpy(&ma_host, &ma_device, sizeof(MedArr), cudaMemcpyDeviceToHost);
  cudaMemcpy(&dif_host, &dif_device, sizeof(Diff), cudaMemcpyDeviceToHost);
  cudaMemcpy(&ran_host, &ran_device, sizeof(Range), cudaMemcpyDeviceToHost);
  cudaMemcpy(&ip_host, &ip_device, sizeof(Inpaluse), cudaMemcpyDeviceToHost);
  for (int t = 0; t < tmax; t++) {
    Vel(&aft_device, &bef_device, ma_device, dif_device, ran_device.vr);
    printf("ok\n");
    Sig(&aft_device, &bef_device, ma_device, dif_device, ran_device.sr, ip_device, t);
    printf("ok\n");
    Tau(&aft_device, &bef_device, ma_device, dif_device, ran_device.tr);
    printf("ok\n");
    // 加速度算出＆書き込み
    for(int i = 0; i < outNum; i++){
      Acc(&A[i],&aft_host, &bef_host, dif_host, out[i]);
      fprintf(fp1, "%le,%le,%le," , A[i].x,A[i].y,A[i].z);
    }
    fprintf(fp1,"\n");

    swapBefAft(&aft_host, &bef_host, ran_host);
    progressBar(t, tmax);
  }
  // fclose(fp1);
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
