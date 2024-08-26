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
  Object steel;
  Object *clack;
  Range ran;
  Pml pml;
  Diff dif;
  MedArr ma;
  //hosu data(cpu)
  BefAft bef_host;
  BefAft aft_host;
  //device data(gpu)
  BefAft bef_device;
  BefAft aft_device;
  Inpaluse ip;
  Coord in,out1,out2,out3,out4;
  FILE *fp1,*fp2,*fp3,*fp4;
  char fn1[256],fn2[256],fn3[256],fn4[256];
  int tmp = 0;
  Coord region;
  Coord con_size;
  Coord con_st;
  Coord clack_size;
  Coord clack_st;
  Coord out[256];
  Coord center;
  int outNum;
  int tmax;

  para_in(&region,&center,&con_st,&con_size,&clack_st,&clack_size,out,&ip,&outNum,&tmax);
  //入力地点出力
  printf("in:%d,%d,%d\n", ip.in.x, ip.in.y, ip.in.z);
  printf("center:%d,%d,%d\n",center.x,center.y,center.z);
  initMedium(med);
  printf("med[E_CON].gamma = %le\n", med[E_CON].gamma);
  printf("med[E_CON].khi = %le\n", med[E_CON].khi);
  initDiff(&dif, med);
  printf("dif.dt = %le\n", dif.dt);
  initPml(&pml, med, dif);//pml 32*2
  printf("pml.fm = %le\n", pml.fm);
  initConcrete(&con, med[E_CON], pml, con_st.x, con_st.y, con_st.z, con_size.x, con_size.y, con_size.z);//x:+5,-5 y:+3,-3 z:+2,-0
  //コンクリ部分出力
  printf("start:%d,%d,%d\n",con.sp.x,con.sp.y,con.sp.z);
  printf("size:%d,%d,%d\n",con.range.x,con.range.y,con.range.z);
  initRange(&ran, region.x, region.y, region.z, pml);//PMLOK
  initBefAft(&bef_host, ran);
  initBefAft(&aft_host, ran);
  
  initMedArr(&ma, ran.sr);
  initInpalse(&ip, ran.sr, pml, ip.mode, ip.in.x, ip.in.y, ip.in.z, ip.freq);//
  printf("in:%d,%d,%d\n", ip.in.x, ip.in.y, ip.in.z);
  for(int outnum = 0; outnum < outNum; outnum++){
    printf("out:%d,%d,%d\n",out[outnum].x,out[outnum].y,out[outnum].z);
}
  if(ip.mode == E_SINE){
    printf("sin:%f\n", ip.freq);
  } else if(ip.mode == E_RCOS){
    printf("cos:%f\n", ip.freq);
  }

  insertAir(&ma, ran.sr, med[E_AIR]);
  insertConcrete(&ma, con);//maに格納
 
  int clackNum = 0;
  if(clackNum != 0){
    clack = malloc(sizeof(Object) * clackNum);
    initClack(clack,med[E_AIR], &pml, clack_st.x, clack_st.y, clack_st.z, clack_size.x, clack_size.y, clack_size.z);
    printf("clack:%d,%d,%d\n", clack->sp.x, clack->sp.y, clack->sp.z);
    insertClack(&ma, clack, clackNum);
  }

  insertPml(&ma, ran.sr, pml);
  zeroPadding(&bef_host, ran);
  zeroPadding(&aft_host, ran);

  if(clackNum != 0){
    sprintf(fn1, "./%d_%d_%d_cos/first/(%d,%d,%d)_clack.csv" , tmax, region.x, con_size.x, clack->sp.x - pml.pl1.x - 1, clack->sp.y - pml.pl1.y - 1, clack->sp.z - pml.pl1.z - 1);
  } else {
    sprintf(fn1, "./%d_%d_%d_cos/first_con.csv",tmax, region.x, con_size.x);
  }

  //ファイル名出力
  printf("%.*s\n", (int) sizeof fn1, fn1);
  fp1 = fopen(fn1, "wb");

  Coord_acc A[256];
  int count,counter;
  for (int t = 0; t < tmax; t++) {
    Vel(&aft_host, &bef_host, ma, dif, ran.vr);
    Sig(&aft_host, &bef_host, ma, dif, ran.sr, ip, t);
    Tau(&aft_host, &bef_host, ma, dif, ran.tr);
    // 加速度算出＆書き込み
    for(int i = 0; i < outNum; i++){
      Acc(&A[i],&aft_host, &bef_host, dif, out[i]);
      fprintf(fp1, "%le,%le,%le," , A[i].x,A[i].y,A[i].z);
    }
    fprintf(fp1,"\n");

    swapBefAft(&aft_host, &bef_host, ran);
    progressBar(t, tmax);
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
