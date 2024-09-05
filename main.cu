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
#include "./header/extradition.h"

void progressBar(int now, int max);

int main(void) {
  // 宣言
  Medium med[E_M_END];
  Object con;
  Object *clack;
  Range ran;
  Pml pml;
  Diff dif;
  MedArr ma;
  //host data(cpu)
  BefAft bef;
  BefAft aft;
  Inpaluse ip;
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
  Coord threads;

  // スレッド数
  initCoord(&threads, 128, 128, 128);
  // 外部入力
  para_in(&region,&center,&con_st,&con_size,&clack_st,&clack_size,out,&ip,&outNum,&tmax);
  printf("region:(%d,%d,%d)\n", region.x, region.y, region.z);
  // 媒質パターン設定
  initMedium(med);
  printf("med[E_CON].gamma = %le\n", med[E_CON].gamma);
  printf("med[E_CON].khi = %le\n", med[E_CON].khi);
  // 差分間隔設定(initで直入力しているためinsertなし)
  initDiff(&dif, med);
  printf("dif.dt = %le\n", dif.dt);
  // pml層設定
  initPml(&pml, med, dif);
  printf("pml.fm = %le\n", pml.fm);
  // コンクリート配置
  initConcrete(&con, med[E_CON], pml, con_st.x, con_st.y, con_st.z, con_size.x, con_size.y, con_size.z);//x:+5,-5 y:+3,-3 z:+2,-0
  // コンクリ部分出力
  printf("start:%d,%d,%d\n",con.sp.x,con.sp.y,con.sp.z);
  printf("size:%d,%d,%d\n",con.range.x,con.range.y,con.range.z);
  // 計算領域設定
  initRange(&ran, region.x, region.y, region.z, pml);
  initRange(&ran, region.x, region.y, region.z, pml);
  // hostメモリ確保
  // bef = (BefAft *)malloc(sizeof(BefAft));
  initHostBefAft(&bef, ran);
  initHostBefAft(&aft, ran);
  printf("%lu\n",sizeof(bef));

  // 媒質メモリ確保
  initHostMedArr(&ma, ran.sr);
  // 入力情報メモリ確保
  initHostInpalse(&ip, ran.sr, pml, ip.mode, ip.in.x, ip.in.y, ip.in.z, ip.freq);//
  // 入力情報出力(複数地点になったら要変更)
  printf("in:%d,%d,%d\n", ip.in.x, ip.in.y, ip.in.z);
  if(ip.mode == E_SINE){
    printf("sin:%f\n", ip.freq);
  } else if(ip.mode == E_RCOS){
    printf("cos:%f\n", ip.freq);
  }
  // 計測地点出力
  for(int outnum = 0; outnum < outNum; outnum++){
    printf("out:%d,%d,%d\n",out[outnum].x,out[outnum].y,out[outnum].z);
  }

  insertAir(&ma, ran.sr, med[E_AIR]);
  printf("airok\n");
  insertConcrete(&ma, con);//maに格納
  printf("conok\n");

  ratio = 10;
  max_Patern = con_size.x * con_size.y * con_size.z;
  // max_ClackPatern = (con_size.x - 2) * (con_size.y - 2) * (con_size.z - 2);
  clack_count = max_Patern * ratio / 100;
  if(ratio != 0){
    clack = (Object *)malloc(sizeof(Object) * clack_count);
    printf("clackhalfok\n");
    initClack(clack,med[E_AIR], &pml, clack_st.x, clack_st.y, clack_st.z, clack_size.x, clack_size.y, clack_size.z);
    

    printf("ratio:%d\n", ratio);
    insertClack(&ma, clack, ratio);
  }
  printf("inclackok\n");
  insertPml(&ma, ran.sr, pml);
  printf("pmlok\n");
  zeroPadding(&bef, ran);
  zeroPadding(&aft, ran);
  printf("paddingok\n");
  if(ratio != 0){
    // model_count++;
    sprintf(fn1, "./clack/ratio%d/clack_%d.csv", ratio, (model_count + 1));
    fp1 = fopen(fn1, "wb");
    fprintf(fp1, "sp.x,sp.y,sp.z,ln.x,ln.y,ln,z\n");
    // for(int i = 0; i < clack_count; i++){
    //   fprintf(fp1, "%d,%d,%d,", clack[i].sp.x, clack[i].sp.y, clack[i].sp.z, clack[i].range.x,clack[i].range.y, clack[i].range.z);
    // }
  }
  // fclose(fp1);


  //ファイル名出力
  printf("%.*s\n", (int) sizeof fn1, fn1);
  // fp1 = fopen(fn1, "wb");

  Coord_acc Acc;
  
  for (int t = 0; t < tmax; t++) {
    insertInpulse(&ip, dif, t);

    Vel(&aft, &bef, ma, dif, ran, threads);
    printf("okVel\n");

    Sig(&aft, &bef, ma, dif, ran, ip, t, threads);
    printf("okSig\n");
    Tau(&aft, &bef, ma, dif, ran, threads);
    printf("okTau\n");

    // 加速度算出＆書き込み
    for(int i = 0; i < outNum; i++){
      Acceleration(&Acc, &aft, &bef, dif, out[i]);
      fprintf(fp1, "%le,%le,%le," , Acc.x,Acc.y,Acc.z);
    }
    fprintf(fp1,"\n");

    swapBefAft(&aft, &bef, ran);
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
