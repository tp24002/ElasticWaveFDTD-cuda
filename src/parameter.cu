#include "../header/struct.h"
#include "../header/init.h"
#include "../header/parameter.h"
#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>

#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))


// 静的変数
void StaticVariable(Medium *med, Pml *pml, Range *ran, Diff *dif, Object *air, Object *con, Object *clack, int *tmax, int *outnum, int *innum) {
  // 入力
  DimI3 region;
  initDimI3(&region, 15, 15, 15);
  DimI3 con_start;
  DimI3 con_ran;
  initDimI3(&con_start, 3, 3, 3);
  initDimI3(&con_ran, 0, 0, 0);
  DimI3 clack_start;
  DimI3 clack_ran;
  initDimI3(&clack_start, 3, 3, 3);
  initDimI3(&clack_ran, 0, 0, 0);
  *tmax = 32768;
  // *tmax = 131072;
  *outnum = 2;
  *innum = 1;

  // med
  initMedium(med);
  // dif
  initDiff(dif, med);
  // pml
  initPml(pml, med, *dif);
  // ran
  initRange(ran, region, *pml);
  // air
  air->med = med[E_AIR];
  initDimI3(&air->sp, 0, 0, 0);
  initDimI3(&air->range, ran->sr.Txx.x, ran->sr.Txx.y, ran->sr.Txx.z);
  // con
  con->med = med[E_CON];
  // initDimI3(&con->sp, con_start.x + pml->pl1.x - 1, con_start.y + pml->pl1.y - 1, con_start.z + pml->pl1.z - 1);
  initDimI3(&con->sp, 0, 0, 0);
  initDimI3(&con->range, ran->sr.Txx.x, ran->sr.Txx.y, ran->sr.Txx.z);
  // initDimI3(&con->range, con_ran.x, con_ran.y, con_ran.z);
  // 複数欠陥要検討
  // clack
  clack->med = med[E_AIR];
  initDimI3(&clack->sp, clack_start.x + pml->pl1.x - 1, clack_start.y + pml->pl1.y - 1, clack_start.z + pml->pl1.z - 1);
  initDimI3(&clack->range, clack_ran.x, clack_ran.y, clack_ran.z);
}

// 動的変数(静的変数によって大きさが変わる変数)
// メモリ確保後に実行する必要がある
void DynamicVariable(DimD3 *acc, MedArr *ma, Impulse *ip, Range ran, Object air, Object con, Object clack, Pml pml, DimI3 *out, int outnum) {
  // out
  initDimI3(&out[0],  7 + pml.pl1.x - 1, 8 + pml.pl1.y - 1, 8 + pml.pl1.z - 1);
  initDimI3(&out[1],  9 + pml.pl1.x - 1, 8 + pml.pl1.y - 1, 8 + pml.pl1.z - 1);
  // acc
  insertDimD3(acc, outnum);
  // ma
  insertAir(ma, air, ran);
  insertConcrete(ma, con, ran);
  insertClack(ma, clack, ran);
  insertPml(ma, pml, ran);
  // ip
  initDimI3(&ip[0].in, 8 + pml.pl1.x - 1, 8 + pml.pl1.y - 1, 8 + pml.pl1.z - 1);//pml
  ip->freq = 2e6;
  ip->mode = E_SINE;
  // ip->mode = E_RCOS;

}

void insertDimD3(DimD3 *dd3, int outnum) {
  // 0で初期化
  for(int i = 0; i < outnum; i++) {
    dd3[i].x = 0;
    dd3[i].y = 0;
    dd3[i].z = 0;
  }
}

void insertAir(MedArr *ma, Object air, Range ran) {
  int i, j, k;
  Medium objmed = air.med;

  int X = air.sp.x + air.range.x;
  int Y = air.sp.y + air.range.y;
  int Z = air.sp.z + air.range.z;

  int nx = ran.sr.Txx.x;
  int ny = ran.sr.Txx.y;

  for (k = air.sp.z; k < Z; k++) {
    for (j = air.sp.y; j < Y; j++) {
      for (i = air.sp.x; i < X; i++) {
        int idx = k * nx * ny + j * nx + i;
        ma[idx].ramda = objmed.ramda;
        ma[idx].mu = objmed.G;
        ma[idx].c11 = objmed.ramda + 2.0 * objmed.G;
        ma[idx].rho = objmed.rho;
        ma[idx].zetaxx = objmed.zeta;
        ma[idx].zetaxy = objmed.zeta;
        ma[idx].zetaxz = objmed.zeta;
        ma[idx].zetayx = objmed.zeta;
        ma[idx].zetayy = objmed.zeta;
        ma[idx].zetayz = objmed.zeta;
        ma[idx].zetazx = objmed.zeta;
        ma[idx].zetazy = objmed.zeta;
        ma[idx].zetazz = objmed.zeta;
        ma[idx].gamma = objmed.gamma;
        ma[idx].khi = objmed.khi;
        ma[idx].xi11 = objmed.khi + 2.0 * objmed.gamma;
      }
    }
  }
}

void insertConcrete(MedArr *ma, Object con, Range ran) {
  int i, j, k;
  Medium objmed = con.med;

  int X = con.sp.x + con.range.x;
  int Y = con.sp.y + con.range.y;
  int Z = con.sp.z + con.range.z;

  int nx = ran.sr.Txx.x;
  int ny = ran.sr.Txx.y;

  for (k = con.sp.z; k < Z; k++) {
    for (j = con.sp.y; j < Y; j++) {
      for (i = con.sp.x; i < X; i++) {
        int idx = k * nx * ny + j * nx + i;
        ma[idx].ramda = objmed.ramda;
        ma[idx].mu = objmed.G;
        ma[idx].c11 = objmed.ramda + 2.0 * objmed.G;
        ma[idx].rho = objmed.rho;
        ma[idx].zetaxx = objmed.zeta;
        ma[idx].zetaxy = objmed.zeta;
        ma[idx].zetaxz = objmed.zeta;
        ma[idx].zetayx = objmed.zeta;
        ma[idx].zetayy = objmed.zeta;
        ma[idx].zetayz = objmed.zeta;
        ma[idx].zetazx = objmed.zeta;
        ma[idx].zetazy = objmed.zeta;
        ma[idx].zetazz = objmed.zeta;
        ma[idx].gamma = objmed.gamma;
        ma[idx].khi = objmed.khi;
        ma[idx].xi11 = objmed.khi + 2.0 * objmed.gamma;
      }
    }
  }
}

void insertClack(MedArr *ma, Object clack, Range ran) {
  int i, j, k;
  Medium objmed = clack.med;

  int X = clack.sp.x + clack.range.x;
  int Y = clack.sp.y + clack.range.y;
  int Z = clack.sp.z + clack.range.z;

  int nx = ran.sr.Txx.x;
  int ny = ran.sr.Txx.y;

  for (k = clack.sp.z; k < Z; k++) {
    for (j = clack.sp.y; j < Y; j++) {
      for (i = clack.sp.x; i < X; i++) {
        int idx = k * nx * ny + j * nx + i;
        ma[idx].ramda = objmed.ramda;
        ma[idx].mu = objmed.G;
        ma[idx].c11 = objmed.ramda + 2.0 * objmed.G;
        ma[idx].rho = objmed.rho;
        ma[idx].zetaxx = objmed.zeta;
        ma[idx].zetaxy = objmed.zeta;
        ma[idx].zetaxz = objmed.zeta;
        ma[idx].zetayx = objmed.zeta;
        ma[idx].zetayy = objmed.zeta;
        ma[idx].zetayz = objmed.zeta;
        ma[idx].zetazx = objmed.zeta;
        ma[idx].zetazy = objmed.zeta;
        ma[idx].zetazz = objmed.zeta;
        ma[idx].gamma = objmed.gamma;
        ma[idx].khi = objmed.khi;
        ma[idx].xi11 = objmed.khi + 2.0 * objmed.gamma;
      }
    }
  }
}

void insertPml(MedArr *ma, Pml pml, Range ran) {
  int plx1 = pml.pl1.x, plx2 = pml.pl2.x;
  int ply1 = pml.pl1.y, ply2 = pml.pl2.y;
  int plz1 = pml.pl1.z, plz2 = pml.pl2.z;
  int Txximax = ran.sr.Txx.x, Txxjmax = ran.sr.Txx.y, Txxkmax = ran.sr.Txx.z;
  double zeta_max = pml.fm, ta = pml.ta;
  int i, j, k;

  // x方向
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < plx1; i++) {
        int idx = (Txximax * Txxjmax * k) + (Txximax * j) + i;
        ma[idx].zetaxx += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma[idx].zetayx += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma[idx].zetazx += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma[idx].zetadx = ma[idx].zetaxx / ma[idx].rho;
      }
    }
  }
  // x方向
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = Txximax - plx2; i < Txximax; i++) {
        int idx = (Txximax * Txxjmax * k) + (Txximax * j) + i;
        ma[idx].zetaxx += zeta_max * pow(((double)i - (double)(Txximax - 1 - plx2)) / (double)plx2, ta);
        ma[idx].zetayx += zeta_max * pow(((double)i - (double)(Txximax - 1 - plx2)) / (double)plx2, ta);
        ma[idx].zetazx += zeta_max * pow(((double)i - (double)(Txximax - 1 - plx2)) / (double)plx2, ta);
        ma[idx].zetadx = ma[idx].zetaxx / ma[idx].rho;
      }
    }
  }

  // y方向 (前方)
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < ply1; j++) {
      for (i = 0; i < Txximax; i++) {
        int idx = (Txximax * Txxjmax * k) + (Txximax * j) + i;
        ma[idx].zetaxy += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma[idx].zetayy += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma[idx].zetazy += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma[idx].zetady = ma[idx].zetaxy / ma[idx].rho;
      }
    }
  }

  // y方向 (後方)
  for (k = 0; k < Txxkmax; k++) {
    for (j = Txxjmax - ply2; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        int idx = (Txximax * Txxjmax * k) + (Txximax * j) + i;
        ma[idx].zetaxy += zeta_max * pow(((double)j - (double)(Txxjmax - 1 - ply2)) / (double)ply2, ta);
        ma[idx].zetayy += zeta_max * pow(((double)j - (double)(Txxjmax - 1 - ply2)) / (double)ply2, ta);
        ma[idx].zetazy += zeta_max * pow(((double)j - (double)(Txxjmax - 1 - ply2)) / (double)ply2, ta);
        ma[idx].zetady = ma[idx].zetaxy / ma[idx].rho;
      }
    }
  }

  // z方向 (前方)
  for (k = 0; k < plz1; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        int idx = (Txximax * Txxjmax * k) + (Txximax * j) + i;
        ma[idx].zetaxz += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma[idx].zetayz += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma[idx].zetazz += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma[idx].zetadz = ma[idx].zetaxz / ma[idx].rho;
      }
    }
  }

  // z方向 (後方)
  for (k = Txxkmax - plz2; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        int idx = (Txximax * Txxjmax * k) + (Txximax * j) + i;
        ma[idx].zetaxz += zeta_max * pow(((double)k - (double)(Txxkmax - 1 - plz2)) / (double)plz2, ta);
        ma[idx].zetayz += zeta_max * pow(((double)k - (double)(Txxkmax - 1 - plz2)) / (double)plz2, ta);
        ma[idx].zetazz += zeta_max * pow(((double)k - (double)(Txxkmax - 1 - plz2)) / (double)plz2, ta);
        ma[idx].zetadz = ma[idx].zetaxz / ma[idx].rho;
      }
    }
  }
}

