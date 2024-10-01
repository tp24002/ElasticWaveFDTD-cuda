#include "../header/struct.h"
#include "../header/init.h"
#include "../header/parameter.h"
#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>

#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))


// 静的変数
void StaticVariable(Medium *med, Pml *pml, Range *ran, Diff *dif, Object *air, Object *con, Object *clack, int *tmax, int *outNum) {
  // 入力
  Coord region;
  initCoord(&region, 10, 10, 10);
  Coord con_start;
  Coord con_ran;
  initCoord(&con_start, 3, 3, 3);
  initCoord(&con_ran, 3, 3, 3);
  Coord clack_start;
  Coord clack_ran;
  initCoord(&clack_start, 3, 3, 3);
  initCoord(&clack_ran, 1, 1, 1);
  *tmax = 32768;
  *outNum = 2;
  

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
  initCoord(&air->sp, 0, 0, 0);
  initCoord(&air->range, ran->sr.Txx.x, ran->sr.Txx.y, ran->sr.Txx.z);
  // con
  con->med = med[E_CON];
  initCoord(&con->sp, con_start.x, con_start.y, con_start.z);
  initCoord(&con->range, con_ran.x, con_ran.y, con_ran.z);
  // 複数欠陥要検討
  // clack
  clack->med = med[E_AIR];
  initCoord(&clack->sp, clack_start.x, clack_start.y, clack_start.z);
  initCoord(&clack->range, clack_ran.x, clack_ran.y, clack_ran.z);
}

// 動的変数(静的変数によって大きさが変わる変数)
// メモリ確保後に実行する必要がある
void DynamicVariable(AccCoord *acc, MedArr *ma, Impulse *ip, Range ran, Object air, Object con, Object clack, Pml pml, Coord *out, int outNum) {
  // out
  initCoord(&out[0], 3, 3, 3);
  initCoord(&out[1], 10, 10, 10);
  // acc
  insertAccCoord(acc, outNum);
  // ma
  insertAir(ma, air, ran);
  insertConcrete(ma, con, ran);
  insertClack(ma, clack, ran);
  insertPml(ma, pml, ran);
  // ip
  ip->freq = 2.0e6;
  ip->mode = E_RCOS;
  initCoord(&ip->in, 10, 10, 10);//pml

}

void insertAccCoord(AccCoord *acc, int outNum) {
  // 0で初期化
  for(int i = 0; i < outNum; i++) {
    acc[i].x = 0;
    acc[i].y = 0;
    acc[i].z = 0;
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

        *(ma->ramda + idx) = objmed.ramda;
        *(ma->mu + idx) = objmed.G;
        *(ma->c11 + idx) = objmed.ramda + 2.0 * objmed.G;
        *(ma->rho + idx) = objmed.rho;
        *(ma->zetaxx + idx) = objmed.zeta;
        *(ma->zetaxy + idx) = objmed.zeta;
        *(ma->zetaxz + idx) = objmed.zeta;
        *(ma->zetayx + idx) = objmed.zeta;
        *(ma->zetayy + idx) = objmed.zeta;
        *(ma->zetayz + idx) = objmed.zeta;
        *(ma->zetazx + idx) = objmed.zeta;
        *(ma->zetazy + idx) = objmed.zeta;
        *(ma->zetazz + idx) = objmed.zeta;
        *(ma->gamma + idx) = objmed.gamma;
        *(ma->khi + idx) = objmed.khi;
        *(ma->xi11 + idx) = objmed.khi + 2.0 * objmed.gamma;
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

        *(ma->ramda + idx) = objmed.ramda;
        *(ma->mu + idx) = objmed.G;
        *(ma->c11 + idx) = objmed.ramda + 2.0 * objmed.G;
        *(ma->rho + idx) = objmed.rho;
        *(ma->zetaxx + idx) = objmed.zeta;
        *(ma->zetaxy + idx) = objmed.zeta;
        *(ma->zetaxz + idx) = objmed.zeta;
        *(ma->zetayx + idx) = objmed.zeta;
        *(ma->zetayy + idx) = objmed.zeta;
        *(ma->zetayz + idx) = objmed.zeta;
        *(ma->zetazx + idx) = objmed.zeta;
        *(ma->zetazy + idx) = objmed.zeta;
        *(ma->zetazz + idx) = objmed.zeta;
        *(ma->gamma + idx) = objmed.gamma;
        *(ma->khi + idx) = objmed.khi;
        *(ma->xi11 + idx) = objmed.khi + 2.0 * objmed.gamma;
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

        *(ma->ramda + idx) = objmed.ramda;
        *(ma->mu + idx) = objmed.G;
        *(ma->c11 + idx) = objmed.ramda + 2.0 * objmed.G;
        *(ma->rho + idx) = objmed.rho;
        *(ma->zetaxx + idx) = objmed.zeta;
        *(ma->zetaxy + idx) = objmed.zeta;
        *(ma->zetaxz + idx) = objmed.zeta;
        *(ma->zetayx + idx) = objmed.zeta;
        *(ma->zetayy + idx) = objmed.zeta;
        *(ma->zetayz + idx) = objmed.zeta;
        *(ma->zetazx + idx) = objmed.zeta;
        *(ma->zetazy + idx) = objmed.zeta;
        *(ma->zetazz + idx) = objmed.zeta;
        *(ma->gamma + idx) = objmed.gamma;
        *(ma->khi + idx) = objmed.khi;
        *(ma->xi11 + idx) = objmed.khi + 2.0 * objmed.gamma;
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
        *(ma->zetaxx + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        *(ma->zetayx + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        *(ma->zetazx + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        *(ma->zetadx + (Txximax * Txxjmax * k) + (Txximax * j) + i) = *(ma->zetaxx + (Txximax * Txxjmax * k) + (Txximax * j) + i) / *(ma->rho + (Txximax * Txxjmax * k) + (Txximax * j) + i);
      }
    }
  }

  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = Txximax - 1; i > Txximax - 1 - plx2; i--) {
        *(ma->zetaxx + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        *(ma->zetayx + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        *(ma->zetazx + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        *(ma->zetadx + (Txximax * Txxjmax * k) + (Txximax * j) + i) = *(ma->zetaxx + (Txximax * Txxjmax * k) + (Txximax * j) + i) / *(ma->rho + (Txximax * Txxjmax * k) + (Txximax * j) + i);
      }
    }
  }

  // y方向
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < ply1; j++) {
      for (i = 0; i < Txximax; i++) {
        *(ma->zetaxy + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        *(ma->zetayy + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        *(ma->zetazy + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        *(ma->zetady + (Txximax * Txxjmax * k) + (Txximax * j) + i) = *(ma->zetaxy + (Txximax * Txxjmax * k) + (Txximax * j) + i) / *(ma->rho + (Txximax * Txxjmax * k) + (Txximax * j) + i);
      }
    }
  }

  for (k = 0; k < Txxkmax; k++) {
    for (j = Txxjmax - 1; j > Txxjmax - 1 - ply2; j--) {
      for (i = 0; i < Txximax; i++) {
        *(ma->zetaxy + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        *(ma->zetayy + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        *(ma->zetazy + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        *(ma->zetady + (Txximax * Txxjmax * k) + (Txximax * j) + i) = *(ma->zetaxy + (Txximax * Txxjmax * k) + (Txximax * j) + i) / *(ma->rho + (Txximax * Txxjmax * k) + (Txximax * j) + i);
      }
    }
  }

  // z方向
  for (k = 0; k < plz1; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        *(ma->zetaxz + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        *(ma->zetayz + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        *(ma->zetazz + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        *(ma->zetadz + (Txximax * Txxjmax * k) + (Txximax * j) + i) = *(ma->zetaxz + (Txximax * Txxjmax * k) + (Txximax * j) + i) / *(ma->rho + (Txximax * Txxjmax * k) + (Txximax * j) + i);
      }
    }
  }

  for (k = Txxkmax - 1; k > Txxkmax - 1 - plz2; k--) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        *(ma->zetaxz + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        *(ma->zetayz + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        *(ma->zetazz + (Txximax * Txxjmax * k) + (Txximax * j) + i) += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        *(ma->zetadz + (Txximax * Txxjmax * k) + (Txximax * j) + i) = *(ma->zetaxz + (Txximax * Txxjmax * k) + (Txximax * j) + i) / *(ma->rho + (Txximax * Txxjmax * k) + (Txximax * j) + i);
      }
    }
  }
}

void insertImpulse(Impulse *ip, Diff dif, int t, Range ran) {
  int nx = ran.sr.Txx.x;
  int ny = ran.sr.Txx.y;

  // 正しいインデックス計算
  int idx = ip->in.z * nx * ny + ip->in.y * nx + ip->in.x;

  if (ip->mode == E_SINE) {
    *(ip->Txx + idx) = 0;
    *(ip->Tyy + idx) = 0;
    *(ip->Tzz + idx) = 0;
  } else if (ip->mode == E_RCOS) {
    if (t < 1. / ip->freq / dif.dt) {
      *(ip->Txx + idx) = 0;
      *(ip->Tyy + idx) = 0;
      *(ip->Tzz + idx) = 8.e3 * 0.5 * (1. - cos(2. * M_PI * ip->freq * (double)t * dif.dt)) / 2.;
    } else {
      *(ip->Txx + idx) = 0;
      *(ip->Tyy + idx) = 0;
      *(ip->Tzz + idx) = 0;
    }
  } else {
    *(ip->Txx + idx) = 0;
    *(ip->Tyy + idx) = 0;
    *(ip->Tzz + idx) = 0;
  }
}
