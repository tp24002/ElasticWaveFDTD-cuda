#include "../header/struct.h"
#include "../header/init.h"
#include "../header/parameter.h"
#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>

#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))


// 静的変数
// 将来的にairの格納もobject型にしたい(汎用的)
void StaticVariable(Range *ran, Pml *pml, Diff *dif, Object *con, Object *clack, Medium *med, MedArr *ma, int *tmax) {
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
    *tmax = 32;


    // med
    initMedium(med);
    // pml
    initPml(pml, med, *dif);
    // ran
    initRange(ran, region, *pml);
    printf("%d\n",ran->sr.Txx.x);
    // dif
    initDiff(dif, med);
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
void DynamicVariable(MedArr *ma, Impulse *ip, Range ran, Medium *med, Object con, Object clack, Pml pml, Diff dif, int tmax) {
    // ma
    printf("s00\n");
    insertAir(ma, ran.sr, med[E_AIR]);
    printf("s00\n");
    insertConcrete(ma, con);
    printf("s00\n");
    insertClack(ma, clack);
    printf("00\n");
    insertPml(ma, ran.sr, pml);

    // ip
    printf("00\n");
    ip->freq = 2.0e6;
    printf("00\n");
    ip->mode = E_RCOS;
    printf("00\n");
    initCoord(&ip->in, 10, 10, 10);
    printf("00\n");
    insertImpulse(ip, dif, tmax);
}

void insertAir(MedArr *ma, SigRan sr, Medium air) {
  int i, j, k;
  for (k = 0; k < sr.Txx.z; k++) {
    for (j = 0; j < sr.Txx.y; j++) {
      for (i = 0; i < sr.Txx.x; i++) {
        ma->ramda[i][j][k] = air.ramda;
        ma->mu[i][j][k] = air.G;
        ma->c11[i][j][k] = air.ramda + 2. * air.G;
        ma->rho[i][j][k] = air.rho;
        ma->zetaxx[i][j][k] = air.zeta;
        ma->zetaxy[i][j][k] = air.zeta;
        ma->zetaxz[i][j][k] = air.zeta;
        ma->zetayx[i][j][k] = air.zeta;
        ma->zetayy[i][j][k] = air.zeta;
        ma->zetayz[i][j][k] = air.zeta;
        ma->zetazx[i][j][k] = air.zeta;
        ma->zetazy[i][j][k] = air.zeta;
        ma->zetazz[i][j][k] = air.zeta;
        ma->gamma[i][j][k] = air.gamma;
        ma->khi[i][j][k] = air.khi;
        ma->xi11[i][j][k] = air.khi + 2. * air.gamma;
        ma->zetadx[i][j][k] = 0.;
        ma->zetady[i][j][k] = 0.;
        ma->zetadz[i][j][k] = 0.;
      }
    }
  }
}

void insertConcrete(MedArr *ma, Object con) {
  int i, j, k;
  Medium objmed = con.med;
  for (k = con.sp.z; k < con.sp.z + con.range.z; k++) {
    for (j = con.sp.y; j < con.sp.y + con.range.y; j++) {
      for (i = con.sp.x; i < con.sp.x + con.range.x; i++) {
        ma->ramda[i][j][k] = objmed.ramda;
        ma->mu[i][j][k] = objmed.G;
        ma->c11[i][j][k] = objmed.ramda + 2. * objmed.G;
        ma->rho[i][j][k] = objmed.rho;
        ma->zetaxx[i][j][k] = objmed.zeta;
        ma->zetaxy[i][j][k] = objmed.zeta;
        ma->zetaxz[i][j][k] = objmed.zeta;
        ma->zetayx[i][j][k] = objmed.zeta;
        ma->zetayy[i][j][k] = objmed.zeta;
        ma->zetayz[i][j][k] = objmed.zeta;
        ma->zetazx[i][j][k] = objmed.zeta;
        ma->zetazy[i][j][k] = objmed.zeta;
        ma->zetazz[i][j][k] = objmed.zeta;
        ma->gamma[i][j][k] = objmed.gamma;
        ma->khi[i][j][k] = objmed.khi;
        ma->xi11[i][j][k] = objmed.khi + 2. * objmed.gamma;
      }
    }
  }
}

void insertClack(MedArr *ma, Object clack) {
    Coord clacksp;
    Medium clackmed;
    clackmed = clack.med;
    clacksp = clack.sp;
    for (int z = 0; z < clack.range.z; z++) {
      for (int y = 0; y < clack.range.y; y++) {
        for (int x = 0; x < clack.range.x; x++) {
          ma->ramda[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.ramda;
          ma->mu[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.G;
          ma->c11[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.ramda + 2. * clackmed.G;
          ma->rho[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.rho;
          ma->zetaxx[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.zeta;
          ma->zetaxy[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.zeta;
          ma->zetaxz[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.zeta;
          ma->zetayx[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.zeta;
          ma->zetayy[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.zeta;
          ma->zetayz[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.zeta;
          ma->zetazx[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.zeta;
          ma->zetazy[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.zeta;
          ma->zetazz[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.zeta;
          ma->gamma[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.gamma;
          ma->khi[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.khi;
          ma->xi11[clacksp.x + x][clacksp.y + y][clacksp.z + z] = clackmed.khi + 2. * clackmed.gamma;
        }
      }
    }
}

void insertPml(MedArr *ma, SigRan sr, Pml pml) {
  int plx1 = pml.pl1.x, plx2 = pml.pl2.x;
  int ply1 = pml.pl1.y, ply2 = pml.pl2.y;
  int plz1 = pml.pl1.z, plz2 = pml.pl2.z;
  int Txximax = sr.Txx.x, Txxjmax = sr.Txx.y, Txxkmax = sr.Txx.z;
  double zeta_max = pml.fm, ta = pml.ta;
  int i, j, k;
  //x方向
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < plx1; i++) {
        ma->zetaxx[i][j][k] += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma->zetayx[i][j][k] += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma->zetazx[i][j][k] += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma->zetadx[i][j][k] = ma->zetaxx[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = Txximax - 1; i > Txximax - 1 - plx2; i--) {
        ma->zetaxx[i][j][k] += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        ma->zetayx[i][j][k] += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        ma->zetazx[i][j][k] += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        ma->zetadx[i][j][k] = ma->zetaxx[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  //y方向
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < ply1; j++) {
      for (i = 0; i < Txximax; i++) {
        ma->zetaxy[i][j][k] += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma->zetayy[i][j][k] += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma->zetazy[i][j][k] += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma->zetady[i][j][k] = ma->zetaxy[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  for (k = 0; k < Txxkmax; k++) {
    for (j = Txxjmax - 1; j > Txxjmax - 1 - ply2; j--) {
      for (i = 0; i < Txximax; i++) {
        ma->zetaxy[i][j][k] += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        ma->zetayy[i][j][k] += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        ma->zetazy[i][j][k] += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        ma->zetady[i][j][k] = ma->zetaxy[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  //z方向
  for (k = 0; k < plz1; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        ma->zetaxz[i][j][k] += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma->zetayz[i][j][k] += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma->zetazz[i][j][k] += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma->zetadz[i][j][k] = ma->zetaxz[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  for (k = Txxkmax - 1; k > Txxkmax - 1 - plz2; k--) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        ma->zetaxz[i][j][k] += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        ma->zetayz[i][j][k] += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        ma->zetazz[i][j][k] += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        ma->zetadz[i][j][k] = ma->zetaxz[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
}

void insertImpulse(Impulse *ip, Diff dif, int tmax) {
    for(int i = 0; i < tmax; i++) {
       if (ip->mode == E_SINE) {
        ip->Txx[ip->in.x][ip->in.y][ip->in.z] = 0;
        ip->Tyy[ip->in.x][ip->in.y][ip->in.z] = 0;
        ip->Tzz[ip->in.x][ip->in.y][ip->in.z] = 0;
        } else if (ip->mode == E_RCOS) {
            if (tmax < 1. / ip->freq / dif.dt) {
            ip->Txx[ip->in.x][ip->in.y][ip->in.z] = 0;///* 8.e3 * 0.5 * */(1. - cos(2. * M_PI * ip.freq * (double)t * dif.dt)) / 2.;
            ip->Tyy[ip->in.x][ip->in.y][ip->in.z] = 0;///* 8.e3 * 0.5 * */(1. - cos(2. * M_PI * ip.freq * (double)t * dif.dt)) / 2.;
            ip->Tzz[ip->in.x][ip->in.y][ip->in.z] = 8.e3 * 0.5 * (1. - cos(2. * M_PI * ip->freq * (double)tmax * dif.dt)) / 2.;
            } else {
            ip->Txx[ip->in.x][ip->in.y][ip->in.z] = 0.;
            ip->Tyy[ip->in.x][ip->in.y][ip->in.z] = 0.;
            ip->Tzz[ip->in.x][ip->in.y][ip->in.z] = 0.;
            }
        } else {
            ip->Txx[ip->in.x][ip->in.y][ip->in.z] = 0.;
            ip->Tyy[ip->in.x][ip->in.y][ip->in.z] = 0.;
            ip->Tzz[ip->in.x][ip->in.y][ip->in.z] = 0.;
        } 
    }
}
