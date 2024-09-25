#include "../header/struct.h"
#include "../header/init.h"
#include "../header/parameter.h"
#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>

#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))


// 静的変数
__global__ void StaticVariable(Medium *med, Pml *pml, Range *ran, Diff *dif, Object *air, Object *con, Object *clack, MedArr *ma, int *tmax, int *outNum) {
    // 入力
    Coord region;
    initDeviceCoord(&region, 10, 10, 10);
    Coord con_start;
    Coord con_ran;
    initDeviceCoord(&con_start, 3, 3, 3);
    initDeviceCoord(&con_ran, 3, 3, 3);
    Coord clack_start;
    Coord clack_ran;
    initDeviceCoord(&clack_start, 3, 3, 3);
    initDeviceCoord(&clack_ran, 1, 1, 1);
    *tmax = 2;
    *outNum = 2;

    // med
    initMedium<<<1,1>>>(med);
    printf("Med:%f\n", med[0].rho);
    printf("Med:%f\n", med[1].rho);
    // pml
    initPml<<<1,1>>>(pml, med, *dif);
    printf("pml.fm:%f\n", pml->fm);
    printf("pml.ta:%f\n", pml->ta);
    // ran
    initRange<<<1,1>>>(ran, region, *pml);
    printf("Range:%d,%d,%d\n", ran->sr.Txx.x, ran->sr.Txx.y, ran->sr.Txx.z);
    // dif
    initDiff<<<1,1>>>(dif, med);
    printf("sp   diff:%f,%f,%f\n", dif->dx, dif->dy, dif->dz);
    printf("time diff:%e\n", dif->dt);
    // air
    air->med = med[E_AIR];
    initDeviceCoord(&air->sp, 0, 0, 0);
    initDeviceCoord(&air->range, ran->sr.Txx.x, ran->sr.Txx.y, ran->sr.Txx.z);
    // con
    con->med = med[E_CON];
    initDeviceCoord(&con->sp, con_start.x, con_start.y, con_start.z);
    initDeviceCoord(&con->range, con_ran.x, con_ran.y, con_ran.z);
    // 複数欠陥要検討
    // clack
    clack->med = med[E_AIR];
    initDeviceCoord(&clack->sp, clack_start.x, clack_start.y, clack_start.z);
    initDeviceCoord(&clack->range, clack_ran.x, clack_ran.y, clack_ran.z);
}

// 動的変数(静的変数によって大きさが変わる変数)
// メモリ確保後に実行する必要がある
__global__ void DynamicVariable(BefAft *bef, BefAft *aft, AccCoord *acc, MedArr *ma, Impulse *ip, Range *ran, Medium *med, Object *air, Object *con, Object *clack, Pml *pml, Diff *dif, Coord *out, int *outNum) {
  
  // befaft
  insertBefAft(bef, ran);
  insertBefAft(aft, ran);
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
  initDeviceCoord(&ip->in, 10, 10, 10);

  initDeviceCoord(&out[0], 3, 3, 3);
  initDeviceCoord(&out[1], 2, 4, 3);
}

__device__ void insertBefAft(BefAft *ba, Range *ran) {
  // SigArr 初期化
  int sx = ran->sr.Txx.x, sy = ran->sr.Txx.y, sz = ran->sr.Txx.z;
  memset(ba->sa.Txx, 0, sizeof(double) * sx * sy * sz);
  memset(ba->sa.Txxx, 0, sizeof(double) * sx * sy * sz);
  memset(ba->sa.Txxy, 0, sizeof(double) * sx * sy * sz);
  memset(ba->sa.Txxz, 0, sizeof(double) * sx * sy * sz);
  
  int syy_x = ran->sr.Tyy.x, syy_y = ran->sr.Tyy.y, syy_z = ran->sr.Tyy.z;
  memset(ba->sa.Tyy, 0, sizeof(double) * syy_x * syy_y * syy_z);
  memset(ba->sa.Tyyx, 0, sizeof(double) * syy_x * syy_y * syy_z);
  memset(ba->sa.Tyyy, 0, sizeof(double) * syy_x * syy_y * syy_z);
  memset(ba->sa.Tyyz, 0, sizeof(double) * syy_x * syy_y * syy_z);

  int szz_x = ran->sr.Tzz.x, szz_y = ran->sr.Tzz.y, szz_z = ran->sr.Tzz.z;
  memset(ba->sa.Tzz, 0, sizeof(double) * szz_x * szz_y * szz_z);
  memset(ba->sa.Tzzx, 0, sizeof(double) * szz_x * szz_y * szz_z);
  memset(ba->sa.Tzzy, 0, sizeof(double) * szz_x * szz_y * szz_z);
  memset(ba->sa.Tzzz, 0, sizeof(double) * szz_x * szz_y * szz_z);

  // TauArr 初期化
  int txy_x = ran->tr.Txy.x, txy_y = ran->tr.Txy.y, txy_z = ran->tr.Txy.z;
  memset(ba->ta.Txy, 0, sizeof(double) * txy_x * txy_y * txy_z);
  memset(ba->ta.Txyx, 0, sizeof(double) * txy_x * txy_y * txy_z);
  memset(ba->ta.Txyy, 0, sizeof(double) * txy_x * txy_y * txy_z);

  int tyz_x = ran->tr.Tyz.x, tyz_y = ran->tr.Tyz.y, tyz_z = ran->tr.Tyz.z;
  memset(ba->ta.Tyz, 0, sizeof(double) * tyz_x * tyz_y * tyz_z);
  memset(ba->ta.Tyzy, 0, sizeof(double) * tyz_x * tyz_y * tyz_z);
  memset(ba->ta.Tyzz, 0, sizeof(double) * tyz_x * tyz_y * tyz_z);

  int tzx_x = ran->tr.Tzx.x, tzx_y = ran->tr.Tzx.y, tzx_z = ran->tr.Tzx.z;
  memset(ba->ta.Tzx, 0, sizeof(double) * tzx_x * tzx_y * tzx_z);
  memset(ba->ta.Tzxz, 0, sizeof(double) * tzx_x * tzx_y * tzx_z);
  memset(ba->ta.Tzxx, 0, sizeof(double) * tzx_x * tzx_y * tzx_z);

  // VelArr 初期化
  int vx_x = ran->vr.Vx.x, vx_y = ran->vr.Vx.y, vx_z = ran->vr.Vx.z;
  memset(ba->va.Vx, 0, sizeof(double) * vx_x * vx_y * vx_z);
  memset(ba->va.Vxx, 0, sizeof(double) * vx_x * vx_y * vx_z);
  memset(ba->va.Vxy, 0, sizeof(double) * vx_x * vx_y * vx_z);
  memset(ba->va.Vxz, 0, sizeof(double) * vx_x * vx_y * vx_z);

  int vy_x = ran->vr.Vy.x, vy_y = ran->vr.Vy.y, vy_z = ran->vr.Vy.z;
  memset(ba->va.Vy, 0, sizeof(double) * vy_x * vy_y * vy_z);
  memset(ba->va.Vyx, 0, sizeof(double) * vy_x * vy_y * vy_z);
  memset(ba->va.Vyy, 0, sizeof(double) * vy_x * vy_y * vy_z);
  memset(ba->va.Vyz, 0, sizeof(double) * vy_x * vy_y * vy_z);

  int vz_x = ran->vr.Vz.x, vz_y = ran->vr.Vz.y, vz_z = ran->vr.Vz.z;
  memset(ba->va.Vz, 0, sizeof(double) * vz_x * vz_y * vz_z);
  memset(ba->va.Vzx, 0, sizeof(double) * vz_x * vz_y * vz_z);
  memset(ba->va.Vzy, 0, sizeof(double) * vz_x * vz_y * vz_z);
  memset(ba->va.Vzz, 0, sizeof(double) * vz_x * vz_y * vz_z);
}

__device__ void insertAccCoord(AccCoord *acc, int *outNum) {
  // 0で初期化
  for(int i = 0; i < *outNum; i++) {
    acc[i].x = 0;
    acc[i].y = 0;
    acc[i].z = 0;
  }
}

__device__ void insertAir(MedArr *ma, Object *air, Range *ran) {
    int i, j, k;
    Medium objmed = air->med;

    int X = air->sp.x + air->range.x;
    int Y = air->sp.y + air->range.y;
    int Z = air->sp.z + air->range.z;


    int ny = ran->sr.Txx.y;
    int nz = ran->sr.Txx.z;

    for (k = air->sp.z; k < Z; k++) {
        for (j = air->sp.y; j < Y; j++) {
            for (i = air->sp.x; i < X; i++) {
                int idx = k * ny * nz + j * nz + i;

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

__device__ void insertConcrete(MedArr *ma, Object *con, Range *ran) {
  int i, j, k;
  Medium objmed = con->med;

  int X = con->sp.x + con->range.x;
  int Y = con->sp.y + con->range.y;
  int Z = con->sp.z + con->range.z;

  int ny = ran->sr.Txx.y;
  int nz = ran->sr.Txx.z;

  for (k = con->sp.z; k < Z; k++) {
    for (j = con->sp.y; j < Y; j++) {
      for (i = con->sp.x; i < X; i++) {
        int idx = k * ny * nz + j * nz + i;

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

__device__ void insertClack(MedArr *ma, Object *clack, Range *ran) {
  Coord clacksp;
  Medium clackmed;
  clackmed = clack->med;
  clacksp = clack->sp;

  int ny = ran->sr.Txx.y;
  int nz = ran->sr.Txx.z;

  for (int z = 0; z < clack->range.z; z++) {
    for (int y = 0; y < clack->range.y; y++) {
      for (int x = 0; x < clack->range.x; x++) {
        int idx = (clacksp.x + x) * ny * nz + (clacksp.y + y) * nz + (clacksp.z + z);

        *(ma->ramda + idx) = clackmed.ramda;
        *(ma->mu + idx) = clackmed.G;
        *(ma->c11 + idx) = clackmed.ramda + 2.0 * clackmed.G;
        *(ma->rho + idx) = clackmed.rho;
        *(ma->zetaxx + idx) = clackmed.zeta;
        *(ma->zetaxy + idx) = clackmed.zeta;
        *(ma->zetaxz + idx) = clackmed.zeta;
        *(ma->zetayx + idx) = clackmed.zeta;
        *(ma->zetayy + idx) = clackmed.zeta;
        *(ma->zetayz + idx) = clackmed.zeta;
        *(ma->zetazx + idx) = clackmed.zeta;
        *(ma->zetazy + idx) = clackmed.zeta;
        *(ma->zetazz + idx) = clackmed.zeta;
        *(ma->gamma + idx) = clackmed.gamma;
        *(ma->khi + idx) = clackmed.khi;
        *(ma->xi11 + idx) = clackmed.khi + 2.0 * clackmed.gamma;
      }
    }
  }
}

__device__ void insertPml(MedArr *ma, Pml *pml, Range *ran) {
  int plx1 = pml->pl1.x, plx2 = pml->pl2.x;
  int ply1 = pml->pl1.y, ply2 = pml->pl2.y;
  int plz1 = pml->pl1.z, plz2 = pml->pl2.z;
  int Txximax = ran->sr.Txx.x, Txxjmax = ran->sr.Txx.y, Txxkmax = ran->sr.Txx.z;
  double zeta_max = pml->fm, ta = pml->ta;
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

__global__ void insertImpulse(Impulse *ip, Diff *dif, int *t, Range *ran) {
  int ny = ran->sr.Txx.y;
  int nz = ran->sr.Txx.z;

  // 正しいインデックス計算
  int idx = ip->in.x * ny * nz + ip->in.y * nz + ip->in.z;

  if (ip->mode == E_SINE) {
    *(ip->Txx + idx) = 0;
    *(ip->Tyy + idx) = 0;
    *(ip->Tzz + idx) = 0;
  } else if (ip->mode == E_RCOS) {
    if (*t < 1. / ip->freq / dif->dt) {
      *(ip->Txx + idx) = 0;
      *(ip->Tyy + idx) = 0;
      *(ip->Tzz + idx) = 8.e3 * 0.5 * (1. - cos(2. * M_PI * ip->freq * (double)*t * dif->dt)) / 2.;
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
