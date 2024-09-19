#define _USE_MATH_DEFINES
#include "../header/insert.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "../header/struct.h"

void zeroPadSig(SigArr *sa, SigRan sr) {
  int i, j, k;
  for (k = 0; k < sr.Txx.z; k++) {
    for (j = 0; j < sr.Txx.y; j++) {
      for (i = 0; i < sr.Txx.x; i++) {
        sa->Txx[i][j][k] = 0.;
        sa->Txxx[i][j][k] = 0.;
        sa->Txxy[i][j][k] = 0.;
        sa->Txxz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < sr.Tyy.z; k++) {
    for (j = 0; j < sr.Tyy.y; j++) {
      for (i = 0; i < sr.Tyy.x; i++) {
        sa->Tyy[i][j][k] = 0.;
        sa->Tyyx[i][j][k] = 0.;
        sa->Tyyy[i][j][k] = 0.;
        sa->Tyyz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < sr.Tzz.z; k++) {
    for (j = 0; j < sr.Tzz.y; j++) {
      for (i = 0; i < sr.Tzz.x; i++) {
        sa->Tzz[i][j][k] = 0.;
        sa->Tzzx[i][j][k] = 0.;
        sa->Tzzy[i][j][k] = 0.;
        sa->Tzzz[i][j][k] = 0.;
      }
    }
  }
}

void zeroPadTau(TauArr *ta, TauRan tr) {
  int i, j, k;
  for (k = 0; k < tr.Txy.z; k++) {
    for (j = 0; j < tr.Txy.y; j++) {
      for (i = 0; i < tr.Txy.x; i++) {
        ta->Txy[i][j][k] = 0.;
        ta->Txyx[i][j][k] = 0.;
        ta->Txyy[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < tr.Tyz.z; k++) {
    for (j = 0; j < tr.Tyz.y; j++) {
      for (i = 0; i < tr.Tyz.x; i++) {
        ta->Tyz[i][j][k] = 0.;
        ta->Tyzy[i][j][k] = 0.;
        ta->Tyzz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < tr.Tzx.z; k++) {
    for (j = 0; j < tr.Tzx.y; j++) {
      for (i = 0; i < tr.Tzx.x; i++) {
        ta->Tzx[i][j][k] = 0.;
        ta->Tzxz[i][j][k] = 0.;
        ta->Tzxx[i][j][k] = 0.;
      }
    }
  }
}

void zeroPadVel(VelArr *va, VelRan vr) {
  int i, j, k;
  for (k = 0; k < vr.Vx.z; k++) {
    for (j = 0; j < vr.Vx.y; j++) {
      for (i = 0; i < vr.Vx.x; i++) {
        va->Vx[i][j][k] = 0.;
        va->Vxx[i][j][k] = 0.;
        va->Vxy[i][j][k] = 0.;
        va->Vxz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < vr.Vy.z; k++) {
    for (j = 0; j < vr.Vy.y; j++) {
      for (i = 0; i < vr.Vy.x; i++) {
        va->Vy[i][j][k] = 0.;
        va->Vyx[i][j][k] = 0.;
        va->Vyy[i][j][k] = 0.;
        va->Vyz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < vr.Vz.z; k++) {
    for (j = 0; j < vr.Vz.y; j++) {
      for (i = 0; i < vr.Vz.x; i++) {
        va->Vz[i][j][k] = 0.;
        va->Vzx[i][j][k] = 0.;
        va->Vzy[i][j][k] = 0.;
        va->Vzz[i][j][k] = 0.;
      }
    }
  }
}

void zeroPadding(BefAft *ba, Range ran) {
  zeroPadSig(&ba->sa, ran.sr);
  zeroPadTau(&ba->ta, ran.tr);
  zeroPadVel(&ba->va, ran.vr);
}

void zeroPadImpulse(Impulse *ip, Range ran) {
  int i, j, k;
  for (k = 0; k < ran.sr.Txx.z; k++) {
    for (j = 0; j < ran.sr.Txx.y; j++) {
      for (i = 0; i < ran.sr.Txx.x; i++) {
        ip->Txx[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < ran.sr.Tyy.z; k++) {
    for (j = 0; j < ran.sr.Tyy.y; j++) {
      for (i = 0; i < ran.sr.Tyy.x; i++) {
        ip->Tyy[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < ran.sr.Tzz.z; k++) {
    for (j = 0; j < ran.sr.Tzz.y; j++) {
      for (i = 0; i < ran.sr.Tzz.x; i++) {
        ip->Tzz[i][j][k] = 0.;
      }
    }
  }
  ip->freq = 0;
  ip->mode = 0;
  ip->in.x = 0;
  ip->in.y = 0;
  ip->in.z = 0;
}