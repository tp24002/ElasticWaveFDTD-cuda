int i, j, k;
  Medium objmed = air.med;
  for (k = air.sp.z; k < air.sp.z + air.range.z; k++) {
    for (j = air.sp.y; j < air.sp.y + air.range.y; j++) {
      for (i = air.sp.x; i < air.sp.x + air.range.x; i++) {
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