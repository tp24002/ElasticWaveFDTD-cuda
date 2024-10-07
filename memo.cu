  int idx;
  for(int i = 0; i < ran_h.sr.Txx.x; i++) {
    for(int j = 0; j < ran_h.sr.Txx.y; j++) {
      idx = 40 * ran_h.sr.Txx.x * ran_h.sr.Txx.y + j * ran_h.sr.Txx.x + i;
      // if(ma_h->rho[idx] == 1.205){
      //   printf("o");
      // } else {
      //   printf("x");
      // }
      printf("%d", i);
    }
    printf("\n");
  }