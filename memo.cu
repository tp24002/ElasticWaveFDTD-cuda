void MemoryBefAftToDevice(BefAft *h, BefAft *d, Range ran) {
    cudaMemcpy(&h, &d, sizeof(BefAft), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h->sa, &d->sa, sizeof(SigArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h->ta, &d->ta, sizeof(TauArr), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h->va, &d->va, sizeof(VelArr), cudaMemcpyDeviceToHost);
    for(int i = 0; i < ran.sr.Txx.x; i++) {
        for(int j = 0; j < ran.sr.Txx.y; j++){
            cudaMemcpy(&h->sa.Txx, &d->sa.Txx, sizeof(double) * ran.sr.Txx.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Txxx, &d->sa.Txxx, sizeof(double) * ran.sr.Txx.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Txxy, &d->sa.Txxy, sizeof(double) * ran.sr.Txx.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Txxz, &d->sa.Txxz, sizeof(double) * ran.sr.Txx.z, cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i < ran.sr.Tyy.x; i++) {
        for(int j = 0; j < ran.sr.Tyy.y; j++){
            cudaMemcpy(&h->sa.Tyy[i][j], &d->sa.Tyy[i][j], sizeof(double) * ran.sr.Tyy.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tyyx[i][j], &d->sa.Tyyx[i][j], sizeof(double) * ran.sr.Tyy.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tyyy[i][j], &d->sa.Tyyy[i][j], sizeof(double) * ran.sr.Tyy.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tyyz[i][j], &d->sa.Tyyz[i][j], sizeof(double) * ran.sr.Tyy.z, cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i < ran.sr.Tzz.x; i++) {
        for(int j = 0; j < ran.sr.Tzz.y; j++){
            cudaMemcpy(&h->sa.Tzz[i][j], &d->sa.Tzz[i][j], sizeof(double) * ran.sr.Tzz.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tzzx[i][j], &d->sa.Tzzx[i][j], sizeof(double) * ran.sr.Tzz.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tzzy[i][j], &d->sa.Tzzy[i][j], sizeof(double) * ran.sr.Tzz.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->sa.Tzzz[i][j], &d->sa.Tzzz[i][j], sizeof(double) * ran.sr.Tzz.z, cudaMemcpyDeviceToHost);
        }
    }
    //////
    for(int i = 0; i < ran.tr.Txy.x; i++) {
        for(int j = 0; j < ran.tr.Txy.y; j++){
            cudaMemcpy(&h->ta.Txy[i][j], &d->ta.Txy[i][j], sizeof(double) * ran.tr.Txy.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Txyx[i][j], &d->ta.Txyx[i][j], sizeof(double) * ran.tr.Txy.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Txyy[i][j], &d->ta.Txyy[i][j], sizeof(double) * ran.tr.Txy.z, cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i < ran.tr.Tyz.x; i++) {
        for(int j = 0; j < ran.tr.Tyz.y; j++){
            cudaMemcpy(&h->ta.Tyz[i][j], &d->ta.Tyz[i][j], sizeof(double) * ran.tr.Tyz.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Tyzy[i][j], &d->ta.Tyzy[i][j], sizeof(double) * ran.tr.Tyz.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Tyzz[i][j], &d->ta.Tyzz[i][j], sizeof(double) * ran.tr.Tyz.z, cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i < ran.tr.Tzx.x; i++) {
        for(int j = 0; j < ran.tr.Tzx.y; j++){
            cudaMemcpy(&h->ta.Tzx[i][j], &d->ta.Tzx[i][j], sizeof(double) * ran.tr.Tzx.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Tzxz[i][j], &d->ta.Tzxz[i][j], sizeof(double) * ran.tr.Tzx.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->ta.Tzxx[i][j], &d->ta.Tzxx[i][j], sizeof(double) * ran.tr.Tzx.z, cudaMemcpyDeviceToHost);
        }
    }
    ///////
    for(int i = 0; i < ran.vr.Vx.x; i++) {
        for(int j = 0; j < ran.vr.Vx.y; j++){
            cudaMemcpy(&h->va.Vx, &d->va.Vx, sizeof(double) * ran.vr.Vx.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vxx, &d->va.Vxx, sizeof(double) * ran.vr.Vx.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vxy, &d->va.Vxy, sizeof(double) * ran.vr.Vx.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vxz, &d->va.Vxz, sizeof(double) * ran.vr.Vx.z, cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i < ran.vr.Vy.x; i++) {
        for(int j = 0; j < ran.vr.Vy.y; j++){
            cudaMemcpy(&h->va.Vy, &d->va.Vy, sizeof(double) * ran.vr.Vy.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vyx, &d->va.Vyx, sizeof(double) * ran.vr.Vy.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vyy, &d->va.Vyy, sizeof(double) * ran.vr.Vy.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vyz, &d->va.Vyz, sizeof(double) * ran.vr.Vy.z, cudaMemcpyDeviceToHost);
        }
    }
    for(int i = 0; i < ran.vr.Vz.x; i++) {
        for(int j = 0; j < ran.vr.Vz.y; j++) {
            cudaMemcpy(&h->va.Vz, &d->va.Vz, sizeof(double) * ran.vr.Vz.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vzx, &d->va.Vzx, sizeof(double) * ran.vr.Vz.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vzy, &d->va.Vzy, sizeof(double) * ran.vr.Vz.z, cudaMemcpyDeviceToHost);
            cudaMemcpy(&h->va.Vzz, &d->va.Vzz, sizeof(double) * ran.vr.Vz.z, cudaMemcpyDeviceToHost);
        }
    }
}