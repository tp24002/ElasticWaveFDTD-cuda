void transferBefAft(BefAft *ba_d, BefAft *ba_h, Range ran) {
    int N = ran.sr.Txx.x * ran.sr.Txx.y * ran.sr.Txx.z;
    // Step 1: デバイスからホストへ構造体をコピー
    cudaMemcpy(ba_h, ba_d, sizeof(BefAft), cudaMemcpyDeviceToHost);

    // Step 2: デバイス側のポインタを保存
    double *device_Txx = ba_h->sa.Txx;
    double *device_Tyy = ba_h->sa.Tyy;
    double *device_Tzz = ba_h->sa.Tzz;

    double *device_Txy = ba_h->ta.Txy;
    double *device_Tyz = ba_h->ta.Tyz;
    double *device_Tzx = ba_h->ta.Tzx;

    double *device_Vx = ba_h->va.Vx;
    double *device_Vy = ba_h->va.Vy;
    double *device_Vz = ba_h->va.Vz;

    // Step 3: ホスト側でメモリを確保
    ba_h->sa.Txx = (double *)malloc(N * sizeof(double));
    ba_h->sa.Tyy = (double *)malloc(N * sizeof(double));
    ba_h->sa.Tzz = (double *)malloc(N * sizeof(double));

    ba_h->ta.Txy = (double *)malloc(N * sizeof(double));
    ba_h->ta.Tyz = (double *)malloc(N * sizeof(double));
    ba_h->ta.Tzx = (double *)malloc(N * sizeof(double));

    ba_h->va.Vx = (double *)malloc(N * sizeof(double));
    ba_h->va.Vy = (double *)malloc(N * sizeof(double));
    ba_h->va.Vz = (double *)malloc(N * sizeof(double));

    // Step 4: デバイスからホストへデータをコピー
    cudaMemcpy(ba_h->sa.Txx, device_Txx, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->sa.Tyy, device_Tyy, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->sa.Tzz, device_Tzz, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaMemcpy(ba_h->ta.Txy, device_Txy, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->ta.Tyz, device_Tyz, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->ta.Tzx, device_Tzx, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaMemcpy(ba_h->va.Vx, device_Vx, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->va.Vy, device_Vy, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ba_h->va.Vz, device_Vz, N * sizeof(double), cudaMemcpyDeviceToHost);
}
