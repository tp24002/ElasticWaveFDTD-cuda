    cudaMemcpy(&d_medarrptr->ramda, &h_medarrptr->ramda, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->mu, &h_medarrptr->mu, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->c11, &h_medarrptr->c11, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->rho, &h_medarrptr->rho, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->zetaxx, &h_medarrptr->zetaxx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetaxy, &h_medarrptr->zetaxy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetaxz, &h_medarrptr->zetaxz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->zetayx, &h_medarrptr->zetayx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetayy, &h_medarrptr->zetayy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetayz, &h_medarrptr->zetayz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->zetazx, &h_medarrptr->zetazx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetazy, &h_medarrptr->zetazy, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetazz, &h_medarrptr->zetazz, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->gamma, &h_medarrptr->gamma, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->khi, &h_medarrptr->khi, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->xi11, &h_medarrptr->xi11, sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(&d_medarrptr->zetadx, &h_medarrptr->zetadx, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetady, &h_medarrptr->zetady, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&d_medarrptr->zetadz, &h_medarrptr->zetadz, sizeof(double*), cudaMemcpyHostToDevice);