// --------------------------------------------------------------------
//
// title                  :gridding.cu
// description            :Gridding process.
// author                 :
//
// --------------------------------------------------------------------

#include "gridding.h"

/* Initialize output spectrals and weights. */
void init_output(){
    uint32_t num = h_zyx[0] * h_zyx[1] * h_zyx[2];
    h_datacube = RALLOC(double, num);
    h_weightscube = RALLOC(double, num);
    for(uint32_t i = 0; i < num; ++i){
        h_datacube[i] = 0.;
        h_weightscube[i] = 0.;
    }
}

/* Sinc function with simple singularity check. */
double sinc(double x){
    if(fabs(x) < 1.e-10)
        return 1.;
    else
        return sin(x) / x;
}

/* Grid-kernel definitions. get weight*/
double kernel_func_ptr(double distance, double bearing){
    if(h_GMaps.kernel_type == GAUSS1D){   // GAUSS1D
        return exp(-distance * distance * h_kernel_params[0]);
    }
    else if(h_GMaps.kernel_type == GAUSS2D){  // GAUSS2D
        double ellarg = (\
                pow(h_kernel_params[0], 2.0)\
                    * pow(sin(bearing - h_kernel_params[2]), 2.0)\
                + pow(h_kernel_params[1], 2.0)\
                    * pow(cos(bearing - h_kernel_params[2]), 2.0));
        double Earg = pow(distance / h_kernel_params[0] /\
                       h_kernel_params[1], 2.0) / 2. * ellarg;
        return exp(-Earg);
    }
    else if(h_GMaps.kernel_type == TAPERED_SINC){ // TAPERED_SINC
        double arg = PI * distance / h_kernel_params[0];
        return sinc(arg / h_kernel_params[2])\
            * exp(pow(-(arg / h_kernel_params[1]), 2.0));
    }
}

void hcgrid (
        double *h_lons,
        double *h_lats,
        double *h_data,
        double *h_weights,
        double *h_xwcs,
        double *h_ywcs,
        double *h_datacube,
        double *h_weightscube,
        uint64_t *h_hpx_idx) {
        // uint32_t warp_id = blockIdx.x * (blockDim.x / 32) + threadIdx.x / 32;
        // uint32_t tid = ((warp_id % h_GMaps.block_warp_num) * 32 + threadIdx.x % 32) * h_GMaps.factor;
        // printf("\nhere\n"); 
        // printf("%f\n", h_GMaps.sphere_radius);
        uint32_t idx;
        uint32_t xcoord = h_zyx[2];
        uint32_t ycoord = h_zyx[1];
        uint32_t ncoords = h_zyx[1] * h_zyx[2];
        for(idx = 0; idx < h_GMaps.data_shape; idx ++){
            // printf("%d\n", idx);
            double alpha = h_lons[idx] * DEG2RAD;
            double beta = h_lats[idx] * DEG2RAD;
            double in_data = h_data[idx];
            double in_weights = h_weights[idx];
            /* find startpoint*/
            double ubound = h_lats[idx] - h_GMaps.sphere_radius, dbound = h_lats[idx] + h_GMaps.sphere_radius;
            double lbound =  h_lons[idx] + h_GMaps.sphere_radius, rbound = h_lons[idx] - h_GMaps.sphere_radius;
            int lx = 0, rx = 0, uy = 0, dy =0;
            int sy = 0, sp, op;
            /*while(h_xwcs[lx] < lbound && lx < ncoords){
                lx ++;
             }
            lx --;
            while(h_xwcs[rx] < rbound && rx < ncoords){
                rx ++;
            }
            while(h_ywcs[uy] > ubound && uy < ncoords){
                uy += xcoord;
            }
            uy--;
            while(h_ywcs[dy] > dbound && dy < ncoords){
                dy += xcoord;
            }
            for(int xx = 0; xx < xcoord * ycoord; xx ++){
                while(h_xwcs[xx])
            }*/
            while(h_ywcs[sy] < ubound && sy < ncoords){ 
                sy += xcoord;
            }
            if(sy > xcoord){
                sy -= xcoord;
            }
            sp = sy;
            while(h_xwcs[sp] > lbound && sp < ncoords){
                sp ++;
            }
            if(sp > sy){
                sp --;
            }
            /* Gridding*/
            for(int k = sp; h_ywcs[k] < dbound && k < ncoords; k += xcoord){
                for(op = k; h_xwcs[op] > rbound && op < ncoords; op ++){
                    double ga = h_xwcs[op] * DEG2RAD;
                    double  gb = h_ywcs[op] * DEG2RAD;
                    double sdist = true_angular_distance(alpha, beta, ga, gb) * RAD2DEG;
                    double sbear = 0.;
                    if (h_GMaps.bearing_needed) {
                        sbear = great_circle_bearing(alpha, beta, ga, gb);
                    }
                    if(sdist < h_GMaps.sphere_radius){
                        double sweight = kernel_func_ptr(sdist, sbear);
                        double tweight = in_weights * sweight;
                        h_datacube[op] += in_data * tweight;
                        h_weightscube[op] += tweight;
                    }
                }
            }
    }
    return; 

}


/* Gridding process. */
void solve_gridding(const char *infile, const char *tarfile, const char *outfile, const char *sortfile, const int& param, const int &bDim) {
    double iTime1 = cpuSecond();
    // Read input points.
    //reah_input_map_hdf5(infile);
    // printf("\nhere\n");
    read_input_map(infile);
    // Read output map.
    read_output_map(tarfile);

    // Set wcs for output pixels.
    set_WCS();

    // Initialize output spectrals and weights.
    init_output();

//    iTime2 = cpuSecond();
    // Block Indirect Sort i nput points by their healpix indexes.
    // if (param == THRUST) { 
    //     init_input_with_thrust(param);
    // } else {
    //     init_input_with_cpu(param);
    // }

    double iTime3 = cpuSecond();
    // Alloc data for GPU.
    // data_alloc();

    double iTime4 = cpuSecond();
    // Send data from CPU to GPU.
    // data_h2d();
    printf("h_zyx[1]=%d, h_zyx[2]=%d, ", h_zyx[1], h_zyx[2]);
    // for(int i = 0; i < h_zyx[1]; i++){
    //     for(int j = 0; j < h_zyx[2]; j++){
    //         printf("%f ", h_ywcs[i*90 + j]);
    //     }
    //     printf("\n");
    // }

    // Set block and thread.
    // dim3 block(bDim);
    // dim3 grid((h_GMaps.block_warp_num * h_zyx[1] - 1) / (block.x / 32) + 1);
    // printf("grid.x=%d, block.x=%d, ", grid.x, block.x);

    // Get start time.
    // cudaEvent_t start, stop;
    // HANDLE_ERROR(cudaEventCreate(&start));
    // HANDLE_ERROR(cudaEventCreate(&stop));
    // HANDLE_ERROR(cudaEventRecord(start, 0));
    
    hcgrid(h_lons, h_lats, h_data, h_weights, h_xwcs, h_ywcs, h_datacube, h_weightscube, h_hpx_idx);
    
    // Get stop time.
    // printf("kernel elapsed time=%f, ", elapsedTime);

    // Send data from GPU to CPU
    // data_d2h();

    // Write output FITS file
    write_output_map(outfile);

    // Write sorted input FITS file
    if (sortfile) {
        write_ordered_map(infile, sortfile);
    }

    // Release data
    // data_free();
    // HANDLE_ERROR( cudaEventDestroy(start) );
    // HANDLE_ERROR( cudaEventDestroy(stop) );
    // HANDLE_ERROR( cudaDeviceReset() );

    double iTime5 = cpuSecond();
    double iElaps = (iTime5 - iTime1) * 1000.;
    printf("solving_gridding time=%f\n", iElaps);
}