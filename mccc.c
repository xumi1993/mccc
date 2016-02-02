/***********************************
 * Function MCCC determines optimum relative delay times for a 
 * set of seismograms based on the VanDecar & Crosson multi-channel
 * cross-correlation algorithm.
 *
 * 	Author:  Mijian Xu
 *
 *  History:
 *      17/01/2016	Initial coding
***********************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sac.h>
#include "sac.h"
#include <fftw3.h>

int main(int argc, char **argv) {
    int i, j, k, m, error, isdetrend, tref, cut, verbose,
        sta_num, npts, nft, nt, order, passes, isfilter,
        nts, nt1, nt2, o;
    float t1, t2, dt, sum1, sum2, t0,
          transition_bandwidth, attenuation, low, high,
          f1, f2;
    float **data, *trace, *tdel, *data_sum, *corrcef;
    char filter[2];
    char **indir, *kname, *proto;
    SACHEAD hd;
    fftw_complex **fft_all, *srcin, *srcout, *ffis, *ffjs;
    fftw_plan p1, p2, p3;
    double *ccf;
    FILE *fp;
    int argmax(double *a, int n);
    float coef(float* a, float* b, int nft);
    
    k = 0;
    cut = 0;
    error = 0;
    isdetrend = 0;
    isfilter = 0;
    verbose = 0;
    transition_bandwidth = 0.0;
    attenuation = 0.0;
    passes = 2;
/* input parameters */
    for (i=1; !error && i < argc; i++) {
        if (argv[i][0] == '-') {
          ++k;
          switch(argv[i][1]) {
            case 'D':
                isdetrend = 1;
                break;
            case 'C':
                cut = 1;
                sscanf(&argv[i][2],"%d/%f/%f",&tref,&t1,&t2);
                break;
            case 'F':
                isfilter = 1;
                j = sscanf(&argv[i][2],"%2s/%d/%f/%f",&filter,&order,&f1,&f2);
                if (j<3) {error=TRUE;}
                if (j==4&&!strcmp(filter,"bp")){
                   proto = strdup("BP");
                   low = f1; high = f2;
                }
                if (j==3&&strcmp(filter,"hp")){
                   proto = strdup("HP");
                   low = f1; high = 0;
                }
                if (j==3&&strcmp(filter,"lp")){
                   proto = strdup("LP");
                   low = 0; high = f1;
                }
                break;
            case 'v':
                verbose = 1;
                break;
            default:
               error = 1;
               break;
          }
        }
    }
    
    if (argc == 1 || error) {
       fprintf(stderr,"Usage: mccc [-D] [-Ffilter/order/f1[/f2]] [-Ctmark/t1/t2] [-v] sac_files\n\
    -D: Remean traces\n\
    -C: window data [tmark+t1,tmark+t2] (off)\n\
	    mark = -5(b), -3(o), -2(a), 0(t0)\n\
    -F: Fiter traces\n\
        filter = [bp|hp|lp]\n\
        order = filter order\n\
        if [hp|lp] specified f1 is the corner freq.\n\
        if [bp] specified f1 f2 are the corner freq.\n\
    -v: Print file name at corr-correlation\n");
       return -1;
    }


/* Read and cut sacfiles */
    sta_num = argc-k-1;
    data = (float **)calloc(sta_num,sizeof(float *));
    indir = (char **)calloc(sta_num,sizeof(char *));
    k = 0;
    for (i=1; i < argc; i++){
        if (argv[i][0] != '-'){
            kname = strdup(&argv[i][0]);
            indir[k] = (char *)calloc(200,sizeof(char));
            indir[k] = strdup(&argv[i][0]); 
            trace = read_sac(kname, &hd);
            switch (tref){
                case -3:
                    t0 = hd.o;
                    break;
                case -5:
                    t0 = hd.b;
                    break;
                case 0:
                    t0 = hd.t0;
                    break;
                case -2:
                    t0 = hd.a;
                    break;
                default:
                    fprintf(stderr,"Error tmark.");
                    return -1;
                    break;
            }
            dt = hd.delta;
//            printf("%f\n",t0);
//            npts = (int) rint((t2-t1)/dt);
            nts = hd.npts;
            o = hd.o;
            if (isdetrend){detrend(trace, nts);}
            if (isfilter){
                xapiir(trace, nts, SAC_BUTTERWORTH,
                transition_bandwidth, attenuation,
                order, proto,
                low, high, 
	            dt, passes);} 
            nt1 = (int) rint((t0-o+t1)/dt);
            nt2 = (int) rint((t0-o+t2)/dt);
            npts = nt2-nt1+1;
            data[k] = (float *)calloc(npts,sizeof(float));
            m = 0;
            for (j=nt1; j<nt2+1; j++){
                data[k][m] = trace[j];
                ++m;
            }
            ++k;
        }
    }

/* Do mccc */
    nft = npts*2;
    int tcc [sta_num][sta_num];
    float tcc_f [sta_num][sta_num];
    for (i=0; i<sta_num; i++){
        for (j=0; j<sta_num; j++){
            tcc[i][j] = 0;
            tcc_f[i][j] = 0.;
        }
    }
    ccf = (double *)calloc(nft,sizeof(double *));
    srcin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nft);
    srcout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nft);
    fft_all = (fftw_complex**) fftw_malloc(sizeof(fftw_complex*)*nft);  
    ffis = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nft);
    ffjs = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nft);

    for (i=0; i<sta_num; i++){
        fft_all[i] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nft);
        for (k=0; k<npts; k++){
            srcin[k][0] = (double)data[i][k];
            srcin[k][1] = 0.;
        }
        for (k=npts; k<nft; k++){
            srcin[k][0] = 0.;
            srcin[k][1] = 0.;
        }
        p1 = fftw_plan_dft_1d(nft,srcin,srcout,FFTW_FORWARD,FFTW_ESTIMATE);
        fftw_execute(p1);
        for (k=0; k<nft; k++){
            fft_all[i][k][0] = srcout[k][0];
            fft_all[i][k][1] = srcout[k][1];
        }
        fftw_destroy_plan(p1);
    }
    for (i=0; i<sta_num-1; i++){
        if (verbose==1){printf("%s\t%d/%d\n",indir[i],i+1,sta_num);}
        ffis = fft_all[i];
        for (k=0; k<nft; k++){
            ffis[k][1] = -ffis[k][1];           
        }
        for (j=i+1; j<sta_num; j++){
            ffjs = fft_all[j];
            for (k=0; k<nft; k++){
                srcin[k][0] = ffis[k][0]*ffjs[k][0]-ffis[k][1]*ffjs[k][1];
                srcin[k][1] = ffis[k][0]*ffjs[k][1]+ffis[k][1]*ffjs[k][0];
            }
            p3 = fftw_plan_dft_c2r_1d(nft,srcin,ccf,FFTW_ESTIMATE);
            fftw_execute(p3);
            for (k=0; k<nft; k++){
                ccf[k] /= nft;
            }
            tcc[i][j] = argmax(ccf,nft);
            fftw_destroy_plan(p3);
        }
    }
    if(verbose==1) printf("%s\t%d/%d\n",indir[sta_num-1],sta_num,sta_num);

    for (i=0; i<sta_num; i++){
        for (j=0; j<sta_num; j++){
            if (tcc[i][j]>npts){
                tcc[i][j] = tcc[i][j]-(nft+1);
            }            
            tcc_f[i][j] = tcc[i][j] * dt;
        }
    }
    tdel = (float *)calloc(sta_num,sizeof(float *));
    for (i=0; i<sta_num; i++){
        sum1 = 0.;
        sum2 = 0.;
        for (j=0; j<i; j++){
            sum1 += tcc_f[j][i];
        }
        for (j=i; j<sta_num; j++){
            sum2 += tcc_f[i][j];
        }
        tdel[i] = (-sum1+sum2)/sta_num;
    }

/* Calculate correlation coefficient */
    data_sum = (float *)calloc(npts,sizeof(float));
    for (i=0; i<npts; i++){
        data_sum[i] = 0.;
        for (j=0; j<sta_num; j++){
            data_sum[i] += data[j][i];
        }
    }
    corrcef = (float *)calloc(sta_num,sizeof(float));
    for (i=0; i<sta_num; i++){
        corrcef[i] = coef(data_sum,data[i],npts);
    }

/* write timeshift and correlation coefficient to tdel.dat */
    if ((fp = fopen("tdel.dat","w+")) == NULL){
        printf("Can't open file.\n");
        exit(1);
    }
    for (i=0; i<sta_num; i++){
        fprintf(fp,"%s\t%f\t%f\n",indir[i], tdel[i],corrcef[i]);
    }

    free(data);
    free(indir);
    free(ccf);
    free(data_sum);
    free(tdel);
    free(corrcef);
    fftw_free(srcin);
    fftw_free(srcout);
    fftw_free(ffis);
    fftw_free(ffjs);
    fftw_free(fft_all);
}


float coef(float* a, float* b, int n){
    double c=0., cef;
    double norm_a=0., norm_b=0.;
    int i;
    for (i=0; i<n; i++){
        c += a[i]*b[i];
        norm_a += a[i]*a[i];
        norm_b += b[i]*b[i];
    }
    cef = c/(sqrt(norm_a)*sqrt(norm_b));
    return((float)cef);
}


int argmax(double *a, int n){
    int i, idx = 0;
    double amp = 0.;
    for(i=0;i<n;i++){
        if (a[i]>amp){
            amp = a[i];
            idx = i;
        }
    }
    return(idx);
}
