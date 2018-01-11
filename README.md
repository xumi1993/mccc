# mccc
Function MCCC determines optimum relative delay times for a set of seismograms based on the VanDecar & Crosson multi-channel cross-correlation algorithm.

## Installation
### Dependence
* SAC
* fftw3

### Modify Makefile
Specify the directories of SAC library in these lines:
```Makefile
SACINC = -I/usr/local/sac/include
SACLIB = -L/usr/local/sac/lib 
```

### Compile Codes 
```Bash
cd mccc
make
```

## Usage
```
mccc [-D] [-Ffilter/order/f1[/f2]] [-Ctmark/t1/t2] [-v] sac_files
    -D: Remean traces
    -C: window data [tmark+t1,tmark+t2]
        mark = -5(b), -3(o), -2(a), 0(t0)
    -F: Fiter traces
        filter = [bp|hp|lp]
        order = filter order
        if [hp|lp] specified f1 is the corner freq.
        if [bp] specified f1 f2 are the corner freq.
    -v: Print file name at corr-correlation
```
The program will creat a "tdel.dat" file with three columes:

    Filename    delay-time  correlation-coefficient
