void	detrend(float *y, int n) {
     int i;
     double a, b, a11, a12, a22, y1, y2;
     y1 = y2 = 0.;
     for(i=0;i<n;i++) {
       y1 += i*y[i];
       y2 += y[i];
     }
     a12 = 0.5*n*(n-1);
     a11 = a12*(2*n-1)/3.;
     a22 = n;
     b = a11*a22-a12*a12;
     a = (a22*y1-a12*y2)/b; 
     b = (a11*y2-a12*y1)/b;
     for(i=0;i<n;i++) {
       y[i] = y[i] - a*i - b;
     }
}
