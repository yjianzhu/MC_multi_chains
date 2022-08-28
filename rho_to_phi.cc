#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double pi = 3.1415926535897932384626433832795028841971693;

int main ()
{
    FILE *fp;
    double Lx; // distance between the walls
    double x; // x coordinate perpendicular to the walls
    int bins; // number of bins in the x direction
    double rhofav, rhosav, rhotav; // distributions of monomer centers

    printf("Enter Lx: ");
    scanf("%lf",&Lx);
    printf("%lf \n",Lx);

    printf("Enter number of bins: ");
    scanf("%d",&bins);
    printf("%d \n",bins);

    double dx=Lx/bins; // width of the bins
    double phifav[bins], phisav[bins]; // volume fractions of flexible and stiff polymers
    int m=0.5/dx;
    double x1, x2;
    double dV[m+1]; // volume of a monomer between positions x1 and x2 
    for(int i=0;i<=m;i++) {
         x1 = (i-0.5)*dx;
         x2 = (i+0.5)*dx;
         if (x2>0.5) x2=0.5; 
         dV[i] = pi*x2*(0.25-x2*x2/3.0)-pi*x1*(0.25-x1*x1/3.0); 
    }

    fp = fopen("rho_ave","r");

    // read in the distributions of the monomer centers and calculate the volume fractions
    for(int i=0;i<bins;i++) phifav[i]=phisav[i]=0.0;
    for(int i=0;i<bins;i++){
         fscanf(fp,"%lf %lf %lf %lf",&x,&rhofav,&rhosav,&rhotav);
         for(int j=-m;j<=m;j++) if(i+j>=0 && i+j<bins) {
              phifav[i+j] += rhofav*dV[abs(j)];
              phisav[i+j] += rhosav*dV[abs(j)];
         }
    }
    fclose(fp);

    // write out the volume fractions
    fp = fopen("phi_ave","w");
    for(int i=0;i<bins;i++){
         fprintf(fp,"%lf %lf %lf %lf\n",(i+0.5)*dx,phifav[i],phisav[i],phifav[i]+phisav[i]);
    }
    fclose(fp);

    return 0;
}

