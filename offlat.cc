#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <io.h>
#include <direct.h>
#include <string>
#include "random.h" //contains random number generator

// GLOBAL VARIABLES

enum { nf = 12}; // number of flexible chains   //临时修改为 长链
enum { ns = 12}; // number of stiff chains
enum { Nf = 60}; // number of beads in flexible chains  //临时修改为长链
enum { Ns = 20}; // number of beads in stiff chains
enum { nbead= ((long long int)ns) * ((long long int)Ns) + ((long long int)nf) * ((long long int)Nf)};
enum { npoly = nf+ns };

// position vectors for beads
double xob[nbead], yob[nbead], zob[nbead];

double Lx = 25, Ly = 8, Lz = 8; // dimensions of simulation box
double V = Lx*Ly*Lz; // volume of simulation box
bool wall=true; // if true, walls are placed at x=0 and x=Lx
double kappa = 1.0; // bending modulus of stiff polymers
double totalBE = 0.0; // total bending energy of stiff polymers

// numbers of MC steps
enum { nrelax_pb = 100000}; // relaxation steps per bead
enum { neval_pb  = 10000000}; // evaluation steps per bead
enum { nskip_pb = 10000}; // number of steps per bead between outputs
const long long int nrelax = ((long long int)nrelax_pb) * ((long long int)nbead);
const long long int neval = ((long long int)neval_pb) * ((long long int)nbead);
const long long int nskip = ((long long int)nskip_pb) * ((long long int)nbead);

// voxel array
const int Lxi=50, Lyi=16, Lzi=16; // ideally, Lxi=2*Lx, etc
double LtoLat[3]; // ratio of Lxi and Lx, etc
int beadhere[Lxi][Lyi][Lzi]; // table of occupied voxels
bool beadhere_of[Lxi][Lyi][Lzi]; // table of voxel overflows (used during equilibration)
bool overlaps=true, overflows=true; // flags for overlaps and overflows 
 
// other lookup tables
int head[npoly], tail[npoly], type[npoly], parent[nbead];

//other constants
const double pi = 3.1415926535897932384626433832795028841971693;

// variables for subroutine conc()
const int bins=1000;
int stats=0;
double rhos[bins],rhof[bins],rhosav[bins],rhofav[bins], BEav;

// variables for subroutine stats()
int long nfchains=0, nschains=0;
double r0f2=0.0,rgf2=0.0,r0s2=0.0,rgs2=0.0;

// variables for subroutine tistr()
char tms[50];
time_t time0, time1;

//==============================================================
// Checks if file exists
//--------------------------------------------------------------
bool fexist(char *filename){
    if (FILE * file = fopen(filename, "r"))
    {
        fclose(file);
        return true;
    }
    return false;
}
//==============================================================
// Checks if directory exists
//--------------------------------------------------------------
void dexist(std::string dirname){
    if(_access(dirname.c_str(),0)==-1)
        _mkdir(dirname.c_str());
    return ;
}
//==============================================================
// finds the time in a string
//--------------------------------------------------------------
void tistr()
{
    long int tisec = time(NULL)-time0;
    
    int secs,mins, hours, days;
    
    int lm=60;
    int lh=60*lm;
    int ld=24*lh;
    
    days = tisec/ld;
    hours=(tisec-(days*ld))/lh;
    mins=(tisec-(days*ld)-(hours*lh))/lm;
    secs=tisec-(days*ld)-(hours*lh)-(mins*lm);
    
    if(days>0){
        sprintf(tms,"%d, %02d:%02d:%02d",days,hours,mins,secs);
    } else if(hours>0){
        sprintf(tms,"%d:%02d:%02d",hours,mins,secs);
    } else if(mins>0){
        sprintf(tms,"%d:%02d",mins,secs);
    } else {
        sprintf(tms,"%ds",secs);
    }
    
}
//==============================================================
// determines voxel coordinates
//--------------------------------------------------------------
void voxel_coord(double x, double y, double z, int *vx, int *vy, int *vz) {

    x -= floor(x/Lx)*Lx;
    y -= floor(y/Ly)*Ly;
    z -= floor(z/Lz)*Lz;

    *vx = x*LtoLat[0];
    *vy = y*LtoLat[1];
    *vz = z*LtoLat[2];

    return;
}
//==============================================================
// construct voxel table
//--------------------------------------------------------------
void voxels(){
    int x,y,z;
    overflows=false;
    for(x=0;x<Lxi;x++) for(y=0;y<Lyi;y++) for(z=0;z<Lzi;z++){
        beadhere[x][y][z]=-1;
        beadhere_of[x][y][z]=false;
    }
    
    for(int i=0;i<nbead;i++){
        voxel_coord(xob[i],yob[i],zob[i],&x,&y,&z);
        if(beadhere[x][y][z]==-1){ //checks for overlap
            beadhere[x][y][z]=i;
        } else {
            beadhere_of[x][y][z]=true;
            overflows=true;
        }
    }
}
//==============================================================
// output all particle positions
//--------------------------------------------------------------
void out_positions (){
    FILE *fp;
    char fname[100];
    sprintf(fname,"positions");
    fp = fopen(fname,"w");
    fprintf(fp,"%d %d %d %d\n",nf,ns,Nf,Ns);
    for(int i=0;i<nbead;i++){
        fprintf(fp,"%.17lf %.17lf %.17lf\n",xob[i],yob[i],zob[i]); 
        // prints a large number of digits to avoid rounding inaccuracies
    }
    fclose(fp);
}
//==============================================================
// input all particle positions
//--------------------------------------------------------------
int in_positions (){
    FILE *fp;
    int nft,nst,Nft,Nst;
    char fname[100]; //stored fname to make it easier to add funtionality
    sprintf(fname,"positions");
    if(!fexist(fname)) {
        printf("\"%s\" does not exist.\n",fname);
        return 1;
    } else {
        printf("Reading particle positions from \"%s\"\n",fname);
    }
    fp = fopen(fname,"r");
    fscanf(fp,"%d %d %d %d\n",&nft,&nst,&Nft,&Nst);
    if(!(nf==nft && ns==nst && Nf==Nft && Ns==Nst)){
        printf("Parameters do not match:\n");
        printf("Current: %d %d %d %d\n",nf,ns,Nf,Ns);
        printf("Old:     %d %d %d %d\n",nft,nst,Nft,Nst);
        return 1;
    }
    for(int i=0;i<nbead;i++){
        fscanf(fp,"%lf %lf %lf",&xob[i],&yob[i],&zob[i]);
    }
    fclose(fp);
    return 0;
}
//==============================================================
// initializes variables
//--------------------------------------------------------------
void initialize(){
    int mono;

    LtoLat[0] = (1.0*Lxi)/Lx;
    LtoLat[1] = (1.0*Lyi)/Ly;
    LtoLat[2] = (1.0*Lzi)/Lz;
    for (int i=0; i<3; i++){
      if(LtoLat[i]>2.0) printf("WARNING: LtoLat[%1d] is too large\n",i);
      if(LtoLat[i]<sqrt(3.0)) printf("WARNING: LtoLat[%1d] is too small\n",i);
    }

    int lastflex=nf*Nf, N;
    for(int i=0;i<npoly;i++){
        if(i<nf){
            type[i] = 0; // flexible polymer
            head[i] = i*Nf; // first bead of polymer
            tail[i] = head[i]+Nf-1; // last bead of polymer
        } else {
            type[i] = 1; // stiff polymer
            head[i] = nf*Nf + (i-nf)*Ns; // first bead of polymer
            tail[i] = head[i]+Ns-1; // last bead of polymer
        }
    }
    mono=0;
    for(int i=0;i<npoly;i++){
        N=tail[i]-head[i]+1;
        for(int j=0;j<N;j++){
            parent[mono]=i;
            mono++;
        }
    }
    
    voxels();
    for(int i=0;i<bins;i++){
        rhosav[i]=0;
        rhofav[i]=0;
    }
}
//==============================================================
// set monomer positions 
//--------------------------------------------------------------
void set_positions(int read=1){
    int N,mono;
    double newx,newy,newz,theta,phi,dx,dy,dz;
    int xint,yint,zint, read_failed=0;
    
    if(read==1){ //read in configuration
        read_failed = in_positions();
    }
    if(read==0 || read_failed==1) { //generate random configuration
        mono=0;
        for(int i=0;i<npoly;i++){
            if(i<nf) N=Nf;
            else N=Ns;
            
            if(wall) xob[mono]=genrand_res53()*(Lx-1.0)+0.5;
            else xob[mono]=genrand_res53()*Lx; 
            yob[mono]=genrand_res53()*Ly;
            zob[mono]=genrand_res53()*Lz;
            
            mono++;
            for(int j=1;j<N;j++){
                do {
                    dx = 2.0*genrand_res53()-1.0;
                } while (wall && (fabs(xob[mono-1]+dx-0.5*Lx)>0.5*Lx-0.5));
                phi  =genrand_res53()*2.0*pi;
                dy = sqrt(1.0-dx*dx)*cos(phi);
                dz = sqrt(1.0-dx*dx)*sin(phi);
                xob[mono] = xob[mono-1]+dx;
                yob[mono] = yob[mono-1]+dy;
                zob[mono] = zob[mono-1]+dz;
                mono++;
            }
        }
    }
}
//==============================================================
// checks if two monomers overlap
//--------------------------------------------------------------
bool tooclose(double x1, double y1, double z1, double x2, double y2, double z2){

    double dx = fabs(x1-x2);
    double dy = fabs(y1-y2);
    double dz = fabs(z1-z2);

    if(!wall) dx = fmod(dx+0.5*Lx,Lx)-0.5*Lx;
    dy = fmod(dy+0.5*Ly,Ly)-0.5*Ly;
    dz = fmod(dz+0.5*Lz,Lz)-0.5*Lz;
    
    if(dx*dx+dy*dy+dz*dz<0.999999) return true;

    return false;
}
//==============================================================
// checks if there is overlap with another monomer
//--------------------------------------------------------------
bool overlap (double xnew, double ynew, double znew, int oldmono=-1){
    int x, y, z; 

    voxel_coord(xnew,ynew,znew,&x,&y,&z);
    int xc,yc,zc,mono; //x,y,z to check
    
    if(beadhere_of[x][y][z]) return true;
    
    if((beadhere[x][y][z] > -1) && (beadhere[x][y][z] != oldmono)) return true;
    
    for(int dx=-2;dx<3;dx++){
        xc=x+dx;
        if(xc<0) xc+=Lxi;
        if(xc>=Lxi) xc-=Lxi;
        for(int dy=-2;dy<3;dy++){
            yc=y+dy;
            if(yc<0) yc+=Lyi;
            if(yc>=Lyi) yc-=Lyi;
            for(int dz=-2;dz<3;dz++){
                zc=z+dz;
                if(zc<0) zc+=Lzi;
                if(zc>=Lzi) zc-=Lzi;

                mono=beadhere[xc][yc][zc];
                if(mono>-1){ //check if bead is present in voxel
                    if(tooclose(xob[mono],yob[mono],zob[mono],xnew,ynew,znew) && mono!=oldmono){
                        return true; //there is a bead in the way....
                    }
                }
            }
        }
    }
    return false;
}
//==============================================================
// counts the number of overlaps in the system
//--------------------------------------------------------------
int olaps (){
    int oops=0;
    int x,y,z,xc,yc,zc,mono;
    double xx,yy,zz;
    for(int i=0;i<nbead;i++){
        xx=xob[i], yy=yob[i], zz=zob[i];
        voxel_coord(xx,yy,zz,&x,&y,&z);
        for(int dx=-2;dx<3;dx++){
            xc=x+dx;
            if(xc<0) xc+=Lxi;
            if(xc>=Lxi) xc=xc%Lxi;
            for(int dy=-2;dy<3;dy++){
                yc=y+dy;
                if(yc<0) yc+=Lyi;
                if(yc>=Lyi) yc=yc%Lyi;
                for(int dz=-2;dz<3;dz++){
                    zc=z+dz;  
                    if(zc<0) zc+=Lzi;
                    if(zc>=Lzi) zc=zc%Lzi;

                    mono=beadhere[xc][yc][zc];
                    if(mono>i){ //check if bead is present in voxel
                        if(tooclose(xx,yy,zz,xob[mono],yob[mono],zob[mono])){
                            oops++;
                        }
                    }
                }
            }
        }
    }
    
    return oops;
}
//==============================================================
// bending energy calculated from the bond unit vectors
//--------------------------------------------------------------
double BE(double ux1, double uy1, double uz1, double ux2, double uy2, double uz2){
    double dot = ux1*ux2 + uy1*uy2 + uz1*uz2;
    
    //return -kappa*dot; // bending energy
    //return kappa*(1-dot)*(1-dot); // aternative bending energy
    return kappa*acos(dot)*acos(dot);   //zhu
}
//==============================================================
// bending energy of entire polymer
//--------------------------------------------------------------
double BEchain(int poly){

    double ux1, uy1, uz1, ux2, uy2, uz2; // bond vectors
    double BEsum=0.0; 
    int mono = head[poly];
    int N = tail[poly]-mono+1;

    ux2=xob[mono]-xob[mono+1];
    uy2=yob[mono]-yob[mono+1]; 
    uz2=zob[mono]-zob[mono+1];
 
    for(int i=1;i<N-1;i++){
      ux1=ux2, uy1=uy2, uz1=uz2;
      ux2=xob[mono+i]-xob[mono+i+1];
      uy2=yob[mono+i]-yob[mono+i+1]; 
      uz2=zob[mono+i]-zob[mono+i+1];
      BEsum += BE(ux1,uy1,uz1,ux2,uy2,uz2);
    }
    
    return BEsum;
}
//==============================================================
// execute the snake move by shifting all the monomers
//--------------------------------------------------------------
void rollsnake (int poly, bool newhead, double x, double y, double z){
    int mono,xi,yi,zi,N=tail[poly]-head[poly]+1;

    if(newhead){
        //remove tail from lattice
        mono=tail[poly];
        voxel_coord(xob[mono],yob[mono],zob[mono],&xi,&yi,&zi);
        beadhere[xi][yi][zi]=-1;
        
        for(int i=N-1;i>0;i--){
            mono=head[poly]+i;
            
            //update positions
            xob[mono]=xob[mono-1];
            yob[mono]=yob[mono-1];
            zob[mono]=zob[mono-1];
            
            //update lattice
            voxel_coord(xob[mono],yob[mono],zob[mono],&xi,&yi,&zi);
            beadhere[xi][yi][zi]=mono;
        }
        mono=head[poly];
        xob[mono]=x;
        yob[mono]=y;
        zob[mono]=z;
        
        //place new head in voxel
        voxel_coord(x,y,z,&xi,&yi,&zi);
        beadhere[xi][yi][zi]=mono;
        
    } else { //new tail (tail is large i end)
        //remove head from lattice
        mono=head[poly];
        voxel_coord(xob[mono],yob[mono],zob[mono],&xi,&yi,&zi);
        beadhere[xi][yi][zi]=-1;
        for(int i=0;i<N-1;i++){
            mono=head[poly]+i;
            
            //update positions
            xob[mono]=xob[mono+1];
            yob[mono]=yob[mono+1];
            zob[mono]=zob[mono+1];
            
            //update lattice
            voxel_coord(xob[mono],yob[mono],zob[mono],&xi,&yi,&zi);
            beadhere[xi][yi][zi]=mono;
        }
        mono=tail[poly];
        xob[mono]=x;
        yob[mono]=y;
        zob[mono]=z;
        
        //place new head in voxel
        voxel_coord(x,y,z,&xi,&yi,&zi);
        beadhere[xi][yi][zi]=mono;
    }
}
//==============================================================
// snakemove move              
//--------------------------------------------------------------
int snake (){
    int poly = genrand_res53()*(nf+ns); //randomly choose polymer
    int newhead=true; // move in direction of head
    if(genrand_res53()>0.5) newhead=false; // move in direction of tail
    double dBE, newx, newy, newz, phi, ux1, uy1, uz1, ux2, uy2, uz2;
    int m1, m2, m3; // monomers of old bond
    int m4, m5; // two monomers of new bond

    if(newhead){
        m1=tail[poly]; // monomer to be deleted
        m2=m1-1;
        m3=m2-1;
        m4=head[poly]; // monomer where new bead is attached
        m5=m4+1;
    } else {
        m1=head[poly]; // monomer to be deleted 
        m2=m1+1;
        m3=m2+1;
        m4=tail[poly]; // monomer where new bead is attached
        m5=m4-1;
    }

    ux1 = xob[m4]-xob[m5]; 
    uy1 = yob[m4]-yob[m5]; 
    uz1 = zob[m4]-zob[m5]; 
   
    do { // choose a random displacement that does not cause backfolding 
        phi = genrand_res53()*2.0*pi;
        ux2 = 2.0*genrand_res53()-1.0;
        uy2 = sqrt(1.0-ux2*ux2)*cos(phi);
        uz2 = sqrt(1.0-ux2*ux2)*sin(phi);
    } while (ux1*ux2+uy1*uy2+uz1*uz2<-0.5);

    newx = xob[m4]+ux2;
    newy = yob[m4]+uy2;
    newz = zob[m4]+uz2;

    if(wall && (newx<0.5 || newx>(Lx-0.5))) return 0;
    if(overlap(newx,newy,newz,m1)) return 0; //reject if there is an overlap

    // get voxel coordinates of the monomer to be removed   
    int xi, yi, zi;
    voxel_coord(xob[m1],yob[m1],zob[m1],&xi,&yi,&zi);

    //if(type[poly]==1 && !beadhere_of[xi][yi][zi]){ //only do if stiff and no overflow   
    if(!beadhere_of[xi][yi][zi]){   //zhu 修改为flexible 也算能量
        // calculate bending energy of new bond
        dBE = BE(ux1,uy1,uz1,ux2,uy2,uz2);

        // subtract bending energy of deleted bond
        ux1=xob[m1]-xob[m2], uy1=yob[m1]-yob[m2], uz1=zob[m1]-zob[m2];
        ux2=xob[m2]-xob[m3], uy2=yob[m2]-yob[m3], uz2=zob[m2]-zob[m3];
        dBE -= BE(ux1,uy1,uz1,ux2,uy2,uz2);

        double P = exp(-dBE);
        if(genrand_res53()>P) return 0; //reject
        totalBE += dBE; // update total bond energy
    }

    rollsnake(poly,newhead,newx,newy,newz); //Accept! i.e. move polymer via snake move

    // check if there was an overflow at position where monomer was removed
    if (beadhere_of[xi][yi][zi]) voxels();
    
    return 1;
}
//==============================================================
// crankshaft move  
//--------------------------------------------------------------
int crank (){
    int poly = genrand_res53()*(nf+ns); //randomly choose polymer
    int N = tail[poly]-head[poly]+1;
    int mono0 = floor(genrand_res53()*(N-2.0)); //choose monomer along chain, neglect head and tail
    int mono = mono0 + head[poly] + 1;
        
    int m1=mono+1,m2=mono-1;
    double x1=xob[m1];
    double y1=yob[m1];
    double z1=zob[m1];
    
    double x2=xob[m2];
    double y2=yob[m2];
    double z2=zob[m2];
    
    double x=xob[mono];
    double y=yob[mono];
    double z=zob[mono];
    
    // generate vector from m1 to m2
    double dx12 = x2-x1, dy12 = y2-y1, dz12 = z2-z1;
    double R = sqrt(1.0-0.25*(dx12*dx12+dy12*dy12+dz12*dz12));

    // generate vector to midpoint between m1 and m2
    double xmid=0.5*(x1+x2), ymid=0.5*(y1+y2), zmid=0.5*(z1+z2) ;

    // generate vector from midpoint to mono
    double xp1=x-xmid, yp1=y-ymid, zp1=z-zmid;

    // this vector should have a length of R, but best to make sure it is accurate
    double scale = R/sqrt(xp1*xp1+yp1*yp1+zp1*zp1);
    xp1 *= scale, yp1 *= scale, zp1 *= scale;

    //generate a perpendicular vector using cross product
    double xp2 = yp1*dz12 - zp1*dy12;
    double yp2 = zp1*dx12 - xp1*dz12;
    double zp2 = xp1*dy12 - yp1*dx12;

    // scale new perpendicular vector to have length R
    scale = R/sqrt(xp2*xp2+yp2*yp2+zp2*zp2);
    xp2 *= scale, yp2 *= scale, zp2 *= scale;

    // rotate mono by a random angle phi
    double phi = genrand_res53()*2.0*pi, cosine=cos(phi), sine=sin(phi);
    double dx = cosine*xp1 + sine*xp2;
    double dy = cosine*yp1 + sine*yp2;
    double dz = cosine*zp1 + sine*zp2;
    
    //new monomer position
    double newx = xmid + dx;
    double newy = ymid + dy;
    double newz = zmid + dz;
    
    //wall
    if(wall && (newx<0.5 || newx>(Lx-0.5))) return 0;
    
    if(overlap(newx, newy, newz, mono)) return 0;
    
    // get voxel coordinates of the monomer to be moved
    int oxi, oyi, ozi;
    voxel_coord(x,y,z,&oxi,&oyi,&ozi);

    //if(type[poly]==1 && !beadhere_of[oxi][oyi][ozi]){ //only do if stiff and no overflow
    if(!beadhere_of[oxi][oyi][ozi]){ //zhu 都计算
        double ux1, uy1, uz1, ux2, uy2, uz2; // bond vectors 
        double dBE=0; // change in total bond energy
        //old middle is x,y,z
        //new middle is newx, newy, newz
        //m1 is higher, m2 is lower

        if(mono0<N-3){ //only if a bond exists at higher position
            ux1=x-x1, uy1=y-y1, uz1=z-z1;
            ux2=x1-xob[m1+1], uy2=y1-yob[m1+1], uz2=z1-zob[m1+1];
            dBE -= BE(ux1,uy1,uz1,ux2,uy2,uz2); // for old bond at mono+1

            ux1=newx-x1, uy1=newy-y1, uz1=newz-z1;
            dBE += BE(ux1,uy1,uz1,ux2,uy2,uz2); // for new bond at mono+1
        }
 
        if(mono0>0){ //only if a bond exists at lower position
            ux1=x-x2, uy1=y-y2, uz1=z-z2;
            ux2=x2-xob[m2-1], uy2=y2-yob[m2-1], uz2=z2-zob[m2-1];
            dBE -= BE(ux1,uy1,uz1,ux2,uy2,uz2); // for old bond at mono-1

            ux1=newx-x2, uy1=newy-y2, uz1=newz-z2;
            dBE += BE(ux1,uy1,uz1,ux2,uy2,uz2); // for new bond at mono-1
        }

        double P = exp(-dBE);
        if(genrand_res53()>P) return 0; //reject
        totalBE += dBE; // update total bond energy
    }
        
    // accept move and update positions
    xob[mono]=newx;
    yob[mono]=newy;
    zob[mono]=newz;
    
    // update lookup tables
    int nxi, nyi, nzi;
    voxel_coord(newx,newy,newz,&nxi,&nyi,&nzi); // voxel of new position
    beadhere[oxi][oyi][ozi]=-1; // remove old lattice position
    beadhere[nxi][nyi][nzi]=mono; // add new lattice position
    if(beadhere_of[oxi][oyi][ozi]) voxels();
    
    return 1;
}
//==============================================================
// rotation move for end monomers              
//--------------------------------------------------------------
int rotation (){
    int poly = genrand_res53()*(nf+ns); //randomly choose polymer
    int movehead=true; // move head monomer
    if(genrand_res53()>0.5) movehead=false; // move tail monomer
    double dBE, newx, newy, newz, phi, ux1, uy1, uz1, ux2, uy2, uz2;
    int m1, m2, m3; // monomers of the altered bond

    if(movehead){
        m1=head[poly]; // monomer to be moved
        m2=m1+1;
        m3=m2+1;
    } else {
        m1=tail[poly]; // monomer to be moved
        m2=m1-1;
        m3=m2-1;
    }

    // direction of fixed bond
    ux1 = xob[m2]-xob[m3]; 
    uy1 = yob[m2]-yob[m3]; 
    uz1 = zob[m2]-zob[m3]; 
   
    do { // choose a random displacement that does not cause backfolding 
        phi = genrand_res53()*2.0*pi;
        ux2 = 2.0*genrand_res53()-1.0;
        uy2 = sqrt(1.0-ux2*ux2)*cos(phi);
        uz2 = sqrt(1.0-ux2*ux2)*sin(phi);
    } while (ux1*ux2+uy1*uy2+uz1*uz2<-0.5);

    newx = xob[m2]+ux2;
    newy = yob[m2]+uy2;
    newz = zob[m2]+uz2;

    if(wall && (newx<0.5 || newx>(Lx-0.5))) return 0;
    if(overlap(newx,newy,newz,m1)) return 0; //reject if there is an overlap

    // get voxel coordinates of the monomer to be moved   
    int oxi, oyi, ozi;
    voxel_coord(xob[m1],yob[m1],zob[m1],&oxi,&oyi,&ozi);

    //if(type[poly]==1 && !beadhere_of[oxi][oyi][ozi]){ //only do if stiff and no overflow
    if(!beadhere_of[oxi][oyi][ozi]){ //zhu 
        // calculate bending energy of new bond
        dBE = BE(ux1,uy1,uz1,ux2,uy2,uz2);

        // subtract bending energy of deleted bond
        ux2=xob[m1]-xob[m2], uy2=yob[m1]-yob[m2], uz2=zob[m1]-zob[m2];
        dBE -= BE(ux1,uy1,uz1,ux2,uy2,uz2);

        double P = exp(-dBE);
        if(genrand_res53()>P) return 0; //reject
        totalBE += dBE; // update total bond energy
    }
    
    // if gotten this far, accept!
    // update positions
    xob[m1]=newx;
    yob[m1]=newy;
    zob[m1]=newz;
    
    //update lattice
    int nxi, nyi, nzi;
    voxel_coord(newx,newy,newz,&nxi,&nyi,&nzi); // voxel of new position
    beadhere[oxi][oyi][ozi]=-1; // remove old lattice position
    beadhere[nxi][nyi][nzi]=m1; // add new lattice position
    if(beadhere_of[oxi][oyi][ozi]) voxels();
    
    return 1;
}
//==============================================================
// swap move              
//--------------------------------------------------------------
int swap_chain (){
    if((int)Ns!=(int)Nf) return 0;
    int polyf = genrand_res53()*(nf); //choose a flexible chain
    int polys = genrand_res53()*(ns)+nf; //choose a stiff chain
    
    double dBE = BEchain(polyf)-BEchain(polys);
    
    double P = exp(-dBE);
    if(genrand_res53()>P) return 0; //reject
    totalBE += dBE; // update total bond energy
    
    double tmp;
    int monos, monof, xis, yis, zis, xif, yif, zif;
    for(int i=0;i<Ns;i++){ //if got this far, accept! swap chains
        monos = head[polys]+i; 
        monof = head[polyf]+i;
 
        //swap monomers in voxels
        voxel_coord(xob[monos],yob[monos],zob[monos],&xis,&yis,&zis);
        voxel_coord(xob[monof],yob[monof],zob[monof],&xif,&yif,&zif);
        beadhere[xis][yis][zis] = monof;
        beadhere[xif][yif][zif] = monos;
        
        //swap positions
        tmp = xob[monos];
        xob[monos] = xob[monof];
        xob[monof] = tmp;
        
        tmp = yob[monos];
        yob[monos] = yob[monof];
        yob[monof] = tmp;
        
        tmp = zob[monos];
        zob[monos] = zob[monof];
        zob[monof] = tmp;
    }
    
    return 1;
}
//==============================================================
// output to dump file, setup for ovito
//--------------------------------------------------------------
void dump (int flag, int tt, double append=0){
    int mol,typee, sitee;
    FILE *fp;
    char fname[30];
    int ptwall=1;
    int natoms=nbead;
    double x,y,z;
    //if(ptwall==0) natoms=nbead;
    sprintf(fname,"dumps/dump%.2lf-%d",append,tt);
    fp = fopen(fname,"w");
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",tt);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,"%d\n",nbead);
    fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(fp,"0 %lf\n",Lx);
    fprintf(fp,"0 %lf\n",Ly);
    fprintf(fp,"0 %lf\n",Lz);
    fprintf(fp,"ITEM: ATOMS id mol type x y z\n");
    for(int i=0;i<nbead;i++){ //nbead or vac
        mol=parent[i];
        typee = type[mol];
        x = xob[i]-floor(xob[i]/Lx)*Lx;
        y = yob[i]-floor(yob[i]/Ly)*Ly;
        z = zob[i]-floor(zob[i]/Lz)*Lz;
        fprintf(fp,"%d %d %d %lf %lf %lf\n",i,mol,typee,x,y,z);
    }
    fclose(fp);

    if (flag) {
    // output concentration profiles    
    sprintf(fname,"rho_ave");
    fp = fopen(fname,"w");
    for(int i=0;i<bins;i++){
        x = (i+0.5)*Lx/bins;
        fprintf(fp,"%lf %lf %lf %lf\n",x,rhofav[i],rhosav[i],rhofav[i]+rhosav[i]);
    }
    fclose(fp);

    // output average chain dimensions
    if(nfchains*nschains>0){
        sprintf(fname,"poly_size");
        fp = fopen(fname,"a");
        fprintf(fp,"%ld %lf %lf  ",nfchains,r0f2,rgf2);
        fprintf(fp,"%ld %lf %lf\n",nschains,r0s2,rgs2);
        fclose(fp);
    }

    // output average bending energy
    sprintf(fname,"bend_energy");
    fp = fopen(fname,"a");
    fprintf(fp,"%lf\n",BEav);
    }

    return;
}
//==============================================================
// seed random number generator
//--------------------------------------------------------------
void seed_number_generator()
{
	unsigned long int iseed;
	FILE *file_ptr;

	iseed = time(NULL);
	iseed += 123456789;

	init_genrand(iseed);

	file_ptr = fopen("random_seed.txt", "w");
	fprintf(file_ptr, "%15lu\n", iseed);
	fclose(file_ptr);
}
//==============================================================
// collect statistics for concentration profiles
//--------------------------------------------------------------
void conc (){
    int N,mono,xi;
    double x,y,z, ADD=bins/V;
    for(int i=0;i<bins;i++){
        rhos[i]=0;
        rhof[i]=0;
    }
    
    for(int poly=0;poly<npoly;poly++){
        N = tail[poly]-head[poly]+1;
        for(int j=0;j<N;j++){
            mono=head[poly]+j;
            x = xob[mono]-floor(xob[mono]/Lx)*Lx;
            xi=bins*x/Lx;
            if(type[poly]==0) rhof[xi]+=ADD;
            else rhos[xi]+=ADD;
        }
    }
    
    stats ++;
    
    for(int i=0;i<bins;i++){
        rhofav[i] += (rhof[i]-rhofav[i])/(1.0*stats);
        rhosav[i] += (rhos[i]-rhosav[i])/(1.0*stats);
    }

    // collect statistics for average bending energy
    BEav += (totalBE/((Ns-2.0)*ns)-BEav)/(1.0*stats);
}
//==============================================================
// collect statistics for R0 and Rg
//--------------------------------------------------------------
void poly_size (){
    double dx,dy,dz,r02,xcm,ycm,zcm,rg2;
    int m1,mN,middle,N;

    // only select polymers with xmin < x< xmax
    double xmin = 0.4*Lx, xmax = 0.6*Lx;

    for(int poly=0;poly<npoly;poly++){
        m1 = head[poly];
        mN = tail[poly];
        N = mN-m1+1;
        middle = (m1+mN)/2;
        if (xob[middle]>xmin && xob[middle]<xmax) {
            dx = xob[mN]-xob[m1];
            dy = yob[mN]-yob[m1];
            dz = zob[mN]-zob[m1];
            r02 = dx*dx+dy*dy+dz*dz;

            xcm = ycm = zcm = 0.0;
            for(int m=m1;m<=mN;m++){
                xcm += xob[m];
                ycm += yob[m];
                zcm += zob[m];
            }
            xcm /= N, ycm /= N, zcm /= N;
            rg2 = 0.0;
            for(int m=m1;m<=mN;m++){
                dx = xob[m]-xcm;
                dy = yob[m]-ycm;
                dz = zob[m]-zcm;
                rg2 += dx*dx+dy*dy+dz*dz;
            }
            rg2 /= N;

            if(type[poly]==0) {
                nfchains++;
                r0f2 +=(r02-r0f2)/(1.0*nfchains);
                rgf2 +=(rg2-rgf2)/(1.0*nfchains);
            } else {
                nschains++;
                r0s2 +=(r02-r0s2)/(1.0*nschains);
                rgs2 +=(rg2-rgs2)/(1.0*nschains);
            }
        }
    }
}
//==============================================================
// testing subroutine 
//--------------------------------------------------------------
void test (int flag){
    int i,mono,poly,N,x,y,z;
    double dx,dy,dz,dr2,sum=0.0;

    if(wall){
        for(i=0;i<nbead;i++) if(xob[i]<0.5||xob[i]>Lx-0.5)
            printf("WARNING: bead %d is out of range\n",i);
    }

    for(poly=0;poly<nf+ns;poly++){
        mono=head[poly];
        N=tail[poly]-head[poly]-1;
        for(i=0;i<N;i++) {
            dx = xob[mono+i+1]-xob[mono+i];
            dy = yob[mono+i+1]-yob[mono+i];
            dz = zob[mono+i+1]-zob[mono+i];
            dr2 = dx*dx+dy*dy+dz*dz;
            if (fabs(dr2-1.0)>0.000001) 
                 printf("WARNING: beads %2d and %2d of polymer %d are separated by %lf\n",i,i+1,poly,sqrt(dr2));
        }
    }

    if (flag){
        //for(poly=nf;poly<ns+nf;poly++) sum += BEchain(poly);
        for(poly=0;poly<ns+nf;poly++) sum += BEchain(poly); //zhu
        if(fabs(sum-totalBE)>0.00001) printf("WARNING: total bending energy is wrong\n");
    }

    int copy_beadhere[Lxi][Lyi][Lzi];
    bool copy_beadhere_of[Lxi][Lyi][Lzi];
    for(x=0;x<Lxi;x++) for(y=0;y<Lyi;y++) for(z=0;z<Lzi;z++){
        copy_beadhere[x][y][z]=beadhere[x][y][z];
        copy_beadhere_of[x][y][z]=beadhere_of[x][y][z];
    }
    voxels();
    bool flag1=false, flag2=false;
    for(x=0;x<Lxi;x++) for(y=0;y<Lyi;y++) for(z=0;z<Lzi;z++){
        if(copy_beadhere[x][y][z]!=beadhere[x][y][z]) flag1=true;
        if(copy_beadhere_of[x][y][z]!=beadhere_of[x][y][z]) flag2=true;
    }
    if(flag1) printf("WARNING: error in beadhere table\n");
    if(flag2) printf("WARNING: error in beadhere_of table\n");
}
//==============================================================
int main ()
{
    printf("Parameters:\n");
    printf("n: %d %d\n",ns,nf);
    printf("N: %d %d\n",Ns,Nf);
    printf("L: %lf %lf %lf\n",Lx,Ly,Lz);
    printf("wall: %d\n",wall);
    printf("kappa: %lf\n\n",kappa);
    
    printf("nrelax_pb: %u\n",nrelax_pb);
    printf("neval_pb: %u\n",neval_pb);
    printf("nskip_pb: %u\n",nskip_pb);
    
    seed_number_generator();
    printf("Seed setup.\n");
    set_positions(0);
    printf("Positions initialized.\n");
    initialize();
    printf("Initalized locations.\n");

    printf("Check if dumps folder exists?\n");
    dexist("dumps");

    dump(0,-1);
    printf("Output for first time\n");
    test(0); // set flag=0 to avoid testing total bending energy

    double dx,dy,dz;
    int mono;
    double cracc=0, snacc=0, rotcc=0, swacc=0;
    
    time0 = time(NULL);
    
    int nover=0, oops=1;
    
    printf("\nStart equilibration\n");
    cracc=0, snacc=0, rotcc=0, swacc=0;
    for(long int i=0;i<=nrelax||oops>1;i++){
        snacc+=snake();
        cracc+=crank();
        rotcc+=rotation();
        //swacc+=swap_chain();  //zhu
                        
        if(i%nskip==0) {
            dump(0,int(i/nbead));
            nover=0;
            for(int ii=0;ii<Lxi;ii++) for(int jj=0;jj<Lyi;jj++) for(int kk=0;kk<Lzi;kk++) 
                if(beadhere_of[ii][jj][kk]) nover++;
            oops=olaps();
            out_positions();
            tistr();
            printf("i=%ld/%u, overflows=%d, overlaps=%d c:%lg s:%lg r:%lg sw:%lg (Time: %s)\n",
                    i/nbead,nrelax_pb,nover,oops,cracc/(1.0+i),snacc/(1.0+i),rotcc/(1.0+i),swacc/(1.0+i),tms);
            fflush(stdout); //flush print buffer
        }
    }

    // calculate total bending energy once there is no overflow 
    totalBE=0.0;    
    for(int poly=0;poly<ns+nf;poly++) 
        totalBE += BEchain(poly);
    
    printf("\nCollect statistics\n");
    cracc=0, snacc=0, rotcc=0, swacc=0;
    for(long int i=1;i<neval;i++){
        snacc+=snake();
        cracc+=crank();
        rotcc+=rotation();
        //swacc+=swap_chain();  //zhu
        
        if(i%5000==0) conc();
        if(i%5000==0) poly_size();

        if(i%nskip==0) {
            dump(1,int((i+nrelax)/nbead));
            out_positions();
            tistr();
            printf("j=%ld/%u, c:%lg s:%lg r:%lg sw:%lg (Time: %s)\n",
                 i/nbead,neval_pb,cracc/(1.0+i),snacc/(1.0+i),rotcc/(1.0+i),swacc/(1.0+i),tms);
            fflush(stdout);
            test(1);
        }
    }
    
    printf("\nDONE!\n");
    return 0;
}

