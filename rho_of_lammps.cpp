#include<iostream>
#include<algorithm>
#include<fstream>
#include<string>
#include<vector>
#include<array>
#include<cmath>
// #include <io.h>
// #include <direct.h>

int N;
double Lx,Ly,Lz;
std::vector<int> number_type;
const int bins=1000;
double rhos[bins],rhof[bins],rhosav[bins]{0},rhofav[bins]{0}; // short and long 
int stats=0;
int save_step=1000;
int save_time=1;

//==============================================================
// Checks if directory exists//linux and windows not same
//--------------------------------------------------------------
// void dexist(std::string dirname){
//     if(_access(dirname.c_str(),0)==-1)
//         _mkdir(dirname.c_str());
//     return ;
// }

void conc(std::vector<std::array<double,3>> &x)
{
    int zi,n_type=1;
    double V=Lx*Ly*Lz;
    double z, ADD=bins/V;
    for(int i=0;i<bins;i++)
    {
        rhos[i]=0;
        rhof[i]=0;
    }

    stats ++;

    for(int i=0;i<x.size();i++)
    {
        z=x[i][2]-std::floor(x[i][2]/Lz)*Lz;
        zi=bins*z/Lz;
        if(i>=number_type[n_type])
            rhof[zi]+=ADD;
        else
            rhos[zi]+=ADD;
    }
    for(int i=0;i<bins;i++){
        rhofav[i] += (rhof[i]-rhofav[i])/(1.0*stats);
        rhosav[i] += (rhos[i]-rhosav[i])/(1.0*stats);
    }

    if(stats==save_step)
    {
        //dexist("dumps");
        std::string fname;
        fname="rho_ave"+std::to_string(stats*save_time);
        save_time++;
        std::fstream write;
        write.open(fname,std::ios::out);
        for(int i=0;i<bins;i++){
            double z = (i+0.5)*Lz/bins;
            write<<z<<'\t'<<rhofav[i]<<'\t'<<rhosav[i]<<'\t'<<rhofav[i]+rhosav[i]<<std::endl;
        }
        write.close();

        for(int i=0;i<bins;i++)
        {
            rhosav[i]=0;
            rhofav[i]=0;
        }
        stats=0;
    }
}

int  get_para(std::fstream &read)
{
    double low_x=0,low_y=0,low_z=0,high_x=0,high_y=0,high_z=0;
    std::string aline;
    if(!std::getline(read,aline))
        return 0;

    read>>N;
    for(int i=0;i<8;i++)
        if(!std::getline(read,aline))
            return 0;

    read>>low_x>>high_x;
    if(!std::getline(read,aline))
        return 0;
    read>>low_y>>high_y;
    if(!std::getline(read,aline))
        return 0;
    read>>low_z>>high_z;
    Lx=high_x-low_x;
    Ly=high_y-low_y;
    Lz=high_z-low_z;

    return 1;
}
//设计成这样不需要new和delete
int read_data(std::vector<std::array<double,3>> &x,std::fstream &read)
{
    int NB,atom_type,old_atom_type=0;
    
    std::string  aline;

    read>>NB;
    if(NB!=N)
    {
        std::cout<<"atom number isn't same!"<<std::endl;
        return 0;
    }

    if(!std::getline(read,aline))   return 0;
    if(!std::getline(read,aline))   return 0;
    
    for(int i=0;i<NB;i++)
    {
        x.push_back(std::array<double,3>{0});
        read>>atom_type>>x[i][0]>>x[i][1]>>x[i][2];
        if(atom_type!=old_atom_type)
        {
            number_type.push_back(i);
            old_atom_type=atom_type;
        }
    }

    return 1;
}

void cal_concertration(std::string readFileName,std::string dataName)
{
    std::fstream read;
    read.open(readFileName,std::ios::in);
    if(!get_para(read))
    {
        std::cout<<"error of reading"<<std::endl;
        exit(-1);
    }

    std::vector<std::array<double,3>> x;
    read.close();

    read.open(dataName,std::ios::in);
    while (read_data(x,read))
    {
        conc(x);
        x.clear();
        number_type.clear();
    }
    
    read.close();

}

int main(int argc, char *argv[])
{
    std::string readfile=argv[argc-1];
    cal_concertration(readfile,"dump.lammpstrj");
}