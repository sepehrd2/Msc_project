#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <ctime>

#define N 200                                   //Numbers of bead in membrane
#define M 4                                     //Numbers of bead in rod
#define N_r 10                                 //Number of rod
#define PI 3.14159265358979323846             //Pi number

//Function for making the force radial to the membrane
void Radial_Force_single(double F[2],double R[2],double FR[2]){
    double cosa,sina;
    double FFR;
    cosa=R[0]/sqrt(R[0]*R[0]+R[1]*R[1]);                //cos of attached bead of membrane
    sina=R[1]/sqrt(R[0]*R[0]+R[1]*R[1]);                //sin of attached bead of membrane
    FFR=F[0]*cosa+F[1]*sina;                            //Radial component of parallel force
    FR[0]=FFR*cosa;                                     //x component
    FR[1]=FFR*sina;                                     //y component
}
//Function for making the force tangential to a membrane
double Tangential_Force_single(double F[2],double R[2]){
    double cosa,sina;
    double FFT;
    cosa=R[0]/sqrt(R[0]*R[0]+R[1]*R[1]);          //cos of attached bead of membrane
    sina=R[1]/sqrt(R[0]*R[0]+R[1]*R[1]);          //ain of attached bead of membrane
    FFT=F[1]*cosa-F[0]*sina;                      //tangential component of parallel force
    return FFT;
}
//Function for upwarding the attached bead
int changing_bead(int bead){
    if (bead>(N-1))
    {
        bead=bead-N;                //avoiding the number of bead to be bigger than N-1
    }
    if (bead<0)
    {
        bead=bead+N;                //avoiding the number of bead to be smaller than zero
    }
    return bead;
}
//Function for finding the distance between two points
double dis(double r1[2],double r2[2]){
    double d;
    d=sqrt((r1[0]-r2[0])*(r1[0]-r2[0])+(r1[1]-r2[1])*(r1[1]-r2[1]));
    return (d);
}
//Muller Box function
double muller(){
    double z0,v1,v2,R;
    do
    {
        v1=2.0*rand()/double(RAND_MAX+1.0)-1.0;
        v2=2.0*rand()/double(RAND_MAX+1.0)-1.0;
        R=v1*v1+v2*v2;
    }
    while (R>1.0 || R==0.0);
    z0=v1*sqrt(-2.0*log(R)/R);
    return z0;
}
//Function for changing the attached bead
int changing_attachement(double A[N][2],int bead_C,double V_r_T,double t){
    double d_r,d=0,remind_d,r1[2],r2[2];
    int bead_number=0;
    d_r=V_r_T*t;                                                  //distance for translation
    if (V_r_T>0)
    {
        do{
        r1[0]=A[changing_bead(bead_C+bead_number)][0];r1[1]=A[changing_bead(bead_C+bead_number)][1];
        r2[0]=A[changing_bead(bead_C+bead_number+1)][0];r2[1]=A[changing_bead(bead_C+bead_number+1)][1];
        d=d+dis(r1,r2);
        bead_number=bead_number+1;}
        while (d<abs(d_r));                                      //checking how many beads should be passed
    }
    if (V_r_T<0)
    {
        do{
        r1[0]=A[changing_bead(bead_C-bead_number)][0];r1[1]=A[changing_bead(bead_C-bead_number)][1];
        r2[0]=A[changing_bead(bead_C-bead_number-1)][0];r2[1]=A[changing_bead(bead_C-bead_number-1)][1];
        d=d+dis(r1,r2);
        bead_number=bead_number-1;}
        while (d<abs(d_r));
    }
    if (V_r_T==0)
    {
        bead_number=1;
        d=0;
    }
    remind_d=(d-abs(d_r))/dis(r1,r2);                               //distance between the last bead of rod and next bead of membrane
    if (remind_d<muller())
    {
        bead_C=bead_C+bead_number;                                  //next bead
    }
    else
    {
        bead_C=bead_C+bead_number-bead_number/abs(bead_number);     //previous bead
    }
    if (bead_C>(N-1))
    {
        bead_C=bead_C-N;
    }
    if (bead_C<0)
    {
        bead_C=bead_C+N;
    }
    return bead_C;
}
//Function for finding the absolute of number
double absolute(double a){
    if (a<0)
    {
        a=-a;
    }
    return a;
}
//Function for finding the positions
void Position(double R,double A[N][2]){
    double theta;
    int i;
    theta=0;
    for (i=0;i<N;i++)
    {
        A[i][0]=R*cos(theta);A[i][1]=R*sin(theta);                 //x and y component
        theta=theta+2.0*PI/N;
    }
}
//Function for finding the angles of beads in membrane
void angle_M(double A[N][2],double Angles[N]){
    int i;
    for (i=0;i<N;i++)
    {
        if (A[i][1]>=0)
            Angles[i]=acos(A[i][0]/sqrt(A[i][0]*A[i][0]+A[i][1]*A[i][1]));
         if (A[i][1]<0)
            Angles[i]=2.0*PI-acos(A[i][0]/sqrt(A[i][0]*A[i][0]+A[i][1]*A[i][1]));
    }
}
//Function for defining rods
void rod(double B[M][2],double R,double angle,double x,double y){
    double r1[2],L;
    int i;
    L=2.0*R*M;                                            //length of the rod
    r1[0]=cos(angle);r1[1]=sin(angle);                    //direction of rod
    B[0][0]=x;B[0][1]=y;                                  //attached bead of rod to membrane
    for (i=1;i<M;i++)
    {
        B[i][0]=B[i-1][0]-2.0*R*r1[0];
        B[i][1]=B[i-1][1]-2.0*R*r1[1];
    }
}
//Function for calculating dot product
double dot(double r1[2],double r2[2]){
    double d;
    d=(r1[0]*r2[0]+r1[1]*r2[1]);
    return (d);
}
//Function for defining the cross product
double cross(double r1[2],double r2[2]){
    double d;
    d=r1[0]*r2[1]-r2[0]*r1[1];
    return (d);
}
//Function for defining the polymer radius
double radius(double A[N][2]){
    int i;
    double r,S=0.0;
    for (i=0;i<N;i++)
    {
        S=S+sqrt(A[i][0]*A[i][0]+A[i][1]*A[i][1]);
    }
    r=S/N;
    return (r);
}
//Function for mirror BC
void mirror1(double A[N][2],double AA[N+2][2]){
    int i,j;
    for (i=1;i<(N+1);i++)
    {
        for (j=0;j<2;j++)
        {
            AA[i][j]=A[i-1][j];
        }
    }
    AA[0][0]=A[N-1][0];AA[0][1]=A[N-1][1];
    AA[N+1][0]=A[0][0];AA[N+1][1]=A[0][1];
}
//Function for mirror BC
void mirror2(double A[N][2],double AA[N+4][2]){
    int i,j;
    for (i=2;i<(N+2);i++)
    {
        for (j=0;j<2;j++)
        {
            AA[i][j]=A[i-2][j];
        }
    }
    AA[1][0]=A[N-1][0];AA[1][1]=A[N-1][1];
    AA[0][0]=A[N-2][0];AA[0][1]=A[N-2][1];
    AA[N+2][0]=A[0][0];AA[N+2][1]=A[0][1];
    AA[N+3][0]=A[1][0];AA[N+3][1]=A[1][1];
}
//Function for omitting the tangential forces
void Radial_Force(double A[N][2],double F[N][2]){
    int i,j;
    double FF;
    for (i=0;i<N;i++)
    {
        FF=F[i][0]*(A[i][0])/sqrt(A[i][0]*A[i][0]+A[i][1]*A[i][1])+F[i][1]*(A[i][1])/sqrt(A[i][0]*A[i][0]+A[i][1]*A[i][1]);
        for (j=0;j<2;j++)
        {
            F[i][j]=FF*(A[i][j])/sqrt(A[i][0]*A[i][0]+A[i][1]*A[i][1]);
        }
    }
}
//Function for finding the distance between all particles
void Pdis(double A[N][2],double lave[N]){
    double AA[N+2][2],r1[2],r2[2],r3[2];
    int i;
    mirror1(A,AA);
    for (i=1;i<(N+1);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        r3[0]=AA[i-1][0];r3[1]=AA[i-1][1];
        lave[i-1]=0.5*(dis(r1,r3)+dis(r1,r2));
    }
}
//Function for SSLJ between membrane_one rod
void membrane_rod_SSLJ(double A[N][2],double B[M][2],double ALPHA_M,double E_M,double FR[3],double FSM[N][2]){
    double FF,d,r1[2],r2[2],S1=0.0,S2=0.0,S3=0.0,midr[2],R[M][2],FSR[M][2],rcut_M,ENERGY_CONST_M;
    int i,j,k;
    rcut_M=sqrt(pow(2.0,1.0/3.0)-ALPHA_M*ALPHA_M);
    ENERGY_CONST_M=(E_M/(4.0*pow(1./ALPHA_M,6)*(pow(1./ALPHA_M,6)-1.0)+1.0));
    for (i=0;i<N;i++)
    {
        for (k=0;k<2;k++)
        {
            FSM[i][k]=0;
        }
    }
    for (i=0;i<M;i++)
    {
        for (k=0;k<2;k++)
        {
            FSR[i][k]=0;
        }
    }
    //midr[0]=(B[0][0]+B[M-1][0])/2.0;midr[1]=(B[0][1]+B[M-1][1])/2.0;
    for (i=0;i<M;i++)
    {
        r1[0]=B[i][0];r1[1]=B[i][1];
        R[i][0]=r1[0]-B[0][0];
        R[i][1]=r1[1]-B[0][1];
    }
    for (i=0;i<N;i++)
    {
        r1[0]=A[i][0];r1[1]=A[i][1];
        for (j=1;j<M;j++)   //j start from 1 because the 0 bead is attached bead
        {
            r2[0]=B[j][0];r2[1]=B[j][1];
            d=dis(r1,r2);
            if (d<rcut_M)
            {
                for (k=0;k<2;k++)
                {
                    FF=24.0*ENERGY_CONST_M*(2.0*pow(ALPHA_M*ALPHA_M+d*d,-7)-pow(ALPHA_M*ALPHA_M+d*d,-4))*(r1[k]-r2[k]);
                    FSM[i][k]=FF+FSM[i][k];
                    FSR[j][k]=FSR[j][k]-FF;
                }
            }
        }
    }
    for (i=1;i<M;i++)
    {
        S1=S1+FSR[i][0];                          //bead level
        S2=S2+FSR[i][1];                          //bead level
        r1[0]=R[i][0];r1[1]=R[i][1];
        r2[0]=FSR[i][0];r2[1]=FSR[i][1];
        S3=S3+cross(r1,r2);                       //bead level
    }
    FR[0]=S1;      //rod level
    FR[1]=S2;      //rod level
    FR[2]=S3;      //rod level
}
//Function for SSLJ between one rod_ one rod
void rod_rod_SSLJ(double B[N_r][M][2],double ALPHA_r,double E_r,double FR[N_r][3]){
    double FF,d,r1[2],r2[2],S1=0.0,S2=0.0,S3=0.0,midr[N_r][2],R[N_r][M][2],rcut_r,ENERGY_CONST_r,FSR[N_r][M][2],BB[M][2];
    int i,j,k,ii,jj;
    for (i=0;i<N_r;i++)
    {
        for (j=0;j<M;j++)
        {
        for (k=0;k<2;k++)
        {
            FSR[i][j][k]=0.0;
        }
        }
    }
    rcut_r=sqrt(pow(2.0,1.0/3.0)-ALPHA_r*ALPHA_r);
    ENERGY_CONST_r=(E_r/(4.0*pow(1./ALPHA_r,6)*(pow(1./ALPHA_r,6)-1.0)+1.0));
    for (ii=0;ii<N_r;ii++)
    {
        for (jj=(ii+1);jj<N_r;jj++)
        {
            for (i=1;i<M;i++)
            {
                r1[0]=B[ii][i][0];r1[1]=B[ii][i][1];
                for (j=1;j<M;j++)
                {
                    r2[0]=B[jj][j][0];r2[1]=B[jj][j][1];
                    d=dis(r1,r2);
                    if (d<rcut_r)
                    {
                        for (k=0;k<2;k++)
                        {
                            FF=24.0*ENERGY_CONST_r*(2.0*pow(ALPHA_r*ALPHA_r+d*d,-7)-pow(ALPHA_r*ALPHA_r+d*d,-4))*(r1[k]-r2[k]);
                            FSR[ii][i][k]=FF+FSR[ii][i][k];
                            FSR[jj][j][k]=-FF+FSR[jj][j][k];
                        }
                    }
                }
            }
        }
    }
    for (i=0;i<N_r;i++)
    {
        for (j=0;j<M;j++)
        {
            for (k=0;k<2;k++)
            {
                BB[j][k]=B[i][j][k];
            }
        }
        for (j=0;j<M;j++)
        {
        r1[0]=BB[j][0];r1[1]=BB[j][1];
        R[i][j][0]=r1[0]-BB[0][0];
        R[i][j][1]=r1[1]-BB[0][1];
        }
    }
    for (j=0;j<N_r;j++)
    {
    for (i=1;i<M;i++)
    {
        S1=S1+FSR[j][i][0];                     //bead level
        S2=S2+FSR[j][i][1];                     //bead level
        r1[0]=R[j][i][0];r1[1]=R[j][i][1];
        r2[0]=FSR[j][i][0];r2[1]=FSR[j][i][1];
        S3=S3+cross(r1,r2);                     //bead level
    }
    FR[j][0]=S1;      //rod level
    FR[j][1]=S2;      //rod level
    FR[j][2]=S3;      //rod level
    S1=0.0;
    S2=0.0;
    S3=0.0;
    }
}
//Function for SSLJ
void SSLJ(double A[N][2],double B[N_r][M][2],double FR[N_r][3],double FM[N][2],double ALPHA_M,double E_M,double ALPHA_r,double E_r){
    double d,r1[2],r2[2],S1=0.0,S2=0.0,S3=0.0,FSM[N][2],FR1[3],FRR1[N_r][3],FR2[N_r][3],BB[M][2];
    int i,j,k;
    for (i=0;i<N_r;i++)
    {
    for (j=0;j<3;j++)
    {
        FRR1[i][j]=0;
    }
    }
    for (j=0;j<N;j++)
    {
        for (k=0;k<2;k++)
        {
            FM[j][k]=0;
        }
    }
    for (i=0;i<N_r;i++)
    {
        for (j=0;j<M;j++)
        {
            for (k=0;k<2;k++)
            {
                BB[j][k]=B[i][j][k];
            }
        }
        membrane_rod_SSLJ(A,BB,ALPHA_M,E_M,FR1,FSM);
        for (j=0;j<3;j++)
        {
            FRR1[i][j]=FR1[j]+FRR1[i][j];
        }
        for (j=0;j<N;j++)
        {
            for (k=0;k<2;k++)
            {
                FM[j][k]=FSM[j][k]+FM[j][k];
            }
        }
    }
    rod_rod_SSLJ(B,ALPHA_r,E_r,FR2);
    for (i=0;i<N_r;i++)
    {
        for (j=0;j<3;j++)
        {
            FR[i][j]=FRR1[i][j]+FR2[i][j];
        }
    }
}
//Function for being sure that all the beads of rod are in the cell
int bead_position(double B[M][2],double R){
    int i;
    double r,S=0.0;
    for (i=0;i<M;i++)
    {
        r=B[i][0]*B[i][0]+B[i][1]*B[i][1];
        if (r>(R*R))
        {
            S=1;
        }
            else
        {
        S=0.0;
        }
    }

    return S;
}
//Function for finding the stretching forces
void SForce(double A[N][2],double F[N][2],double ks,double ls){
    double r1[2],r2[2],r3[2],d1,d2,AA[N+2][2];
    int i,j;
    mirror1(A,AA);
    for (i=1;i<(N+1);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        r3[0]=AA[i-1][0];r3[1]=AA[i-1][1];
        d1=dis(r1,r2);
        d2=dis(r1,r3);
        for (j=0;j<2;j++)
        {
            F[i-1][j]=-ks*((r1[j]-r3[j])*(d2-ls)/d2+(r1[j]-r2[j])*(d1-ls)/d1);
        }
    }
   Radial_Force(A,F);
}
//Function for finding the bending forces
void BForce(double A[N][2],double F[N][2],double k) {
    double r1[2],r2[2],r3[2],r4[2],r5[2],r6[2],r7[2],r8[2],d1,d2,d3,d4,AA[N+4][2];
    int i,j;
    mirror2(A,AA);
        for (i=2;i<(N+2);i++)
        {
            r1[0]=AA[i][0];r1[1]=AA[i][1];
            r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
            r3[0]=AA[i-1][0];r3[1]=AA[i-1][1];
            r4[0]=AA[i-2][0];r4[1]=AA[i-2][1];
            r5[0]=AA[i+2][0];r5[1]=AA[i+2][1];
            r6[0]=r5[0]-2*r2[0]+r1[0];r6[1]=r5[1]-2*r2[1]+r1[1];
            r7[0]=r2[0]-2*r1[0]+r3[0];r7[1]=r2[1]-2*r1[1]+r3[1];
            r8[0]=r1[0]-2*r3[0]+r4[0];r8[1]=r1[1]-2*r3[1]+r4[1];
            d1=dis(r1,r2);
            d2=dis(r1,r3);
            d3=dis(r2,r5);
            d4=dis(r3,r4);
            for (j=0;j<2;j++)
            {
                F[i-2][j]=(12.0*k/(pow((d1+d2),4))*((r1[j]-r3[j])/d1-(r2[j]-r1[j])/d3)*(dot(r6,r6)+dot(r7,r7)+dot(r8,r8))-8.0*k/pow((d1+d3),3)*(r5[j]-4.0*r2[j]+6.0*r1[j]-4.0*r3[j]+r4[j]))/2.5;
            }
        }
        Radial_Force(A,F);
}
//Function for finding the Length Force
void LForce(double A[N][2],double F[N][2],double kl,double L0){
    double AA[N+2][2],r1[2],r2[2],r3[2],d1,d2,S=0;
    int i,j;
    mirror1(A,AA);
    for (i=1;i<(N+1);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        d1=dis(r1,r2);
        S=S+d1;
    }
    for (i=1;i<(N+1);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        r3[0]=AA[i-1][0];r3[1]=AA[i-1][1];
        d1=dis(r1,r3);
        d2=dis(r1,r2);
        for (j=0;j<2;j++)
        {
            F[i-1][j]=-kl*(S-L0)*((r1[j]-r3[j])/d1-(r2[j]-r1[j])/d2);
        }
    }
    Radial_Force(A,F);
}
//Function for finding the Area Force
void AForce(double A[N][2],double F[N][2],double ka,double A0){
    double AA[N+2][2],flag,sign,r1[2],r2[2],r3[2],r4[2],r5[2],S=0.0,Area;
    int i,j;
    mirror1(A,AA);
    for (i=1;i<(N+1);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        Area=(cross(r1,r2))/2.0;
        S=Area+S;
    }
    for (i=1;i<(N+1);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        r3[0]=AA[i-1][0];r3[1]=AA[i-1][1];
        r4[0]=r2[0]-r1[0];r4[1]=r2[1]-r1[1];
        r5[0]=r1[0]-r3[0];r5[1]=r1[1]-r3[1];
        flag=cross(r1,r4);
        if (flag>0)
        {
            sign=1;
        }
        else
        {
            sign=-1;
        }
        for (j=0;j<2;j++)
        {
            F[i-1][j]=-ka*(S-A0)*(sign*(dot(r3,r3)*r1[j]-dot(r1,r3)*r3[j])/(2*absolute(cross(r3,r5)))+sign*(dot(r2,r2)*r1[j]-dot(r1,r2)*r2[j])/(2.0*absolute(cross(r1,r4))));
        }
    }
    Radial_Force(A,F);
}
//Function for finding the stretching energy
double SEnergy(double A[N][2],double ks,double ls){
    double AA[N+2][2],r1[2],r2[2],d1,S=0;
    int i;
    mirror1(A,AA);
    for (i=1;i<(N+1);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        d1=dis(r1,r2);
        S=S+0.5*ks*(d1-ls)*(d1-ls);
    }
    return (S);
}
//Function for finding the bending energy
double BEnergy(double A[N][2],double kk){
    double kb,AA[N+4][2],r1[2],r2[2],r3[2],r4[2],r5[2],d1,d2,Lave,S=0;
    int i;
    mirror2(A,AA);
    for (i=2;i<(N+2);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        r3[0]=AA[i-1][0];r3[1]=AA[i-1][1];
        r4[0]=AA[i+2][0];r4[1]=AA[i+2][1];
        d1=dis(r1,r2);
        d2=dis(r1,r3);
        r5[0]=r4[0]-2*r2[0]+r1[0];r5[1]=r4[1]-2*r2[1]+r1[1];
        Lave=(d1+d2)/2;
        kb=kk/pow(Lave,3);
        S=S+kb/2*dot(r5,r5);
    }
    return (S);
}
//Function for defining the area energy
double AreaEnergy(double A[N][2],double ka,double A0){
    double Area,r1[2],r2[2],AA[N+2][2],S=0,E;
    int i;
    mirror1(A,AA);
    for (i=1;i<(N+1);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        Area=0.5*(cross(r1,r2));
        S=S+(Area);
    }
    E=0.5*ka*(S-A0)*(S-A0);
    return (E);
}
//Function for defining the length energy
double LengthEnergy(double A[N][2],double kl,double L0){
    double AA[N+2][2],E,S=0,r1[2],r2[2],d1;
    int i;
    mirror1(A,AA);
    for (i=1;i<(N+1);i++)
    {
        r1[0]=AA[i][0];r1[1]=AA[i][1];
        r2[0]=AA[i+1][0];r2[1]=AA[i+1][1];
        d1=dis(r2,r1);
        S=S+d1;
    }
    E=0.5*kl*(S-L0)*(S-L0);
    return (E);
}
//Function for finding the analytical stretching energy                   //For checking
double ASEnergy(double A[N][2],double ks,double ls){
    double E,theta,R;
    R=radius(A);
    theta=2*PI/N;
    E=ks*(N/2)*(2*R*sin(theta/2)-ls)*(2*R*sin(theta/2)-ls);
    return (E);
}
//Function for finding the analytical bending energy                      //For checking
double ABEnergy(double A[N][2],double kb){
    double E,theta,R;
    R=radius(A);
    theta=2*PI/N;
    E=kb*N*(1-cos(theta))/(2*R*sin(theta/2));
    return E;
}
//Function for calculating the analitical Area energy                     //For checking
double AAreaEnergy(double A[N][2],double ka,double A0){
    double theta,E,R;
    theta=2*PI/N;
    R=radius(A);
    E=0.5*ka*(N*R*R*sin(theta/2)*cos(theta/2)-A0)*(N*R*R*sin(theta/2)*cos(theta/2)-A0);
    return E;
}
//Function for calculating the analitical Length energy                   //For checking
double ALengthEnergy(double A[N][2],double kl,double L0){
    double R,theta,E;
    theta=2*PI/N;
    R=radius(A);
    E=0.5*kl*(2*N*R*sin(theta/2)-L0)*(2*N*R*sin(theta/2)-L0);
    return E;
}
//Function for calculating the analitical Stretch Force                   //For checking
double ASForce(double A[N][2],double ks,double ls){
    double F,theta,R;
    R=radius(A);
    theta=2*PI/N;
    F=-2*ks*sin(theta/2)*(2*R*sin(theta/2)-ls);
    return F;
}
//Function for calculating the analitical Bending Force                   //For checking
double ABForce(double A[N][2],double kb){
    double F,theta,R;
    R=radius(A);
    theta=2*PI/N;
    F=kb*(1-cos(theta))/(2*R*R*sin(theta/2));
    return F;
}
//Function for caculating the analytical Area Force                       //For checking
double AAForce(double A[N][2],double ka,double A0){
    double F,R,theta;
    R=radius(A);
    theta=2*PI/N;
    F=-ka*(2*R*sin(theta/2)*cos(theta/2))*(R*R*N*sin(theta/2)*cos(theta/2)-A0);
    return F;
}
//Function for calculating the analytical Length Force                    //For checking
double ALForce(double A[N][2],double kl,double L0){
    double F,R,theta;
    R=radius(A);
    theta=2*PI/N;
    F=-kl*(2*sin(theta/2))*(2*R*N*sin(theta/2)-L0);
    return F;
}

using namespace std;

int main(){

    /* ******************************************Parameters**************************************************** */
    double Rs,ks,kk,ls,theta0,R,gama,t,ka,kl,L0,kb,A0,sigma02,noise_var,sina[N_r],cosa[N_r],rcut,E_M,gama0,alpha_M;
    double R_A_M[N_r][2],F_A_M_R[2],F_A_M_T,FFLUC,V_r_xy[2],V_r_T,r_d,bead_d,r_bead1[2],r_bead2[2],remain_d,V_r_pa;
    double V_r_pe,A[N][2],V[N][2],F[N][2],FS[N][2],FB[N][2],FL[N][2],FA[N][2],lave[N],FFLUC_M[N][2],FFLUC_r[N_r][3];
    double H[N],B[N_r][M][2],noise_var_r[3],V_r[N_r][3],r_angle[N_r],F1,F2,F3,FSSLJ_r[N_r][3],FSSLJ_M[N][2],Angles[N];
    double F_A_M[N_r][2],F_r_xy[2],angle,L,x,y,gamar[3],PForce,R_r,alpha_r,E_r,BB[M][2];
    int i,j,bead_C[N_r],K,k,ii,iii,jj,fluc,z;
    /* ******************************************Parameters**************************************************** */

    /* ************************************************Files*************************************************** */
    ofstream sepehr;
    sepehr.open("Positions.XYZ");
    ofstream sepehr2;
    sepehr2.open("Fluctuation.txt");
    ofstream sepehr3;
    sepehr3.open("rod.txt");
    /* ************************************************Files*************************************************** */

    /* *******************************************Initial Values*********************************************** */

    srand(time(0));                                                   //Seeding the rand function every second
    ks=1.0;                                                           //Spring constant
    kk=1.0;                                                           //Bending constant
    ka=1.0;                                                           //Area stifness
    kl=1.0;                                                           //Length stifness
    kb=kk;                                                            //Bending stifness (analytical formula)
    gama0=1.00;                                                       //Friction coefficent
    theta0=(2.0*PI)/N;                                                //Equilbrium angle
    alpha_M=0.5;                                                      //alpha in Energy SSLJ membrane_rod
    E_M=70.0;                                                         //Energy barrier membrane_rod
    alpha_r=1.0;//alpha_M;                                             //alpha in Energy SSLJ rod_rod
    R_r=sqrt(pow(2.0,1.0/3.0)-alpha_r*alpha_r);                       //Radius of beads in rod
    E_r=20.0;                                                          //energy barrier SSLJ rod_rod
    //ls=2.0*R_r;                                                       //bond length
    //Rs=ls/(2.0*sin(theta0/2.0));                                      //membrane radius
    Rs=10;
    R=10.000014472401796;                                                             //Radious of ring
    ls=2.0*Rs*sin(theta0/2.0);                                       //Equilibrium length
    L0=N*ls;                                                          //equilinrium length of membrane
    A0=N*Rs*Rs*sin(theta0/2.0)*cos(theta0/2.0);                       //equilibrium area
    PForce=100.00;                                                     //Prpulsion Force
    L=M*2.0*R_r;                                                      //Length of rod
    /* *******************************************Initial Values*********************************************** */

    /* ******************************************Initial Positions********************************************* */
    Position(R,A);                                                     //Membrane positioning
    /*bead_C[0] = rand()%N;                                              //attached bead of membrane
    do{
    r_angle[0]=(rand()/double(RAND_MAX+1.0))*360.0*PI/180.0;           //angle of rod between 0 and 360 degree
    x=A[bead_C[0]][0];                                                 //position of attached bead of rod x
    y=A[bead_C[0]][1];                                                 //position of attached bead of rod y
    rod(BB,R_r,r_angle[0],x,y);
    }
    while(bead_position(BB,R)==1);
    for (j=0;j<M;j++)
    {
        for (k=0;k<2;k++)
        {
            B[0][j][k]=BB[j][k];                                        //different rods
        }
    }*/
    for (i=0;i<N_r;i++)
    {
    //do
    //{
    bead_C[i] = rand()%N;                                              //attached bead of membrane
    do{
    r_angle[i]=(rand()/double(RAND_MAX+1.0))*360.0*PI/180.0;           //angle of rod between 0 and 360 degree
    x=A[bead_C[i]][0];                                                 //position of attached bead of rod x
    y=A[bead_C[i]][1];                                                 //position of attached bead of rod y
    rod(BB,R_r,r_angle[i],x,y);
    }
    while(bead_position(BB,R)==1);
    //}
    //while (abs(bead_C[i]-bead_C[i-1])<M);
    for (j=0;j<M;j++)
    {
        for (k=0;k<2;k++)
        {
            B[i][j][k]=BB[j][k];                                        //different rods
        }
    }
    }

    /* ******************************************Initial Positions********************************************* */

    /* *****************************************Initial rod parameters***************************************** */
    t=0.00001;                                                              //Time step
    gamar[0]=L*gama0;                                                   //friction coeficient parallel
    gamar[1]=2.0*gamar[0];                                              //friction coeficient prependuclar
    gamar[2]=gamar[0]*L*L/6.0;                                          //friction coeficient angular
    sigma02=2.0*0.02/(t*gama0);
    noise_var_r[0]=sqrt(sigma02*L);
    noise_var_r[1]=sqrt(2.0*sigma02*L);
    noise_var_r[2]=sqrt(2.0*sigma02*L*L*L/12.0);
    /* *****************************************Initial rod parameters***************************************** */

    /* *******************************************Numerical solution******************************************* */
    K     =     10000;                                                    //Iterations
    /* *******************************************Numerical solution******************************************* */

    /* *******************************************Start of iteration******************************************* */
    for (ii=0;ii<K;ii++)
    {
    cout<<ii<<"\n";                                                    //Iteration

    /* *************************************************results1*********************************************** */
    if (ii%10==0)                                                      //number of time step for getting the frame
    {
       // sepehr<<"title = 'membrane_rod'"<<"\n"<<"variables = 'x', 'y'"<<"\n"<<"zone "<<"i="<<(N)<<", f=point"<<"\n";
        sepehr<<(N+N_r*M)<<"\n"<<"Atoms. Timestep: "<<ii<<"\n";
    for (z=0;z<N_r;z++)
    {
    for (i=0;i<(M);i++)
    {
        sepehr<<"R "<<B[z][i][0]<<" "<<B[z][i][1]<<" "<<"0"<<"\n";//<<B[i][1]<<"\n";                   //rod x
    }
    }
    //for (z=0;z<N_r;z++)
    //{
    //for (i=0;i<(M);i++)
    //{
    //    sepehr3<<B[z][i][1]<<"\n";                                     //rod y
    //}
    //}
    for (i=0;i<N;i++)
    {
      sepehr<<"M "<<A[i][0]<<" "<<A[i][1]<<" "<<"0"<<"\n";                        //membrane x
    }
    //for (i=0;i<N;i++)
   // {
   //    sepehr<<A[i][1]<<"\n";                                         //membrane y
    //}
    }
    /* *************************************************results1*********************************************** */

    /* ********************************************Membrane forces********************************************* */
        SForce(A,FS,ks,ls);                                             //stretching forces
        BForce(A,FB,kk);                                                //bending forces
        LForce(A,FL,kl,L0);                                             //length constrained forces
        AForce(A,FA,ka,A0);                                             //area constrained forces
    /* ********************************************Membrane forces********************************************* */

        Pdis(A,lave);                                                   //bond length

    /* *******************************************Membrane noise*********************************************** */
        for (fluc=0;fluc<N;fluc++)
        {
            noise_var=sqrt(2.0*sigma02*lave[fluc]);                     //noise variances
            FFLUC=noise_var*muller();                                   //This force must be radial
            FFLUC_M[fluc][0]=FFLUC*A[fluc][0]/sqrt(A[fluc][0]*A[fluc][0]+A[fluc][1]*A[fluc][1]);
            FFLUC_M[fluc][1]=FFLUC*A[fluc][1]/sqrt(A[fluc][0]*A[fluc][0]+A[fluc][1]*A[fluc][1]);
        }
    /* *******************************************Membrane noise*********************************************** */

    /* *********************************************rod noise************************************************** */
        for (i=0;i<N_r;i++)
        {
        F1=muller()*noise_var_r[0]+PForce;                            //Noise and propolsion in parallel
        F2=muller()*noise_var_r[1];                                   //Noise in prepandicular
        F3=muller()*noise_var_r[2];                                   //Noise in theta
        cosa[i]=(B[i][0][0]-B[i][M-1][0])/sqrt((B[i][M-1][0]-B[i][0][0])*(B[i][M-1][0]-B[i][0][0])+(B[i][M-1][1]-B[i][0][1])*(B[i][M-1][1]-B[i][0][1]));
        sina[i]=(B[i][0][1]-B[i][M-1][1])/sqrt((B[i][M-1][0]-B[i][0][0])*(B[i][M-1][0]-B[i][0][0])+(B[i][M-1][1]-B[i][0][1])*(B[i][M-1][1]-B[i][0][1]));
        FFLUC_r[i][0]=F1*cosa[i]-F2*sina[i];                                 //Total force in x direction from the rod
        FFLUC_r[i][1]=F1*sina[i]+F2*cosa[i];                                 //Total force in y direction from the rod
        FFLUC_r[i][2]=F3;                                              //Total force in theta direction from the rod
        }
        SSLJ(A,B,FSSLJ_r,FSSLJ_M,alpha_M,E_M,alpha_r,E_r);                           //SSLJ force for both membrane and the rod
    /* *********************************************rod noise************************************************** */

    /* ****************************************sum of membrane forces****************************************** */
        for (iii=0;iii<N;iii++)
        {
            for (jj=0;jj<2;jj++)
            {
                F[iii][jj]=FFLUC_M[iii][jj]+FS[iii][jj]+FB[iii][jj]+FL[iii][jj]+FA[iii][jj]+FSSLJ_M[iii][jj];   //All the forces that are applied to the membrane
            }
        }
    /* ****************************************sum of membrane forces****************************************** */

    /* **********************************************Interactoion********************************************** */
        for (i=0;i<N_r;i++)
        {
        F_A_M[i][0]=FFLUC_r[i][0]+FSSLJ_r[i][0];                                           //Parallel rod force x
        F_A_M[i][1]=FFLUC_r[i][1]+FSSLJ_r[i][1];                                           //Parallel rod force y
        R_A_M[i][0]=A[bead_C[i]][0];R_A_M[i][1]=A[bead_C[i]][1];                             //Attached bead
        Radial_Force_single(F_A_M[i],R_A_M[i],F_A_M_R);                                    //Making the parallel force of rod radial
        F[bead_C[i]][0]=F[bead_C[i]][0]+F_A_M_R[0];                                        //Adding the radial parallel force of rod to the attached bead x
        F[bead_C[i]][1]=F[bead_C[i]][1]+F_A_M_R[1];                                        //Adding the radial parallel force of rod to the attached bead y
        }
    /* **********************************************Interactoion********************************************** */

    /* *************************************************Sliding************************************************ */
        for (i=0;i<N_r;i++)
        {
        F_r_xy[0]=(FFLUC_r[i][0]+FSSLJ_r[i][0]);                             //All forces of rod x
        F_r_xy[1]=(FFLUC_r[i][1]+FSSLJ_r[i][1]);                             //All forces of rod y
        V_r_pa=(F_r_xy[0]*cosa[i]+F_r_xy[1]*sina[i])/gamar[0];               //velocity of the rod parallel component
        V_r_pe=(F_r_xy[1]*cosa[i]-F_r_xy[0]*sina[i])/gamar[1];               //velocity of the rod prependicular component
        V_r[i][0]=V_r_pa*cosa[i]-V_r_pe*sina[i];                             //velocity of the rod x
        V_r[i][1]=V_r_pa*sina[i]+V_r_pe*cosa[i];                             //velocity of the rod y
        V_r_T=Tangential_Force_single(V_r[i],R_A_M[i]);                      //tangential parallel force of rod
        bead_C[i]=changing_attachement(A,bead_C[i],V_r_T,t);                 //changing the attached bead
        }
    /* *************************************************Sliding************************************************ */

    /* **********************************************repositioning********************************************* */
        for (i=0;i<N_r;i++)
        {
        V_r[i][2]=(FFLUC_r[i][2]+FSSLJ_r[i][2])/gamar[2];                    //angular velocity of the rod
        r_angle[i]=V_r[i][2]*t+r_angle[i];                                   //new orientation of the rod
        for (j=0;j<M;j++)
        {
        for (k=0;k<2;k++)
        {
            BB[j][k]=B[i][j][k];
        }
        }
        rod(BB,R_r,r_angle[i],A[bead_C[i]][0],A[bead_C[i]][1]);              //new position of the rod
        for (j=0;j<M;j++)
        {
        for (k=0;k<2;k++)
        {
            B[i][j][k]=BB[j][k];
        }
        }
        }
       for (j=0;j<N;j++)
        {
            gama=2*lave[j]*gama0;                                   //friction coefficient mebranes bead
            for (k=0;k<2;k++)
            {
                V[j][k]=F[j][k]/gama;                               //velocity of membrane beads
                A[j][k]=A[j][k]+t*V[j][k];                          //new position of membrane beads
            }
            H[j]=sqrt(A[j][0]*A[j][0]+A[j][1]*A[j][1])-R;           //height of beads
        }
    /* **********************************************repositioning********************************************* */

    /* *************************************************checking*********************************************** */
        for (i=0;i<N_r;i++)
        {
            for (j=0;j<M;j++)
            {
                for (k=0;k<2;k++)
                {
                    BB[j][k]=B[i][j][k];
                }
            }
                if (bead_position(BB,R)==1)
                {
                cout<<"Rod crossed the membrane :-)";                 //rod is crossing the membrane or not
                return 0;
                }
        }
        if (A[0][0]>(1.5*R))
        {
            cout<<"It is going to be unsteady :-( "<<"\n";           //unsteady answer
            return 0;
        }
    /* *************************************************checking*********************************************** */

    /* *************************************************results2*********************************************** */
        if (ii%50==0 && ii>10)
        {
        for (i=0;i<N;i++)
        {
        sepehr2<<H[i]<<"\n";                                       //Fluctuation
        }
        }
    /* *************************************************results2*********************************************** */
}   /* ********************************************end of iteration******************************************* */
}

