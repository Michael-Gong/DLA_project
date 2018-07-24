#include <iostream>
//#include <mpi.h>
#include "/usr/include/mpich-x86_64/mpi.h"
#include <fstream>
#include <math.h>
#include <ctime>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <limits>
#include "ovector.h"

// ########constant defined here ########
#define    pi      3.1415926535897932384626
#define    pi2d    57.296
#define    q0      1.6e-19 // C
#define    m0      9.1e-31 // kg
#define    v0      2.99792458e8 // m/s
#define    kb      1.38065e-23  // J/K
#define    mu0     1.257e-6     // N/A^2
#define    epsilon0  8.854e-12  // F/m
#define    h_planck  6.626e-34  // J*s
#define    wavelength  1.0e-6   // m

using namespace std;

void getdata1d(string name, double data[], int size) {
    ifstream file(name);
    if(file.is_open()){
        for(int i = 0; i < size; ++i){
		    file >> data[i];
			//cout << i << "th :"<< data[i] <<endl;
		}
    }
	else
	    cout<<"shit!"<<endl;
}

void linspace(double data[], double start, double end, int number){
    if(number == 1){
	    data[0] = start;
	    return;
	}
    double dstep = (end-start)/(double)(number-1);
    for(int i=0; i<number; ++i){
	    data[i] = dstep*i+start;
	}
	return;
}


void log10space(double data[], double start, double end, int number){ // start and end is the number of exponential 
    double dstep = (end-start)/(double)(number-1);
    for(int i=0; i<number; ++i){
	    data[i] = pow(10,dstep*i+start);
	}
	return;
}

int main(int argc, char **argv)
{
    double frequency = v0*2.0*pi/wavelength;
	double exunit    = m0*v0*frequency/q0;
	double bxunit    = m0*frequency/q0;
    double denunit   = frequency*frequency*epsilon0*m0/(q0*q0);
	cout<<"electric field unit: "<< setprecision(24) << exunit <<endl;
	cout<<"magnetic field unit: "<< setprecision(24) << bxunit <<endl;
	cout<<"density unit nc: "    << setprecision(24) << denunit<<endl;
    
	int       data_size = 20005;
    string    dataname ("./Data/");
    double    px[data_size];    getdata1d(dataname+"px_0.txt", px, data_size);    
    double    py[data_size];    getdata1d(dataname+"py_0.txt", py, data_size);    
    double    pz[data_size];    getdata1d(dataname+"pz_0.txt", pz, data_size);    
    double    grid_x[data_size];    getdata1d(dataname+"x_0.txt", grid_x, data_size); 
    double    grid_y[data_size];    getdata1d(dataname+"y_0.txt", grid_y, data_size); 
    double    grid_z[data_size];    getdata1d(dataname+"z_0.txt", grid_z, data_size); 
    double    timett[data_size];    getdata1d(dataname+"t_0.txt", timett, data_size); 
    double    part_ex[data_size];   getdata1d(dataname+"ex_part_0.txt", part_ex, data_size);
    double    part_ey[data_size];   getdata1d(dataname+"ey_part_0.txt", part_ey, data_size);
    double    part_ez[data_size];   getdata1d(dataname+"ez_part_0.txt", part_ez, data_size);
    double    part_bx[data_size];   getdata1d(dataname+"bx_part_0.txt", part_bx, data_size);
    double    part_by[data_size];   getdata1d(dataname+"by_part_0.txt", part_by, data_size);
    double    part_bz[data_size];   getdata1d(dataname+"bz_part_0.txt", part_bz, data_size);
    double    gg[data_size];
	double    e_vx[data_size],  e_vy[data_size],  e_vz[data_size];
	double    e_ax[data_size],  e_ay[data_size],  e_az[data_size];
    cout<<"here"<<endl;   
	for(int i = 0; i < data_size; ++i){ 
	    gg[i]    = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+1.0);
		e_vx[i]  = px[i]/gg[i];
		e_vy[i]  = py[i]/gg[i];
		e_vz[i]  = pz[i]/gg[i];
        //e_ax[i]  = (-part_ex[i]-(e_vy[i]*part_bz[i]-e_vz[i]*part_by[i])-e_vx[i]*(-e_vx[i]*part_ex[i]-e_vy[i]*part_ey[i]-e_vz[i]*part_ez[i]))/gg[i];
		//e_ay[i]  = (-part_ey[i]-(e_vz[i]*part_bx[i]-e_vx[i]*part_bz[i])-e_vy[i]*(-e_vx[i]*part_ex[i]-e_vy[i]*part_ey[i]-e_vz[i]*part_ez[i]))/gg[i];
		//e_az[i]  = (-part_ez[i]-(e_vx[i]*part_by[i]-e_vy[i]*part_bx[i])-e_vz[i]*(-e_vx[i]*part_ex[i]-e_vy[i]*part_ey[i]-e_vz[i]*part_ez[i]))/gg[i];
                e_ax[i]  = (-part_ex[i]-(e_vy[i]*part_bz[i]-e_vz[i]*part_by[i])-e_vx[i]*(-e_vx[i]*part_ex[i]-e_vy[i]*part_ey[i]-e_vz[i]*part_ez[i]))/gg[i];
		e_ay[i]  = (-part_ey[i]-(e_vz[i]*part_bx[i]-e_vx[i]*part_bz[i])-e_vy[i]*(-e_vx[i]*part_ex[i]-e_vy[i]*part_ey[i]-e_vz[i]*part_ez[i]))/gg[i];
		e_az[i]  = (-part_ez[i]-(e_vx[i]*part_by[i]-e_vy[i]*part_bx[i])-e_vz[i]*(-e_vx[i]*part_ex[i]-e_vy[i]*part_ey[i]-e_vz[i]*part_ez[i]))/gg[i];
	}
    cout<<"here"<<endl;   
    
    int       x_size=60001,  y_size=1,  z_size=1;
    double    x_omega[x_size], y_theta[y_size], z_phi[z_size];
	double    data_I[x_size][y_size][z_size], data_I_t[x_size][y_size][z_size];
    double    norm_fac=1.0/4.0/3.14/epsilon0*q0*q0/4.0/(3.14*3.14)/v0*frequency/(m0*v0*v0);
    cout<<"here"<<endl;   
    linspace(x_omega, 0., 6000., x_size);
    linspace(y_theta, 0./pi2d, 0./pi2d, y_size);
    linspace(z_phi,   0./pi2d, 0./pi2d, z_size);
    cout<<"here"<<endl;   
//############ To calculate integral of the dI/domega##########
    Vector3 amplitude1, data_1, data_2;
	double dt, e_tim_t, phase2;
    for(int i_phi=0; i_phi<z_size; ++i_phi){
      for(int i_theta=0; i_theta<y_size; ++i_theta){
        for(int i_omega=0; i_omega<x_size; ++i_omega){
    //     cout<<"i_phi"<<z_phi[i_phi]<<endl;
	       Vector3 n_dirc( cos(y_theta[i_theta]), sin(y_theta[i_theta])*cos(z_phi[i_phi]), sin(y_theta[i_theta])*sin(z_phi[i_phi]) );
	       //cout<<"n_dirc: x:"<<n_dirc.e1()<<"; y:"<<n_dirc.e2()<<"; z:"<<n_dirc.e3()<<endl;
	       dt=timett[1]-timett[0]; data_1=0; data_2=0;
	       for(int i_time=0; i_time<data_size; ++i_time){
	         Vector3 e_vel_t(e_vx[i_time], e_vy[i_time], e_vz[i_time]);
			 if(isnan(e_vel_t.e1())) cout<<"nan emerges at e_vel_t!"<<endl;
		     //cout<<"vx:"<<e_vel_t.e1()<<"; vy:"<<e_vel_t.e2()<<"; vz:"<<e_vel_t.e3()<<endl;
	         Vector3 e_acc_t(e_ax[i_time], e_ay[i_time], e_az[i_time]);
			 if(isnan(e_acc_t.e1())) cout<<"nan emerges at e_acc_t!"<<endl;
		     //cout<<"ax:"<<e_acc_t.e1()<<"; ay:"<<e_acc_t.e2()<<"; az:"<<e_acc_t.e3()<<endl;
	         Vector3 e_pos_t(grid_x[i_time], grid_y[i_time], grid_z[i_time]);
			 if(isnan(e_pos_t.e1())) cout<<"nan emerges at e_pos_t!"<<endl;
		     //cout<<"xx:"<<e_pos_t.e1()<<"; yy:"<<e_pos_t.e2()<<"; zz:"<<e_pos_t.e3()<<endl;
             e_tim_t = timett[i_time];
		     //cout<<"time:"<<e_tim_t<<endl;
			 //Vector3 temp_vec = (n_dirc-e_vel_t).cross(e_acc_t);
			 Vector3 temp_vec1 = (n_dirc-e_vel_t);
			 Vector3 temp_vec2 = temp_vec1.cross(e_acc_t);
			 //if(isnan(temp_vec2.e1())) cout<<"nan emerges!"<<endl;
	         amplitude1 = n_dirc.cross((n_dirc-e_vel_t).cross(e_acc_t))/pow(1.0-e_vel_t*n_dirc,2);
			 //if(isnan(amplitude1.e1())) cout<<"nan emerges!"<<endl;
			// cout<<"amplitude_x:"<<amplitude1.e1()<<"; amplitude_y:"<<amplitude1.e2()<<endl;
		     phase2     = (e_tim_t-n_dirc*e_pos_t)*x_omega[i_omega];
		 	 //cout<<"phase2:"<<phase2<<endl;
			 data_1     = data_1+amplitude1*cos(phase2)*dt;
			 data_2     = data_2+amplitude1*sin(phase2)*dt;
	//	 	 cout<<"time: "<<e_tim_t<<"; data_1^2:"<<data_1*data_1<<"; data_2^2:"<<data_2*data_2<<endl;
	//	 	 cout<<"time: "<<e_tim_t<<"; vel_1:"<<e_vel_t.e1()<<"; acc_1:"<<e_acc_t.e1()<<"; pos_1:"<<e_pos_t.e1()<<endl;
		   }
		   data_I[i_omega][i_theta][i_phi]  = (data_1*data_1+data_2*data_2);
	           if(isnan(data_I[i_omega][i_theta][i_phi])) cout<<"nan emerges at data_I!"<<endl;
		  // data_I_t[i_phi][i_theta][i_omega]= data_I[i_phi][i_theta][i_omega]*sin(y_theta[i_theta])*norm_fac*(y_theta[1]-y_theta[0])*(z_phi[1]-z_phi[0]);
		  // cout<<"i_phi:"<<z_phi[i_phi]<<"i_theta:"<<y_theta[i_theta]<<"i_omega:"<<x_omega[i_omega]<<"; dI/dw:"<<data_I[i_omega][i_theta][i_phi]<<endl;
		   cout<<"finished "<<double(i_omega+1+i_theta*x_size+i_phi*y_size*x_size)/double(x_size*y_size*z_size)*100.0<<" %"<<endl;
		}
	  }
    }
    //#############################################################

    char filename[100]; //

    sprintf(filename,"./C-Rad/grid_omega_x.txt");
    FILE *fp_omega = fopen(filename,"a");
    for(int i_omega=0; i_omega<x_size; ++i_omega){
        fprintf(fp_omega,"%lf\n",x_omega[i_omega]);  
    }
    fclose(fp_omega);
    fprintf(stdout, "grid_omega_x is already written into ./C-Rad/grid_omega_x.txt\n");

    sprintf(filename,"./C-Rad/grid_theta_y.txt");
    FILE *fp_theta = fopen(filename,"a");
    for(int i_theta=0; i_theta<y_size; ++i_theta){
        fprintf(fp_theta,"%lf\n",y_theta[i_theta]);  
    }
    fclose(fp_theta);
    fprintf(stdout, "grid_theta_y is already written into ./C-Rad/grid_theta_y.txt\n");

    sprintf(filename,"./C-Rad/grid_phi_z.txt");
    FILE *fp_phi   = fopen(filename,"a");
    for(int i_phi=0; i_phi<z_size; ++i_phi){
        fprintf(fp_phi,"%lf\n",z_phi[i_phi]);  
    }
    fclose(fp_phi);
    fprintf(stdout, "grid_phi_z  is already written into ./C-Rad/grid_phi_z.txt\n");

    sprintf(filename,"./C-Rad/data.txt");
    FILE *fp_data = fopen(filename,"a");
    for(int i_omega=0; i_omega<x_size; ++i_omega){
        for(int i_theta=0; i_theta<y_size; ++i_theta){
            for(int i_phi=0; i_phi<z_size; ++i_phi){
                //fprintf(fp_data,"%lf\n",data_I[i_omega][i_theta][i_phi]*sin(y_theta[i_theta])*norm_fac*(y_theta[1]-y_theta[0])*(z_phi[1]-z_phi[0]));      
                fprintf(fp_data,"%lf\n",data_I[i_omega][i_theta][i_phi]);      
            }
        }
    }
    fclose(fp_data);
    fprintf(stdout, "data  is already written into ./C-Rad/data.txt\n");
//        for(int i = 0; i < 4002; ++i)
//	    {
//			cout << i << "th :"<< e_field[i] <<endl;
//		}
    return 0;
}
