#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include "planet.h"

using namespace std;

//We wish to simulate the solar system, using Runge-Kutta 4.
//Applying Newton's second law of motion, we get a set of differential equations.
//These we rewrite into two separate first order

void f(double* A, double *dA, double* M, int m, int n);

int main(int argc, char* argv[]) {

    //Takes in a vector A, containing position and velocity of each planet
    //r which is the distance between each planet, and M, the masses.

    /*fstream inFile;
    inFile.open("distance_from_sun_1.dat", ios::out);
    inFile.close();*/


    int n = atoi(argv[1]);          //Number of planets
    double t_max = atof(argv[2]);   //Number of years
    int k = atoi(argv[3]);          //Case
    int m = 4*n;
    int N = 100000;
    double dt;
    dt = t_max/N;
    cout << dt << endl;
    double pi;
    pi = atan(1)*4;
    double t_0 = 0;
    double h = (t_max - t_0)/N;
    double *A = new double[4*n];
    double *M = new double[n];
    Planet Earth = Planet(1.0, 0.0, 0.0, 2*pi, 3.0e-6);
    Planet Sun = Planet(0.0, 0.0, 0.0, -0.003332, 1.0);  // Jupiter: -0.0026585, All: -0.003332
    Planet Moon = Planet(0.995, 0.0, 0.0, 2.1*pi, 3.69e-8);
    Planet Venus = Planet(0.727, 0.0, 0.0, 7.385, 2.45e-6);
    Planet Mercury = Planet(0.3548, 0.0, 0.0, 10.107, 1.66e-7);
    Planet Mars = Planet(1.642, 0.0, 0.0, 5.085, 3.23e-7);
    Planet Jupiter = Planet(5.168, 0.0, 0.0, 2.764, 9.55e-4); // 10: 9.55e-3, 1000: 9.55e-1
    Planet Saturn = Planet(9.883, 0.0, 0.0, 2.047, 2.86e-4);
    Planet Uranus = Planet(20.06, 0.0, 0.0, 1.435, 4.36e-5);
    Planet Neptune = Planet(30.04, 0.0, 0.0, 1.139, 5.17e-5);
    Planet Pluto = Planet(52.36, 0.0, 0.0, 0.992, 6.61e-9);
    Planet B[n];

    switch(k) {

    case 1:
        B[0] = Earth;
        B[1] = Sun;
        break;

    case 2:
        B[0] = Earth;
        B[1] = Sun;
        B[2] = Moon;
        break;

    case 3:
        B[0] = Earth;
        B[1] = Sun;
        B[2] = Venus;
        B[3] = Mercury;
        B[4] = Mars;
        B[5] = Jupiter;
        B[6] = Saturn;
        B[7] = Uranus;
        B[8] = Neptune;
        B[9] = Pluto;
        break;

    case 4:
        B[0] = Earth;
        B[1] = Sun;
        B[2] = Jupiter;
        break;

    }

    // Initial conditions
    for (int i=0; i < n; i++){

        A[4*i] = B[i].x0;     // xi
        A[4*i+1] = B[i].y0;   // yi
        A[4*i+2] = B[i].vx0;  // vx
        A[4*i+3] = B[i].vy0;  // vy
        M[i] = B[i].M;
        //cout << M[i] << endl;
    }


    double K_1[4*n];
    double K_2[4*n];
    double K_3[4*n];
    double K_4[4*n];
    double K_1_new[4*n];
    double K_2_new[4*n];
    double K_3_new[4*n];
    //double* R = new double;

    fstream newFile;
    newFile.open("data.dat", ios::out);
    newFile << n << endl;
    newFile << t_max << endl;
    newFile << k << endl;

    for (int t=0; t < N; t++) {

        for(int i=0; i < n; i++){
            newFile << A[4*i] << " " << A[4*i+1] << " ";

        }
        newFile << endl;

        f(A, K_1, M, m, n);

        for(int i=0; i < m; i++){
            K_1_new[i] = A[i] + K_1[i]*(h/2);
        }
        f(K_1_new, K_2, M, m, n);

        for(int i=0; i < m; i++){
            K_2_new[i] = A[i] + K_2[i]*(h/2);
        }
        f(K_2_new, K_3, M, m, n);

        for(int i=0; i < m; i++){
            K_3_new[i] = A[i] + K_2[i]*h;
        }
        f(K_3_new, K_4, M, m, n);

        for(int i=0; i < m; i++){
            A[i] = A[i] + (h/6)*(K_1[i] + 2*K_2[i] + 2*K_3[i] + K_4[i]);


        }
    }

    newFile.close();
    return 0;

}

void f(double* A, double* dA, double* M, int m, int n) {

    // Defining constants
    double pi;
    pi = atan(1)*4;
    double G = 4*pi*pi; // 6.67428e-11;
    double r, L;
    double x, y, L_x, L_y;

    for(int i=0; i < m; i++){

        dA[i] = 0;
    }

    //ofstream outFile;
    //outFile.open("distance_from_sun_1.dat", ios::app);

    for(int i=0; i < n; i++){

        dA[4*i] = A[4*i+2];    // dxi/dt = vx_old
        dA[4*i+1] = A[4*i+3];  // dyi/dt = vy_old


        for (int j=0; j < n; j++){
            if(j!=i){

                L_x = A[4*i];
                L_y = A[4*i+1];
                x = A[4*i] - A[4*j];
                y = A[4*i+1] - A[4*j+1];
                r = sqrt(x*x + y*y);
                L = sqrt(L_x*L_x + L_y*L_y);
                //outFile << L << endl;
                //cout << r << endl;
                dA[4*i+2] += -(G*M[j]/(r*r*r))*(x);     // vx_new
                dA[4*i+3] += -(G*M[j]/(r*r*r))*(y);     // vy_new
                //cout << dA[4*i+2] << endl;
                //cout << M[j] << endl;
            }
        }
    }

    //outFile.close();

}




