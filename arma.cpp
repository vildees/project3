#include <iostream>
#include <math.h>
#include <armadillo>
#include "planet.h"
#include "constants.h"

using namespace std;
using namespace arma;

void f(vec & A, vec & dAdt, planet * B, int n);
void dvdt(vec & A, vec & dvxdt, vec & dvydt, planet *B, int n);


int main() {

    int i, j;
    int n = 2;                                     // number of planets
    int N = 1000;
    double t_i = 0.0 ;                             // initial time
    double t_f = 15.;                             // final time
    double dt = (t_f-t_i)/N;
    vec k1, k2, k3, k4;

    vec A = zeros(4*n);
    vec dAdt = zeros(4*n);           // dA/dt
    vec C = zeros(4*n);

    planet Earth = planet(1.0, 0.0, 0.0, 2.*pi-2.3, M_earth);
    planet Sun = planet(0.0, 0.0, 0.0, 0.0, M_sun);
    planet B[n];
    B[0] = Earth;
    B[1] = Sun;

    for (i=0; i<n; i++) {
        j = 4*i;
        A(j) = B[i].x0;                       // x - position
        A(j+1) = B[i].y0;                     // y - position
        A(j+2) = B[i].vx;                     // vx - velocity
        A(j+3) = B[i].vy;                     // vy - velocity
    }

    fstream outFile;
    outFile.open("data.dat", ios::out);

    outFile << n << endl;

    for(i=0; i<n; i++)
        outFile << A(4*i) << " " << A(4*i+1) << " ";
    outFile << endl;

    for(i=0; i<N; i++) {
        f(A, dAdt, B, n);
        k1 = dAdt;
        C = A + k1*(dt/2);
        f(C, dAdt, B, n);
        k2 = dAdt;
        C = A + k2*(dt/2);
        f(C, dAdt, B, n);
        k3 = dAdt;
        C = A + k3*dt;
        f(C, dAdt, B, n);
        k4 = dAdt;
        A += (dt/6.)*(k1 + 2*(k2+k3) + k4);


        for(j=0; j<n; j++)
            outFile << A(4*j) << " " << A(4*j+1) << " "; // << endl;
        outFile << endl;
    }

    outFile.close();
}


void dvdt(vec & A, vec & dvxdt, vec & dvydt, planet * B, int n) {

    int i, j;
    double xi, xj, yi, yj;
    double x, y, r;
    double x_dvdt, y_dvdt;
    double M1, M2;



    for (i=0; i<n; i++) {
        xi = A(4*i);                           // x-position
        yi = A(4*i+1);                         // y-position
        M1 = B[i].M;

        for (j=0; j<n; j++) {
            xj = A(4*j);
            yj = A(4*j+1);
            M2 = B[j].M;

            if (i != j) {
                x = xi-xj;
                y = yi-yj;

                r = sqrt(x*x+y*y);
                x_dvdt += -G*M2*x / (r*r*r);
                y_dvdt += (-G*M2*y) / (r*r*r);
            }
        }

        dvxdt(i) = x_dvdt;
        dvydt(i) = y_dvdt;

        x_dvdt = 0;
        y_dvdt = 0;
    }

}

void f(vec & A, vec & dAdt, planet *B, int n) {

    int i, j;

    vec  dvxdt = zeros(n);
    vec  dvydt = zeros(n);
    dvdt(A, dvxdt, dvydt, B, n);

    for(i=0; i<n; i++) {
        j = 4*i;
        dAdt(j) = A(j+2);                     // dx/dt = vx
        dAdt(j+1) = A(j+3);                   // dy/dt = vy
        dAdt(j+2) = dvxdt(i);                  // dvx/dt
        dAdt(j+3) = dvydt(i);               // dvy/dt

    }

}
