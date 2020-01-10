#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>


using namespace std;
const int n = 100;
int j = 2;

class Jump{


    public:


    double x[n+1];
    double y[n+1];
    double vx[n+1];
    double vy[n+1];
    double ax[n+1];
    double ay[n+1];

    //Assign all the constants involved

    double angle = 42.5*M_PI/180.0; //Only for constant
    double mass = 250; //Only for constant
    double area = 0.93; // ONly for constant
    double density = 1.2; //Only for constant
    double g = 9.80;
    double speed = 67;
    double k = area*density/ (2*mass);
    double dt = 2*speed*sin(angle)/(g*n);
    double d = dt*dt/2;


    double v = 0;
    double p = 0;

    double d2 = 0;
    double d3 = 0;

    //Assign the quantities for the first two points 
    void Initialise_values();

    void run_predictions();

    void outputResults();


};


void Jump::Initialise_values(){
    x[0] = y[0] = 0;
    vx[0] = speed*cos(angle);
    vy[0] = speed*sin(angle);

    v = sqrt(vx[0]*vx[0] + vy[0]*vy[0]);

    ax[0] = -k* pow(v,1.0/3) *vx[0];
    ay[0] = -g-k* pow(v,1.0/3) *vy[0];
    p = vx[0]*ax[0]+vy[0]*ay[0]; //dot product
    
    x[1] = x[0]+dt*vx[0]+d*ax[0];
    y[1] = y[0]+dt*vy[0]+d*ay[0];

    vx[1] = vx[0] + dt*ax[0] -d*k* (v*ax[0] + p*vx[0] / (3 * pow(v,1.0/3)));
    vy[1] = vy[0] + dt*ay[0] -d*k* (v*ay[0] + p*vy[0] / (3 * pow(v,1.0/3)));


    v = sqrt(vx[1]*vx[1]+vy[1]*vy[1]);

    ax[1] = -k*v*vx[1];
    ay[1] = -g-k*v*vy[1];
}

void Jump::run_predictions(){


    // Calculate other position and velocity recursively
    d2 = 2*dt;
    d3 = dt/3;

    for (int i=0; i<n-1; ++i) {

        // Predict the next position and velocity
        x[i+2] = x[i] + d2 * vx[i+1];
        y[i+2] = y[i] + d2 * vy[i+1];
        vx[i+2] = vx[i] + d2 * ax[i+1];
        vy[i+2] = vy[i] + d2 * ay[i+1];
        v = sqrt( vx[i+2] * vx[i+2] + vy[i+2] * vy[i+2]);
        ax[i+2] = -k * pow(v,1.0/3) * vx[i+2];
        ay[i+2] = -g - k * pow(v,1.0/3) * vy[i+2];

        // Correct the new position and velocity
        x[i+2] = x[i]+d3*(vx[i+2]+4*vx[i+1]+vx[i]);
        y[i+2] = y[i]+d3*(vy[i+2]+4*vy[i+1]+vy[i]);
        vx[i+2] = vx[i]+d3*(ax[i+2]+4*ax[i+1]+ax[i]);
        vy[i+2] = vy[i]+d3*(ay[i+2]+4*ay[i+1]+ay[i]);
    }
}


void Jump::outputResults(){
    for(int i=0; i != n; i+=j){
        cout << x[i] << "\t" << y[i] << endl;
    }
}





int main(){


    Jump test;

    test.Initialise_values();
    test.run_predictions();
    test.outputResults();



    return 0;
}

