#include <iostream>
#include <iomanip>


using std::cout; using std::endl;


/* this is the inverse
[[-0.59375    -0.3828125   0.5703125   0.4921875 ]
 [ 0.296875    0.31640625 -0.41015625 -0.12109375]
 [-0.09375    -0.0078125   0.1953125  -0.1328125 ]
 [ 0.625       0.21875    -0.46875    -0.28125   ]]
*/




int main(){


    const size_t n = 4;
    double matA[n][n] = { 
        {5,2,4,6},
        {5,8,7,2},
        {6,4,8,5},
        {5,4,1,3}
    };

    

      
    double c[n];
    double d[n][n];

    double s[n][n][n];






    return 0;
}




void fl(double s[4][4][4], double matA[4][4], double c[], size_t n){   

    
    for(int i=0; i != n; ++i)
        s[0][i][i] = 1;

    for(int i = 1; i != n; ++i){
       matMult(s[i],matA,s[i-1],n);
       c[n-i] = -tr(s[i],n)/i;

       for(int j=0; j !=n; ++j)
        s[i][j][j] += c[n-i];
    }

    double tempMat[4][4];
    matMult(tempMat, matA, s[n-1], n);
    c[0] = -tr(tempMat,n)/n;
}



void matMult(double z[4][4], double x[4][4], double y[4][4], size_t n){

    for (int i = 0; i != n; i++)
        for(int j=0; j != n; j++)
            for(int k=0; k != n; k++)
            z[i][j] += x[i][k] * y[k][j];
}



double tr(double x[4][4], size_t n){
    double sum = 0;
    for(int i = 0; i!= n; i++)
        sum += x[i][i];
}



