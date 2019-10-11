#include <iostream>
#include <iomanip>
#include <string>


using std::cout; using std::endl;
using std::cin; using std::istream;
using std::streamsize; using std::setprecision;
using std::string; using std::fixed;



/* this is the inverse
[[-0.59375    -0.3828125   0.5703125   0.4921875 ]
 [ 0.296875    0.31640625 -0.41015625 -0.12109375]
 [-0.09375    -0.0078125   0.1953125  -0.1328125 ]
 [ 0.625       0.21875    -0.46875    -0.28125   ]]
*/


void fl(double s[4][4][4], double matA[4][4], double c[], size_t);
void matMult(double z[4][4], double x[4][4], double y[4][4], size_t);
double tr(double x[4][4], size_t n);


void printMatrix(double matForPrint[4][4]){
         
    size_t n = 4;
    streamsize prec = cout.precision();
    streamsize desiredPrec = 7;

    cout << setprecision(desiredPrec) << "[";

    for(int i=0; i < n; ++i){
        cout << (i == 0 ? "[" : " [");
        for(int j=0; j < n; ++j)
            cout << fixed <<matForPrint[i][j] << " ";
        cout << "]" << (i==n-1 ? "" : "\n");
    }
    cout << "]"<< endl;
    
}





/*
A function that implements fadeev-laverrier method to put the inverse of matA in the 2-d array passed in as "s".
The size of the matrix must be passed in as "n".
*/
void fl(double s[4][4][4], double matA[4][4], double c[], size_t n){   

    //Initialise s to have 1 on diagonal and zeros elsewhere 
    for (int i = 0; i < n; i++){
        for(int j=0; j < n; j++){
            for(int k=0; k < n; k++){
                if(j == k){
                    s[i][j][k] = 1;
                } else {
                    s[i][j][k] = 0;
                }
    }}}



    for(int i = 1; i < n; ++i){
       matMult(s[i],matA,s[i-1],n);
       c[n-i] = -tr(s[i],n)/i;

       for(int j=0; j !=n; ++j)
        s[i][j][j] += c[n-i];
    }

    //Create a temporary array to put matrix multiplication result in
    double tempMat[4][4];
    matMult(tempMat, matA, s[n-1], n);
    //Calculate c as required by Fadeev-Leverrier method
    c[0] = -tr(tempMat,n)/n;
}


/*
Matrix multiplication function: pass in any two matrices "x" & "y" and the function will put the result of x*y into "z"
Must also supply the dimensions of the matrices as "n".
Will only work for square matrices as this program intends to find inverses.
*/
void matMult(double z[4][4], double x[4][4], double y[4][4], size_t n){

    //Initialise the matrix to get rid of random values
    for (int i = 0; i < n; i++)
        for(int j=0; j < n; j++)
            z[i][j] = 0;

    //Implementation of matrix multiplication looping over all entries
    for (int i = 0; i < n; i++)
        for(int j=0; j < n; j++)
            for(int k=0; k < n; k++)
            z[i][j] += x[i][k] * y[k][j];
}


/*
Calculates the trace of a matrix "x", i.e. the sum of its diagonal entries.
Must also pass in the size of the array as "n"
*/
double tr(double x[4][4], size_t n){

    double sum = 0;

    for(int i = 0; i < n; i++)
        sum += x[i][i];
    return sum;
}




int main(){

    //Set the size of the matrix here
    const size_t n = 4;

    //Matrix for inverse entered here
    double matA[n][n] = { 
        {5,2,4,6},
        {5,8,7,2},
        {6,4,8,5},
        {5,4,1,3}
    };
    


    //Create arrays to hold relevant matrices for Fadeev-Leverrier
    double d[n][n];
    double c[n];
    double s[n][n][n];

    //Call Fadeev-Leverrier method
    fl(s,matA,c,n);

    //Populate the inverse matrix elementwise
    for(int i=0; i <n; ++i)
        for(int j=0; j < n; ++j)
            d[i][j] = -s[n-1][i][j]/c[0];
    
    
    printMatrix(matA);
    printMatrix(d);



    cout << "End of program press any key to finish";
    cin.get();



    return 0;
}






