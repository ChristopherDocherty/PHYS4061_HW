#include <iostream>
#include <iomanip>
#include <string>


using std::cout; using std::endl;
using std::cin; using std::istream;
using std::streamsize; using std::setprecision;
using std::string; using std::fixed;




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
    cout << "] \n" << setprecision(prec) <<endl;

    
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

    //The matrices used for testing inverse function
    double matA[n][n] = { 
        {5,2,4,6},
        {5,8,7,2},
        {6,4,8,5},
        {5,4,1,3}
    };

    double matB[n][n] = {
        {55,23,1,45},
        {5,65,14,26},
        {25,14,88,2},
        {25,25,18,99}
    };

    //Inverse matrices as calculated in numpy - See extra python file for proof of results

    string matAinverseText =
        "[[-0.59375    -0.3828125   0.5703125   0.4921875 ]"
         "\n [ 0.296875    0.31640625 -0.41015625 -0.12109375]"
         "\n [-0.09375    -0.0078125   0.1953125  -0.1328125 ]"
         "\n [ 0.625       0.21875    -0.46875    -0.28125   ]]";
    
    string matBinverseText = 
        "[[ 0.02147166 -0.00479942  0.00226743 -0.0085452 ]"
         "\n [ 0.00156592  0.01719185 -0.00169073 -0.00519266]"
         "\n [-0.0062426  -0.00130587  0.01103738  0.00295753]"
         "\n [-0.00468255 -0.00289197 -0.00215243  0.01303243]]";





    //Create arrays to hold relevant matrices for Fadeev-Leverrier
    double matAinverse[n][n];
    double matBinverse[n][n];
    double c[n];
    double s[n][n][n];

    //Call Fadeev-Leverrier method
    fl(s,matA,c,n);

    //Populate the inverse matrix elementwise
    for(int i=0; i <n; ++i)
        for(int j=0; j < n; ++j)
            matAinverse[i][j] = -s[n-1][i][j]/c[0];

    //Repeated for the second matrix
    fl(s,matB,c,n);

    for(int i=0; i <n; ++i)
        for(int j=0; j < n; ++j)
            matBinverse[i][j] = -s[n-1][i][j]/c[0];
    
    

    //printing the results for comparison

    cout << "For the matrix A given by:" << endl;
    printMatrix(matA);
    cout << "The expected inverse is: \n" << matAinverseText << endl;
    cout << "\n The inverse was found to be: " << endl;
    printMatrix(matAinverse);


    cout << "\n" << "For the matrix B given by: " << endl;
    printMatrix(matB);
    cout << "The expected inverse is: \n" << matBinverseText << endl;
    cout << "\n The inverse was found to be: " << endl;
    printMatrix(matBinverse);





    cout << " \n End of program press any key to terminate";
    cin.get();



    return 0;
}