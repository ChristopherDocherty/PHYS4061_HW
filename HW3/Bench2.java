// A program to solve the problem of a person sitting
// on a bench with the relaxation scheme.
import java.lang.*;

public class Bench2 {

final static int n = 100, m = 2;

public static void main(String argv[]) {

    double u[] = new double[n+1];
    double d[] = new double[n+1];
    double s[] = new double[n+1];
    
    double l = 3, l2 = l/2, h = l/n, h2 = h*h;
    double x0 = 0.25, x2 = x0*x0, e0 = 1/Math.E;
    double x = 0, rho = 3, g = 9.8, f0 = 200;
    double y = 1e9*Math.pow(0.03,3)*0.2/3;
    double u0 = 0.032, del = 1e-7;
    int nmax = 10000;

    for (double p = 0.1; p <= 1.0; p+=0.1 ){

       
        // Evaluate the source in the equation
        for (int i=0; i<=n; ++i) {
            s[i] = rho*g;
            x = h*i-l2;

            if (Math.abs(x) < x0)
                s[i] += f0*(Math.exp(-x*x/x2)-e0);

            s[i] *= h2/y;
        }

        for (int i=1; i<n; ++i) {

            x = Math.PI*h*i/l;
            u[i] = u0*Math.sin(x);
            d[i] = 1;
        }

        d[0] = d[n] = 1;

        relax(u, d, s, p, del, nmax);

        // Output the result in every m time step
        x = 0;

            //Commented out to only leave information about covergence speed
            /*
            double mh = m*h;

            for (int i=0; i<n; i+=m) {
                System.out.println(x + " " + 100*u[i]);
                x += mh;
            }
            */
            
        }
}

    // Method to complete one step of relaxation.
    public static void relax(double u[], double d[],
    double s[], double p, double del, int nmax) {

        int n = u.length-1;
        double q = 1-p, fi = 0;
        double du = 2*del;
        int k = 0;

        while ((du>del) && (k<nmax)) {

            du = 0;

            for (int i=1; i<n; ++i) {
                fi = u[i];
                u[i] = p * u[i]
                    + q * ((d[i+1] + d[i]) * u[i+1] + (d[i] + d[i-1]) * u[i-1] + 2 * s[i]) / (4 * d[i]);
                fi = u[i]-fi;
                du += fi*fi;
            }
            du = Math.sqrt(du/n);
            k++;
        }

        System.out.println("Convergence" +
        " found after " + k + " iterations");

        if (k==nmax) System.out.println("Convergence not" +
        " found after " + nmax + " iterations");
    }
}