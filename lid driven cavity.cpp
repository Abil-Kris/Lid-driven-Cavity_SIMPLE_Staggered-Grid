#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

void SIMPLE(double** u_vel, double** v_vel, double** p_1, int n1, double nu_f, double rho_f, double d1) {

    double* x_grid = new double [n1+2];
    double* y_grid = new double [n1+2];

    double* x_face = new double [n1+1];
    double* y_face = new double [n1+1];
    
    for(int i=1; i < n1+1; ++i) {
        x_grid [i] = (d1/(2.0 * n1)) + (d1 / n1) * (i-1);
        y_grid [i] = (d1/(2.0 * n1)) + (d1 / n1) * (i-1);
    }
    x_grid[0] = 0.0;
    y_grid[0] = 0.0;
    x_grid[n1+1] = d1;
    y_grid[n1+1] = d1;

    for(int i=0; i < n1 + 1; ++i) {
        x_face [i] = (d1/n1) * i;
        y_face [i] = (d1/n1) * i;
    }

    double delta = d1 / n1;
    double residue = 1.0;

    double de, dn, ds, dw;
    double fe, fn, fs, fw;
    double b, alpha = 0.5;
    int m = (n1 + 1) / 2;
    int iter = 0;

    //initialization

    vector<vector<double>> ae_u(n1 + 2, vector<double>(n1 + 1));
    vector<vector<double>> an_u(n1 + 2, vector<double>(n1 + 1));
    vector<vector<double>> as_u(n1 + 2, vector<double>(n1 + 1));
    vector<vector<double>> aw_u(n1 + 2, vector<double>(n1 + 1));
    vector<vector<double>> ap_u(n1 + 2, vector<double>(n1 + 1));
    vector<vector<double>> u_vel_c(n1 + 2, vector<double>(n1 + 1));
    
    vector<vector<double>> ae_v(n1 + 1, vector<double>(n1 + 2));
    vector<vector<double>> an_v(n1 + 1, vector<double>(n1 + 2));
    vector<vector<double>> as_v(n1 + 1, vector<double>(n1 + 2));
    vector<vector<double>> aw_v(n1 + 1, vector<double>(n1 + 2));
    vector<vector<double>> ap_v(n1 + 1, vector<double>(n1 + 2));
    vector<vector<double>> v_vel_c(n1 + 1, vector<double>(n1 + 2));

    vector<vector<double>> ae_p(n1, vector<double>(n1));
    vector<vector<double>> an_p(n1, vector<double>(n1));
    vector<vector<double>> as_p(n1, vector<double>(n1));
    vector<vector<double>> aw_p(n1, vector<double>(n1));
    vector<vector<double>> ap_p(n1, vector<double>(n1));
    vector<vector<double>> p_c (n1, vector<double>(n1));

    for (int i=0;i < n1+2; ++i) 
        for (int j =0;j < n1+1; ++j) 
            ap_u[i][j] = 0.0;
    
    for (int i=0;i < n1+1; ++i)
        for (int j =0;j < n1+2; ++j) 
            ap_v[i][j] = 0.0;

    //SIMPLE starts here
    
    while (residue > 0.001 && iter < 100000) {
        
        iter += 1;
        if (iter % 100 == 0)
            cout << "Running Iteration: " << iter << " Residue: " << residue << endl;
        residue = 0.0;
    
        //u velocity

        for ( int i =1; i < n1+1; ++i) {
            for(int j =1; j < n1; ++j) {

            de = nu_f * delta / ( x_face[j+1] - x_face[j] );
            dw = nu_f * delta / ( x_face[j] - x_face[j-1] );
            dn = nu_f * delta / ( y_grid[i+1] - y_grid[i] );
            ds = nu_f * delta / ( y_grid[i] - y_grid[i-1] );

            fe = rho_f * delta * ( u_vel[i][j+1] + u_vel[i][j] ) * 0.5;
            fw = rho_f * delta * ( u_vel[i][j] + u_vel[i][j-1] ) * 0.5;
            fn = rho_f * delta * ( v_vel[i][j+1] + v_vel[i][j] ) * 0.5;
            fs = rho_f * delta * ( v_vel[i-1][j+1] + v_vel[i-1][j] ) * 0.5;

            ae_u[i][j] = de - 0.5 * fe;
            aw_u[i][j] = dw + 0.5 * fw;
            an_u[i][j] = dn - 0.5 * fn;
            as_u[i][j] = ds + 0.5 * fs;
            ap_u[i][j] = ae_u[i][j] + aw_u[i][j] + as_u[i][j] + an_u[i][j] + (fe - fw + fn - fs); // Source term is zero here

            b = ( p_1[i-1][j-1] - p_1[i-1][j] ) * delta;
            u_vel_c[i][j] = ( ae_u[i][j] * u_vel[i][j+1] + aw_u[i][j] * u_vel[i][j-1] + as_u[i][j] * u_vel[i-1][j] + an_u[i][j] * u_vel[i+1][j] + b) / ap_u[i][j];

            }
        }

        //v velocity

        for ( int i =1; i < n1; ++i) {
            for ( int j =1; j <n1+1; ++j) {

            de = nu_f * delta / ( x_grid[j+1] - x_grid[j] );
            dw = nu_f * delta / ( x_grid[j] - x_grid[j-1] );
            dn = nu_f * delta / ( y_face[i+1] - y_face[i] );
            ds = nu_f * delta / ( y_face[i] - y_face[i-1] );

            fe = rho_f * delta * ( u_vel[i+1][j] + u_vel[i][j] ) * 0.5;
            fw = rho_f * delta * ( u_vel[i+1][j-1] + u_vel[i][j-1] ) * 0.5;
            fn = rho_f * delta * ( v_vel[i+1][j] + v_vel[i][j] ) * 0.5;
            fs = rho_f * delta * ( v_vel[i-1][j] + v_vel[i][j] ) * 0.5;
            
            ae_v[i][j] = de - 0.5 * fe;
            aw_v[i][j] = dw + 0.5 * fw;
            an_v[i][j] = dn - 0.5 * fn;
            as_v[i][j] = ds + 0.5 * fs;
            ap_v[i][j] = ae_v[i][j] + aw_v[i][j] + as_v[i][j] + an_v[i][j] + (fe - fw + fn - fs); // Source term is zero here

            b = ( p_1[i-1][j-1] - p_1[i][j-1] ) * delta;
            v_vel_c[i][j] = ( ae_v[i][j] * v_vel[i][j+1] + aw_v[i][j] * v_vel[i][j-1] + as_v[i][j] * v_vel[i-1][j] + an_v[i][j] * v_vel[i+1][j] + b) / ap_v[i][j];

            }
        }

        //pressure correction
        for (int i=0;i < n1; ++i) {   // initializing p_c to zero
            for (int j =0;j < n1; ++j) {
                p_c[i][j] = 0.0;
            }
        }

        for (int i=0;i < n1; ++i) {
            for (int j =0;j < n1; ++j) {
                if (i==0)
                    as_p[i][j] = 0.0;
                else {
                    as_p[i][j] = rho_f * (delta / ap_v[i][j+1]) * delta;
                    p_c[i][j] = p_c[i][j] + as_p[i][j] * p_c[i-1][j];
                }
                if (j==0)
                    aw_p[i][j] = 0.0;
                else {
                    aw_p[i][j] = rho_f * (delta / ap_u[i+1][j]) * delta;
                    p_c[i][j] = p_c[i][j] + aw_p[i][j] * p_c[i][j-1];
                }
                if (i==n1-1)
                    an_p[i][j] = 0.0;
                else{
                    an_p[i][j] = rho_f * (delta / ap_v[i+1][j+1]) * delta;
                    p_c[i][j] = p_c[i][j] + an_p[i][j] * p_c[i+1][j];
                }
                if (j==n1-1) 
                    ae_p[i][j] = 0.0;
                else {
                    ae_p[i][j] = rho_f * (delta / ap_u[i+1][j+1]) * delta;
                    p_c[i][j] = p_c[i][j] + ae_p[i][j] * p_c[i][j+1];
                }
                b = rho_f * delta * (u_vel_c[i+1][j] - u_vel_c[i+1][j+1] + v_vel_c[i][j+1] - v_vel_c[i+1][j+1]);
                residue += abs(b);
                p_c[i][j] = p_c[i][j] + b;
                ap_p[i][j] = ae_p[i][j] + aw_p[i][j] + as_p[i][j] + an_p[i][j]; 

                if ( i == m && j == m ) {   // reference pressure point
                    p_1[i][j] = 0.0;
                    p_c[i][j] = 0.0;
                    continue;
                }

                if (ap_p[i][j] == 0.0)
                    cout << "division by zero in "<<i<<"  "<<j;
                p_c[i][j] = p_c[i][j] / ap_p[i][j];
            }
        }

        //pressure, velocity update

        for (int i=0; i<n1; ++i) {
            for (int j=0; j<n1; ++j) {
                p_1[i][j] = p_1[i][j] + alpha * (p_c[i][j]);
            }
        }

        for (int i=1; i<n1+1; ++i) {
            for (int j=1; j<n1; ++j) {
                u_vel[i][j] = u_vel_c[i][j] + alpha * (delta / ap_u[i][j]) * ( p_c[i-1][j-1] - p_c[i-1][j] );
            }
        }

        for (int i=1; i<n1; ++i) {
            for (int j=1; j<n1+1; ++j) {
                v_vel[i][j] = v_vel_c[i][j] + alpha * (delta / ap_v[i][j]) * ( p_c[i-1][j-1] - p_c[i][j-1] );
            }
        }

    }
    cout << "Iteration End: " << iter << " Residue: " << residue << endl;

}

void result (double** u_vel, double** v_vel, double** p_1, int n1) {

    double u_visual[n1][n1], v_visual[n1][n1];

    ofstream fout_u1("u_vel.txt");
    for (int i = 0; i < n1 + 2; ++i) {
        for (int j = 0; j < n1 + 1; ++j) {
            fout_u1 << u_vel[i][j] << " ";
        }
        fout_u1 << "\n";
    }
    fout_u1.close();

    // Write v_vel to file
    ofstream fout_v1("v_vel.txt");
    for (int i = 0; i < n1 + 1; ++i) {
        for (int j = 0; j < n1 + 2; ++j) {
            fout_v1 << v_vel[i][j] << " ";
        }
        fout_v1 << "\n";
    }
    fout_v1.close();

    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            u_visual[i][j] = 0.5 * (u_vel[i+1][j] + u_vel[i+1][j+1]);
            v_visual[i][j] = 0.5 * (v_vel[i][j+1] + v_vel[i+1][j+1]);
        }
    }

    ofstream fout_u2("u_visual.txt");
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            fout_u2 << u_visual[i][j] << " ";
        }
        fout_u2 << "\n";
    }
    fout_u2.close();

    // Write v_vel to file
    ofstream fout_v2("v_visual.txt");
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            fout_v2 << v_visual[i][j] << " ";
        }
        fout_v2 << "\n";
    }
    fout_v2.close();

    // Write p_1 to file
    ofstream fout_p("p_1.txt");
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            fout_p << p_1[i][j] << " ";
        }
        fout_p << "\n";
    }
    fout_p.close();
}

int main() {
    int n1;
    double nu_f, rho_f, d1;
    cout << "Enter the number of grids\n";
    cin >> n1;
    cout << "Enter the dynamic visc and density\n";
    cin >> nu_f >> rho_f;
    cout << "Enter the length\n";
    cin >> d1;  

    double** u_vel = new double*[n1+2];  // n1+2 * n1+1
    double** v_vel = new double*[n1+1];  // n1+1 * n1+2
    double** p_1 = new double*[n1];   //n1 * n1

    for (int i = 0; i < n1+2; i++) {
        u_vel[i] = new double[n1+1];
    }

    for (int i = 0; i < n1+1; i++) {
        v_vel[i] = new double[n1+2];
    }

     for (int i = 0; i < n1; i++) {
        p_1[i] = new double[n1];
    }

    for ( int i = 0; i < n1+2; ++i) {
        for(int j = 0; j < n1+1; ++j){
            u_vel[i][j] = 0.0;
            if (i == n1+1)
                u_vel[i][j] = 1.0;      // Top lid moving with velocity 1.0
        }
    }

    for ( int i = 0; i < n1+1; ++i)
        for(int j = 0; j < n1+2; ++j)
            v_vel[i][j] = 0.0;
     
    for ( int i =0; i < n1; ++i)
        for(int j =0; j < n1; ++j)
            p_1[i][j] = 0.01;

    SIMPLE(u_vel, v_vel, p_1, n1, nu_f, rho_f, d1);
    result (u_vel, v_vel, p_1, n1);
    return 0;
}
    

    
 