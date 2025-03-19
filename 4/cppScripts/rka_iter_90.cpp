#include <iostream>
#include <Eigen>
#include <cfloat>
#include <cmath>
 
using namespace std;
using namespace Eigen;

Matrix3f f(const Vector3f& r_V, double sepX, double sepY) {
    Eigen::Matrix3f e01;
    Eigen::Matrix3f e02;
    Eigen::Matrix3f e03;
    Eigen::Matrix3f e04;
    Eigen::Matrix3f e11;
    Eigen::Matrix3f e12;
    Eigen::Matrix3f e13;
    Eigen::Matrix3f e14;    

    double x00 = r_V(0) + sepX;
    double x01 = r_V(1) + sepY;
    double x02 = r_V(2);
    
    double x10 = r_V(0) - sepX;
    double x11 = r_V(1) - sepY;
    double x12 = r_V(2);
    double r02 = (x00 * x00 + x01 * x01);
    double r12 = (x10 * x10 + x11 * x11);
    
    double aI = 1.0;
    double IpqXpXq0 = ((x00 * x00 - x01 * x01));


    e01 << 1.0, 0.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, 0.0, 0.0;
    
    e02 << (x00 * x00), 0.0, 0.0,
          0.0, -1.0 * (x01 * x01), 0.0,
          0.0, 0.0, 0.0;

    e03 << IpqXpXq0, 0.0, 0.0,
          0.0, IpqXpXq0, 0.0,
          0.0, 0.0, IpqXpXq0;

    e04 << x00 * x00 * IpqXpXq0, x00 * x01 * IpqXpXq0, 0.0,
           x00 * x01 * IpqXpXq0, x01 * x01 * IpqXpXq0, 0.0,
           0.0, 0.0, 0.0;
    

    double IpqXpXq1 = (-1.0 * x10*x10) + (x11 * x11);
    e11 << -1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 0.0;

    e12 << -1.0 * (x10 * x10), 0.0, 0.0,
            0.0, (x11 * x11), 0.0,
            0.0, 0.0, 0.0;

    // Const
    e13 << IpqXpXq1, 0.0, 0.0,
            0.0, IpqXpXq1, 0.0,
            0.0, 0.0, IpqXpXq1;

    e14 << x10 * x10 * IpqXpXq1, x10 * x11 * IpqXpXq1, 0.0,
            x10 * x11 * IpqXpXq1, x11 * x11 * IpqXpXq1, 0.0,
            0.0, 0.0, 0.0;
    

    return ((-6 * e01) + (60 * e02 + 15 * e03) / (r02) + (-105 * e04 / (r02*r02))) + 
           ((-6 * e11) + (60 * e12 + 15 * e13) / (r12) + (-105 * e14 / (r12*r12)));

    
}

double eigen_solve_val(Matrix3f E_temp, int icity, Vector3f r) {
    // Compute eigenvalues and eigenvectors
    EigenSolver<Matrix3f> solver(E_temp);
    double x = r(0);
    
    Vector3cf eigenvalues = solver.eigenvalues();   // Complex eigenvalues

    double maxEigenvalue = -1e6; // Large negative value for initialization
    int maxIndex = -1;

    if (x < 0) {
        icity *= -1.0;
    }

    // Loop through eigenvalues to find the largest with the same sign as `icity`
    for (int i = 0; i < 3; ++i) {
        double realVal = eigenvalues[i].real(); // Extract real part of eigenvalue
        
        if (icity * realVal > 0 && realVal * icity > maxEigenvalue) {
            maxEigenvalue = realVal * icity * -1;
            maxIndex = i;
        
    }
    }
    double result;
    if (maxIndex < 0) {
        result = 0.0;
    } else {
        result = eigenvalues(maxIndex).real();
    }

    //std::cout << solver.eigenvectors() << "\n";
    //std::cout << eigenvalues << "\n";
    //std::cout << maxIndex << "\n";
    //std::cout << "\n";

    return result;
    }

// Function to get the eigenvector corresponding to the largest eigenvalue with the given sign
Vector3f eigen_solve(Matrix3f E_temp, int icity, Vector3f r) {
    // Compute eigenvalues and eigenvectors
    EigenSolver<Matrix3f> solver(E_temp);
    Vector3cf eigenvalues = solver.eigenvalues();   // Complex eigenvalues
    Matrix3cf eigenvectors = solver.eigenvectors(); // Corresponding eigenvectors
    Vector3f result;
    double x = r(0);

    double maxEigenvalue = -1e6; // Large negative value for initialization
    int maxIndex = -1;
    if (x < 0) {
        icity *= -1.0;
    }

    // Loop through eigenvalues to find the largest with the same sign as `icity`
    for (int i = 0; i < 3; ++i) {
        double realVal = eigenvalues[i].real(); // Extract real part of eigenvalue
        
        if (icity * realVal > 0 && realVal * icity > maxEigenvalue) {
            maxEigenvalue = realVal * icity * -1;
            maxIndex = i;
           
    }
    }
    // Convert complex eigenvector to real vector (assuming it's real-valued)
    if (maxIndex < 0) {
        result << 0.0,0.0,0.0;
        return result;
    } 
    
    result = eigenvectors.col(maxIndex).real();
    
    return result / result.norm();
}


extern "C" {

struct vect {
    double x[10000];
    double y[10000];
    double z[10000];
    double m[10000];
    int its;
};

vect rka_iter_double(double sepX, double sepY, double seed_x, double seed_y, double seed_z, int num_its, int icity, double ending_tolerance, double delta_0, double safety, double h0) {    
    // Butcher tableau for RK45
    static const double a2 = 1.0 / 5.0, a3 = 3.0 / 10.0, a4 = 3.0 / 5.0, a5 = 1.0, a6 = 7.0 / 8.0;
    static const double b21 = 1.0 / 5.0;
    static const double b31 = 1.0 / 40.0, b32 = 9.0 / 40.0;
    static const double b41 = 3.0 / 10.0, b42 = -9.0 / 10.0, b43 = 6.0 / 5.0;
    static const double b51 = -11.0 / 54.0, b52 = -5.0/2.0, b53 = -70.0 / 27.0, b54 = -35.0 / 27.0;
    static const double b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0;
    static const double c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 621.0, c6 = 512.0 / 1771.0;
    static const double cd1 = 37.0 / 378.0 - 2825.0 / 27648.0, cd3 = 250.0 / 621.0 - 18575.0 / 48384.0, cd4 = 125.0 / 621.0 - 13525.0 / 55296.0, cd5 = -277.0 / 14336.0, cd6 = 512.0 / 1771.0 - 1.0 / 4.0;
    bool sizing = true;
    vect r_change_vect;

    Vector3f r = {seed_x, seed_y, seed_z};
    Vector3f r_past = {0.0, 0.0, 0.0};
    Vector3f r_change;
    int i;
    double h = h0;
    double deltNorm;
    for (i = 0; i < num_its; i++) {
        sizing = true;

        while (sizing) {// Butcher tableau for RK45
            // Compute Runge-Kutta stages
            Matrix3f E1 = f(r, sepX, sepY);
            Vector3f k1 = h * eigen_solve(E1, icity, r);

            Matrix3f E2 = f(r + b21 * k1, sepX, sepY);
            Vector3f k2 = h * eigen_solve(E2, icity, r);

            Matrix3f E3 = f(r + b31 * k1 + b32 * k2, sepX, sepY);
            Vector3f k3 = h * eigen_solve(E3, icity, r);
            
            Matrix3f E4 = f(r + b41 * k1 + b42 * k2 + b43 * k3, sepX, sepY);
            Vector3f k4 = h * eigen_solve(E4, icity, r);

            Matrix3f E5 = f(r + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4, sepX, sepY);
            Vector3f k5 = h * eigen_solve(E5, icity, r);
            
            Matrix3f E6 = f(r + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5, sepX, sepY);
            Vector3f k6 = h * eigen_solve(E6, icity, r);

            Vector3f delta = cd1 * k1 + cd3 * k3 + cd4 * k4 + cd5 * k5 + cd6 * k6;

            deltNorm = delta.norm();

            if (deltNorm == 0.0) {
                r_change = delta;
                break;
                
            } else if ((delta.array().abs() > delta_0).any()) {
                h *= safety * std::pow(std::abs(delta_0 / delta.norm()), 0.25f);
            } else {
                h *= safety * std::pow(std::abs(delta_0 / delta.norm()), 0.2f);

                // Assuming the calculations for r_change involve variables like k1, k3, k4, k6
                r_change = c1 * k1 + c3 * k3 + c4 * k4 + c6 * k6;
                
                // Set step_sizing to false
                sizing = false;
            }
        }
        if (deltNorm == 0.0) {
            break;
        }

        // Account for sign missmatch
        double r_mag = r_change.norm();
        double dot = r_past.dot(r_change/r_mag);

        if (dot < -.9) {
            r_change *= -1;
        }

        double val = eigen_solve_val(f(r, sepX, sepY), icity, r);
        r += r_change;
        r_past = r_change / r_mag;
        
        //std::cout << r << "\n";

        r_change_vect.x[i] = r(0);
        r_change_vect.y[i] = r(1);
        r_change_vect.z[i] = r(2);
        r_change_vect.m[i] = val;

        if (abs(r.norm()) < ending_tolerance) {
            break;
        }
    
    }

    r_change_vect.its = i;


    return r_change_vect;
}

};

int main() {
    //rka_iter(1.0, 1.0, 1.0, 2000, 1, .5);

    return 0;
}
