#include <iostream>
#include <Eigen>
#include <cfloat>
#include <cmath>
 
using namespace std;
using namespace Eigen;

// Define hyperbolic functions
double sech(double x) {
    return 1.0 / cosh(x);
}
double coth(double x) {
    return 1.0 / tanh(x);
}

// Function to compute the matrix
Matrix3f f(const Vector3f& r_V, double R, double S, double vX) {
    Eigen::Matrix3f B;

    // Precompute common
    double x = r_V(0);
    double y = r_V(1);

    double x2 = x * x;
    double y2 = y * y;
    double rS = R * S;
    double rS2 = rS * rS;
    double r = std::sqrt(x2 + y2);
    double r2 = std::pow(r, 2);
    double r3 = std::pow(r, 3);
    double S2 = S*S;
    double cothRS = coth(rS);
    double sechRS = sech(rS);
    double tanhRS = std::tanh(rS);

    double sechRPlus = sech(S * (R + r));
    double sechRMinus = sech(S * (-R + r));
    double sechRPlus2 = sechRPlus * sechRPlus;
    double sechRMinus2 = sechRMinus * sechRMinus;
    double tanhRPlus = std::tanh(S * (R + r));
    double tanhRMinus = std::tanh(S * (-R + r));

    // Matrix components (1st row)
    B(0, 0) =  .25 * cothRS * ((S*x*y * (sechRMinus2 - sechRPlus2) / r3) + (2 * S2 * x * y * ((sechRMinus2 * tanhRMinus) - (sechRPlus2 * tanhRPlus)) / r2));
    B(0, 1) = -.25 * cothRS * ((S*x2 * (sechRMinus2 - sechRPlus2) / r3)  + (2 * S2 * x2    * ((sechRMinus2 * tanhRMinus) - (sechRPlus2 * tanhRPlus)) / r2) + (S * (sechRPlus2 - sechRMinus2) / r));
    B(0, 2) = -1.0 * B(0, 0) + .25 * cothRS * ((S*x2 * (sechRMinus2 - sechRPlus2) / r3) + (S * (sechRPlus2 - sechRMinus2) / r) + (2 * S2 * x2 * (sechRMinus2*tanhRMinus - sechRPlus2*tanhRPlus) / r2));
    // 2nd row
    B(1, 0) =  .25 * cothRS * ((S*y2 * (sechRMinus2 - sechRPlus2) / r3)  + (2 * S2 * y2    * ((sechRMinus2 * tanhRMinus) - (sechRPlus2 * tanhRPlus)) / r2) + (S * (sechRPlus2 - sechRMinus2) / r));
    B(1, 1) = -1.0 * B(0, 0);
    B(1, 2) = -1.0 * B(1, 1) - .25 * cothRS * ((S*y2 * (sechRMinus2 - sechRPlus2) / r3) + S* (sechRPlus2 - sechRMinus2) / r + 2 * S2 * y2 * ((sechRMinus2 * tanhRMinus) - (sechRPlus2 * tanhRPlus)) / r2);
    // 3rd row
    B(2,0) = -.25 * cothRS * (S * (sechRPlus2 - sechRMinus2) / r);
    B(2,1) = -1.0 * B(2,0);
    B(2,2) = 0;

    return B;
}

double eigen_solve_val(Matrix3f E_temp, int icity) {
    // Compute eigenvalues and eigenvectors
    //std::cout << E_temp << "\n";
    EigenSolver<Matrix3f> solver(E_temp);
    
    Vector3cf eigenvalues = solver.eigenvalues();   // Complex eigenvalues

    double maxEigenvalue = -1e6; // Large negative value for initialization
    int maxIndex = -1;

    // Loop through eigenvalues to find the largest with the same sign as `icity`
    for (int i = 0; i < 3; ++i) {
        double realVal = eigenvalues[i].real(); // Extract real part of eigenvalue
        
        if (icity * realVal > 0 && realVal * icity > maxEigenvalue) {
            
            maxEigenvalue = realVal * icity;
            maxIndex = i;
            
        }
    }
    double result;
    if (maxIndex < 0) {
        result = 0.0;
    } else {
        result = eigenvalues(maxIndex).real();
    }

    // Test for trace free
    if ((eigenvalues(0) + eigenvalues(1) + eigenvalues(2)).real() > 0.00001){
        std::cout << eigenvalues(0) + eigenvalues(1) + eigenvalues(2) << '\n';
    }

    return result;
    }

// Function to get the eigenvector corresponding to the largest eigenvalue with the given sign
Vector3f eigen_solve(Matrix3f E_temp, int icity) {
    // Compute eigenvalues and eigenvectors
    EigenSolver<Matrix3f> solver(E_temp);
    Vector3cf eigenvalues = solver.eigenvalues();   // Complex eigenvalues
    Matrix3cf eigenvectors = solver.eigenvectors(); // Corresponding eigenvectors
    Vector3f result;

    double maxEigenvalue = -1e6; // Large negative value for initialization
    int maxIndex = -1;

    // Loop through eigenvalues to find the largest with the same sign as `icity`
    for (int i = 0; i < 3; ++i) {
        double realVal = eigenvalues[i].real(); // Extract real part of eigenvalue
        
        if (icity * realVal > 0 && realVal * icity > maxEigenvalue) {
            maxEigenvalue = realVal * icity;
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

vect rka_iter(double R, double sigma, double vX, double seed_x, double seed_y, double seed_z, int num_its, int icity, double ending_tolerance, double delta_0, double safety, double h0) {    
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
    double h;

    for (i = 0; i < num_its; i++) {
        h = h0;
        sizing = true;

        while (sizing) {
            // Compute Runge-Kutta stages
            Matrix3f E1 = f(r, R, sigma, vX);
            Vector3f k1 = h * eigen_solve(E1, icity);

            Matrix3f E2 = f(r + b21 * k1, R, sigma, vX);
            Vector3f k2 = h * eigen_solve(E2, icity);

            Matrix3f E3 = f(r + b31 * k1 + b32 * k2, R, sigma, vX);
            Vector3f k3 = h * eigen_solve(E3, icity);

            Matrix3f E4 = f(r + b41 * k1 + b42 * k2 + b43 * k3, R, sigma, vX);
            Vector3f k4 = h * eigen_solve(E4, icity);

            Matrix3f E5 = f(r + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4, R, sigma, vX);
            Vector3f k5 = h * eigen_solve(E5, icity);

            Matrix3f E6 = f(r + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5, R, sigma, vX);
            Vector3f k6 = h * eigen_solve(E6, icity);

            Vector3f delta = cd1 * k1 + cd3 * k3 + cd4 * k4 + cd5 * k5 + cd6 * k6;
            
            if ((delta.array().abs() > delta_0).any()) {
                h *= safety * std::pow(std::abs(delta_0 / delta.norm()), 0.25f);
            } else {
                h *= safety * std::pow(std::abs(delta_0 / delta.norm()), 0.2f);

                r_change = c1 * k1 + c3 * k3 + c4 * k4 + c6 * k6;
                
                // Set step_sizing to false
                sizing = false;
            }
        }

        //std::cout << r_change(0) << ' ' << r_change(1) << " " << '\n';
        
        // Account for sign missmatch
        double r_mag = r_change.norm();
        double dot = r_past.dot(r_change/r_mag);

        if (dot < -.95) {
            r_change *= -1;
        }

        r += r_change;
        r_past = r_change / r_mag;
        double val = eigen_solve_val(f(r, R, sigma, vX), icity);

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
}

int main() {
    //rka_iter(1.0, 1.0, 1.0, 2000, 1, .5);

    return 0;
}
