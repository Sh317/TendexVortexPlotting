#include <iostream>
#include <Eigen>
#include <cfloat>
#include <cmath>
 

// Code to test my implementation of the B field since its a lot more complicated than the E field
using namespace std;
using namespace Eigen;

extern "C" {

struct mat {
    double m[3][3];
};

// Define hyperbolic functions
double sech(double x) {
    return 1.0 / cosh(x);
}
double coth(double x) {
    return 1.0 / tanh(x);
}

// Function to compute the matrix
mat f(double x, double y, double z, double R, double S, double vX) {
        Eigen::Matrix3f B;
        mat BM;
    
        // Precompute common
        double x2 = x * x;
        double y2 = y * y;
        double z2 = z * z;
        double rS = R * S;
        double rS2 = rS * rS;
        double r = std::sqrt(x2 + y2 + z2);
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

        // Process for output
        BM.m[0][0] = B(0,0);
        BM.m[0][1] = B(0,1);
        BM.m[0][2] = B(0,2);
        BM.m[1][0] = B(1,0);
        BM.m[1][1] = B(1,1);
        BM.m[1][2] = B(1,2);
        BM.m[2][0] = B(2,0);
        BM.m[2][1] = B(2,1);
        BM.m[2][2] = B(2,2);
        return BM;
}

}


int main() {
    return 0;
}