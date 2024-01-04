#include "Matrix.h"
#include "simple_gyro.h"
#include <cmath>
#include <fstream>

int main() {

    // In standard basis
    std::vector<double> position = {0.0, 1.0*std::sqrt(3)/2, 1.0/2};
    std::vector<double> Fg = {0.0, 0.0, -0.4};
    std::vector<double> zAxis = {0.0, 0.0, 1.0};

    //In rotating body basis
    std::vector<double> angularVelocity = {0.0, 1.0, 50.0};


    // simulation parameters
    double dt = 0.00001;
    double t = 0.0;
    double tmax = 200.0;
    double I1 = 10.0;
    double I2 = 10.0;
    double I3 = 20.0;

    //Initialize csv file
    std::ofstream output;
    int STEPS = 1000;
    output.open("simple_gyro_v2.csv");

    int index = 0;

    while(t < tmax) {
        std::vector<double> torque_standard_basis = crossProduct(position, Fg);

        //convert angular velocity to standard basis for printing
        std::vector<double> angularVelocity_standard_basis = {0.0, 0.0, 0.0};

        std::vector<double> e3 = normalize(position);
        std::vector<double> e2 = normalize(crossProduct(position, zAxis));
        std::vector<double> e1 = normalize(crossProduct(e2, e3));

        // change of basis matrix
        Matrix B(e1[0], e2[0], e3[0],
                 e1[1], e2[1], e3[1],
                 e1[2], e2[2], e3[2]);

        //std::cout << B_.inverse() << std::endl;

        angularVelocity_standard_basis = B * angularVelocity;

        std::vector<double> velocity = crossProduct(angularVelocity_standard_basis, position);

        if(index % STEPS == 0) {
            //print time and position to csv file
            output << t << "," << position[0] << "," << position[1] << "," << position[2] << "," << angularVelocity_standard_basis[0] << "," << angularVelocity_standard_basis[1] << "," << angularVelocity_standard_basis[2] << "," << torque_standard_basis[0] << "," << torque_standard_basis[1] << "," << torque_standard_basis[2] << "," << velocity[0] << "," << velocity[1] << "," << velocity[2] <<std::endl;
        }
        //Find angular acceleration using Euler equations
        // (dL/dt)space = (dL/dt)body + w x L
        std::vector<double> angular_acceleration(3);
        angular_acceleration[0] = torque_standard_basis[0]/I1 + ((I2-I3)/I1)*angularVelocity[1]*angularVelocity[2];
        angular_acceleration[1] = torque_standard_basis[1]/I2 + ((I3-I1)/I2)*angularVelocity[2]*angularVelocity[0];
        angular_acceleration[2] = torque_standard_basis[2]/I3 + ((I1-I2)/I3)*angularVelocity[0]*angularVelocity[1];

        // std::cout << torque_standard_basis[0]/I1 + ((I2-I3)/I1)*angularVelocity[1]*angularVelocity[2] << std::endl;
        // std::cout << angular_acceleration[0] << std::endl;
        // std::cout << angular_acceleration[1] << std::endl;
        // std::cout << angular_acceleration[2] << std::endl;
        // exit(1);
        angularVelocity[0] += dt * angular_acceleration[0];
        angularVelocity[1] += dt * angular_acceleration[1];
        angularVelocity[2] += dt * angular_acceleration[2];

        //Rotate Position Vector

        angularVelocity_standard_basis = B * angularVelocity;

        // //Step 1: Create Change of basis matrix so that we can define the new z' axis
        // // as the direction of the angualr velocity. This way we can then apply a
        // // rotation about the z axis for the given amount of degrees depending on
        // // angular velocity = change in angle / change in time
        // // change in angle = angular velocity * change in time
        // // This change in angle is a rotation in the basis where the z axis is the
        // // angular velocity

        // std::vector<double> e3p = normalize(angularVelocity);
        // std::vector<double> e2p = normalize(crossProduct(angularVelocity, zAxis));
        // std::vector<double> e1p = normalize(crossProduct(e2, e3));

        // // change of basis matrix
        // Matrix Bp(e1p[0], e2p[0], e3p[0],
        //          e1p[1], e2p[1], e3p[1],
        //          e1p[2], e2p[2], e3p[2]);

        // //Step 2: Create rotation matrix. We want to rotate an angle = angular velocity * change in time
        // Matrix Rz(std::cos(mag(angularVelocity) * dt), -std::sin(mag(angularVelocity) * dt), 0.0,
        //           std::sin(mag(angularVelocity) * dt), std::cos(mag(angularVelocity) * dt), 0.0,
        //           0.0, 0.0, 1.0);

        // //Step 3: Apply change of basis matrix to position vector
        // std::vector<double> transformed_position = Bp.inverse() * position;

        // //Step 4: Apply rotation to position in new basis
        // std::vector<double> rotated_transformed_position = Rz * transformed_position;

        // //Step 5: Put the new rotated position back into standard coordinates
        // position = Bp * rotated_transformed_position;

        //std::cout << t << std::endl;

        velocity = crossProduct(angularVelocity_standard_basis, position);
        position[0] += velocity[0] * dt;
        position[1] += velocity[1] * dt;
        position[2] += velocity[2] * dt;
        position = normalize(position);

        //Total Energy should remain constant
        //E0 = 1/2(I3)(w3_0)^2 + Fg*z_0
        //E = 1/2(I1)(w1)^2 + 1/2(I2)(w2)^2 + 1/2(I3)(w3)^2 + Fg*z
        //w3 will remain constant
        //It will be difficult to readjust z

        //We can try normalizing w1 and w2 instead so error does not cause these values to explode
        //Constant = I1*w1*w1 + I2*w2*w2
        //We can calculate what 2*(E-1/2(I3)(w3)^2-Fg*z) should be and then compare to I1*w1*w1 + I2*w2*w2

        //First Option
        //Solve equation Expected = I1*(w1+stuff)^2 + I2*(w2+stuff)^2
        //I1 * (w1^2 + 2*w1*stuff + stuff^2) + I2 * (w2^2 + 2*w2*stuff + stuff^2)
        //(I1 * w1^2 + I2 * w2^2 - Constant) + (I1*2*w1 + I2*2*w2) * stuff + (I1 + I2) * stuff^2 = 0
        //May not work if no real roots

        //Second Option
        //

        t+=dt;
        index++;
    }

    return 0;
}

std::vector<double> crossProduct(std::vector<double> a, std::vector<double> b) {
    std::vector<double> result;
    result.push_back(a[1] * b[2] - a[2] * b[1]);
    result.push_back(a[2] * b[0] - a[0] * b[2]);
    result.push_back(a[0] * b[1] - a[1] * b[0]);
    return result;
}

double mag(std::vector<double> vec) {
    return std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

std::vector<double> normalize(std::vector<double> vec) {
    if(mag(vec) == 0.0) {
        return vec;
    }
    double magnitude = mag(vec);
    std::vector<double> result;
    result.push_back(vec[0]/magnitude);
    result.push_back(vec[1]/magnitude);
    result.push_back(vec[2]/magnitude);
    return result;
}