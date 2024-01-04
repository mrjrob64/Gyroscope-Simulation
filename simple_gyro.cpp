#include "Matrix.h"
#include "simple_gyro.h"
#include <cmath>
#include <fstream>

int main() {

    // In standard basis
    std::vector<double> position = {0.0, 0.1/std::sqrt(2), 0.1/std::sqrt(2)};
    std::vector<double> Fg = {0.0, 0.0, 0.0};
    std::vector<double> angularVelocity = {0.0, 0.0, 100.0};
    std::vector<double> zAxis = {0.0, 0.0, 1.0};

    // simulation parameters
    double dt = 0.0000001;
    double t = 0.0;
    double tmax = 1.0;
    double I1 = 0.1;
    double I2 = 0.1;
    double I3 = 0.1;

    //Initialize csv file
    std::ofstream output;
    output.open("simple_gyro.csv");

    while(t < tmax) {
        std::vector<double> torque_standard_basis = crossProduct(position, Fg);

        //print time and position to csv file
        output << t << "," << position[0] << "," << position[1] << "," << position[2] << "," << angularVelocity[0] << "," << angularVelocity[1] << "," << angularVelocity[2] << "," << torque_standard_basis[0] << "," << torque_standard_basis[1] << "," << torque_standard_basis[2] << std::endl;

        
        std::vector<double> torque_body_basis;

        // use change of basis to convert torque from standard basis to body basis where z-axis is aligned with position vector
        // first get orthonormal basis for body basis
        std::vector<double> e3 = normalize(position);
        std::vector<double> e2 = normalize(crossProduct(position, zAxis));
        std::vector<double> e1 = normalize(crossProduct(e2, e3));

        // change of basis matrix
        Matrix B(e1[0], e2[0], e3[0],
                 e1[1], e2[1], e3[1],
                 e1[2], e2[2], e3[2]);
        
        // convert torque from standard basis to body basis
        torque_body_basis = B * torque_standard_basis;

        // calculate angular acceleration
        std::vector<double> angularAcceleration = {torque_body_basis[0]/I1, torque_body_basis[1]/I2, torque_body_basis[2]/I3};
        //convert angular acceleration back to standard basis
        angularAcceleration = B.inverse() * angularAcceleration;

        // update angular velocity
        angularVelocity[0] += angularAcceleration[0] * dt;
        angularVelocity[1] += angularAcceleration[1] * dt;
        angularVelocity[2] += angularAcceleration[2] * dt;

        // update position
        // first define rotation matrix about z-axis
        Matrix Rz(std::cos(angularVelocity[2] * dt), -std::sin(angularVelocity[2] * dt), 0.0,
                  std::sin(angularVelocity[2] * dt), std::cos(angularVelocity[2] * dt), 0.0,
                  0.0, 0.0, 1.0);

        // get position in body basis
        std::vector<double> position_body_basis = B * position;
        
        // rotate position in body basis
        position_body_basis = Rz * position_body_basis;

        // convert position back to standard basis
        position = B.inverse() * position_body_basis;

        // std::vector<double> velocity = crossProduct(angularVelocity, position);

        // position[0] += velocity[0] * dt;
        // position[1] += velocity[1] * dt;
        // position[2] += velocity[2] * dt;

        // // renormalize position
        // position = normalize(position);
        // position[0] *= 0.1;
        // position[1] *= 0.1;
        // position[2] *= 0.1;

        //angular velocity needs to be rotated by amount position was rotated
        
        // update time
        t += dt;

        //if position is below ground, get rid of velocity in that direction
        // if (position[2] < 0.0) {
        //     position[2] = 0.0;
        // }
        // friction with ground
        // apply friction to angular velocity
        // double frictionCoeff = 0.001;
        // for (int i = 0; i < 3; ++i) {
        //     angularVelocity[i] *= std::exp(-frictionCoeff * dt);
        // }
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
    double magnitude = mag(vec);
    std::vector<double> result;
    result.push_back(vec[0]/magnitude);
    result.push_back(vec[1]/magnitude);
    result.push_back(vec[2]/magnitude);
    return result;
}