// Matrix.h

#include <iostream>
#include <vector>

class Matrix {
public:

    // Default constructor
    Matrix() {
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                data[i][j] = 0.0;
            }
        }
        data[0][0] = 1.0;
        data[1][1] = 1.0;
        data[2][2] = 1.0;
    }

    // Constructor with elements
    Matrix(double a11, double a12, double a13,
           double a21, double a22, double a23,
           double a31, double a32, double a33) {
        data[0][0] = a11; data[0][1] = a12; data[0][2] = a13;
        data[1][0] = a21; data[1][1] = a22; data[1][2] = a23;
        data[2][0] = a31; data[2][1] = a32; data[2][2] = a33;
    }

    // Matrix inversion
    Matrix inverse() {
        //Is matrix invertible?
        double det = data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
                     data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
                     data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
        if(det == 0) {
            std::cout << "Matrix is not invertible" << std::endl;
            //end program
            exit(1);
            return Matrix();
        }

        //Calculate inverse
        Matrix result;
        result.data[0][0] = (1/det) * (data[1][1] * data[2][2] - data[1][2] * data[2][1]);
        result.data[0][1] = (1/det) * (data[0][2] * data[2][1] - data[0][1] * data[2][2]);
        result.data[0][2] = (1/det) * (data[0][1] * data[1][2] - data[0][2] * data[1][1]);
        result.data[1][0] = (1/det) * (data[1][2] * data[2][0] - data[1][0] * data[2][2]);
        result.data[1][1] = (1/det) * (data[0][0] * data[2][2] - data[0][2] * data[2][0]);
        result.data[1][2] = (1/det) * (data[0][2] * data[1][0] - data[0][0] * data[1][2]);
        result.data[2][0] = (1/det) * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
        result.data[2][1] = (1/det) * (data[0][1] * data[2][0] - data[0][0] * data[2][1]);
        result.data[2][2] = (1/det) * (data[0][0] * data[1][1] - data[0][1] * data[1][0]);

        return result;

    }

    // Overloading the * operator for matrix multiplication
    Matrix operator*(const Matrix& other) const {
        Matrix result;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.data[i][j] = 0;
                for (int k = 0; k < 3; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }

        return result;
    }

    // Overloading the * operator for matrix-vector multiplication
    std::vector<double> operator*(const std::vector<double>& vec) const {
        std::vector<double> result;

        for (int i = 0; i < 3; ++i) {
            result.push_back(0);
            for (int j = 0; j < 3; ++j) {
                result[i] += data[i][j] * vec[j];
            }
        }

        return result;
    }

    // Overloading the << operator for easy printing
    friend std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                os << mat.data[i][j] << " ";
            }
            os << std::endl;
        }
        return os;
    }

private:
    double data[3][3];
};
