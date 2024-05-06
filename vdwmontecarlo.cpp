#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <random>
#include <ctime>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <regex>
#include <tuple>
#include <utility>

using namespace std;

double LJPotentialBetweenTwoAtoms(double Sigma, double Epsilon, double distance) {
    double SigmaDividedByDistance = Sigma / distance;
    double Energy = 4 * Epsilon * (pow(SigmaDividedByDistance, 12) - pow(SigmaDividedByDistance, 6));
    return Energy;
}

double LJPotentialBetweenTwoCubes(double Sigma, double Epsilon, double LJCutOff,
                                   const double atomList1[][3], int size1,
                                   const double atomList2[][3], int size2,
                                   double BoxLength) {
    double CubeCubeEnergy = 0;
    int row1 = size1 * size1 * size1;
    int row2 = size2 * size2 * size2;

    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < row2; ++j) {
            double Atom1X = atomList1[i][0];
            double Atom1Y = atomList1[i][1];
            double Atom1Z = atomList1[i][2];
            double Atom2X = atomList2[j][0];
            double Atom2Y = atomList2[j][1];
            double Atom2Z = atomList2[j][2];

            double dx = Atom1X - Atom2X;
            dx = dx - BoxLength * round(dx / BoxLength);
            double dy = Atom1Y - Atom2Y;
            dy = dy - BoxLength * round(dy / BoxLength);
            double dz = Atom1Z - Atom2Z;
            dz = dz - BoxLength * round(dz / BoxLength);
            double distance = sqrt(dx * dx + dy * dy + dz * dz);

            if (distance <= LJCutOff) {
                double energy = LJPotentialBetweenTwoAtoms(Sigma, Epsilon, distance);
                CubeCubeEnergy += energy;
            }
        }
    }

    return CubeCubeEnergy;
}

typedef double Vector3[3];
typedef double Matrix3x3[3][3];

// Function to multiply a 3x3 matrix with a 3D vector
void MultiplyMatrixVector(const Matrix3x3& matrix, const Vector3 input, Vector3 result) {
    for (int i = 0; i < 3; ++i) {
        result[i] = 0;
        for (int j = 0; j < 3; ++j) {
            result[i] += input[j] * matrix[j][i];
        }
    }
}

void vectorSubtraction(const double a[3], const double b[3], double result[3]) {
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}

void normalizeVector(double vec[3]) {
    double norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    vec[0] /= norm;
    vec[1] /= norm;
    vec[2] /= norm;
}

void vectorAddition(const double a[3], const double b[3], double result[3]) {
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
}

void vectorScalarProduct(const double vec[3], double scalar, double result[3]) {
    result[0] = scalar * vec[0];
    result[1] = scalar * vec[1];
    result[2] = scalar * vec[2];
}

void CubeRotation(const Vector3& CubeCentroid, const Vector3 VectorX, const Vector3 VectorY, const Vector3 VectorZ, const Vector3& RotationCenter,
                  const Vector3& RotationAngles, Vector3& CubeCentroidResult, Vector3 CubeVectorResult[3]) {

    Vector3 RotVectorX, RotVectorY, RotVectorZ;

    double Alpha = RotationAngles[0] * M_PI / 180.0;
    double Beta = RotationAngles[1] * M_PI / 180.0;
    double Gamma = RotationAngles[2] * M_PI / 180.0;

    Matrix3x3 RotationX = {{1, 0, 0}, {0, cos(Alpha), -sin(Alpha)}, {0, sin(Alpha), cos(Alpha)}};
    Matrix3x3 RotationY = {{cos(Beta), 0, sin(Beta)}, {0, 1, 0}, {-sin(Beta), 0, cos(Beta)}};
    Matrix3x3 RotationZ = {{cos(Gamma), -sin(Gamma), 0}, {sin(Gamma), cos(Gamma), 0}, {0, 0, 1}};

    Vector3 RotationCenter_To_CubeCentroid;
    for (int i = 0; i < 3; ++i) {
        RotationCenter_To_CubeCentroid[i] = CubeCentroid[i] - RotationCenter[i];
    }

    Vector3 tempresult1;
    Vector3 tempresult2;
    MultiplyMatrixVector(RotationX, VectorX, tempresult1);
    MultiplyMatrixVector(RotationY, tempresult1, tempresult2);
    MultiplyMatrixVector(RotationZ, tempresult2, RotVectorX);

    MultiplyMatrixVector(RotationX, VectorY, tempresult1);
    MultiplyMatrixVector(RotationY, tempresult1, tempresult2);
    MultiplyMatrixVector(RotationZ, tempresult2, RotVectorY);

    MultiplyMatrixVector(RotationX, VectorZ, tempresult1);
    MultiplyMatrixVector(RotationY, tempresult1, tempresult2);
    MultiplyMatrixVector(RotationZ, tempresult2, RotVectorZ);

    MultiplyMatrixVector(RotationX, RotationCenter_To_CubeCentroid, tempresult1);
    MultiplyMatrixVector(RotationY, tempresult1, tempresult2);
    MultiplyMatrixVector(RotationZ, tempresult2, RotationCenter_To_CubeCentroid);

    double lengthX = sqrt(RotVectorX[0] * RotVectorX[0] + RotVectorX[1] * RotVectorX[1] + RotVectorX[2] * RotVectorX[2]);
    double lengthY = sqrt(RotVectorY[0] * RotVectorY[0] + RotVectorY[1] * RotVectorY[1] + RotVectorY[2] * RotVectorY[2]);
    double lengthZ = sqrt(RotVectorZ[0] * RotVectorZ[0] + RotVectorZ[1] * RotVectorZ[1] + RotVectorZ[2] * RotVectorZ[2]);

    RotVectorX[0] /= lengthX;
    RotVectorX[1] /= lengthX;
    RotVectorX[2] /= lengthX;

    RotVectorY[0] /= lengthY;
    RotVectorY[1] /= lengthY;
    RotVectorY[2] /= lengthY;

    RotVectorZ[0] /= lengthZ;
    RotVectorZ[1] /= lengthZ;
    RotVectorZ[2] /= lengthZ;

    for (int i = 0; i < 3; ++i) {
        CubeCentroidResult[i] = RotationCenter_To_CubeCentroid[i] + RotationCenter[i];
        CubeVectorResult[0][i] = RotVectorX[i];
        CubeVectorResult[1][i] = RotVectorY[i];
        CubeVectorResult[2][i] = RotVectorZ[i];
    }
}

void matrixMultiply(const double A[3][3], const double B[3][3], double result[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = 0;
            for (int k = 0; k < 3; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void findVerticesOfCube(double CubeSideLength, double COMX, double COMY, double COMZ, const Vector3 VectorX, const Vector3 VectorY, const Vector3 VectorZ, double Vertex[][3]){

    double matrix[8][3] = {
        {-1, +1, +1},
        {-1, -1, +1},
        {-1, -1, -1},
        {-1, +1, -1},
        {+1, +1, +1},
        {+1, -1, +1},
        {+1, -1, -1},
        {+1, +1, -1}
    };

    double x,y,z;
    for (int i = 0; i < 8; ++i) {
        x = matrix[i][0];
        y = matrix[i][1];
        z = matrix[i][2];
        Vertex[i][0] = COMX + (CubeSideLength / 2.0) * (x * VectorX[0] + y * VectorY[0] + z * VectorZ[0]);
        Vertex[i][1] = COMY + (CubeSideLength / 2.0) * (x * VectorX[1] + y * VectorY[1] + z * VectorZ[1]);
        Vertex[i][2] = COMZ + (CubeSideLength / 2.0) * (x * VectorX[2] + y * VectorY[2] + z * VectorZ[2]);
    }
}

void Reorientation(int i, int j, double Cube2SideLength, double BoxLength,
                   Vector3 Particle1Centroid, Vector3 Particle2Centroid,
                   Vector3 Particle1VectorsX, Vector3 Particle1VectorsY, Vector3 Particle1VectorsZ,
                   Vector3 Particle2VectorsX, Vector3 Particle2VectorsY, Vector3 Particle2VectorsZ, double Vertex[][3]) {

    // Create temporary variables to store intermediate values
    double Size1Temp, Size2Temp;
    Vector3 Particle1CentroidTemp, Particle2CentroidTemp;
    Vector3 Particle1VectorsXTemp, Particle1VectorsYTemp, Particle1VectorsZTemp;
    Vector3 Particle2VectorsXTemp, Particle2VectorsYTemp, Particle2VectorsZTemp;

    Vector3 P1P2;
    vectorSubtraction(Particle2Centroid, Particle1Centroid, P1P2);
    for (int k = 0; k < 3; ++k) {
       P1P2[k] = P1P2[k]-BoxLength*round(P1P2[k]/BoxLength);
    }

    double P1P2Norm = sqrt(P1P2[0] * P1P2[0] + P1P2[1] * P1P2[1] + P1P2[2] * P1P2[2]);
    //cout << "P1P2 " << P1P2[0] << " " << P1P2[1] << " " << P1P2[2] << endl;

    double AngleBetweenP1P2AndX = std::acos(std::max(-1.0, std::min(1.0, P1P2[0] / P1P2Norm)));
    double AngleBetweenP1P2AndXNegative = std::acos(std::max(-1.0, std::min(1.0, -P1P2[0] / P1P2Norm)));
    double AngleBetweenP1P2AndY = std::acos(std::max(-1.0, std::min(1.0, P1P2[1] / P1P2Norm)));
    double AngleBetweenP1P2AndYNegative = std::acos(std::max(-1.0, std::min(1.0, -P1P2[1] / P1P2Norm)));
    double AngleBetweenP1P2AndZ = std::acos(std::max(-1.0, std::min(1.0, P1P2[2] / P1P2Norm)));
    double AngleBetweenP1P2AndZNegative = std::acos(std::max(-1.0, std::min(1.0, -P1P2[2] / P1P2Norm)));

    // Find out which axis (and direction) has the smallest angle with the P1P2 vector
    double MinimumAngle = min({AngleBetweenP1P2AndX, AngleBetweenP1P2AndXNegative, AngleBetweenP1P2AndY,
                                    AngleBetweenP1P2AndYNegative, AngleBetweenP1P2AndZ, AngleBetweenP1P2AndZNegative});

    //cout << "Minimum angle " << MinimumAngle << endl;
    //cout << AngleBetweenP1P2AndX << " " << AngleBetweenP1P2AndXNegative << " " << AngleBetweenP1P2AndY << " " << AngleBetweenP1P2AndYNegative << " " << AngleBetweenP1P2AndZ << " " << AngleBetweenP1P2AndZNegative << endl;

    double temp1, temp2;
    bool swapx = false;
    bool swapy = false;
    bool swapz = false;

    if (AngleBetweenP1P2AndXNegative == MinimumAngle) {
       swapx = true;
       for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
    } else if (AngleBetweenP1P2AndYNegative == MinimumAngle) {
       swapy = true;
       for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
    } else if (AngleBetweenP1P2AndZNegative == MinimumAngle) {
       swapz = true;
       for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
    }

    // Define the angles to rotate cube 2 about cube 1
    double VirtualMoveAlpha, VirtualMoveBeta, VirtualMoveGamma;

    /*double sinBeta1 = Particle1VectorsX[2];
    double cosBeta1 = sqrt(1-sinBeta1*sinBeta1);
    double sinAlpha1 = -Particle1VectorsY[2]/cosBeta1;
    double cosAlpha1 = Particle1VectorsZ[2]/cosBeta1;
    double sinGamma1 = -Particle1VectorsX[1]/cosBeta1;
    double cosGamma1 = Particle1VectorsX[0]/cosBeta1;

    double sinBeta2 = Particle2VectorsX[2];
    double cosBeta2 = sqrt(1-sinBeta2*sinBeta2);
    double sinAlpha2 = -Particle2VectorsY[2]/cosBeta2;
    double cosAlpha2 = Particle2VectorsZ[2]/cosBeta2;
    double sinGamma2 = -Particle2VectorsX[1]/cosBeta2;
    double cosGamma2 = Particle2VectorsX[0]/cosBeta2;
    */
    /*if (abs(Particle1VectorsY[0]) > abs(Particle1VectorsX[0])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsX[k];
          temp2 = Particle1VectorsY[k];
          Particle1VectorsX[k] = temp2;
          Particle1VectorsY[k] = temp1;
       }
    } else if (abs(Particle1VectorsZ[0]) > abs(Particle1VectorsX[0])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsX[k];
          temp2 = Particle1VectorsZ[k];
          Particle1VectorsX[k] = temp2;
          Particle1VectorsZ[k] = temp1;
       }
    }

    if (abs(Particle1VectorsX[1]) > abs(Particle1VectorsY[1])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsX[k];
          temp2 = Particle1VectorsY[k];
          Particle1VectorsX[k] = temp2;
          Particle1VectorsY[k] = temp1;
       }
    } else if (abs(Particle1VectorsZ[1]) > abs(Particle1VectorsY[1])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsY[k];
          temp2 = Particle1VectorsZ[k];
          Particle1VectorsY[k] = temp2;
          Particle1VectorsZ[k] = temp1;
       }
    }

    if (abs(Particle1VectorsX[2]) > abs(Particle1VectorsZ[2])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsX[k];
          temp2 = Particle1VectorsZ[k];
          Particle1VectorsX[k] = temp2;
          Particle1VectorsZ[k] = temp1;
       }
    } else if (abs(Particle1VectorsY[2]) > abs(Particle1VectorsZ[2])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsY[k];
          temp2 = Particle1VectorsZ[k];
          Particle1VectorsY[k] = temp2;
          Particle1VectorsZ[k] = temp1;
       }
    }

   if (abs(Particle2VectorsY[0]) > abs(Particle2VectorsX[0])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsX[k];
          temp2 = Particle2VectorsY[k];
          Particle2VectorsX[k] = temp2;
          Particle2VectorsY[k] = temp1;
       }
    } else if (abs(Particle2VectorsZ[0]) > abs(Particle2VectorsX[0])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsX[k];
          temp2 = Particle2VectorsZ[k];
          Particle2VectorsX[k] = temp2;
          Particle2VectorsZ[k] = temp1;
       }
    }

    if (abs(Particle2VectorsX[1]) > abs(Particle2VectorsY[1])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsX[k];
          temp2 = Particle2VectorsY[k];
          Particle2VectorsX[k] = temp2;
          Particle2VectorsY[k] = temp1;
       }
    } else if (abs(Particle2VectorsZ[1]) > abs(Particle2VectorsY[1])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsY[k];
          temp2 = Particle2VectorsZ[k];
          Particle2VectorsY[k] = temp2;
          Particle2VectorsZ[k] = temp1;
       }
    }

    if (abs(Particle2VectorsX[2]) > abs(Particle2VectorsZ[2])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsX[k];
          temp2 = Particle2VectorsZ[k];
          Particle2VectorsX[k] = temp2;
          Particle2VectorsZ[k] = temp1;
       }
    } else if (abs(Particle2VectorsY[2]) > abs(Particle2VectorsZ[2])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsY[k];
          temp2 = Particle2VectorsZ[k];
          Particle2VectorsY[k] = temp2;
          Particle2VectorsZ[k] = temp1;
       }
    }*/

    // Generate rotation matrix for cube 1
    double R1[3][3] = {{Particle1VectorsX[0], Particle1VectorsX[1], Particle1VectorsX[2]},
                       {Particle1VectorsY[0], Particle1VectorsY[1], Particle1VectorsY[2]},
                       {Particle1VectorsZ[0], Particle1VectorsZ[1], Particle1VectorsZ[2]}};

    /*double inv_1x[3][3] = {{1.0, 0.0, 0.0},
                          {0.0, cosAlpha1, sinAlpha1},
                          {0.0, -sinAlpha1, cosAlpha1}};

    double inv_1y[3][3] = {{cosBeta1, 0.0, -sinBeta1},
                           {0.0, 1.0, 0.0},
                           {sinBeta1, 0.0, cosBeta1}};

    double inv_1z[3][3] = {{cosGamma1, sinGamma1, 0.0},
                           {-sinGamma1, cosGamma1, 0.0},
                           {0.0, 0.0, 1.0}};

    cout << "Inverse rotation x matrix" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << inv_1x[k][l] << " " ;
        }
        cout << endl;
    }

    cout << "Inverse rotation y matrix" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << inv_1y[k][l] << " " ;
        }
        cout << endl;
    }
    cout << "Inverse rotation z matrix" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << inv_1z[k][l] << " " ;
        }
        cout << endl;
    }

    double R2_x[3][3] = {{1.0, 0.0, 0.0},
                          {0.0, cosAlpha2, -sinAlpha2},
                          {0.0, sinAlpha2, cosAlpha2}};

    double R2_y[3][3] = {{cosBeta2, 0.0, sinBeta2},
                           {0.0, 1.0, 0.0},
                           {-sinBeta2, 0.0, cosBeta2}};

    double R2_z[3][3] = {{cosGamma2, -sinGamma2, 0.0},
                           {sinGamma2, cosGamma2, 0.0},
                           {0.0, 0.0, 1.0}};

    double resultx[3][3], resulty[3][3], resultz[3][3];
    matrixMultiply(R2_x, inv_1x, resultx);
    matrixMultiply(R2_y, inv_1y, resulty);
    matrixMultiply(R2_z, inv_1z, resultz);

    cout << "Resultx" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << resultx[k][l] << " " ;
        }
        cout << endl;
    }

    cout << "Resulty" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << resulty[k][l] << " " ;
        }
        cout << endl;
    }

    cout << "Resultz" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << resultz[k][l] << " " ;
        }
        cout << endl;
    }

    double ryz[3][3], inv_1[3][3];
    matrixMultiply(resulty, resultz, ryz);
    matrixMultiply(resultx, ryz, inv_1);
    */
    // Calculate the inverse matrix of cube 1, rotation matrix is orthogonal matrix, so inverse of it is equal to its transpose
    double inv_1[3][3];
    //cout << "Inverse rotation matrix" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            inv_1[k][l] = R1[l][k];
            //cout << inv_1[k][l] << " " ;
        }
        //cout << endl;
    }

    Vector3 RotatedParticle2Centroid, RotatedParticle2VectorsX, RotatedParticle2VectorsY, RotatedParticle2VectorsZ;
    //cout << "Particle2VectorsX[0] " << Particle2VectorsX[0] << " Particle2VectorsX[1] " << Particle2VectorsX[1] << " Particle2VectorsX[2] " << Particle2VectorsX[2] << endl;
    MultiplyMatrixVector(inv_1, P1P2, RotatedParticle2Centroid);
    MultiplyMatrixVector(inv_1, Particle2VectorsX, RotatedParticle2VectorsX);
    //cout << "RotatedParticle2VectorsX[0] " << RotatedParticle2VectorsX[0] << " RotatedParticle2VectorsX[1] " << RotatedParticle2VectorsX[1] << " RotatedParticle2VectorsX[2] " << RotatedParticle2VectorsX[2] << endl;
    MultiplyMatrixVector(inv_1, Particle2VectorsY, RotatedParticle2VectorsY);
    MultiplyMatrixVector(inv_1, Particle2VectorsZ, RotatedParticle2VectorsZ);

    double CubeVertex[8][3];
    findVerticesOfCube(Cube2SideLength,RotatedParticle2Centroid[0],RotatedParticle2Centroid[1],RotatedParticle2Centroid[2],RotatedParticle2VectorsX,RotatedParticle2VectorsY,RotatedParticle2VectorsZ,CubeVertex);

    double Norm = sqrt(RotatedParticle2Centroid[0]*RotatedParticle2Centroid[0]+RotatedParticle2Centroid[1]*RotatedParticle2Centroid[1]+RotatedParticle2Centroid[2]*RotatedParticle2Centroid[2]);
    AngleBetweenP1P2AndX = std::acos(std::max(-1.0, std::min(1.0, RotatedParticle2Centroid[0] / Norm)));
    AngleBetweenP1P2AndXNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticle2Centroid[0] / Norm)));
    AngleBetweenP1P2AndY = std::acos(std::max(-1.0, std::min(1.0, RotatedParticle2Centroid[1] / Norm)));
    AngleBetweenP1P2AndYNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticle2Centroid[1] / Norm)));
    AngleBetweenP1P2AndZ = std::acos(std::max(-1.0, std::min(1.0, RotatedParticle2Centroid[2] / Norm)));
    AngleBetweenP1P2AndZNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticle2Centroid[2] / Norm)));
    // Find out which axis (and direction) has the smallest angle with the P1P2 vector
    MinimumAngle = min({AngleBetweenP1P2AndX, AngleBetweenP1P2AndXNegative, AngleBetweenP1P2AndY,
                                    AngleBetweenP1P2AndYNegative, AngleBetweenP1P2AndZ, AngleBetweenP1P2AndZNegative});

    if (AngleBetweenP1P2AndX == MinimumAngle) {
        Particle2Centroid[0] = RotatedParticle2Centroid[0];
        Particle2Centroid[1] = RotatedParticle2Centroid[1];
        Particle2Centroid[2] = RotatedParticle2Centroid[2];
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = CubeVertex[k][0];
           Vertex[k][1] = CubeVertex[k][1];
           Vertex[k][2] = CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndY == MinimumAngle) {
        Particle2Centroid[0] = RotatedParticle2Centroid[1];
        Particle2Centroid[1] = -RotatedParticle2Centroid[0];
        Particle2Centroid[2] = RotatedParticle2Centroid[2];
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = CubeVertex[k][1];
           Vertex[k][1] = -CubeVertex[k][0];
           Vertex[k][2] = CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndZ == MinimumAngle) {
        Particle2Centroid[0] = RotatedParticle2Centroid[2];
        Particle2Centroid[1] = RotatedParticle2Centroid[1];
        Particle2Centroid[2] = -RotatedParticle2Centroid[0];
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = CubeVertex[k][2];
           Vertex[k][1] = CubeVertex[k][1];
           Vertex[k][2] = -CubeVertex[k][0];
        }
    } else if (AngleBetweenP1P2AndXNegative == MinimumAngle) { // in case reorienting centroid again results in -x axis
        Particle2Centroid[0] = -RotatedParticle2Centroid[0];
        Particle2Centroid[1] = -RotatedParticle2Centroid[1];
        Particle2Centroid[2] = -RotatedParticle2Centroid[2];
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = -CubeVertex[k][0];
           Vertex[k][1] = -CubeVertex[k][1];
           Vertex[k][2] = -CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndYNegative == MinimumAngle) {
        Particle2Centroid[0] = -RotatedParticle2Centroid[1];
        Particle2Centroid[1] = RotatedParticle2Centroid[0];
        Particle2Centroid[2] = -RotatedParticle2Centroid[2];
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = -CubeVertex[k][1];
           Vertex[k][1] = CubeVertex[k][0];
           Vertex[k][2] = -CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndZNegative == MinimumAngle) {
        Particle2Centroid[0] = -RotatedParticle2Centroid[2];
        Particle2Centroid[1] = -RotatedParticle2Centroid[1];
        Particle2Centroid[2] = RotatedParticle2Centroid[0];
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = -CubeVertex[k][2];
           Vertex[k][1] = -CubeVertex[k][1];
           Vertex[k][2] = CubeVertex[k][0];
        }
    }
    // Reset cube 1 orientation to standard basis vectors
    Particle1Centroid[0] = 0.0;
    Particle1Centroid[1] = 0.0;
    Particle1Centroid[2] = 0.0;
    Particle1VectorsX[0] = 1.0;
    Particle1VectorsX[1] = 0.0;
    Particle1VectorsX[2] = 0.0;
    Particle1VectorsY[0] = 0.0;
    Particle1VectorsY[1] = 1.0;
    Particle1VectorsY[2] = 0.0;
    Particle1VectorsZ[0] = 0.0;
    Particle1VectorsZ[1] = 0.0;
    Particle1VectorsZ[2] = 1.0;

   for (int k = 0; k < 3; ++k) {
      Particle2VectorsX[k] = RotatedParticle2VectorsX[k];
      Particle2VectorsY[k] = RotatedParticle2VectorsY[k];
      Particle2VectorsZ[k] = RotatedParticle2VectorsZ[k];
   }

}



/////////////////////////////////////  vdW energy calculations ///////////////////////////////////////////

// Computes the analytical expressions for the energies when dx=dx(y,z)
double ecalcgen(double a1, double a2, double b1, double b2, double c1, double c2, double c3,
                 double ym, double yp, double sig, double cene, double n, double c3now, int caseind) {

    double u1, u2, energy;
    sig = 1.0; // all distances are written in terms of sigma

    if (caseind == 1) {
        u1 = (yp - ym) * pow((a1 * c2 + c3now), 1 - n);
        u2 = (ym - yp) * pow((a2 * c2 + c3now), 1 - n);
        energy = cene * (u1 + u2) * pow(sig, n - 2) / (c2 * (n - 1) * (n - 2));
    } else if (caseind == 2) {
        u1 = (yp - ym) * pow((a1 * c2 + c3now), 1 - n);
        u2 = pow((a2 * c2 + c3now + (c1 + b2 * c2) * yp), 2 - n) - pow((a2 * c2 + c3now + (c1 + b2 * c2) * ym), 2 - n);
        u2 /= (c1 + b2 * c2);
        energy = cene * (u1 + u2) * pow(sig, n - 2) / (c2 * (n - 1) * (n - 2));
    } else if (caseind == 3) {
        u1 = pow((a1 * c2 + c3now + (c1 + b1 * c2) * ym), 2 - n) -
             pow((a1 * c2 + c3now + (c1 + b1 * c2) * yp), 2 - n);
        u1 /= (c1 + b1 * c2);
        u2 = (ym - yp) * pow((a2 * c2 + c3now), 1 - n);
        energy = cene * (u1 + u2) * pow(sig, n - 2) / (c2 * (n - 1) * (n - 2));
    } else if (caseind == 4) {
        u1 = pow((a1 * c2 + c3now + (c1 + b1 * c2) * ym), 2 - n) -
             pow((a1 * c2 + c3now + (c1 + b1 * c2) * yp), 2 - n);
        u1 /= (c1 + b1 * c2);
        u2 = pow((a2 * c2 + c3now + (c1 + b2 * c2) * yp), 2 - n) -
             pow((a2 * c2 + c3now + (c1 + b2 * c2) * ym), 2 - n);
        u2 /= (c1 + b2 * c2);
        energy = cene * (u1 + u2) * pow(sig, n - 2) / (c2 * (n - 1) * (n - 2));
    }

    return energy;
}

// Computes the analytical expressions for the energies when dx=dx(y)
double ecalcc2(double a1, double a2, double b1, double b2, double c1, double c2, double c3,
               double ym, double yp, double sig, double cene, double n, double c3now, int caseind) {

    double energy;
    sig = 1.0; // all distances are written in terms of sigma

    double u1 = pow((c1 * yp + c3now), 1 - n) * ((a1 * c1 - a2 * c1) * (n - 2) + (b1 - b2) * (c3now + c1 * (n - 1) * yp));
    double u2 = -pow((c1 * ym + c3now), 1 - n) * ((a1 * c1 - a2 * c1) * (n - 2) + (b1 - b2) * (c3now + c1 * (n - 1) * ym));

    energy = cene * pow(sig, n - 2) * (u1 + u2) / (pow(c1, 2) * (n - 2) * (n - 1));
    return energy;
}

// Computes the analytical expressions for the energies when dx=dx(z)
double ecalcc1(double a1, double a2, double b1, double b2, double c1, double c2, double c3,
               double ym, double yp, double sig, double cene, double n, double c3now, int caseind) {

    double u1, u2, energy;
    sig = 1.0; // all distances are written in terms of sigma

    if (caseind == 1) {
        energy = cene * pow(sig, n - 2) * ((yp - ym) * pow((a2 * c2 + c3now), 1 - n) + (ym - yp) * pow((a1 * c2 + c3now), 1 - n)) / (c2 * (1 - n));
    } else if (caseind == 2) {
        u1 = (yp - ym) * pow((a2 * c2 + c3now), 1 - n);
        u2 = pow((a1 * c2 + c3now + b1 * c2 * ym), 2 - n) - pow((a1 * c2 + c3now + b1 * c2 * yp), 2 - n);
        u2 /= (b1 * c2 * (2 - n));
        energy = cene * pow(sig, n - 2) * (u1 + u2) / (c2 * (1 - n));
    } else if (caseind == 3) {
        u1 = (ym - yp) * pow((a1 * c2 + c3now), 1 - n);
        u2 = pow((a2 * c2 + c3now + b2 * c2 * yp), 2 - n) - pow((a2 * c2 + c3now + b2 * c2 * ym), 2 - n);
        u2 /= (b2 * c2 * (2 - n));
        energy = cene * pow(sig, n - 2) * (u1 + u2) / (c2 * (1 - n));
    } else if (caseind == 4) {
        u1 = -pow((a2 * c2 + c3now + b2 * c2 * ym), 2 - n) + pow((a2 * c2 + c3now + b2 * c2 * yp), 2 - n);
        u1 /= (2 * b2 * c2 - b2 * c2 * n);
        u2 = pow((a1 * c2 + c3now + b1 * c2 * ym), 2 - n) - pow((a1 * c2 + c3now + b1 * c2 * yp), 2 - n);
        u2 /= (2 * b1 * c2 - b1 * c2 * n);
        energy = cene * (u1 + u2) * pow(sig, n - 2) / (c2 * (1 - n));
    }

    return energy;
}

// Comparator function to sort indices based on the desired column
bool compareIndices(const int& i, const int& j, int column, const double arr[4][3]) {
    return arr[i][column] < arr[j][column];
}

vector<int> sortIndicesByDesiredColumn(const double arr[4][3], int column) {
    // Create an array of indices
    vector<int> indices(4);
    iota(indices.begin(), indices.end(), 0);

    // Sort indices based on the second column
    sort(indices.begin(), indices.end(),
              [&arr, column](const int& i, const int& j) {
                  return compareIndices(i, j, column, arr);
              });

    return indices;
}

double calc_gen(double sig, double lc, double atre1, double atre2, double repe, double catr1, double catr2,double crep, double facein[][3], double region[][2], double Epsilon) {
    double uatr, urep, utot;
    utot = 0;
    sig = sig/10; // make A nm since this formula is written in terms of nm
    double dslim = 2.0;
    double slopetol1 = pow(10,-2);
    double slopetol2 = pow(10,-6);
    double yp = region[1][0];
    double ym = region[0][0];
    double b2 = (region[1][1] - region[0][1]) / (region[1][0] - region[0][0]);
    double a2 = region[0][1] - b2 * region[0][0];
    double b1 = (region[3][1] - region[2][1]) / (region[3][0] - region[2][0]);
    double a1 = region[2][1] - b1 * region[2][0];

    double temper[4][3];
    int column = 1; // second column
    vector<int> sortedIndices = sortIndicesByDesiredColumn(facein,column);
    // Assign rows based on the sorted indices
    for (int i = 0; i < 4; ++i) {
        int originalIndex = sortedIndices[i];
        copy(begin(facein[originalIndex]), end(facein[originalIndex]), begin(temper[i]));
    }

    double face[4][3];
    face[0][0] = temper[0][0];
    face[0][1] = temper[0][1];
    face[0][2] = temper[0][2];

    face[3][0] = temper[3][0];
    face[3][1] = temper[3][1];
    face[3][2] = temper[3][2];

    if (temper[1][2] > temper[2][2]) {
        face[1][0] = temper[1][0];
        face[1][1] = temper[1][1];
        face[1][2] = temper[1][2];

        face[2][0] = temper[2][0];
        face[2][1] = temper[2][1];
        face[2][2] = temper[2][2];
    } else {
        face[1][0] = temper[2][0];
        face[1][1] = temper[2][1];
        face[1][2] = temper[2][2];

        face[2][0] = temper[1][0];
        face[2][1] = temper[1][1];
        face[2][2] = temper[1][2];
    }

    double AA = (face[1][1] - face[0][1]) * (face[2][2] - face[0][2]) - (face[1][2] - face[0][2]) * (face[2][1] - face[0][1]);
    double BB = (face[1][2] - face[0][2]) * (face[2][0] - face[0][0]) - (face[1][0] - face[0][0]) * (face[2][2] - face[0][2]);
    double CC = (face[1][0] - face[0][0]) * (face[2][1] - face[0][1]) - (face[1][1] - face[0][1]) * (face[2][0] - face[0][0]);

    double c1 = -BB / AA;
    double c2 = -CC / AA;
    double c3 = face[0][0] - lc - c1 * face[0][1] - c2 * face[0][2];
    double n;

    if (c3 > dslim) {
        if (abs(c2) < slopetol1 && abs(c1) < slopetol1) {
            double area = ((region[0][1] - region[2][1]) + (region[1][1] - region[3][1])) * (region[1][0] - region[0][0]) / 2.0;
            n = atre2;
            uatr = catr2 * area * pow(c3, -n);
            urep = 0.0;
        } else if (abs(c1) < slopetol1) {
            int caseindc1;
            if (abs(b1) < slopetol2 && abs(b2) < slopetol2) {
                caseindc1 = 1;
            } else if (abs(b2) < slopetol2) {
                caseindc1 = 2;
            } else if (abs(b1) < slopetol2) {
                caseindc1 = 3;
            } else {
                caseindc1 = 4;
            }
            //cout << "Parameters" << endl;
            //cout << a1 << a2 << b1 << b2 << c1 << c2 << c3 << ym << yp << sig << catr2 << atre2 << c3 << caseindc1 << endl;
            uatr = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3, caseindc1);
            urep = 0.0;
        } else if (abs(c2) < slopetol1) {
            n = atre2;
            uatr = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3, 1);
            urep = 0.0;
        } else {
            double check1 = abs(c1 + b1 * c2);
            double check2 = abs(c1 + b2 * c2);
            int caseindgen;
            if (check1 < slopetol2 && check2 < slopetol2) {
                caseindgen = 1;
            } else if (check1 < slopetol2) {
                caseindgen = 2;
            } else if (check2 < slopetol2) {
                caseindgen = 3;
            } else {
                caseindgen = 4;
            }
            uatr = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3, caseindgen);
            urep = 0.0;
        }
    } else {
        double c3t = dslim;
        if (abs(c2) < slopetol1 && abs(c1) < slopetol1) {
            double area = ((region[0][1] - region[2][1]) + (region[1][1] - region[3][1])) * (region[1][0] - region[0][0]) / 2.0;
            n = atre1;
            uatr = catr1 * area * pow(c3, -n);
            double nt2 = atre2;
            double nt1 = atre1;
            double uatrt2 = catr2 * area * pow(c3t, -nt2);
            double uatrt1 = catr1 * area * pow(c3t, -nt1);
            uatr = uatr + uatrt2 - uatrt1;
            n = repe;
            urep = crep * area * pow(c3, -n);
            double ureptemp = crep * area * pow(c3t, -n);
            urep = urep - ureptemp;
        } else if (abs(c1) < slopetol1) {
            int caseindc1;
            if (abs(b1) < slopetol2 && abs(b2) < slopetol2) {
                caseindc1 = 1;
            } else if (abs(b2) < slopetol2) {
                caseindc1 = 2;
            } else if (abs(b1) < slopetol2) {
                caseindc1 = 3;
            } else {
                caseindc1 = 4;
            }
            uatr = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3, caseindc1);
            double uatrt1 = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3t, caseindc1);
            double uatrt2 = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3t, caseindc1);
            uatr = uatr + uatrt2 - uatrt1;
            n = repe;
            urep = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3, caseindc1);
            double ureptemp = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3t, caseindc1);
            urep = urep - ureptemp;
        } else if (abs(c2) < slopetol1) {
            uatr = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3, 1);
            double uatrt1 = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3t, 1);
            double uatrt2 = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3t, 1);
            uatr = uatr + uatrt2 - uatrt1;
            n = repe;
            urep = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3, 1);
            double ureptemp = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3t, 1);
            urep = urep - ureptemp;
        } else {
            double check1 = abs(c1 + b1 * c2);
            double check2 = abs(c1 + b2 * c2);
            int caseindgen;
            if (check1 < slopetol2 && check2 < slopetol2) {
                caseindgen = 1;
            } else if (check1 < slopetol2) {
                caseindgen = 2;
            } else if (check2 < slopetol2) {
                caseindgen = 3;
            } else {
                caseindgen = 4;
            }
            uatr = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3, caseindgen);
            double uatrt1 = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3t, caseindgen);
            double uatrt2 = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3t, caseindgen);
            uatr = uatr + uatrt2 - uatrt1;
            n = repe;
            urep = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3, caseindgen);
            double ureptemp = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3t, caseindgen);
            urep = urep - ureptemp;
        }
    }

    uatr = uatr * Epsilon;
    urep = urep * Epsilon;
    utot = uatr + urep;
    return utot;
}

struct Region {
    double points[4][2];
};

vector<Region> region_calc(double facein[][3], double lc) {

    const double reg_tol = lc + pow(10, -9);
    double temper[4][3];
    vector<int> sortedIndices = sortIndicesByDesiredColumn(facein,1);
    // Assign rows based on the sorted indices
    for (int i = 0; i < 4; ++i) {
        int originalIndex = sortedIndices[i];
        copy(begin(facein[originalIndex]), end(facein[originalIndex]), begin(temper[i]));
    }

    double face2[4][3];
    face2[0][0] = temper[0][0];
    face2[0][1] = temper[0][1];
    face2[0][2] = temper[0][2];

    face2[3][0] = temper[3][0];
    face2[3][1] = temper[3][1];
    face2[3][2] = temper[3][2];

    if (temper[1][2] > temper[2][2]) {
        face2[1][0] = temper[1][0];
        face2[1][1] = temper[1][1];
        face2[1][2] = temper[1][2];

        face2[2][0] = temper[2][0];
        face2[2][1] = temper[2][1];
        face2[2][2] = temper[2][2];
    } else {
        face2[1][0] = temper[2][0];
        face2[1][1] = temper[2][1];
        face2[1][2] = temper[2][2];

        face2[2][0] = temper[1][0];
        face2[2][1] = temper[1][1];
        face2[2][2] = temper[1][2];
    }

    int xcount=0;
    vector<double> xsave;

    /*cout <<"face2" << endl;
    for (int i = 0; i < 4; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << face2[i][j] << " ";
       }
       std::cout << std::endl;
    }*/

    // Line1
    double L1_slope = (face2[1][1] - face2[0][1]) / (face2[1][2] - face2[0][2]);
    double L1y1 = (-lc - face2[1][2]) * L1_slope + face2[1][1];
    double L1y2 = (lc - face2[1][2]) * L1_slope + face2[1][1];

    //cout << "L1_slope " << L1_slope << " L1y1 " << L1y1 << " L1y2 " << L1y2 << endl;
    if (L1y1 > face2[0][1] && L1y1 < face2[1][1] && L1y1 <= lc && L1y1 >= -lc) {
       xcount += 1;
       xsave.push_back(L1y1);
    }

    if (L1y2 > face2[0][1] && L1y2 < face2[1][1] && L1y2 <= lc && L1y2 >= -lc) {
       xcount += 1;
       xsave.push_back(L1y2);
    }

    //cout << "Line 1 xcount " << xcount << endl;
    // Line2
    double L2_slope = (face2[3][1] - face2[1][1]) / (face2[3][2] - face2[1][2]);
    double L2y1 = (-lc - face2[1][2]) * L2_slope + face2[1][1];
    double L2y2 = (lc - face2[1][2]) * L2_slope + face2[1][1];

    //cout << "L2_slope " << L2_slope << " L2y1 " << L2y1 << " L2y2 " << L2y2 << endl;
    if (L2y1 < face2[3][1] && L2y1 > face2[1][1] && L2y1 <= lc && L2y1 >= -lc) {
        xcount += 1;
        xsave.push_back(L2y1);
    }

    if (L2y2 < face2[3][1] && L2y2 > face2[1][1] && L2y2 <= lc && L2y2 >= -lc) {
        xcount += 1;
        xsave.push_back(L2y2);
    }
    //cout << "Line 2 xcount " << xcount << endl;
    // Line3
    double L3_slope = (face2[2][1] - face2[0][1]) / (face2[2][2] - face2[0][2]);
    double L3y1 = (-lc - face2[2][2]) * L3_slope + face2[2][1];
    double L3y2 = (lc - face2[2][2]) * L3_slope + face2[2][1];

    //cout << "L3_slope " << L3_slope << " L3y1 " << L3y1 << " L3y2 " << L3y2 << endl;
    if (L3y1 > face2[0][1] && L3y1 < face2[2][1] && L3y1 <= lc && L3y1 >= -lc) {
        xcount += 1;
        xsave.push_back(L3y1);
    }

    if (L3y2 > face2[0][1] && L3y2 < face2[2][1] && L3y2 <= lc && L3y2 >= -lc) {
        xcount += 1;
        xsave.push_back(L3y2);
    }
    //cout << "Line 3 xcount " << xcount << endl;
    // Line4
    double L4_slope = (face2[3][1] - face2[2][1]) / (face2[3][2] - face2[2][2]);
    double L4y1 = (-lc - face2[2][2]) * L4_slope + face2[2][1];
    double L4y2 = (lc - face2[2][2]) * L4_slope + face2[2][1];

    //cout << "L4_slope " << L4_slope << " L4y1 " << L4y1 << " L4y2 " << L4y2 << endl;
    if (L4y1 < face2[3][1] && L4y1 > face2[2][1] && L4y1 <= lc && L4y1 >= -lc) {
        xcount += 1;
        xsave.push_back(L4y1);
    }
    //cout << "Part 1 xcount " << xcount << endl;
    if (L4y2 < face2[3][1] && L4y2 > face2[2][1] && L4y2 <= lc && L4y2 >= -lc) {
        xcount += 1;
        xsave.push_back(L4y2);
    }
    //cout << "Part 2 xcount " << xcount << endl;
    for (int i = 0; i < 4; ++i) {
        if (face2[i][1] >= -lc && face2[i][1] <= lc && face2[i][2] >= -lc && face2[i][2] <= lc) {
            xcount += 1;
            xsave.push_back(face2[i][1]);
        }
    }
    //cout << "Part 3 xcount " << xcount << endl;
    if (face2[0][1] < -lc) {
        xcount += 1;
        xsave.push_back(-lc);
    }
    //cout << "Part 4 xcount " << xcount << endl;

    if (face2[3][1] > lc) {
        xcount++;
        xsave.push_back(lc);
    }
    //cout << "Part 5 xcount " << xcount << endl;
    vector<double> xsort = xsave;
    // Print all rows in the sorted vector
    //std::cout << "Sorted vector (xsort):" << std::endl;
    //for (const auto& value : xsort) {
    //    std::cout << value << " ";
    //}
    //std::cout << std::endl;

    sort(xsort.begin(), xsort.end());
    auto it = unique(xsort.begin(), xsort.end());
    xsort.resize(std::distance(xsort.begin(), it));

    vector<Region> region;
    int regioncount = 0;

    for (size_t i1 = 1; i1 < xsort.size(); ++i1) {
        double RYLB = xsort[i1 - 1];
        double RYUB = xsort[i1];

        double slope, RZUB1, RZUB2, RZLB1, RZLB2;

        if (RYLB < face2[1][1]) {
            slope = (face2[1][2] - face2[0][2]) / (face2[1][1] - face2[0][1]);
            if (isinf(slope)) {
                RZUB1 = face2[1][2];
            } else {
                RZUB1 = slope * (RYLB - face2[1][1]) + face2[1][2];
            }
        } else {
            slope = (face2[3][2] - face2[1][2]) / (face2[3][1] - face2[1][1]);
            if (isinf(slope)) {
                RZUB1 = face2[1][2];
            } else {
                RZUB1 = slope * (RYLB - face2[1][1]) + face2[1][2];
            }
        }

        if (RYUB < face2[1][1]) {
            slope = (face2[1][2] - face2[0][2]) / (face2[1][1] - face2[0][1]);
            if (isinf(slope)) {
                RZUB2 = face2[1][2];
            } else {
                RZUB2 = slope * (RYUB - face2[1][1]) + face2[1][2];
            }
        } else {
            slope = (face2[3][2] - face2[1][2]) / (face2[3][1] - face2[1][1]);
            if (isinf(slope)) {
                RZUB2 = face2[1][2];
            } else {
                RZUB2 = slope * (RYUB - face2[1][1]) + face2[1][2];
            }
        }

        if (RYLB < face2[2][1]) {
            slope = (face2[2][2] - face2[0][2]) / (face2[2][1] - face2[0][1]);
            if (isinf(slope)) {
                RZLB1 = face2[2][2];
            } else {
                RZLB1 = slope * (RYLB - face2[2][1]) + face2[2][2];
            }
        } else {
            slope = (face2[3][2] - face2[2][2]) / (face2[3][1] - face2[2][1]);
            if (isinf(slope)) {
                RZLB1 = face2[2][2];
            } else {
                RZLB1 = slope * (RYLB - face2[2][1]) + face2[2][2];
            }
        }

        if (RYUB < face2[2][1]) {
            slope = (face2[2][2] - face2[0][2]) / (face2[2][1] - face2[0][1]);
            if (isinf(slope)) {
                RZLB2 = face2[2][2];
            } else {
                RZLB2 = slope * (RYUB - face2[2][1]) + face2[2][2];
            }
        } else {
            slope = (face2[3][2] - face2[2][2]) / (face2[3][1] - face2[2][1]);
            if (isinf(slope)) {
                RZLB2 = face2[2][2];
            } else {
                RZLB2 = slope * (RYUB - face2[2][1]) + face2[2][2];
            }
        }

        if (RZUB1 > lc) {
           RZUB1 = lc;
        }
        if (RZUB2 > lc) {
           RZUB2 = lc;
        }
        if (RZLB1 < -lc) {
           RZLB1 = -lc;
        }
        if (RZLB2 < -lc) {
           RZLB2 = -lc;
        }

        if (RZUB1 >= -reg_tol && RZUB2 >= -reg_tol && RZLB1 <= reg_tol && RZLB2 <= reg_tol) {
            Region newRegion;
            newRegion.points[0][0] = RYLB;
            newRegion.points[0][1] = RZUB1;
            newRegion.points[1][0] = RYUB;
            newRegion.points[1][1] = RZUB2;
            newRegion.points[2][0] = RYLB;
            newRegion.points[2][1] = RZLB1;
            newRegion.points[3][0] = RYUB;
            newRegion.points[3][1] = RZLB2;
            region.push_back(newRegion);
            regioncount++;
        }
    }

    //cout << "regioncount " << regioncount << endl;
    return region;
}

double calc_ver(double Particle2VerticesEdge[][3],double AASigma, double catr1, double atre1, double catr2, double atre2, double crep, double repe, double Epsilon, double lc) {

    vector<vector<int>> facever = {{0, 3, 2, 1}, {2, 3, 7, 6}, {4, 7, 6, 5}, {0, 4, 5, 1}, {1, 2, 6, 5}, {0, 3, 7, 4}};
    vector<vector<int>> fid2d = {{0, 3}, {0, 5}, {3, 5}, {0, 4}, {3, 4}, {0, 1}, {1, 4}, {1, 5}, {2, 3}, {2, 5}, {2, 4}, {1, 2}};
    vector<vector<int>> fid3d = {{0, 3, 5}, {0, 3, 4}, {0, 1, 4}, {0, 1, 5}, {2, 3, 5}, {2, 3, 4}, {1, 2, 4}, {1, 2, 5}};

    // Find the minimum value along the first column
    auto minRow = min_element(Particle2VerticesEdge, Particle2VerticesEdge+8,
        [](const auto& a, const auto& b) {
            return a[0] < b[0];
        });

    // Get the row index of the minimum value
    int minindex = distance(Particle2VerticesEdge, minRow);
    //cout << "minindex " << minindex << endl;
    int mincount = 0;
    double minval = Particle2VerticesEdge[minindex][0];
    //cout << "minval " << minval << endl;
    double minedgeY = Particle2VerticesEdge[minindex][1];
    double minedgeZ = Particle2VerticesEdge[minindex][2];
    vector<int> minind;

    for (int i = 0; i < 8; ++i) {
       if (abs(Particle2VerticesEdge[i][0]-minval) < 0.01) {
          mincount++;
          minind.push_back(i);
       }
    }
    //cout << "comx" << comx << endl;
    //cout << "mincount " << mincount << endl;
    //cout << "mindind" << endl;
    //for (int i = 0; i < minind.size(); ++i) {
    //    std::cout << minind[i] << " ";
    //}
    cout << endl;

    double comx = 0.0;
    double comy = 0.0;
    double comz = 0.0;

    double ver[8][3];
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          ver[i][j] = Particle2VerticesEdge[i][j] / AASigma;
          if (j == 0) {
             ver[i][j] = ver[i][j] - minval/AASigma+lc;
          }
       }
    }

    double utot, uval;
    if (mincount == 4) {
       int rcount = 0;
       double tempv[4][3];

       for (int i = 0; i < mincount; ++i) {
          copy(ver[minind[i]], ver[minind[i]] + 3, tempv[i]);
       }

       /*cout <<"tempv" << endl;
       for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 3; ++j) {
             std::cout << tempv[i][j] << " ";
          }
          std::cout << std::endl;
       }*/

       double facevs[4][3];
       int column = 2; // third column
       vector<int> sortedIndices = sortIndicesByDesiredColumn(tempv,column);
       // Assign rows based on the sorted indices
       for (int i = 0; i < mincount; ++i) {
          int originalIndex = sortedIndices[i];
          copy(begin(tempv[originalIndex]), end(tempv[originalIndex]), begin(facevs[i]));
       }
       // Assign values to facev based on the sorted tempv
       double facev[4][3];

       if (facevs[2][1] > facevs[3][1]) {
          copy(facevs[3], facevs[3] + 3, facev[0]);
          copy(facevs[2], facevs[2] + 3, facev[1]);
       } else {
          copy(facevs[2], facevs[2] + 3, facev[0]);
          copy(facevs[3], facevs[3] + 3, facev[1]);
       }
       if (facevs[0][1] > facevs[1][1]) {
          copy(facevs[0], facevs[0] + 3, facev[2]);
          copy(facevs[1], facevs[1] + 3, facev[3]);
       } else {
          copy(facevs[1], facevs[1] + 3, facev[2]);
          copy(facevs[0], facevs[0] + 3, facev[3]);
       }

       /*cout <<"facev" << endl;
       for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 3; ++j) {
             std::cout << facev[i][j] << " ";
          }
          std::cout << std::endl;
       }*/

       utot = 0;
       double tempface[4][3];
       double tempregion[4][2];

       vector<Region> region = region_calc(facev, lc);
       for (size_t i = 0; i < region.size(); ++i) {
          for (int j = 0; j < 4; ++j) {
             for (int k = 0; k < 3; ++k) {
                tempface[j][k] = facev[j][k];
                if (k == 0) {
                   tempface[j][k] += 1 + minval/AASigma - lc;
                }
             }
             for (int m = 0; m < 2; ++m) {
                tempregion[j][m] = region[i].points[j][m];
             }
          }


          /*std::cout << "tempface:" << std::endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 3; ++j) {
                std::cout << tempface[i][j] << " ";
             }
             std::cout << std::endl;
          }

          // Printing tempregion
          std::cout << "tempregion:" << std::endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 2; ++j) {
                std::cout << tempregion[i][j] << " ";
             }
             std::cout << std::endl;
          }*/

          uval = calc_gen(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,tempface,tempregion,Epsilon);
          cout << "uval " << uval << endl;
          utot += uval;
       }
    } else if (mincount == 2) {
       double rcount = 0;
       utot = 0;
       int fid;
       if (minind[0] == 0) {
          if (minind[1] == 1) {
             fid = 0;
          } else if (minind[1] == 3) {
             fid = 1;
          } else if (minind[1] == 4) {
             fid = 2;
          }
       } else if (minind[0] == 1) {
          if (minind[1] == 2) {
             fid = 3;
          } else if (minind[1] == 5) {
             fid = 4;
          }
       } else if (minind[0] == 2) {
          if (minind[1] == 3) {
             fid = 5;
          } else if (minind[1] == 6) {
             fid = 6;
          }
       } else if (minind[0] == 3) {
          if (minind[1] == 7) {
             fid = 7;
          }
       } else if (minind[0] == 4) {
          if (minind[1] == 5) {
             fid = 8;
          } else if (minind[1] == 7) {
             fid = 9;
          }
       } else if (minind[0] == 5) {
          if (minind[1] == 6) {
             fid = 10;
          }
       } else if (minind[0] == 6) {
          if (minind[1] == 7) {
             fid = 11;
          }
       }
       double tempv[4][3];
       for (int i1 = 0; i1 < 2; ++i1) {
          for (int i2 = 0; i2 < 4; ++i2) {
             int row = fid2d[fid][i1];
             int faceIndex = facever[row][i2];
             for (int i3 = 0; i3 < 3; ++i3) {
                tempv[i2][i3] = ver[faceIndex][i3];
             }
          }
          /*cout <<"tempv" << endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 3; ++j) {
                std::cout << tempv[i][j] << " ";
             }
             std::cout << std::endl;
          }*/
          double facevs[4][3];
          int column = 2; // third column
          vector<int> sortedIndices = sortIndicesByDesiredColumn(tempv,column);
          for (int i = 0; i < 4; ++i) {
             int originalIndex = sortedIndices[i];
             copy(begin(tempv[originalIndex]), end(tempv[originalIndex]), begin(facevs[i]));
          }

          double facev[4][3];
          if (facevs[2][1] > facevs[3][1]) {
             facev[0][0] = facevs[3][0];
             facev[0][1] = facevs[3][1];
             facev[0][2] = facevs[3][2];

             facev[1][0] = facevs[2][0];
             facev[1][1] = facevs[2][1];
             facev[1][2] = facevs[2][2];
          } else {
             facev[0][0] = facevs[2][0];
             facev[0][1] = facevs[2][1];
             facev[0][2] = facevs[2][2];

             facev[1][0] = facevs[3][0];
             facev[1][1] = facevs[3][1];
             facev[1][2] = facevs[3][2];
          }

          if (facevs[0][1] > facevs[1][1]) {
             facev[2][0] = facevs[0][0];
             facev[2][1] = facevs[0][1];
             facev[2][2] = facevs[0][2];

             facev[3][0] = facevs[1][0];
             facev[3][1] = facevs[1][1];
             facev[3][2] = facevs[1][2];
          } else {
             facev[2][0] = facevs[1][0];
             facev[2][1] = facevs[1][1];
             facev[2][2] = facevs[1][2];

             facev[3][0] = facevs[0][0];
             facev[3][1] = facevs[0][1];
             facev[3][2] = facevs[0][2];
          }
          /*cout <<"facev" << endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 3; ++j) {
                std::cout << facev[i][j] << " ";
             }
             std::cout << std::endl;
          }*/
          vector<Region> region = region_calc(facev, lc);
          double tempface[4][3];
          double tempregion[4][2];

          for (size_t i = 0; i < region.size(); ++i) {
             for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 3; ++k) {
                   tempface[j][k] = facev[j][k];
                   if (k == 0) {
                      tempface[j][k] += 1 + minval/AASigma - lc;
                   }
                }
                for (int m = 0; m < 2; ++m) {
                   tempregion[j][m] = region[i].points[j][m];
                }
             }

             /*std::cout << "tempface:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 3; ++j) {
                   std::cout << tempface[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

          // Printing tempregion
             /*std::cout << "tempregion:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 2; ++j) {
                   std::cout << tempregion[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

             uval = calc_gen(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,tempface,tempregion,Epsilon);
             cout << "uval " << uval << endl;
             utot += uval;
          }
       }
    }  else if (mincount == 1) {
       double rcount = 0;
       utot = 0;
       double tempv[4][3];

       for (int i1 = 0; i1 < 3; ++i1) {
          for (int i2 = 0; i2 < 4; ++i2) {
             for (int j = 0; j < 3; ++j) {
                int row = fid3d[minindex][i1];
                int faceindex = facever[row][i2];
                tempv[i2][j] = ver[faceindex][j];
             }
          }
          /*cout <<"tempv" << endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 3; ++j) {
                std::cout << tempv[i][j] << " ";
             }
             std::cout << std::endl;
          }*/
          double facevs[4][3];
          int column = 2; // third column
          vector<int> sortedIndices = sortIndicesByDesiredColumn(tempv,column);
          for (int i = 0; i < 4; ++i) {
             int originalIndex = sortedIndices[i];
             copy(begin(tempv[originalIndex]), end(tempv[originalIndex]), begin(facevs[i]));
          }

          double facev[4][3];
          if (facevs[2][1] > facevs[3][1]) {
             facev[0][0] = facevs[3][0];
             facev[0][1] = facevs[3][1];
             facev[0][2] = facevs[3][2];

             facev[1][0] = facevs[2][0];
             facev[1][1] = facevs[2][1];
             facev[1][2] = facevs[2][2];
          } else {
             facev[0][0] = facevs[2][0];
             facev[0][1] = facevs[2][1];
             facev[0][2] = facevs[2][2];

             facev[1][0] = facevs[3][0];
             facev[1][1] = facevs[3][1];
             facev[1][2] = facevs[3][2];
          }

          if (facevs[0][1] > facevs[1][1]) {
             facev[2][0] = facevs[0][0];
             facev[2][1] = facevs[0][1];
             facev[2][2] = facevs[0][2];

             facev[3][0] = facevs[1][0];
             facev[3][1] = facevs[1][1];
             facev[3][2] = facevs[1][2];
          } else {
             facev[2][0] = facevs[1][0];
             facev[2][1] = facevs[1][1];
             facev[2][2] = facevs[1][2];

             facev[3][0] = facevs[0][0];
             facev[3][1] = facevs[0][1];
             facev[3][2] = facevs[0][2];
          }
          /*cout <<"facev" << endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 3; ++j) {
                std::cout << facev[i][j] << " ";
             }
             std::cout << std::endl;
          }*/
          vector<Region> region = region_calc(facev, lc);
          double tempface[4][3];
          double tempregion[4][2];

          for (size_t i = 0; i < region.size(); ++i) {
             for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 3; ++k) {
                   tempface[j][k] = facev[j][k];
                   if (k == 0) {
                      tempface[j][k] += 1 + minval/AASigma - lc;
                   }
                }
                for (int m = 0; m < 2; ++m) {
                   tempregion[j][m] = region[i].points[j][m];
                }
             }
             /*std::cout << "tempface:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 3; ++j) {
                   std::cout << tempface[i][j] << " ";
                }
                std::cout << std::endl;
             }

             // Printing tempregion
             std::cout << "tempregion:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 2; ++j) {
                   std::cout << tempregion[i][j] << " ";
                }
                std::cout << std::endl;
             }*/
             uval = calc_gen(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,tempface,tempregion,Epsilon);
             cout << "uval " << uval << endl;
             utot += uval;
          }
       }
    } else {
       utot = 0;
    }

    return utot;
}


void InitialConfiguration(int NumberOfParticles, double BoxLength, double AASigma,
                          unordered_map<double, double>& SizeRatioDictionary,
                          vector<double>& Sizes, vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,
                          double (*VectorX)[3], double (*VectorY)[3], double (*VectorZ)[3],
                          int& restartStep, ofstream& LAMMPSTrajectoryFile, ofstream& EnergyFile, ofstream& TimeFile) {

    srand(11);

    vector<double> keys_list;
    vector<double> values_list;

    for (const auto& entry : SizeRatioDictionary) {
        keys_list.push_back(entry.first);
        values_list.push_back(entry.second);
    }

    double valuesSum = accumulate(values_list.begin(), values_list.end(), 0.0);
    int ParticleID = 0;
    vector<int> CumulativeNumberOfParticlesList;

    for (double size : keys_list) {
        int NumberOfCurrentSize = round(SizeRatioDictionary[size] / valuesSum * NumberOfParticles);
        CumulativeNumberOfParticlesList.push_back(NumberOfCurrentSize);
    }

    for (size_t i = 0; i < CumulativeNumberOfParticlesList.size(); ++i) {
        if (ParticleID < NumberOfParticles) {
            int NumberOfCurrentSize = CumulativeNumberOfParticlesList[i];
            for (int j = 0; j < NumberOfCurrentSize; ++j) {
                Sizes[ParticleID] = keys_list[i];
                ParticleID++;
            }
        } else {
            break; // Exit the loop if ParticleID exceeds NumberOfParticles
        }
    }

    double BoxLengthHalf = BoxLength / 2.0;
    double SigmaCut = 2 * AASigma;
    double SigmaCutSquare = pow(SigmaCut, 2);

    double FirstPositionX = ((double)rand()/(RAND_MAX));
    double FirstPositionY = ((double)rand()/(RAND_MAX));
    double FirstPositionZ = ((double)rand()/(RAND_MAX));
    double FirstCubeAlpha = ((double)rand()/(RAND_MAX)) * 90;
    double FirstCubeBeta = ((double)rand()/(RAND_MAX)) * 90;
    double FirstCubeGamma = ((double)rand()/(RAND_MAX)) * 90;

    for (int i = 0; i < NumberOfParticles; ++i) {
        VectorX[i][0] = 1.0;
        VectorX[i][1] = 0.0;
        VectorX[i][2] = 0.0;
        VectorY[i][0] = 0.0;
        VectorY[i][1] = 1.0;
        VectorY[i][2] = 0.0;
        VectorZ[i][0] = 0.0;
        VectorZ[i][1] = 0.0;
        VectorZ[i][2] = 1.0;
    }

    Vector3 RotationCenter_To_CubeCentroid;
    Vector3 FirstCubeVectors[3];
    Vector3 CubeVectors[3];

    CubeRotation({FirstPositionX, FirstPositionY, FirstPositionZ},
             VectorX[0],VectorY[0],VectorZ[0],
             {FirstPositionX, FirstPositionY, FirstPositionZ},
             {FirstCubeAlpha, FirstCubeBeta, FirstCubeGamma},
             RotationCenter_To_CubeCentroid,
             FirstCubeVectors);

    for (int i = 0; i < 3; ++i) {
        VectorX[0][i] = FirstCubeVectors[0][i];
        VectorY[0][i] = FirstCubeVectors[1][i];
        VectorZ[0][i] = FirstCubeVectors[2][i];
    }

    Rx[0] = BoxLength * FirstPositionX - BoxLengthHalf;
    Ry[0] = BoxLength * FirstPositionY - BoxLengthHalf;
    Rz[0] = BoxLength * FirstPositionZ - BoxLengthHalf;

    for (int i = 1; i < NumberOfParticles; ++i) {
        std::cout << "Insertion Of Molecule " << i + 1 << " Successful\n";
        bool Repeat = true;
        while (Repeat) {
            Repeat = false;
            double PositionX = BoxLength * ((double)rand()/(RAND_MAX)) - BoxLengthHalf;
            double PositionY = BoxLength * ((double)rand()/(RAND_MAX)) - BoxLengthHalf;
            double PositionZ = BoxLength * ((double)rand()/(RAND_MAX)) - BoxLengthHalf;
            double Alpha = ((double)rand()/(RAND_MAX)) * 90;
            double Beta = ((double)rand()/(RAND_MAX)) * 90;
            double Gamma = ((double)rand()/(RAND_MAX)) * 90;

            CubeRotation({PositionX, PositionY, PositionZ},
                 VectorX[i], VectorY[i], VectorZ[i],
                 {PositionX, PositionY, PositionZ},
                 {Alpha, Beta, Gamma},
                 RotationCenter_To_CubeCentroid,
                 CubeVectors);

            double CubeSideLength1 = Sizes[i] * AASigma;

            for (int K = 0; K < i; ++K) {
                double CubeSideLength2 = Sizes[K] * AASigma;
                double DistanceFromAnOldAtomX = PositionX - Rx[K];
                double DistanceFromAnOldAtomY = PositionY - Ry[K];
                double DistanceFromAnOldAtomZ = PositionZ - Rz[K];

                DistanceFromAnOldAtomX -= BoxLength * round(DistanceFromAnOldAtomX / BoxLength);
                DistanceFromAnOldAtomY -= BoxLength * round(DistanceFromAnOldAtomY / BoxLength);
                DistanceFromAnOldAtomZ -= BoxLength * round(DistanceFromAnOldAtomZ / BoxLength);

                double DistanceSquare = pow(DistanceFromAnOldAtomX, 2)+pow(DistanceFromAnOldAtomY, 2)+pow(DistanceFromAnOldAtomZ, 2);

                if (DistanceSquare < pow(0.866 * (CubeSideLength1 + CubeSideLength2), 2)) {
                    Repeat = true;
                    break;
                }
            }

            for (int j = 0; j < 3; ++j) {
                VectorX[i][j] = CubeVectors[0][j];
                VectorY[i][j] = CubeVectors[1][j];
                VectorZ[i][j] = CubeVectors[2][j];
            }

            Rx[i] = PositionX;
            Ry[i] = PositionY;
            Rz[i] = PositionZ;
        }
    }
    std::cout << "Start from random configuration\n";
}


void findAtomPositionsInCube(double CubeSize, double Sigma, double COMX, double COMY, double COMZ, const Vector3 VectorX, const Vector3 VectorY, const Vector3 VectorZ, double atomList[][3]) {

    // Determine coordinates based on CubeSize
    double numCoordinates = (CubeSize - 1) / 2;

    // Fill the 2D array with coordinates
    int index = 0;
    for (double x = -numCoordinates; x < CubeSize/2; x += 1) {
        for (double y = -numCoordinates; y < CubeSize/2; y += 1) {
            for (double z = -numCoordinates; z < CubeSize/2; z += 1) {
                atomList[index][0] = x * Sigma;
                atomList[index][1] = y * Sigma;
                atomList[index][2] = z * Sigma;
                ++index;
            }
        }
    }

    double x, y, z;
    // Translate and rotate the coordinates
    for (int i = 0; i < index; ++i) {
        x = atomList[i][0];
        y = atomList[i][1];
        z = atomList[i][2];
        atomList[i][0] = COMX + 1.0 * (x * VectorX[0] + y * VectorY[0] + z * VectorZ[0]);
        atomList[i][1] = COMY + 1.0 * (x * VectorX[1] + y * VectorY[1] + z * VectorZ[1]);
        atomList[i][2] = COMZ + 1.0 * (x * VectorX[2] + y * VectorY[2] + z * VectorZ[2]);
    }
}

double calculate_angle_between_vectors(Vector3 vector_a, Vector3 vector_b){
    double dot_product = vector_a[0] * vector_b[0] + vector_a[1] * vector_b[1] + vector_a[2] * vector_b[2];
    double norm_a = sqrt(vector_a[0] * vector_a[0] + vector_a[1] * vector_a[1] + vector_a[2] * vector_a[2]);
    double norm_b = sqrt(vector_b[0] * vector_b[0] + vector_b[1] * vector_b[1] + vector_b[2] * vector_b[2]);
    double angle = acos(dot_product / (norm_a * norm_b));
    double angleDegrees = fmod(abs(round(angle * (180.0 / M_PI))), 360.0);
    return angleDegrees;
}

double wrapAngle(double angle) {
    if (angle <= 90 && angle >= 0) {
        // angle = angle;  // No change needed in this case
    } else if (angle > 90 && angle <= 180) {
        angle -= 90;
    } else if (angle > 180 && angle <= 270) {
        angle -= 180;
    } else if (angle > 270 && angle <= 360) {
        angle -= 270;
    } else if (angle > -90 && angle < 0) {
        angle = 90 + angle;
    } else if (angle > -180 && angle <= -90) {
        angle = 180 + angle;
    } else if (angle > -270 && angle <= -180) {
        angle = 270 + angle;
    } else if (angle > -360 && angle <= -270) {
        angle = 360 + angle;
    }
    return angle;
}

bool CheckOverlapForTwoCubes(const int i, const int j, double BoxLength, double AASigma, double CubeSideLength1, double CubeSideLength2,
                              const Vector3 Particle1Centroid, const Vector3 Particle2Centroid,
                              const Vector3 Particle1VectorX, const Vector3 Particle1VectorY, const Vector3 Particle1VectorZ,
                              const Vector3 Particle2VectorX, const Vector3 Particle2VectorY, const Vector3 Particle2VectorZ) {

    double Size1 = CubeSideLength1 / AASigma;
    double Size2 = CubeSideLength2 / AASigma;
    Vector3 ParticleC1, ParticleC2, VectorX_1, VectorY_1, VectorZ_1, VectorX_2, VectorY_2, VectorZ_2;
    for (int k = 0; k < 3; ++k) {
       ParticleC1[k] = Particle1Centroid[k];
       ParticleC2[k] = Particle2Centroid[k];
       VectorX_1[k] = Particle1VectorX[k];
       VectorY_1[k] = Particle1VectorY[k];
       VectorZ_1[k] = Particle1VectorZ[k];
       VectorX_2[k] = Particle2VectorX[k];
       VectorY_2[k] = Particle2VectorY[k];
       VectorZ_2[k] = Particle2VectorZ[k];
    }


    double cube1initial[8][3];
    findVerticesOfCube(CubeSideLength1,ParticleC1[0],ParticleC1[1],ParticleC1[2],VectorX_1,VectorY_1,VectorZ_1,cube1initial);
    double cube2initial[8][3];
    findVerticesOfCube(CubeSideLength2,ParticleC2[0],ParticleC2[1],ParticleC2[2],VectorX_2,VectorY_2,VectorZ_2,cube2initial);

    /*cout << "Vertex 1 initial" << endl;
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube1initial[i][j] << " ";
       }
       std::cout << std::endl;
    }
    cout << "Vertex 2 initial" << endl;
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube2initial[i][j] << " ";
       }
       std::cout << std::endl;
    }*/

    double cube2Vertex[8][3];
    Reorientation(i, j, CubeSideLength2, BoxLength, ParticleC1, ParticleC2, VectorX_1, VectorY_1, VectorZ_1, VectorX_2, VectorY_2, VectorZ_2, cube2Vertex);

    double Cube1X = ParticleC1[0];
    double Cube1Y = ParticleC1[1];
    double Cube1Z = ParticleC1[2];
    double Cube2X = ParticleC2[0];
    double Cube2Y = ParticleC2[1];
    double Cube2Z = ParticleC2[2];

    if (max(Size1, Size2) == 1) {
        // Single particle
        double center_distance_sq = pow((Cube1X - Cube2X), 2) + pow((Cube1Y - Cube2Y), 2) + pow((Cube1Z - Cube2Z), 2);
        double max_val_sq = pow((CubeSideLength1 + CubeSideLength2), 2) / 4.0;

        if (center_distance_sq > max_val_sq) {
            return false;
        }
    } else {
        double center_distance_sq = pow((Cube1X - Cube2X), 2) + pow((Cube1Y - Cube2Y), 2) + pow((Cube1Z - Cube2Z), 2);
        double max_val_sq = (pow((CubeSideLength1 + CubeSideLength2), 2) * 3.0) / 4.0;

        if (center_distance_sq > max_val_sq) {
            return false;
        }
        double cube1Vertex[8][3];
        findVerticesOfCube(CubeSideLength1,Cube1X,Cube1Y,Cube1Z,VectorX_1,VectorY_1,VectorZ_1,cube1Vertex);

        // Check if any vertex is inside the other cube
        double CubeSideLengthHalf1 = CubeSideLength1/2.0;
        for (int i = 0; i < 8; ++i) {
           double vertex_x = cube2Vertex[i][0];
           double vertex_y = cube2Vertex[i][1];
           double vertex_z = cube2Vertex[i][2];
           vertex_x = vertex_x-BoxLength*round(vertex_x/BoxLength);
           vertex_y = vertex_y-BoxLength*round(vertex_y/BoxLength);
           vertex_z = vertex_z-BoxLength*round(vertex_z/BoxLength);
           //cout << "vertex " << vertex_x << " " << vertex_y << " " << vertex_z << endl;
           if (-CubeSideLengthHalf1 <= vertex_x && vertex_x <= CubeSideLengthHalf1 && -CubeSideLengthHalf1 <= vertex_y && vertex_y <= CubeSideLengthHalf1 && -CubeSideLengthHalf1 <= vertex_z && vertex_z <= CubeSideLengthHalf1) {
              cout << "overlap" << endl;
              return true;
           }
        }
        // Check if any edge intersects with other cube interior
        vector<int> sortedIndices = sortIndicesByDesiredColumn(cube2Vertex,0);
        for (int i = 1; i < 4; ++i) {
           int min_ind = sortedIndices[0];
           int next_ind = sortedIndices[i];
           double slope_x = (cube2Vertex[next_ind][0]-cube2Vertex[min_ind][0])/Size2;
           double slope_y = (cube2Vertex[next_ind][1]-cube2Vertex[min_ind][1])/Size2;
           double slope_z = (cube2Vertex[next_ind][2]-cube2Vertex[min_ind][2])/Size2;
           for (int j = 0; j < Size2; ++j) {
              double line_x = slope_x*j+cube2Vertex[min_ind][0];
              double line_y = slope_y*j+cube2Vertex[min_ind][1];
              double line_z = slope_z*j+cube2Vertex[min_ind][2];
              if (-CubeSideLengthHalf1 <= line_x && line_x <= CubeSideLengthHalf1 && -CubeSideLengthHalf1 <= line_y && line_y <= CubeSideLengthHalf1 && -CubeSideLengthHalf1 <= line_z && line_z <= CubeSideLengthHalf1) {
                 cout << "overlap" << endl;
                 return true;
              }
           }
        }


        /*cout << "Vertex 1 final" << endl;
        for (int i = 0; i < 8; ++i) {
           for (int j = 0; j < 3; ++j) {
              std::cout << cube1Vertex[i][j] << " ";
           }
           std::cout << std::endl;
        }
        cout << "Vertex 2 final" << endl;
        for (int i = 0; i < 8; ++i) {
           for (int j = 0; j < 3; ++j) {
              std::cout << cube2Vertex[i][j] << " ";
           }
           std::cout << std::endl;
        }*/


        for (int i = 0; i < 3; ++i) {
            Vector3 axis_vec = {0.0, 0.0, 0.0};
            axis_vec[i] = 1.0;
            double max_proj1 = 0.0;
            double max_proj2 = 0.0;
            double min_proj1 = 1000000000.0; //very big number to initialize
            double min_proj2 = 1000000000.0; //very big number to initialize
            for (int j = 0; j < 8; ++j) {
                double proj1 = cube1Vertex[j][0]*axis_vec[0]+cube1Vertex[j][1]*axis_vec[1]+cube1Vertex[j][2]*axis_vec[2];
                double proj2 = cube2Vertex[j][0]*axis_vec[0]+cube2Vertex[j][1]*axis_vec[1]+cube2Vertex[j][2]*axis_vec[2];
                if (proj1 < min_proj1) {
                   min_proj1 = proj1;
                }
                if (proj2 < min_proj2) {
                   min_proj2 = proj2;
                }
                if (proj1 > max_proj1) {
                   max_proj1 = proj1;
                }
                if (proj2 > max_proj2) {
                   max_proj2 = proj2;
                }
            }
            //cout << axis_vec[0] << " " << axis_vec[1] << " " << axis_vec[2] << endl;
            //cout << min_proj1 << " " << min_proj2 << " " << max_proj1 << " " << max_proj2 << endl;
            if (max_proj1 < min_proj2 || max_proj2 < min_proj1) {
               return false;
            }
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (i != j) {
                   Vector3 axis_vec = {0.0, 0.0, 0.0};
                   axis_vec[i] = 1.0;
                   axis_vec[j] = 1.0;
                   double max_proj1 = 0.0;
                   double max_proj2 = 0.0;
                   double min_proj1 = 1000000000.0; //very big number to initialize
                   double min_proj2 = 1000000000.0; //very big number to initialize
                   for (int j = 0; j < 8; ++j) {
                       double proj1 = cube1Vertex[j][0]*axis_vec[0]+cube1Vertex[j][1]*axis_vec[1]+cube1Vertex[j][2]*axis_vec[2];
                       double proj2 = cube2Vertex[j][0]*axis_vec[0]+cube2Vertex[j][1]*axis_vec[1]+cube2Vertex[j][2]*axis_vec[2];
                       if (proj1 < min_proj1) {
                          min_proj1 = proj1;
                       }
                       if (proj2 < min_proj2) {
                          min_proj2 = proj2;
                       }
                       if (proj1 > max_proj1) {
                          max_proj1 = proj1;
                       }
                       if (proj2 > max_proj2) {
                          max_proj2 = proj2;
                       }
                   }
                   //cout << axis_vec[0] << " " << axis_vec[1] << " " << axis_vec[2] << endl;
                   //cout << min_proj1 << " " << min_proj2 << " " << max_proj1 << " " << max_proj2 << endl;
                   if (max_proj1 < min_proj2 || max_proj2 < min_proj1) {
                      return false;
                   }
                }
            }
        }
    }
    return true;
}


double EnergyBetweenTwoParticles(string model, double AASigma, double CGSigma, double CGEpsilon, double cutoffCG, double Hamaker, double AtomDensity, int Size_i, int Size_j, double BoxLength, const Vector3 ParticleiCentroid, const Vector3 ParticlejCentroid, const Vector3 VectorX_i, const Vector3 VectorY_i, const Vector3 VectorZ_i, const Vector3 VectorX_j, const Vector3 VectorY_j, const Vector3 VectorZ_j) {
    double CubeCubeEnergy = 0.0;
    if (model == "AACG") {
        // AA/CG LJ energy calculation
        double COMX, COMY, COMZ;
        COMX = ParticleiCentroid[0];
        COMY = ParticleiCentroid[1];
        COMZ = ParticleiCentroid[2];
        int NumberofRows = pow(Size_i * AASigma / CGSigma, 3);
        int CGCubeSize1 = Size_i * AASigma / CGSigma;
        double atomList1[NumberofRows][3];
        findAtomPositionsInCube(CGCubeSize1, CGSigma, COMX, COMY, COMZ, VectorX_i, VectorY_i, VectorZ_i, atomList1);

        COMX = ParticlejCentroid[0];
        COMY = ParticlejCentroid[1];
        COMZ = ParticlejCentroid[2];
        NumberofRows = pow(Size_j * AASigma / CGSigma, 3);
        int CGCubeSize2 = Size_j * AASigma / CGSigma;
        double atomList2[NumberofRows][3];
        findAtomPositionsInCube(CGCubeSize2, CGSigma, COMX, COMY, COMZ, VectorX_j, VectorY_j, VectorZ_j, atomList2);

        // Energy
        CubeCubeEnergy = LJPotentialBetweenTwoCubes(CGSigma, CGEpsilon, cutoffCG, atomList1, CGCubeSize1, atomList2, CGCubeSize2, BoxLength);
    } else if (model == "vdW") {
        // vdW energy calculation
        int i = 0;
        int j = 1;

        double Particle1Centroid[3] = {ParticleiCentroid[0], ParticleiCentroid[1], ParticleiCentroid[2]};
        double Particle2Centroid[3] = {ParticlejCentroid[0], ParticlejCentroid[1], ParticlejCentroid[2]};
        Vector3 VectorX_1,VectorX_2,VectorY_1,VectorY_2,VectorZ_1,VectorZ_2;

        for (int k = 0; k < 3; ++k) {
           VectorX_1[k] = VectorX_i[k];
           VectorY_1[k] = VectorY_i[k];
           VectorZ_1[k] = VectorZ_i[k];
           VectorX_2[k] = VectorX_j[k];
           VectorY_2[k] = VectorY_j[k];
           VectorZ_2[k] = VectorZ_j[k];
        }
        //cout << Particle1Centroid[0] << " " << Particle1Centroid[1] << " " << Particle1Centroid[2] << endl;
        //cout << Particle2Centroid[0] << " " << Particle2Centroid[1] << " " << Particle2Centroid[2] << endl;
        double Particle2VerticesEdge[8][3];
        double Cube2SideLength = Size_j * AASigma;

        Reorientation(i,j,Cube2SideLength,BoxLength,Particle1Centroid,Particle2Centroid,VectorX_1,VectorY_1,VectorZ_1,VectorX_2,VectorY_2,VectorZ_2,Particle2VerticesEdge);
        /*
        cout << "Reorientation" << endl;
        cout << "Particle2Centroid[0] " << Particle2Centroid[0] << " Particle2Centroid[1] " << Particle2Centroid[1] << " Particle2Centroid[2] " << Particle2Centroid[2] << endl;
        cout << "VectorX_1[0] " << VectorX_1[0] << " VectorX_1[1] " << VectorX_1[1] << " VectorX_1[2] " << VectorX_1[2] << endl;
        cout << "VectorY_1[0] " << VectorY_1[0] << " VectorY_1[1] " << VectorY_1[1] << " VectorY_1[2] " << VectorY_1[2] << endl;
        cout << "VectorZ_1[0] " << VectorZ_1[0] << " VectorZ_1[1] " << VectorZ_1[1] << " VectorZ_1[2] " << VectorZ_1[2] << endl;
        cout << "VectorX_2[0] " << VectorX_2[0] << " VectorX_2[1] " << VectorX_2[1] << " VectorX_2[2] " << VectorX_2[2] << endl;
        cout << "VectorY_2[0] " << VectorY_2[0] << " VectorY_2[1] " << VectorY_2[1] << " VectorY_2[2] " << VectorY_2[2] << endl;
        cout << "VectorZ_2[0] " << VectorZ_2[0] << " VectorZ_2[1] " << VectorZ_2[1] << " VectorZ_2[2] " << VectorZ_2[2] << endl;
        */
        /*for (int i = 0; i < 8; ++i) {
           for (int j = 0; j < 3; ++j) {
              std::cout << Particle2VerticesEdge[i][j] << " ";
           }
           std::cout << std::endl;
        }*/

        // Parameters
        double catr1 = -7.75596;
        double atre1 = 3.4339;
        double catr2 = -4.65093;
        double atre2 = 2.7167;
        double crep = 3.85898;
        double repe = 10.66482;
        double lc = Size_j / 2;

        double Epsilon = (Hamaker * (pow(10, -19)) / (4 * M_PI * M_PI * pow(AtomDensity, 2) * pow(AASigma, 6) * 6.9477 * 2 * (pow(10, -21))));
        // cout << "Epsilon " << Epsilon << endl;
        // Call function to calculate potential energy between two particles
        CubeCubeEnergy = calc_ver(Particle2VerticesEdge, AASigma, catr1, atre1, catr2, atre2, crep, repe, Epsilon, lc);
        if (std::isinf(std::abs(CubeCubeEnergy))) {
            CubeCubeEnergy = 0.0;
        }
    }

    return CubeCubeEnergy;
}

pair<double, bool> OneParticleEnergy(string model, int i, double Size_i, double cutoffCG, double Hamaker, double AtomDensity, const double Rx_i, const double Ry_i, const double Rz_i, 
                                     vector<double>& Sizes, vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,const double VectorX[][3], const double VectorY[][3], const double VectorZ[][3],
                         int NumberOfParticles, double BoxLength, double CGEpsilon, double CGSigma, double AASigma, vector<int>& clusterParticleList){
   // Calculate energy of one cube
   bool overLapFlag = false;
   double CurrentAtomTotalPotentialEnergy = 0;
   double CurrentPairPotentialEnergy = 0;
   Vector3 VectorX_i,VectorY_i,VectorZ_i, VectorX_j, VectorY_j, VectorZ_j;

   for (int j = 0; j < NumberOfParticles; ++j) {
       if (j != i && find(clusterParticleList.begin(), clusterParticleList.end(), j) == clusterParticleList.end()) {
          Size_i = Sizes[i];
          cout << i << " and " << j << " investigation" << endl;
          double Size_j = Sizes[j];
          for (int k = 0; k < 3; ++k) {
             VectorX_i[k] = VectorX[i][k];
             VectorY_i[k] = VectorY[i][k];
             VectorZ_i[k] = VectorZ[i][k];
          }
          VectorX_j[0] = VectorX[j][0];VectorX_j[1] = VectorX[j][1];VectorX_j[2] = VectorX[j][2];
          VectorY_j[0] = VectorY[j][0];VectorY_j[1] = VectorY[j][1];VectorY_j[2] = VectorY[j][2];
          VectorZ_j[0] = VectorZ[j][0];VectorZ_j[1] = VectorZ[j][1];VectorZ_j[2] = VectorZ[j][2];
          double Rx_j = Rx[j];
          double Rx_ij = Rx_i - Rx_j;
          double Ry_j = Ry[j];
          double Ry_ij = Ry_i - Ry_j;
          double Rz_j = Rz[j];
          double Rz_ij = Rz_i - Rz_j;
          /*cout << "i " << i << " and j " << j << endl;
          cout << "VectorX_i[0] " << VectorX_i[0] << " VectorX_i[1] " << VectorX_i[1] << " VectorX_i[2] " << VectorX_i[2] << endl;
          cout << "VectorY_i[0] " << VectorY_i[0] << " VectorY_i[1] " << VectorY_i[1] << " VectorY_i[2] " << VectorY_i[2] << endl;
          cout << "VectorZ_i[0] " << VectorZ_i[0] << " VectorZ_i[1] " << VectorZ_i[1] << " VectorZ_i[2] " << VectorZ_i[2] << endl;
          cout << "VectorX_j[0] " << VectorX_j[0] << " VectorX_j[1] " << VectorX_j[1] << " VectorX_j[2] " << VectorX_j[2] << endl;
          cout << "VectorY_j[0] " << VectorY_j[0] << " VectorY_j[1] " << VectorY_j[1] << " VectorY_j[2] " << VectorY_j[2] << endl;
          cout <<"VectorZ_j[0] " << VectorZ_j[0] << " VectorZ_j[1] " << VectorZ_j[1] << " VectorZ_j[2] " << VectorZ_j[2] << endl;
          cout << "Rxi " << Rx_i << " Ryi " << Ry_i << " Rzi " << Rz_i << endl;
          cout << "Rxj " << Rx_j << " Ryj " << Ry_j << " Rzj " << Rz_j << endl;*/
          Rx_ij = Rx_ij - BoxLength * round(Rx_ij / BoxLength);
          Ry_ij = Ry_ij - BoxLength * round(Ry_ij / BoxLength);
          Rz_ij = Rz_ij - BoxLength * round(Rz_ij / BoxLength);
          double RijSquare = Rx_ij * Rx_ij + Ry_ij * Ry_ij + Rz_ij * Rz_ij;
          Vector3 ParticleiCentroid = {0.0, 0.0, 0.0};
          Vector3 ParticlejCentroid = {-Rx_ij, -Ry_ij, -Rz_ij};
          double MaxSize = max(Size_i,Size_j);
          double CutOff =  2.5 * MaxSize * CGSigma;
          double CutOffSquare = CutOff * CutOff;
          if (RijSquare < CutOffSquare) {
             double CubeSideLength1 = Size_i*AASigma;
             double CubeSideLength2 = Size_j*AASigma;
             if (CheckOverlapForTwoCubes(i,j,BoxLength,AASigma,CubeSideLength1,CubeSideLength2,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j)) {
                CurrentPairPotentialEnergy = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j);
                overLapFlag = true;
                cout << i << " and " << j << " overlap" << endl;
             }else {
                CurrentPairPotentialEnergy = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j);
                if (CurrentPairPotentialEnergy > 0 || isnan(CurrentPairPotentialEnergy)) { // as a measure in case our overlapping algorithm did not detect the overlap
                   overLapFlag = true;
                }
             }
             CurrentAtomTotalPotentialEnergy += CurrentPairPotentialEnergy;
             cout << "CurrentPairPotentialEnergy " << CurrentPairPotentialEnergy << endl;
          }
       }
   }
   // cout << "CurrentAtomTotalPotentialEnergy " << CurrentAtomTotalPotentialEnergy << endl;
   return make_pair(CurrentAtomTotalPotentialEnergy, overLapFlag);
}

double TotalEnergy(vector<double> Sizes, vector<double> Rx, vector<double> Ry, vector<double> Rz, double VectorX[][3],double VectorY[][3],double VectorZ[][3], int NumberOfParticles, double BoxLength, double CGSigma, double CGEpsilon, double AASigma, string model, double cutoffCG, double Hamaker, double AtomDensity){
   //Calculate total energy of the system
   Vector3 VectorX_i,VectorX_j,VectorY_i,VectorY_j,VectorZ_i,VectorZ_j;
   double CurrentPairPotentialEnergy = 0;
   double CurrentAtomTotalPotentialEnergy = 0;
   bool overLapFlag = false;

   for (int i = 0; i < NumberOfParticles; ++i) {
       for (int j = i+1; j < NumberOfParticles; ++j) {
           double Size_i = Sizes[i];
           double Size_j = Sizes[j];
           cout << i << " and " << j << " investigated" << endl;
           VectorX_i[0] = VectorX[i][0];VectorX_i[1] = VectorX[i][1];VectorX_i[2] = VectorX[i][2];
           VectorY_i[0] = VectorY[i][0];VectorY_i[1] = VectorY[i][1];VectorY_i[2] = VectorY[i][2];
           VectorZ_i[0] = VectorZ[i][0];VectorZ_i[1] = VectorZ[i][1];VectorZ_i[2] = VectorZ[i][2];

           VectorX_j[0] = VectorX[j][0];VectorX_j[1] = VectorX[j][1];VectorX_j[2] = VectorX[j][2];
           VectorY_j[0] = VectorY[j][0];VectorY_j[1] = VectorY[j][1];VectorY_j[2] = VectorY[j][2];
           VectorZ_j[0] = VectorZ[j][0];VectorZ_j[1] = VectorZ[j][1];VectorZ_j[2] = VectorZ[j][2];

           double Rx_ij = Rx[i] - Rx[j];
           double Ry_ij = Ry[i] - Ry[j];
           double Rz_ij = Rz[i] - Rz[j];
           Rx_ij = Rx_ij - BoxLength * round(Rx_ij / BoxLength);
           Ry_ij = Ry_ij - BoxLength * round(Ry_ij / BoxLength);
           Rz_ij = Rz_ij - BoxLength * round(Rz_ij / BoxLength);
           //cout << "Total Energy" << endl;
           //cout << "i " << i << " and j " << j << endl;
           //cout << "VectorX_i[0] " << VectorX_i[0] << " VectorX_i[1] " << VectorX_i[1] << " VectorX_i[2] " << VectorX_i[2] << endl;
           //cout << "VectorY_i[0] " << VectorY_i[0] << " VectorY_i[1] " << VectorY_i[1] << " VectorY_i[2] " << VectorY_i[2] << endl;
           //cout << "VectorZ_i[0] " << VectorZ_i[0] << " VectorZ_i[1] " << VectorZ_i[1] << " VectorZ_i[2] " << VectorZ_i[2] << endl;
           //cout << "VectorX_j[0] " << VectorX_j[0] << " VectorX_j[1] " << VectorX_j[1] << " VectorX_j[2] " << VectorX_j[2] << endl;
           //cout << "VectorY_j[0] " << VectorY_j[0] << " VectorY_j[1] " << VectorY_j[1] << " VectorY_j[2] " << VectorY_j[2] << endl;
           //cout << "VectorZ_j[0] " << VectorZ_j[0] << " VectorZ_j[1] " << VectorZ_j[1] << " VectorZ_j[2] " << VectorZ_j[2] << endl;
           //cout << "Rxj " << Rx[j] << " Ryj " << Ry[j] << " Rzj " << Rz[j] << endl;
           //cout << "Rxi " << Rx[i] << " Ryi " << Ry[i] << " Rzi " << Rz[i] << endl;
           //cout << "Rxj " << Rx[j] << " Ryj " << Ry[j] << " Rzj " << Rz[j] << endl;
           double RijSquare = Rx_ij * Rx_ij + Ry_ij * Ry_ij + Rz_ij * Rz_ij;
           Vector3 ParticleiCentroid = {0.0, 0.0, 0.0};
           Vector3 ParticlejCentroid = {-Rx_ij, -Ry_ij, -Rz_ij};
           //cout << "Here totatl energy part" << endl;
           //cout << ParticlejCentroid[0] << " " << ParticlejCentroid[1] << " " << ParticlejCentroid[2] << endl;
           double MaxSize = max(Size_i,Size_j);
           double CutOff =  2.5 * MaxSize * CGSigma;
           double CutOffSquare = CutOff * CutOff;
           if (RijSquare < CutOffSquare) {
              double CubeSideLength1 = Size_i*AASigma;
              double CubeSideLength2 = Size_j*AASigma;
              if (CheckOverlapForTwoCubes(i,j,BoxLength,AASigma,CubeSideLength1,CubeSideLength2,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j)) {
                 CurrentPairPotentialEnergy = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j);
                 overLapFlag = true;
                 cout << i << " and " << j << " overlap" << endl;
              }else {
                 cout << ParticlejCentroid[0] << " " << ParticlejCentroid[1] << " " << ParticlejCentroid[2] << endl;
                 CurrentPairPotentialEnergy = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j);
              }
              CurrentAtomTotalPotentialEnergy += CurrentPairPotentialEnergy;
              cout << "CurrentPairPotentialEnergy " << CurrentPairPotentialEnergy << endl;
           }
       }
   }
   //cout << "CurrentAtomTotalPotentialEnergy" << endl;

   return CurrentAtomTotalPotentialEnergy;
}

void apply_periodic_boundary(Vector3 position, double BoxLength, Vector3 result){
   for (int i = 0; i < 3; ++i) {
       result[i] = position[i] - BoxLength * round(position[i]/ BoxLength);
   }
}


void calculate_cluster_centroid(vector<array<double, 3>> particles, double BoxLength, Vector3 centroid, double ClusterTemp[][3]){
   size_t numRows = particles.size();
   Vector3 result;
   Vector3 convertedElement;
   array<double, 3> element;

   for (size_t i = 0; i < numRows; ++i) {
       element = particles[i];
       copy(element.begin(), element.end(), convertedElement);
       apply_periodic_boundary(convertedElement, BoxLength, result);
       ClusterTemp[i][0] = result[0];
       ClusterTemp[i][1] = result[1];
       ClusterTemp[i][2] = result[2];
   }

   double X_sum = ClusterTemp[0][0];
   double Y_sum = ClusterTemp[0][1];
   double Z_sum = ClusterTemp[0][2];
   centroid[0] = X_sum;
   centroid[1] = Y_sum;
   centroid[2] = Z_sum;
   double X_dist, Y_dist, Z_dist;

   for (int i = 1; i < numRows; ++i) {
       X_dist = ClusterTemp[i][0]-centroid[0];
       Y_dist = ClusterTemp[i][1]-centroid[1];
       Z_dist = ClusterTemp[i][2]-centroid[2];
       X_dist = X_dist - BoxLength * round(X_dist / BoxLength);
       Y_dist = Y_dist - BoxLength * round(Y_dist / BoxLength);
       Z_dist = Z_dist - BoxLength * round(Z_dist / BoxLength);
       X_sum += centroid[0]+X_dist;
       Y_sum += centroid[1]+Y_dist;
       Z_sum += centroid[2]+Z_dist;
       ClusterTemp[i][0] = centroid[0]+X_dist;
       ClusterTemp[i][1] = centroid[1]+Y_dist;
       ClusterTemp[i][2] = centroid[2]+Z_dist;
   }

   double Xaverage = X_sum / numRows;
   double Yaverage = Y_sum / numRows;
   double Zaverage = Z_sum / numRows;

   centroid[0] = Xaverage;
   centroid[1] = Yaverage;
   centroid[2] = Zaverage;
   apply_periodic_boundary(centroid, BoxLength, centroid);
}


void writeTrajectory(int TrajectoryInterval, string& Style, int Step, int NumberOfParticles, double AASigma, double CGCubeSize, double CGSigma, double BoxLengthHalf,
                     vector<double>& Sizes, vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,
                     double VectorX[][3], double VectorY[][3], double VectorZ[][3]) {
    if (Step % TrajectoryInterval == 0) {
        /*if (Style == "Virtual") {
            ofstream LAMMPSTrajectoryFile("LAMMPSTrajectory.lammpstrj", std::ios::app);
            LAMMPSTrajectoryFile << "ITEM: TIMESTEP\n" << Step << "\n";
            LAMMPSTrajectoryFile << "ITEM: NUMBER OF ATOMS\n" << NumberOfParticles << "\n";
            LAMMPSTrajectoryFile << "ITEM: BOX BOUNDS pp pp pp\n" << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << "ITEM: ATOMS mol id type x y z radius quatw quati quatj quatk\n";

            for (int i = 0; i < NumberOfParticles; ++i) {
                double CubeSideLength = Sizes[i] * AASigma;
                Rotation rotation(VectorX[i], VectorY[i], VectorZ[i]);
            /*    LAMMPSTrajectoryFile << i << " " << i << " 1 " << Rx[i] << " " << Ry[i] << " " << Rz[i] << " " 
                                     << CubeSideLength/2 << " " << Quaterion[1] << " " << Quaterion[2] << " " 
                                     << Quaterion[3] << " " << -Quaterion[0] << "\n";
            }

        }*/ if (Style == "AA") {
            int numAtoms = CGCubeSize * CGCubeSize * CGCubeSize * NumberOfParticles;

            ofstream LAMMPSTrajectoryFile("LAMMPSTrajectory.lammpstrj", std::ios::app);
            LAMMPSTrajectoryFile << "ITEM: TIMESTEP\n" << Step << "\n";
            LAMMPSTrajectoryFile << "ITEM: NUMBER OF ATOMS\n" << numAtoms << "\n";
            LAMMPSTrajectoryFile << "ITEM: BOX BOUNDS pp pp pp\n" << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << "ITEM: ATOMS mol id type x y z radius\n";

            double COMX,COMY,COMZ;
            int NumAtomsPerCube = pow(CGCubeSize, 3);
            double atomList[NumAtomsPerCube][3];

            for (int i = 0; i < NumberOfParticles; ++i) {
                COMX = Rx[i];
                COMY = Ry[i];
                COMZ = Rz[i];
                findAtomPositionsInCube(CGCubeSize, CGSigma, COMX, COMY, COMZ, VectorX[i], VectorY[i], VectorZ[i], atomList);

                for (int j = 0; j < NumAtomsPerCube; ++j) {
                    double atomX = atomList[j][0];
                    double atomY = atomList[j][1];
                    double atomZ = atomList[j][2];

                    LAMMPSTrajectoryFile << i << " " << i * NumAtomsPerCube + j + 1 << " 2 " << atomX << " " << atomY << " " << atomZ << " " << CGSigma / 2 << "\n";
                }
            }
        } else if (Style == "Vertex") {
            std::ofstream LAMMPSTrajectoryFile("LAMMPSTrajectory.lammpstrj", std::ios::app);
            LAMMPSTrajectoryFile << "ITEM: TIMESTEP\n" << Step << "\n";
            LAMMPSTrajectoryFile << "ITEM: NUMBER OF ATOMS\n" << 8 * NumberOfParticles << "\n";
            LAMMPSTrajectoryFile << "ITEM: BOX BOUNDS pp pp pp\n" << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << "ITEM: ATOMS mol id type x y z radius\n";
            double Vertex[8][3];
            double COMX,COMY,COMZ;

            for (int i = 0; i < NumberOfParticles; ++i) {
                double CubeSideLength = Sizes[i] * AASigma;
                COMX = Rx[i];
                COMY = Ry[i];
                COMZ = Rz[i];
                findVerticesOfCube(CubeSideLength, COMX, COMY, COMZ, VectorX[i], VectorY[i], VectorZ[i], Vertex);
                for (int j = 0; j < 8; ++j) {
                    double atomX = Vertex[j][0];
                    double atomY = Vertex[j][1];
                    double atomZ = Vertex[j][2];

                    LAMMPSTrajectoryFile << i << " " << i * 8 + j + 1 << " 2 " << atomX << " " << atomY << " " << atomZ << " " << AASigma / 2 << "\n";
                }
            }
        }
    }
}

void readRestart(const std::string& restartFile, int NumberOfParticles,
                 std::vector<double>& Sizes, std::vector<double>& Rx,
                 std::vector<double>& Ry, std::vector<double>& Rz,
                 double VectorX[][3], double VectorY[][3], double VectorZ[][3],
                 int& restartStep) {
    std::ifstream inputFile(restartFile);

    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open file " << restartFile << std::endl;
        return;
    }

    std::string line;
    std::getline(inputFile, line);  // Read and discard the header line
    int num = 0;

    // Iterate through each line in the CSV file
    while (std::getline(inputFile, line)) {

        // Use a string stream to parse the values
        std::stringstream ss(line);
        std::string value;
        std::vector<double> values;

        // Iterate through each value in the line
        while (std::getline(ss, value, ',')) {
            values.push_back(std::stod(value));
        }

        // Now you have the values in the 'values' vector
        // Example usage:
        double size = values[0];
        double rx = values[1];
        double ry = values[2];
        double rz = values[3];
        double vxx = values[4];
        double vxy = values[5];
        double vxz = values[6];
        double vyx = values[7];
        double vyy = values[8];
        double vyz = values[9];
        double vzx = values[10];
        double vzy = values[11];
        double vzz = values[12];
        Sizes[num] = size;
        Rx[num] = rx;
        Ry[num] = ry;
        Rz[num] = rz;
        VectorX[num][0] = vxx;
        VectorX[num][1] = vxy;
        VectorX[num][2] = vxz;
        VectorY[num][0] = vyx;
        VectorY[num][1] = vyy;
        VectorY[num][2] = vyz;
        VectorZ[num][0] = vzx;
        VectorZ[num][1] = vzy;
        VectorZ[num][2] = vzz;

        num += 1;
    }

    // Close the input file
    inputFile.close();
}


int writeRestart(int RestartFileInterval, int lastRestartStep, vector<double>& Sizes,
                 vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,
                 double VectorX[][3], double VectorY[][3], double VectorZ[][3], int Step) {
    // write restart file
    if (Step % RestartFileInterval == 0) {
        string lastRestartFile = "restart_" + std::to_string(lastRestartStep) + ".csv";
        ofstream outputFile("restart_" + std::to_string(Step) + ".csv");

        if (outputFile.is_open()) {
            outputFile << "Sizes,Rx,Ry,Rz,VectorXx,VectorXy,VectorXz,VectorYx,VectorYy,VectorYz,VectorZx,VectorZy,VectorZz\n";

            for (size_t i = 0; i < Sizes.size(); ++i) {
                outputFile << Sizes[i] << "," << Rx[i] << "," << Ry[i] << "," << Rz[i] << ","
                           << VectorX[i][0] << "," << VectorX[i][1] << "," << VectorX[i][2] << ","
                           << VectorY[i][0] << "," << VectorY[i][1] << "," << VectorY[i][2] << ","
                           << VectorZ[i][0] << "," << VectorZ[i][1] << "," << VectorZ[i][2] << "\n";
            }

            outputFile.close();

            if (std::filesystem::exists(lastRestartFile)) {
                std::filesystem::remove(lastRestartFile);
            }

            lastRestartStep = Step;
        }
    }

    return lastRestartStep;
}

void writeEnergy(int EnergyOutputInterval, int Step, double SystemPotentialEnergy, ofstream& EnergyFile) {
    // Write out energies every EnergyOutputInterval steps
    if (Step % EnergyOutputInterval == 0) {
       EnergyFile << Step << " " << SystemPotentialEnergy << "\n" << flush;
    }
}

int main() {
    string Style = "AA"; //  AA, Virutal, Vertex
    string model = "vdW"; //AACG, vdW

    int NumberOfParticles = 8;
    double BoxLength = 125.0;
    double BoxLengthHalf = BoxLength/2.0;

    double kB = 0.0019872;  // Boltzmann constant, kcal/(mol*K)
    double Temperature = 300.0;  // Temperature (K)
    double kBT = kB * Temperature;  // kB * Temperature
    double AASigma = 2.63;
    double MaxRotation = 30 * M_PI / 180; //maximum allowed rotation for a cube in each trial step
    double MaxMove = AASigma  * 4; // maximum allowed translation for a cube in each trial step
    int EquilibrationSteps = 0; // Equilibration time steps
    int ProductionSteps = 100000000; // Production time steps
    int EnergyOutputInterval = 1000;
    int TrajectoryInterval = 1000;
    int RestartFileInterval = 10000;
    int lastRestartStep = 0;
    double CubeSideAASize = 12;
    unordered_map<double, double> SizeRatioDictionary = {{CubeSideAASize, 1.0}};
    double CGSigma, CGCubeSize, cutoffCG, AAEpsilon, CGEpsilon, Hamaker, AtomDensity;

    if (model == "AACG") {
       CGSigma = AASigma * 6;
       cutoffCG = 2.5 * CGSigma;
       CGCubeSize = CubeSideAASize * AASigma / CGSigma;
       AAEpsilon = 0.24 * 20 /219.71115;
       // Dictionary for scale factors
       unordered_map<double, double> ScaleFactorDictionary = {
        {1 * AASigma, 1.0},
        {2 * AASigma, 4.55724609191523},
        {3 * AASigma, 11.650078714945},
        {4 * AASigma, 23.4823671146823},
        {6 * AASigma, 67.683078067836547},
        {12 * AASigma, 438.9791464846652}
        };

        double ScaleFactor;
        if (ScaleFactorDictionary.find(CGSigma) != ScaleFactorDictionary.end()) {
        ScaleFactor = ScaleFactorDictionary[CGSigma];
        }
        // Scale epsilon of AA to match minimum of CG energy with minimum of AA energy
        CGEpsilon = AAEpsilon * ScaleFactor;
    } else if (model == "vdW") {
        CGSigma = 1 * AASigma; // vdW is not CG model, this CG size is just for trajectory file outputing AA style purposes
        CGCubeSize = CubeSideAASize; // vdW is not CG model, this CG size is just for trajectory file outputing AA style purposes
        double ScaleFactor = 0.884470464305972;
        Hamaker = 1.5 * ScaleFactor * 20.0 / 219.71115; // Hamaker constant, (10E^19 J)
        AtomDensity = 0.0585648358351751; // atom density in a cube, atoms/Angstrom^3
    }

    vector<double> Sizes(NumberOfParticles, 0.0);
    vector<double> Rx(NumberOfParticles, 0.0);
    vector<double> Ry(NumberOfParticles, 0.0);
    vector<double> Rz(NumberOfParticles, 0.0);
    double VectorX[NumberOfParticles][3];
    double VectorY[NumberOfParticles][3];
    double VectorZ[NumberOfParticles][3];

    int restartStep = 0;
    int Step = 0;

    ofstream LAMMPSTrajectoryFile("LAMMPSTrajectory.lammpstrj");
    ofstream EnergyFile("Energies.out");
    ofstream TimeFile("Time.out");
    //ofstream RDFFile("rdf.out");
    //ofstream MeanRDFFile("mean_rdf.out");
    //ofstream MeanAggFile("meanagg.out");
    //RDFFile.flush();
    //MeanRDFFile.flush();
    //MeanAggFile.flush();
    // Call the function
    InitialConfiguration(NumberOfParticles, BoxLength, AASigma, SizeRatioDictionary,
                         Sizes, Rx, Ry, Rz, VectorX, VectorY, VectorZ,
                         restartStep, LAMMPSTrajectoryFile, EnergyFile, TimeFile);

    //string restartFile = "restart_29900.csv";
    //readRestart(restartFile,NumberOfParticles,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,restartStep);

    writeTrajectory(TrajectoryInterval, Style, Step, NumberOfParticles,AASigma,CGCubeSize,CGSigma,BoxLengthHalf,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ);
    lastRestartStep = writeRestart(RestartFileInterval,lastRestartStep,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,Step);
    ofstream ProbFile("prob.out");
    ProbFile.flush();
    double SystemPotentialEnergy = TotalEnergy(Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,NumberOfParticles,BoxLength,CGSigma,CGEpsilon,AASigma,model,cutoffCG,Hamaker,AtomDensity);

    double Attempt = 0.0;
    double Accept = 0.0;
    int CountBoltzmann = 0;
    double SingleMoveAccept = 0.0;
    double SingleMoveNumber = 0.0;
    double ClusterMoveAccept = 0.0;
    double ClusterMoveNumber = 0.0;

    // For rdf calculation
    //int num_bins = 1000;
    //double normalized_bin_counts[num_bins];
    //double mean_bin_counts[num_bins];
    //for (int k = 0; k < num_bins; ++k) {
      // mean_bin_counts[k] = 0.0;
    //}
    //double bin_width = 0.1;
    //double normalization_constant = NumberOfParticles * (NumberOfParticles / pow(BoxLength,3.0)) * 4.0 * M_PI * bin_width;

    for (int Step = 1 + restartStep; Step <= restartStep + EquilibrationSteps + ProductionSteps; ++Step) {
        //double mean_agg = 0.0;
        //for (int k = 0; k < num_bins; ++k) {
          // normalized_bin_counts[k] = 0.0;
        //}
        bool EarlyTermination = false;
        bool overLapFlag = false;
        ProbFile << " " << endl;
        ProbFile << "Step " << Step << endl;
        cout << "Step " << Step << endl;
        cout << endl;
        if (Step == 1) {
           cout << "\nStarting Equilibration\n";
        }
        if (Step == EquilibrationSteps + 1) {
           cout << "\nStarting Production\n";
        }

        if (Step % 100 == 0) {
           if (Step <= EquilibrationSteps) {
              cout << "Equilibration Step " << Step << "\n";
           } else {
              cout << "Production Step " << Step << "\n";
           }
        }

        // Choose a random particle i as seed
        int i = static_cast<int>(NumberOfParticles * ((double)rand()/(RAND_MAX)));
        ProbFile << "Seed particle: " << i << endl;

        // Position and orientation of the seed
        double seedX = Rx[i]; // center of mass in x, Angstrom
        double seedY = Ry[i]; // center of mass in y, Angstrom
        double seedZ = Rz[i]; // center of mass in z, Angstrom
        Vector3 seedCentroid = {seedX, seedY, seedZ};
        Vector3 seedVectorX, seedVectorY, seedVectorZ;

        for (int k = 0; k < 3; ++k) {
             seedVectorX[k] = VectorX[i][k]; // x component of orientation vector
             seedVectorY[k] = VectorY[i][k]; // y component of orientation vector
             seedVectorZ[k] = VectorZ[i][k]; // z component of orientation vector
        }

        // Create a list to recruit particles into a cluster
        vector<int> clusterParticleList;
        // Recruit the seed particle as the first member of the cluster
        clusterParticleList.push_back(i);
        // Create a list for particles haven't been recruited to the cluster
        vector<int> ParticlesNotInClusterList(NumberOfParticles);
        for (int ParticleID = 0; ParticleID < NumberOfParticles; ++ParticleID) {
            ParticlesNotInClusterList[ParticleID] = ParticleID;
        }
        // Particles that for sure will be recruited to cluster but due to multiple of such particles bonding with a linker at the same time
        // They will be revisited later after going through one branch
        // Each of these particles will begin a new branch of recruitment
        vector<int> ToBeRecruitedParticlesList;
        // Temporary List for the j(s) linked with current i
        vector<int> TemperaryjList;

        // Removing the seed particle from ParticlesNotInClusterList
        ParticlesNotInClusterList.erase(remove(ParticlesNotInClusterList.begin(), ParticlesNotInClusterList.end(), i), ParticlesNotInClusterList.end());

        vector<double> reverseMoveProbabilityList, forwardMoveProbabilityList, unacceptedMoveProbabilityList, unacceptedReverseMoveProbabilityList;
        vector<int> linkConnectionList;

        double MaxClusterSize;
        double singleOrCluster = ((double)rand()/(RAND_MAX));
        if (singleOrCluster < 0.5) {
           MaxClusterSize = 1; // single move with 50% probability
           MaxMove = AASigma * 4;
           MaxRotation = 10 * M_PI / 180;
           ProbFile << "Single move selected" << endl;
        } else {
           MaxClusterSize = NumberOfParticles;
           MaxMove = AASigma * 48;
        }

        double tranOrRot = ((double)rand()/(RAND_MAX));
        bool isTranslation = false;
        bool isRotation = false;
        double VirtualMoveX, VirtualMoveY, VirtualMoveZ, VirtualMoveAlpha, VirtualMoveBeta, VirtualMoveGamma;

        if (tranOrRot <= 0.5) {
           // Translation
           isTranslation = true;
           VirtualMoveX = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxMove;
           VirtualMoveY = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxMove;
           VirtualMoveZ = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxMove;
           VirtualMoveAlpha = 0.0;
           VirtualMoveBeta = 0.0;
           VirtualMoveGamma = 0.0;
        }else {
           // Rotation
           isRotation = true;
           VirtualMoveX = 0;
           VirtualMoveY = 0;
           VirtualMoveZ = 0;
           double randomAxis = ((double)rand()/(RAND_MAX));
           // randomly choose an axis from x,y,z to rotate
           if (randomAxis <= 0.333) { // pick x-axis to rotate, the other two rotation angles are zero
              VirtualMoveAlpha = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxRotation;
              VirtualMoveBeta = 0.0;
              VirtualMoveGamma = 0.0;
           } else if (randomAxis <= 0.666){
              VirtualMoveAlpha = 0.0;
              VirtualMoveBeta = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxRotation;
              VirtualMoveGamma = 0.0;
           } else {
              VirtualMoveAlpha = 0.0;
              VirtualMoveBeta = 0.0;
              VirtualMoveGamma = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxRotation;
           }
        }

        Matrix3x3 VirtualMoveRotationX = {{1, 0, 0}, {0, cos(VirtualMoveAlpha), -sin(VirtualMoveAlpha)}, {0, sin(VirtualMoveAlpha), cos(VirtualMoveAlpha)}};
        Matrix3x3 VirtualMoveRotationY = {{cos(VirtualMoveBeta), 0, sin(VirtualMoveBeta)}, {0, 1, 0}, {-sin(VirtualMoveBeta), 0, cos(VirtualMoveBeta)}};
        Matrix3x3 VirtualMoveRotationZ = {{cos(VirtualMoveGamma), -sin(VirtualMoveGamma), 0}, {sin(VirtualMoveGamma), cos(VirtualMoveGamma), 0}, {0, 0, 1}};

        Matrix3x3 ReverseMoveRotationX = {{1, 0, 0}, {0, cos(-VirtualMoveAlpha), -sin(-VirtualMoveAlpha)}, {0, sin(-VirtualMoveAlpha), cos(-VirtualMoveAlpha)}};
        Matrix3x3 ReverseMoveRotationY = {{cos(-VirtualMoveBeta), 0, sin(-VirtualMoveBeta)}, {0, 1, 0}, {-sin(-VirtualMoveBeta), 0, cos(-VirtualMoveBeta)}};
        Matrix3x3 ReverseMoveRotationZ = {{cos(-VirtualMoveGamma), -sin(-VirtualMoveGamma), 0}, {sin(-VirtualMoveGamma), cos(-VirtualMoveGamma), 0}, {0, 0, 1}};

        ProbFile << "VirtualMoveX " << VirtualMoveX << endl;
        ProbFile << "VirtualMoveY " << VirtualMoveY << endl;
        ProbFile << "VirtualMoveZ " << VirtualMoveZ << endl;
        ProbFile << "VirtualMoveAlpha " << VirtualMoveAlpha << " VirtualMoveBeta " << VirtualMoveBeta << " VirtualMoveGamma " << VirtualMoveGamma << endl;

        // if maximum cluster size is 1, then there is no need to build cluster
        // because the particle itself is considered a cluster
        bool Loop;
        double overallReverseProbability, overallUnacceptedReverseProbability, overallMoveProbability, overallUnacceptedMoveProbability;
        if (MaxClusterSize > 1) {
           Loop = true; // continue to loop if Loop is True
        } else {
           Loop = false; // continue to loop if Loop is True
           overallReverseProbability = 1.0;
           overallUnacceptedReverseProbability = 1.0;
           overallMoveProbability = 1.0;
           overallUnacceptedMoveProbability = 1.0;
        }
        // to regulate cutoff for cluster size
        double Nc = ((double)rand()/(RAND_MAX));
        ProbFile << "Nc " << (1/Nc) << endl;

        ProbFile << "clusterParticleList: " << endl;
        for (int particleID : clusterParticleList) {
           ProbFile << particleID << " ";
        }
        ProbFile << endl;
        while (Loop) {
           // check if size of cluster exceeds maximum
           if (clusterParticleList.size() >= MaxClusterSize) {
              Loop = false;
              ProbFile << "Cluster exceeds maximum cluster size" << endl;
              break;
           }
           // to abort the link formation procedure if the cluster size exceeds Nc
           if (1.0 / clusterParticleList.size() < Nc) {
              ProbFile << "Early termination" << endl;
              //ProbFile << "clusterParticleList: " << endl;
              //for (int particleID : clusterParticleList) {
              //   ProbFile << particleID << " ";
              //}
              //ProbFile << endl;
              Loop = false; // End recruiting particles to cluster
              EarlyTermination = true;
              break;
           }
           double Size_i = Sizes[i];
           double Rx_i_old = Rx[i];
           double Ry_i_old = Ry[i];
           double Rz_i_old = Rz[i];
           //ProbFile << "Rx_i_old " << Rx_i_old << " Ry_i_old " << Ry_i_old << " Rz_i_old " << Rz_i_old << endl;
           Vector3 ParticleiCentroid_old = {Rx_i_old, Ry_i_old, Rz_i_old};
           Vector3 VectorX_i_old, VectorY_i_old, VectorZ_i_old, VectorX_j_old, VectorY_j_old, VectorZ_j_old;
           for (int k = 0; k < 3; ++k) {
              VectorX_i_old[k] = VectorX[i][k]; // x component of orientation vector
              VectorY_i_old[k] = VectorY[i][k]; // y component of orientation vector
              VectorZ_i_old[k] = VectorZ[i][k]; // z component of orientation vector
           }
           // loop over all j(s) in ParticlesNotInClusterList
           // to find all j(s) that link to the current i
           for (int ParticleNotInCluster : ParticlesNotInClusterList) {
              // Go through each particle in ParticlesNotInClusterList
              // to check i-j interaction to determine if j should be recruited to cluster
              int j = ParticleNotInCluster;
              ProbFile << i << " and " << j << " investigation" << endl;

              ProbFile << "Rx_i_old " << Rx_i_old << " Ry_i_old " << Ry_i_old << " Rz_i_old " << Rz_i_old << endl;
              double Size_i = Sizes[i];
              double Size_j = Sizes[j];
              double Rx_j_old = Rx[j];  // center of mass in x, Angstrom
              double Ry_j_old = Ry[j];  // center of mass in y, Angstrom
              double Rz_j_old = Rz[j];  // center of mass in z, Angstrom
              ProbFile << "Rx_j_old " << Rx_j_old << " Ry_j_old " << Ry_j_old << " Rz_j_old " << Rz_j_old << endl;
              for (int k = 0; k < 3; ++k) {
                 VectorX_j_old[k] = VectorX[j][k]; // x component of orientation vector
                 VectorY_j_old[k] = VectorY[j][k]; // y component of orientation vector
                 VectorZ_j_old[k] = VectorZ[j][k]; // z component of orientation vector
              }

              ProbFile << "VectorX_i_old[0] " << VectorX_i_old[0] << " VectorX_i_old[1] " << VectorX_i_old[1] << " VectorX_i_old[2] " << VectorX_i_old[2] <<endl;
              ProbFile << "VectorY_i_old[0] " << VectorY_i_old[0] << " VectorY_i_old[1] " << VectorY_i_old[1] << " VectorY_i_old[2] " << VectorY_i_old[2] <<endl;
              ProbFile << "VectorZ_i_old[0] " << VectorZ_i_old[0] << " VectorZ_i_old[1] " << VectorZ_i_old[1] << " VectorZ_i_old[2] " << VectorZ_i_old[2] <<endl;
              ProbFile << "VectorX_j_old[0] " << VectorX_j_old[0] << " VectorX_j_old[1] " << VectorX_j_old[1] << " VectorX_j_old[2] " << VectorX_j_old[2] <<endl;
              ProbFile << "VectorY_j_old[0] " << VectorY_j_old[0] << " VectorY_j_old[1] " << VectorY_j_old[1] << " VectorY_j_old[2] " << VectorY_j_old[2] <<endl;
              ProbFile << "VectorZ_j_old[0] " << VectorZ_j_old[0] << " VectorZ_j_old[1] " << VectorZ_j_old[1] << " VectorZ_j_old[2] " << VectorZ_j_old[2] <<endl;
              // relative distance between i and j
              double Rx_ij_old = Rx_i_old - Rx_j_old;
              double Ry_ij_old = Ry_i_old - Ry_j_old;
              double Rz_ij_old = Rz_i_old - Rz_j_old;
              // smallest relative distance between i and j considering periodic boundary conditions
              Rx_ij_old = Rx_ij_old - BoxLength * round(Rx_ij_old/BoxLength);
              Ry_ij_old = Ry_ij_old - BoxLength * round(Ry_ij_old/BoxLength);
              Rz_ij_old = Rz_ij_old - BoxLength * round(Rz_ij_old/BoxLength);
              // absolute distance between i and j
              double ijDistanceSquare_old = Rx_ij_old * Rx_ij_old + Ry_ij_old * Ry_ij_old + Rz_ij_old * Rz_ij_old;
              double ijDistance_old = pow(ijDistanceSquare_old, 0.5);
              // i as reference, move j to the image where ij distance is minimal
              Rx_j_old = Rx_i_old - Rx_ij_old;
              Ry_j_old = Ry_i_old - Ry_ij_old;
              Rz_j_old = Rz_i_old - Rz_ij_old;
              Vector3 ParticlejCentroid_old = {Rx_j_old, Ry_j_old, Rz_j_old};
              ProbFile << "ParticleiCentroid_old[0] " << ParticleiCentroid_old[0] << " ParticleiCentroid_old[1] " << ParticleiCentroid_old[1] << " ParticleiCentroid_old[2] " << ParticleiCentroid_old[2] << endl;
              ProbFile << "ParticlejCentroid_old[0] " << ParticlejCentroid_old[0] << " ParticlejCentroid_old[1] " << ParticlejCentroid_old[1] << " ParticlejCentroid_old[2] " << ParticlejCentroid_old[2] << endl;

              // larger size
              double MaxSize = max(Size_i, Size_j);
              double CutOff = 1.5 * MaxSize * AASigma;
              double CutOffSquare = CutOff*CutOff;
              double ijEnergy_old, ijEnergy_new, CubeSideLength1, CubeSideLength2, ijEnergy_reverse;
              double virtualMoveProbability,reverseMoveProbability;

              if (ijDistanceSquare_old <= CutOffSquare) {
                 ProbFile << i << " and " << j << " investigation " << endl;
                 // check energy between i and j before vitual move
                 ijEnergy_old = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid_old,ParticlejCentroid_old,VectorX_i_old, VectorY_i_old, VectorZ_i_old, VectorX_j_old, VectorY_j_old, VectorZ_j_old);
                 CubeSideLength1 = Size_i * AASigma;
                 CubeSideLength2 = Size_j * AASigma;

                 // particle centroid in new coordinate system, with seed particle as the centroid
                 Vector3 ParticleCentroid = {Rx[i]-seedX,Ry[i]-seedY,Rz[i]-seedZ};
                 apply_periodic_boundary(ParticleCentroid,BoxLength,ParticleCentroid);
                 Vector3 ParticleVectorX, ParticleVectorY, ParticleVectorZ;
                 for (int k = 0; k < 3; ++k) {
                    ParticleVectorX[k] = VectorX[i][k]; // x component of orientation vector
                    ParticleVectorY[k] = VectorY[i][k]; // y component of orientation vector
                    ParticleVectorZ[k] = VectorZ[i][k]; // z component of orientation vector
                    //ProbFile << "ParticleVectorX[" << k << "] " << ParticleVectorX[k] << "ParticleVectorY[" << k << "] " << ParticleVectorY[k]<< "ParticleVectorZ[" << k << "] " << ParticleVectorZ[k] << endl;
                 }
                 // rotate each particle centroid and vectors about the centroid of the seed
                 Vector3 tempresult1, tempresult2;
                 Vector3 RotatedParticleCentroid, VectorX_i_new, VectorY_i_new, VectorZ_i_new;
                 MultiplyMatrixVector(VirtualMoveRotationX, ParticleCentroid, tempresult1);
                 MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
                 MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, RotatedParticleCentroid);

                 MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorX, tempresult1);
                 MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
                 MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorX_i_new);

                 MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorY, tempresult1);
                 MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
                 MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorY_i_new);

                 MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorZ, tempresult1);
                 MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
                 MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorZ_i_new);

                 // shift the particle coordinates back to the original coordinate system
                 double Rx_i_new = RotatedParticleCentroid[0] + seedX;
                 double Ry_i_new = RotatedParticleCentroid[1] + seedY;
                 double Rz_i_new = RotatedParticleCentroid[2] + seedZ;

                 // new position of i after virtual move
                 Rx_i_new = Rx_i_new + VirtualMoveX;
                 Ry_i_new = Ry_i_new + VirtualMoveY;
                 Rz_i_new = Rz_i_new + VirtualMoveZ;

                 // apply periodic boundary conditions to positions
                 Vector3 R_i_new = {Rx_i_new,Ry_i_new,Rz_i_new};
                 apply_periodic_boundary(R_i_new, BoxLength,R_i_new);
                 Rx_i_new = R_i_new[0];
                 Ry_i_new = R_i_new[1];
                 Rz_i_new = R_i_new[2];

                 // relatvie distance between new i and old j
                 double Rx_ij_new = Rx_i_new - Rx_j_old;
                 double Ry_ij_new = Ry_i_new - Ry_j_old;
                 double Rz_ij_new = Rz_i_new - Rz_j_old;

                 // smallest relative distance between new i and old j considering periodic boundary conditions
                 Rx_ij_new = Rx_ij_new - BoxLength * round(Rx_ij_new/BoxLength);
                 Ry_ij_new = Ry_ij_new - BoxLength * round(Ry_ij_new/BoxLength);
                 Rz_ij_new = Rz_ij_new - BoxLength * round(Rz_ij_new/BoxLength);

                 // new i as reference, move j to the image where ij distance is minimal
                 Rx_j_old = Rx_i_new - Rx_ij_new;
                 Ry_j_old = Ry_i_new - Ry_ij_new;
                 Rz_j_old = Rz_i_new - Rz_ij_new;

                 Vector3 ParticleiCentroid_new = {Rx_i_new, Ry_i_new, Rz_i_new};
                 Vector3 ParticlejCentroid_minimal = {Rx_j_old,Ry_j_old,Rz_j_old};

                 // check energy between i and j after virtual move
                 if (CheckOverlapForTwoCubes(i,j,BoxLength,AASigma,CubeSideLength1,CubeSideLength2,ParticleiCentroid_new,ParticlejCentroid_minimal,VectorX_i_new,VectorY_i_new,VectorZ_i_new,VectorX_j_old,VectorY_j_old,VectorZ_j_old)) {
                    ijEnergy_new = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid_new, ParticlejCentroid_old,VectorX_i_new, VectorY_i_new, VectorZ_i_new, VectorX_j_old, VectorY_j_old, VectorZ_j_old);
                    overLapFlag = true;
                 } else {
                    ijEnergy_new = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid_new, ParticlejCentroid_old,VectorX_i_new, VectorY_i_new, VectorZ_i_new, VectorX_j_old, VectorY_j_old, VectorZ_j_old);
                 }
                 try {
                    // Calculate virtualMoveProbability, considering the possibility of OverflowError
                    virtualMoveProbability = max(0.0, 1.0 - exp((ijEnergy_old - ijEnergy_new) / kBT));
                 } catch (const std::overflow_error& e) {
                    virtualMoveProbability = 0.0;
                 }
                 ProbFile << "ijEnergy_old " << ijEnergy_old << endl;
                 ProbFile << "ijEnergy_new " << ijEnergy_new << endl;
                 ProbFile << "virtualMoveProbability " << virtualMoveProbability << endl;
                 // check acceptance probability of this virtual move
                 double rand_num = ((double)rand()/(RAND_MAX));
                 if (virtualMoveProbability >= rand_num) {
                    ProbFile << "acceptedVirtualMoveProbability " << virtualMoveProbability << endl;
                    // add j to TemperaryjList for later use
                    TemperaryjList.push_back(j);
                    ProbFile << "TemperaryjList j " << j << endl;
                    // calculate reverse move acceptance probability
                    // particle centroid in new coordinate system, with seed particle as the centroid
                    ParticleCentroid[0] = Rx[i] - seedX;
                    ParticleCentroid[1] = Ry[i] - seedY;
                    ParticleCentroid[2] = Rz[i] - seedZ;
                    apply_periodic_boundary(ParticleCentroid, BoxLength, ParticleCentroid);
                    for (int k = 0; k < 3; ++k) {
                       ParticleVectorX[k] = VectorX[i][k]; // x component of orientation vector
                       ParticleVectorY[k] = VectorY[i][k]; // y component of orientation vector
                       ParticleVectorZ[k] = VectorZ[i][k]; // z component of orientation vector
                    }
                    MultiplyMatrixVector(ReverseMoveRotationX, ParticleCentroid, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, RotatedParticleCentroid);

                    // shift the particle coordinates back to the original coordinate system
                    double Rx_i_reverse = RotatedParticleCentroid[0] + seedX;
                    double Ry_i_reverse = RotatedParticleCentroid[1] + seedY;
                    double Rz_i_reverse = RotatedParticleCentroid[2] + seedZ;

                    // new position of i after reverse move
                    Rx_i_reverse = Rx_i_reverse - VirtualMoveX;
                    Ry_i_reverse = Ry_i_reverse - VirtualMoveY;
                    Rz_i_reverse = Rz_i_reverse - VirtualMoveZ;

                    Vector3 VectorX_i_reverse, VectorY_i_reverse, VectorZ_i_reverse;
                    MultiplyMatrixVector(ReverseMoveRotationX, VectorX_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorX_i_reverse);

                    MultiplyMatrixVector(ReverseMoveRotationX, VectorY_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorY_i_reverse);

                    MultiplyMatrixVector(ReverseMoveRotationX, VectorZ_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorZ_i_reverse);

                    // pickup the central image for the new location of i after reverse move
                    Rx_i_reverse = Rx_i_reverse - BoxLength*round(Rx_i_reverse/BoxLength);
                    Ry_i_reverse = Ry_i_reverse - BoxLength*round(Ry_i_reverse/BoxLength);
                    Rz_i_reverse = Rz_i_reverse - BoxLength*round(Rz_i_reverse/BoxLength);

                    // relative distance between reversed i and old j
                    double Rx_ij_reverse = Rx_i_reverse - Rx_j_old;
                    double Ry_ij_reverse = Ry_i_reverse - Ry_j_old;
                    double Rz_ij_reverse = Rz_i_reverse - Rz_j_old;

                    Rx_ij_reverse = Rx_ij_reverse - BoxLength * round(Rx_ij_reverse/BoxLength);
                    Ry_ij_reverse = Ry_ij_reverse - BoxLength * round(Ry_ij_reverse/BoxLength);
                    Rz_ij_reverse = Rz_ij_reverse - BoxLength * round(Rz_ij_reverse/BoxLength);

                    // reversed i as reference, move j to the image where ij distance is minimal
                    Rx_j_old = Rx_i_reverse - Rx_ij_reverse;
                    Ry_j_old = Ry_i_reverse - Ry_ij_reverse;
                    Rz_j_old = Rz_i_reverse - Rz_ij_reverse;

                    CubeSideLength1 = Sizes[i] * AASigma;
                    CubeSideLength2 = Sizes[j] * AASigma;
                    Vector3 ParticleiCentroid_reverse = {Rx_i_reverse,Ry_i_reverse,Rz_i_reverse};
                    // check energy between i and j after reverse move
                    if (CheckOverlapForTwoCubes(i,j,BoxLength,AASigma,CubeSideLength1,CubeSideLength2,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old)) {
                       ijEnergy_reverse = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old);
                       ProbFile << i << " and " << j << " overlap in reverse move" << endl;
                    } else {
                       ijEnergy_reverse = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old);
                    }
                    // calculate reverse move acceptance probability
                    try {
                       reverseMoveProbability = max(0.0, 1.0 - exp((ijEnergy_old - ijEnergy_reverse) / kBT));
                    } catch (const std::overflow_error& e) {
                       reverseMoveProbability = 0.0;
                    }
                    reverseMoveProbabilityList.push_back(reverseMoveProbability);
                    forwardMoveProbabilityList.push_back(virtualMoveProbability);
                    ProbFile << "virtualMoveProbability " << virtualMoveProbability << endl;
                    ProbFile << "reverseMoveProbability " << reverseMoveProbability << endl;
                 } else {
                    double qij = 1-virtualMoveProbability;
                    unacceptedMoveProbabilityList.push_back(qij);
                    ProbFile << "unacceptedMoveProbability " << qij << endl;
                    // calculate reverse move acceptance probability
                    ParticleCentroid[0] = Rx[i] - seedX;
                    ParticleCentroid[1] = Ry[i] - seedY;
                    ParticleCentroid[2] = Rz[i] - seedZ;
                    apply_periodic_boundary(ParticleCentroid, BoxLength, ParticleCentroid);
                    MultiplyMatrixVector(ReverseMoveRotationX, ParticleCentroid, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, RotatedParticleCentroid);

                    // shift the particle coordinates back to the original coordinate system
                    double Rx_i_reverse = RotatedParticleCentroid[0] + seedX;
                    double Ry_i_reverse = RotatedParticleCentroid[1] + seedY;
                    double Rz_i_reverse = RotatedParticleCentroid[2] + seedZ;

                    // new position of i after reverse move
                    Rx_i_reverse = Rx_i_reverse - VirtualMoveX;
                    Ry_i_reverse = Ry_i_reverse - VirtualMoveY;
                    Rz_i_reverse = Rz_i_reverse - VirtualMoveZ;

                    Vector3 VectorX_i_reverse, VectorY_i_reverse, VectorZ_i_reverse;
                    MultiplyMatrixVector(ReverseMoveRotationX, VectorX_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorX_i_reverse);

                    MultiplyMatrixVector(ReverseMoveRotationX, VectorY_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorY_i_reverse);

                    MultiplyMatrixVector(ReverseMoveRotationX, VectorZ_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorZ_i_reverse);

                    // pickup the central image for the new location of i after reverse move
                    Rx_i_reverse = Rx_i_reverse - BoxLength*round(Rx_i_reverse/BoxLength);
                    Ry_i_reverse = Ry_i_reverse - BoxLength*round(Ry_i_reverse/BoxLength);
                    Rz_i_reverse = Rz_i_reverse - BoxLength*round(Rz_i_reverse/BoxLength);

                    // relative distance between reversed i and old j
                    double Rx_ij_reverse = Rx_i_reverse - Rx_j_old;
                    double Ry_ij_reverse = Ry_i_reverse - Ry_j_old;
                    double Rz_ij_reverse = Rz_i_reverse - Rz_j_old;

                    Rx_ij_reverse = Rx_ij_reverse - BoxLength * round(Rx_ij_reverse/BoxLength);
                    Ry_ij_reverse = Ry_ij_reverse - BoxLength * round(Ry_ij_reverse/BoxLength);
                    Rz_ij_reverse = Rz_ij_reverse - BoxLength * round(Rz_ij_reverse/BoxLength);

                    // reversed i as reference, move j to the image where ij distance is minimal
                    Rx_j_old = Rx_i_reverse - Rx_ij_reverse;
                    Ry_j_old = Ry_i_reverse - Ry_ij_reverse;
                    Rz_j_old = Rz_i_reverse - Rz_ij_reverse;

                    CubeSideLength1 = Sizes[i] * AASigma;
                    CubeSideLength2 = Sizes[j] * AASigma;
                    Vector3 ParticleiCentroid_reverse = {Rx_i_reverse,Ry_i_reverse,Rz_i_reverse};
                    // check energy between i and j after reverse move
                    if (CheckOverlapForTwoCubes(i,j,BoxLength,AASigma,CubeSideLength1,CubeSideLength2,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old)) {
                       ijEnergy_reverse = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old);
                       ProbFile << i << " and " << j << " overlap in reverse move" << endl;
                    } else {
                       ijEnergy_reverse = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,BoxLength,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old);
                    }
                    // calculate reverse move acceptance probability
                    try {
                       reverseMoveProbability = max(0.0, 1.0 - exp((ijEnergy_old - ijEnergy_reverse) / kBT));
                    } catch (const std::overflow_error& e) {
                       reverseMoveProbability = 0.0;
                    }
                    double qji = 1-reverseMoveProbability;
                    unacceptedReverseMoveProbabilityList.push_back(qji);
                    linkConnectionList.push_back(j);
                    ProbFile << "unacceptedReverseMoveProbability " << qji << endl;
                 }
              }
           }
           overallReverseProbability = 1.0;
           overallMoveProbability = 1.0;
           overallUnacceptedMoveProbability = 1.0;
           overallUnacceptedReverseProbability = 1.0;

           if (reverseMoveProbabilityList.size() > 0) {
              for (double reverseProbability : reverseMoveProbabilityList) {
                 overallReverseProbability *= reverseProbability;
                 ProbFile << "overallReverseProbability " << overallReverseProbability << endl;
              }
           }

           if (unacceptedReverseMoveProbabilityList.size() > 0) {
              for (size_t index = 0; index < unacceptedReverseMoveProbabilityList.size(); ++index) {
                 if (find(clusterParticleList.begin(), clusterParticleList.end(), linkConnectionList[index]) != clusterParticleList.end()) {
                    // Unformed link internal to the cluster
                    double unacceptedReverseProbability = unacceptedReverseMoveProbabilityList[index];
                    overallUnacceptedReverseProbability *= unacceptedReverseProbability;
                    ProbFile << "overallUnacceptedReverseProbability " << overallUnacceptedReverseProbability << endl;
                 }
              }
           }

           if (forwardMoveProbabilityList.size() > 0) {
              for (double forwardMoveProbability : forwardMoveProbabilityList) {
                 overallMoveProbability *= forwardMoveProbability;
                 ProbFile << "overallMoveProbability " << overallMoveProbability << endl;
              }
           }

           if (unacceptedMoveProbabilityList.size() > 0) {
              for (size_t unacceptedindex = 0; unacceptedindex < unacceptedMoveProbabilityList.size(); ++unacceptedindex) {
                  if (find(clusterParticleList.begin(), clusterParticleList.end(), linkConnectionList[unacceptedindex]) != clusterParticleList.end()) {
                     // Unformed link internal to the cluster
                     double unacceptedMoveProbability = unacceptedMoveProbabilityList[unacceptedindex];
                     overallUnacceptedMoveProbability *= unacceptedMoveProbability;
                     ProbFile << "overallUnacceptedMoveProbability " << overallUnacceptedMoveProbability << endl;
                  }
              }
           }

           if (TemperaryjList.size() == 0) {
              if (ToBeRecruitedParticlesList.size() == 0) {
                 Loop = false; // End recruiting particles to cluster
              } else {
                 // Pick a particle from ToBeRecruitedParticlesList as the new i
                 int rand_index = ((double)rand()/(RAND_MAX)) * (ToBeRecruitedParticlesList.size() - 1);
                 i = ToBeRecruitedParticlesList[rand_index];

                 // Add the new i to clusterParticleList
                 if (find(clusterParticleList.begin(), clusterParticleList.end(), i) == clusterParticleList.end()) {
                    clusterParticleList.push_back(i);
                    ProbFile << "clusterParticleList: ";
                    for (int particleID : clusterParticleList) {
                       ProbFile << particleID << " ";
                    }
                    ProbFile << " " << endl;
                 }

                 // Delete the new i from ParticlesNotInClusterList
                 ParticlesNotInClusterList.erase(remove(ParticlesNotInClusterList.begin(), ParticlesNotInClusterList.end(), i), ParticlesNotInClusterList.end());
                 // Delete the new i from ToBeRecruitedParticlesList
                 ToBeRecruitedParticlesList.erase(remove(ToBeRecruitedParticlesList.begin(), ToBeRecruitedParticlesList.end(), i), ToBeRecruitedParticlesList.end());
              }
           } else { // There are some j(s) that link to i
              // Pick a particle from TemperaryjList as the new i
              int rand_ind = ((double)rand()/(RAND_MAX)) * (TemperaryjList.size() - 1);
              i = TemperaryjList[rand_ind];

              // Delete the new i from TemperaryjList
              TemperaryjList.erase(remove(TemperaryjList.begin(), TemperaryjList.end(), i), TemperaryjList.end());

              // Add the new i to clusterParticleList
             if (find(clusterParticleList.begin(), clusterParticleList.end(), i) == clusterParticleList.end()) {
                clusterParticleList.push_back(i);
                ProbFile << "clusterParticleList: ";
                for (int particleID : clusterParticleList) {
                   ProbFile << particleID << " ";
                }
                ProbFile << " " << endl;
             }

             // Delete the new i from ParticlesNotInClusterList
             ParticlesNotInClusterList.erase(remove(ParticlesNotInClusterList.begin(), ParticlesNotInClusterList.end(), i), ParticlesNotInClusterList.end());

             // Remove the new i from ToBeRecruitedParticlesList if it exists there
             ToBeRecruitedParticlesList.erase(remove(ToBeRecruitedParticlesList.begin(), ToBeRecruitedParticlesList.end(), i), ToBeRecruitedParticlesList.end());

             // Add the rest particles (if any) in TemperaryjList to ToBeRecruitedParticlesList
             if (TemperaryjList.size() > 0) {
                for (int Temperaryj : TemperaryjList) {
                   if (find(ToBeRecruitedParticlesList.begin(), ToBeRecruitedParticlesList.end(), Temperaryj) == ToBeRecruitedParticlesList.end()) {
                      ToBeRecruitedParticlesList.push_back(Temperaryj);
                   }
                }
             }
             // Make TemperaryjList empty
             TemperaryjList.clear();
           }
        }
        //////////////////////////  Building Cluster Ends  /////////////////
        // now the cluster is built
        if (EarlyTermination == false) {
           double clusterEnergy_old = 0.0;
           double clusterEnergy_new = 0.0;
           vector<array<double, 3>> particlesPositionList;
           array<double, 3> Ri;
           vector<int> pairList;
           double currentParticleOldEnergy,currentParticleNewEnergy;;
           bool overLapFlag, overLapFlagNew, overLapFlagOld;
           int clusterSize = clusterParticleList.size();
           for (int particleID : clusterParticleList) {
              i = particleID;
              pairList.push_back(i);
              //ProbFile << "pairList values: ";
              //for (int m : pairList) {
              //   ProbFile << m << " ";
              //}
              //ProbFile << std::endl;
              double Size_i = Sizes[i];

              //cout << "Rx old: ";
              //for (double m : Rx) {
              //   cout << m << endl;
              //}
              //cout << "Ry old: ";
              //for (double m : Ry) {
              //   cout << m << endl;
              //}
              //cout << "Rz old: ";
              //for (double m : Rz) {
              //   cout << m << endl;
              //}

              // old position and orientation of this particle before virtual move
              double Rx_i_old = Rx[i];
              double Ry_i_old = Ry[i];
              double Rz_i_old = Rz[i];
              //ProbFile <<"i " << i << " Rx_i_old " << Rx_i_old << " Ry_i_old " << Ry_i_old << " Rz_i_old " << Rz_i_old << endl;
              Ri[0] = Rx_i_old;
              Ri[1] = Ry_i_old;
              Ri[2] = Rz_i_old;
              Vector3 VectorX_i_old,VectorY_i_old,VectorZ_i_old;
              //for (int k = 0; k < 3; ++k) {
              //   VectorX_i_old[k] = VectorX[i][k]; // x component of orientation vector
              //   VectorY_i_old[k] = VectorY[i][k]; // y component of orientation vector
              //   VectorZ_i_old[k] = VectorZ[i][k]; // z component of orientation vector
              //}
              // energy of cluster before cluster move
              pair<double, bool> result_initial = OneParticleEnergy(model,i,Size_i,cutoffCG,Hamaker,AtomDensity,Rx_i_old,Ry_i_old,Rz_i_old,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,NumberOfParticles,BoxLength,CGEpsilon,CGSigma,AASigma,pairList);
              currentParticleOldEnergy = result_initial.first;
              overLapFlagOld = result_initial.second;
              clusterEnergy_old = clusterEnergy_old + currentParticleOldEnergy;
              ProbFile << "currentParticleOldEnergy " << currentParticleOldEnergy << endl;
              particlesPositionList.push_back(Ri);
           }

           // find centroid of the cluster
           Vector3 centroid;
           double ClusterTemp[clusterSize][3];
           calculate_cluster_centroid(particlesPositionList,BoxLength,centroid,ClusterTemp);
           Vector3 n_vector;
           double D_value;
           // damping
           if (isTranslation == true) {
              n_vector[0] = VirtualMoveX;
              n_vector[1] = VirtualMoveY;
              n_vector[2] = VirtualMoveZ;
              double n_vector_norm = sqrt(n_vector[0]*n_vector[0] + n_vector[1]*n_vector[1] + n_vector[2]*n_vector[2]);
              n_vector[0] /= n_vector_norm;
              n_vector[1] /= n_vector_norm;
              n_vector[2] /= n_vector_norm;
              vector<double> RiMinusRcCrossN_vectorNormSquareList;
              vector<double> RadiiOfCubeInCluster;

              for (int i : clusterParticleList) {
                 RadiiOfCubeInCluster.push_back(Sizes[i]*AASigma);
                 auto it = find(clusterParticleList.begin(), clusterParticleList.end(), i);
                 if (it != clusterParticleList.end()) {
                    size_t index = std::distance(clusterParticleList.begin(), it);
                    double Rx_i_old = ClusterTemp[index][0];
                    double Ry_i_old = ClusterTemp[index][1];
                    double Rz_i_old = ClusterTemp[index][2];
                    Vector3 Ri = {Rx_i_old, Ry_i_old, Rz_i_old};
                    Vector3 Rc = {centroid[0], centroid[1], centroid[2]};
                    Vector3 RiminusRc = {Ri[0] - Rc[0], Ri[1] - Rc[1], Ri[2] - Rc[2]};
                    Vector3 RiMinusRcCrossN_vector = {RiminusRc[1]*n_vector[2]-RiminusRc[2]*n_vector[1],RiminusRc[2]*n_vector[0]-RiminusRc[0]*n_vector[2],RiminusRc[0]*n_vector[1]-RiminusRc[1]*n_vector[0]};
                    double RiMinusRcCrossN_vectorNorm = sqrt(RiMinusRcCrossN_vector[0]*RiMinusRcCrossN_vector[0]+RiMinusRcCrossN_vector[1]*RiMinusRcCrossN_vector[1]+RiMinusRcCrossN_vector[2]*RiMinusRcCrossN_vector[2]);
                    double RiMinusRcCrossN_vectorNormSquare = RiMinusRcCrossN_vectorNorm * RiMinusRcCrossN_vectorNorm;
                    RiMinusRcCrossN_vectorNormSquareList.push_back(RiMinusRcCrossN_vectorNormSquare);
                 }
              }
              double sumRadiiOfCubeInCluster = accumulate(RadiiOfCubeInCluster.begin(), RadiiOfCubeInCluster.end(), 0.0);
              double AverageRadiusOfCubeInCluster = sumRadiiOfCubeInCluster / RadiiOfCubeInCluster.size();
              double sumRiMinusRcCrossN_vectorNormSquare = accumulate(RiMinusRcCrossN_vectorNormSquareList.begin(),RiMinusRcCrossN_vectorNormSquareList.end(),0.0);
              double RiMinusRcCrossN_vectorNormSquareAverage = sumRiMinusRcCrossN_vectorNormSquare / RiMinusRcCrossN_vectorNormSquareList.size();
              double EffectiveHydrodynamicRadius = sqrt(RiMinusRcCrossN_vectorNormSquareAverage) + AverageRadiusOfCubeInCluster;
              double Dt = AverageRadiusOfCubeInCluster / EffectiveHydrodynamicRadius;
              D_value = Dt;
           } else if (isRotation == true) {
              if (VirtualMoveAlpha != 0) {
                 n_vector[0] = 1.0;
                 n_vector[1] = 0.0;
                 n_vector[2] = 0.0;
              } else if (VirtualMoveBeta != 0) {
                 n_vector[0] = 0.0;
                 n_vector[1] = 1.0;
                 n_vector[2] = 0.0;
              } else if (VirtualMoveGamma != 0) {
                 n_vector[0] = 0.0;
                 n_vector[1] = 0.0;
                 n_vector[2] = 1.0;
              } else { // if the rotation is set to 0
                 n_vector[0] = 0.0;
                 n_vector[1] = 0.0;
                 n_vector[2] = 1.0;
              }
              vector<double> RiMinusRcCrossN_vectorNormSquareList;
              vector<double> RadiiOfCubeInCluster;

              for (int i : clusterParticleList) {
                 RadiiOfCubeInCluster.push_back(Sizes[i]*AASigma);
                 auto it = find(clusterParticleList.begin(), clusterParticleList.end(), i);
                 if (it != clusterParticleList.end()) {
                    size_t index = std::distance(clusterParticleList.begin(), it);
                    double Rx_i_old = ClusterTemp[index][0];
                    double Ry_i_old = ClusterTemp[index][1];
                    double Rz_i_old = ClusterTemp[index][2];
                    Vector3 Ri = {Rx_i_old, Ry_i_old, Rz_i_old};
                    Vector3 Rc = {centroid[0], centroid[1], centroid[2]};
                    Vector3 RiminusRc = {Ri[0] - Rc[0], Ri[1] - Rc[1], Ri[2] - Rc[2]};
                    Vector3 RiMinusRcCrossN_vector = {RiminusRc[1]*n_vector[2]-RiminusRc[2]*n_vector[1],RiminusRc[2]*n_vector[0]-RiminusRc[0]*n_vector[2],RiminusRc[0]*n_vector[1]-RiminusRc[1]*n_vector[0]};
                    double RiMinusRcCrossN_vectorNorm = sqrt(RiMinusRcCrossN_vector[0]*RiMinusRcCrossN_vector[0]+RiMinusRcCrossN_vector[1]*RiMinusRcCrossN_vector[1]+RiMinusRcCrossN_vector[2]*RiMinusRcCrossN_vector[2]);
                    double RiMinusRcCrossN_vectorNormSquare = RiMinusRcCrossN_vectorNorm * RiMinusRcCrossN_vectorNorm;
                    RiMinusRcCrossN_vectorNormSquareList.push_back(RiMinusRcCrossN_vectorNormSquare);
                 }
              }
              double sumRadiiOfCubeInCluster = accumulate(RadiiOfCubeInCluster.begin(), RadiiOfCubeInCluster.end(), 0.0);
              double AverageRadiusOfCubeInCluster = sumRadiiOfCubeInCluster / RadiiOfCubeInCluster.size();
              double sumRiMinusRcCrossN_vectorNormSquare = accumulate(RiMinusRcCrossN_vectorNormSquareList.begin(),RiMinusRcCrossN_vectorNormSquareList.end(),0.0);
              double RiMinusRcCrossN_vectorNormSquareAverage = sumRiMinusRcCrossN_vectorNormSquare / RiMinusRcCrossN_vectorNormSquareList.size();
              double EffectiveHydrodynamicRadius = sqrt(RiMinusRcCrossN_vectorNormSquareAverage) + AverageRadiusOfCubeInCluster;
              double Dr = pow(AverageRadiusOfCubeInCluster / EffectiveHydrodynamicRadius,3); // rotation damping
              D_value = Dr;
           }
           vector<double> Rx_All_temp;
           vector<double> Ry_All_temp;
           vector<double> Rz_All_temp;
           double VectorX_All_temp[NumberOfParticles][3];
           double VectorY_All_temp[NumberOfParticles][3];
           double VectorZ_All_temp[NumberOfParticles][3];

           for (int k = 0; k < NumberOfParticles; ++k) {
              for (int m = 0; m < 3; ++m) {
                 VectorX_All_temp[k][m] = VectorX[k][m];
                 VectorY_All_temp[k][m] = VectorY[k][m];
                 VectorZ_All_temp[k][m] = VectorZ[k][m];
              }
              Rx_All_temp.push_back(Rx[k]);
              Ry_All_temp.push_back(Ry[k]);
              Rz_All_temp.push_back(Rz[k]);
           }
           Vector3 ParticleCentroid;
           for (int particleID : clusterParticleList) {
              i = particleID;
              ParticleCentroid[0] = Rx[i] - seedX;
              ParticleCentroid[1] = Ry[i] - seedY;
              ParticleCentroid[2] = Rz[i] - seedZ;
              apply_periodic_boundary(ParticleCentroid, BoxLength, ParticleCentroid);

              Vector3 ParticleVectorX, ParticleVectorY, ParticleVectorZ;
              for (int k = 0; k < 3; ++k) {
                 ParticleVectorX[k] = VectorX[i][k]; // x component of orientation vector
                 ParticleVectorY[k] = VectorY[i][k]; // y component of orientation vector
                 ParticleVectorZ[k] = VectorZ[i][k]; // z component of orientation vector
              }

              Matrix3x3 VirtualMoveRotationX = {{1, 0, 0}, {0, cos(VirtualMoveAlpha), -sin(VirtualMoveAlpha)}, {0, sin(VirtualMoveAlpha), cos(VirtualMoveAlpha)}};
              Matrix3x3 VirtualMoveRotationY = {{cos(VirtualMoveBeta), 0, sin(VirtualMoveBeta)}, {0, 1, 0}, {-sin(VirtualMoveBeta), 0, cos(VirtualMoveBeta)}};
              Matrix3x3 VirtualMoveRotationZ = {{cos(VirtualMoveGamma), -sin(VirtualMoveGamma), 0}, {sin(VirtualMoveGamma), cos(VirtualMoveGamma), 0}, {0, 0, 1}};

              Vector3 tempresult1, tempresult2;
              Vector3 RotatedParticleCentroid, VectorX_i_new, VectorY_i_new, VectorZ_i_new;
              MultiplyMatrixVector(VirtualMoveRotationX, ParticleCentroid, tempresult1);
              MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
              MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, RotatedParticleCentroid);

              MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorX, tempresult1);
              MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
              MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorX_i_new);

              MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorY, tempresult1);
              MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
              MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorY_i_new);

              MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorZ, tempresult1);
              MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
              MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorZ_i_new);

              double Rx_i_new = RotatedParticleCentroid[0] + seedX;
              double Ry_i_new = RotatedParticleCentroid[1] + seedY;
              double Rz_i_new = RotatedParticleCentroid[2] + seedZ;

              Rx_i_new = Rx_i_new + VirtualMoveX;
              Ry_i_new = Ry_i_new + VirtualMoveY;
              Rz_i_new = Rz_i_new + VirtualMoveZ;
              Vector3 R_i_new = {Rx_i_new,Ry_i_new,Rz_i_new};
              apply_periodic_boundary(R_i_new, BoxLength,R_i_new);
              Rx_i_new = R_i_new[0];
              Ry_i_new = R_i_new[1];
              Rz_i_new = R_i_new[2];

              Rx_All_temp[i] = Rx_i_new;
              Ry_All_temp[i] = Ry_i_new;
              Rz_All_temp[i] = Rz_i_new;
              for (int k = 0; k < 3; ++k) {
                 //cout << "Here" << endl;
                 //cout << VectorX_i_new[0] << " " << VectorX_i_new[1] << " " << VectorX_i_new[2] << endl;
                 VectorX_All_temp[i][k] = VectorX_i_new[k]; // x component of orientation vector
                 VectorY_All_temp[i][k] = VectorY_i_new[k]; // y component of orientation vector
                 VectorZ_All_temp[i][k] = VectorZ_i_new[k]; // z component of orientation vector
              }
           }
           //cout << "Rx temp: ";
           //for (int m : Rx_All_temp) {
           //    cout << m << endl;
           //}
           //cout << "Ry temp: ";
           //for (int m : Ry_All_temp) {
           //   cout << m << endl;
           //}
           //cout << "Rz temp: ";
           //for (int m : Rz_All_temp) {
           //   cout << m << endl;
           //}
           bool clusterOverlapAfterMove = false;
           pairList.clear();
           // energy of cluster after cluster move
           for (int particleID : clusterParticleList) {
              i = particleID;
              pairList.push_back(i);
              double Size_i = Sizes[i];
              double Rx_i_new = Rx_All_temp[i];
              double Ry_i_new = Ry_All_temp[i];
              double Rz_i_new = Rz_All_temp[i];
              Vector3 VectorX_i_new, VectorY_i_new, VectorZ_i_new;
              for (int k = 0; k < 3; ++k) {
                 VectorX_i_new[k] = VectorX_All_temp[i][k]; // x component of orientation vector
                 VectorY_i_new[k] = VectorY_All_temp[i][k]; // y component of orientation vector
                 VectorZ_i_new[k] = VectorZ_All_temp[i][k]; // z component of orientation vector
              }
              pair<double,bool> result_second = OneParticleEnergy(model,i,Size_i,cutoffCG,Hamaker,AtomDensity,Rx_i_new,Ry_i_new,Rz_i_new,Sizes,Rx_All_temp,Ry_All_temp,Rz_All_temp,VectorX_All_temp,VectorY_All_temp,VectorZ_All_temp,NumberOfParticles,BoxLength,CGEpsilon,CGSigma,AASigma,pairList);
              currentParticleNewEnergy = result_second.first;
              overLapFlagNew = result_second.second;
              if (overLapFlagNew || overLapFlagOld) {
                 clusterOverlapAfterMove = true;
              }
              clusterEnergy_new = clusterEnergy_new + currentParticleNewEnergy;
              ProbFile << "clusterEnergy_new " << clusterEnergy_new << endl;
           }
           double ClusterEnergyChange = clusterEnergy_new - clusterEnergy_old; // delta_V, difference between new and old total energy
           ProbFile << "clusterEnergy_new " << clusterEnergy_new << endl;
           ProbFile << "ClusterEnergyChange " << ClusterEnergyChange << endl;
           //cout << "ClusterEnergyChange " << ClusterEnergyChange << endl;
           //double PotentialEnergy = TotalEnergy(Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,NumberOfParticles,BoxLength,CGSigma,CGEpsilon,AASigma,model,cutoffCG, Hamaker, AtomDensity);
           //ProbFile << "RealEnergy Old " << PotentialEnergy << endl;
           //cout << "RealEnergy Old" << PotentialEnergy << endl;
           //cout << "OldEnergy " << SystemPotentialEnergy << endl;
           Attempt += 1;
           clusterSize = clusterParticleList.size();
           double RandomNumber = ((double)rand()/(RAND_MAX));
           ProbFile << "RandomNumber " << RandomNumber << endl;
           double BoltzmannFactor = exp(-ClusterEnergyChange / kBT);
           ProbFile << "BoltzmannFactor " << BoltzmannFactor << endl;
           if (MaxClusterSize == 1) {
              SingleMoveNumber += 1;
           } else {
              ClusterMoveNumber += 1;
           }
           ProbFile << "overallReverseProbability " << overallReverseProbability << endl;
           ProbFile << "overallUnacceptedReverseProbability " << overallUnacceptedReverseProbability << endl;
           ProbFile << "overallMoveProbability " << overallMoveProbability << endl;
           ProbFile << "overallUnacceptedMoveProbability " << overallUnacceptedMoveProbability << endl;
           ProbFile << "D_value " << D_value << endl;

           double overallReverse = overallReverseProbability * overallUnacceptedReverseProbability;
           double overallMove = overallMoveProbability * overallUnacceptedMoveProbability;
           double W_accept = D_value * min(1.0, overallReverse / overallMove);
           if (MaxClusterSize == 1) {
              W_accept = min(1.0,BoltzmannFactor);
              ProbFile << "W_accept " << W_accept << endl;
              ProbFile << "Single move case" << endl;
           }
           ProbFile << "W_accept " << W_accept << endl;
           ProbFile << "clusterOverlapAfterMove " << clusterOverlapAfterMove << endl;

           if (RandomNumber <= W_accept && !clusterOverlapAfterMove) {  // accept this move if it satisfies Boltzmann factor criteria, make sure no overlap
              SystemPotentialEnergy += ClusterEnergyChange;
              for (int particleID : clusterParticleList) {
                 i = particleID;
                 Rx[i] = Rx_All_temp[i];
                 Ry[i] = Ry_All_temp[i];
                 Rz[i] = Rz_All_temp[i];
                 for (int k = 0; k < 3; ++k) {
                    VectorX[i][k] = VectorX_All_temp[i][k];
                    VectorY[i][k] = VectorY_All_temp[i][k];
                    VectorZ[i][k] = VectorZ_All_temp[i][k];
                 }
              }
              Accept += 1;
              if (MaxClusterSize == 1) {
                 SingleMoveAccept += 1;
                 ProbFile << "Single move accepted!!!" << endl;
              } else {
                 ClusterMoveAccept += 1;
                 ProbFile << "Cluster move accepted!!!" << endl;
              }
              ProbFile << "Accept Move!!!" << endl;
           } else {
              ProbFile << "Reject cluster Move~~~" << endl;
           }
        }
        //cout << "Rx new: ";
        //for (double m : Rx) {
        //   cout << m << endl;
        //}
        //cout << "Ry new: ";
        //for (double m : Ry) {
        //   cout << m << endl;
        //}
        //cout << "Rz new: ";
        //for (double m : Rz) {
        //  cout << m << endl;
        //}
        lastRestartStep = writeRestart(RestartFileInterval,lastRestartStep,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,Step);
        writeEnergy(EnergyOutputInterval, Step, SystemPotentialEnergy, EnergyFile);
        writeTrajectory(TrajectoryInterval,Style,Step,NumberOfParticles,AASigma,CGCubeSize,CGSigma,BoxLengthHalf,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ);
        ProbFile << "New Energy " << SystemPotentialEnergy << endl;
        // cout <<  "New Energy " << SystemPotentialEnergy << endl;
        //double PotentialEnergy = TotalEnergy(Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,NumberOfParticles,BoxLength,CGSigma,CGEpsilon,AASigma,model,cutoffCG, Hamaker,AtomDensity);
        //ProbFile << "Real Energy New " << PotentialEnergy << endl;
        // cout <<  "Real Energy New" << PotentialEnergy << endl;
        double accept_ratio = Accept/Attempt;
        double cluster_accept_ratio = ClusterMoveAccept / ClusterMoveNumber;
        double single_accept_ratio = SingleMoveAccept / SingleMoveNumber;
        ProbFile << "Fraction of accepting moves " << accept_ratio << endl;
        if (ClusterMoveNumber > 0) {
           ProbFile << "Fraction of accepting cluster move " << cluster_accept_ratio << endl;
        }
        if (SingleMoveNumber>0) {
           ProbFile << "Fraction of accepting single move " << single_accept_ratio << endl;
        }

        /*
        vector<bool> inCluster(NumberOfParticles, false);
        for (int k = 0; k < NumberOfParticles; ++k) {
           for (int m = k+1; m < NumberOfParticles; ++m) {
              double pos_x = Rx[k]-Rx[m];
              double pos_y = Ry[k]-Ry[m];
              double pos_z = Rz[k]-Rz[m];
              pos_x = pos_x-BoxLength*round(pos_x/BoxLength);
              pos_y = pos_y-BoxLength*round(pos_y/BoxLength);
              pos_z = pos_z-BoxLength*round(pos_z/BoxLength);
              double dist = sqrt(pos_x*pos_x+pos_y*pos_y+pos_z*pos_z);
              //cout << "dist " << dist << endl;
              int bin_index = static_cast<int>(dist / bin_width);
              if (bin_index < num_bins) {
                 double bin_center = (bin_index+0.5)*bin_width;
                 //cout << "bin_center " << bin_center << endl;
                 //cout << "test" << normalized_bin_counts[bin_index] << endl;
                 normalized_bin_counts[bin_index] += 2.0/(normalization_constant*bin_center*bin_center);
              }
              if (dist < (CubeSideAASize*AASigma*sqrt(3))) {
                 mean_agg += 1;
              }
           }
        }*/

        /*MeanAggFile << mean_agg << endl;
        for (int k = 0; k < num_bins; ++k) {
           mean_bin_counts[k] += normalized_bin_counts[k];
           mean_bin_counts[k] = mean_bin_counts[k]/Step;
           RDFFile << normalized_bin_counts[k] << " " ;
           MeanRDFFile << mean_bin_counts[k] << " ";
        }
        RDFFile << endl;
        MeanRDFFile << endl;*/

    }
    LAMMPSTrajectoryFile.close();
    EnergyFile.close();
    return 0;
}



