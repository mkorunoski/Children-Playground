#include "utils_shadowing.h"
#include <math.h>

void ReduceToUnit(float vector[3])
{
    float length;
    length = (float)sqrt((vector[0]*vector[0]) +
                         (vector[1]*vector[1]) +
                         (vector[2]*vector[2]));
    if(length == 0.0f)
    length = 1.0f;
    vector[0] /= length;
    vector[1] /= length;
    vector[2] /= length;
}

void calcNormal(float v[3][3], float out[3])
{
    float v1[3],v2[3];
    static const int x = 0;
    static const int y = 1;
    static const int z = 2;
    v1[x] = v[0][x] - v[1][x];
    v1[y] = v[0][y] - v[1][y];
    v1[z] = v[0][z] - v[1][z];
    v2[x] = v[1][x] - v[2][x];
    v2[y] = v[1][y] - v[2][y];
    v2[z] = v[1][z] - v[2][z];
    out[x] = v1[y]*v2[z] - v1[z]*v2[y];
    out[y] = v1[z]*v2[x] - v1[x]*v2[z];
    out[z] = v1[x]*v2[y] - v1[y]*v2[x];
    ReduceToUnit(out);
}

void MakeShadowMatrix(GLfloat points[3][3],
                      GLfloat lightPos[4],
                      GLfloat destMat[4][4])
{
    GLfloat planeCoeff[4];
    GLfloat dot;
    calcNormal(points,planeCoeff);

    planeCoeff[3] = - ((planeCoeff[0]*points[2][0]) +
                       (planeCoeff[1]*points[2][1]) +
                       (planeCoeff[2]*points[2][2]));

    dot = planeCoeff[0] * lightPos[0] +
                        planeCoeff[1] * lightPos[1] +
                        planeCoeff[2] * lightPos[2] +
                        planeCoeff[3] * lightPos[3];

    destMat[0][0] = dot - lightPos[0] * planeCoeff[0];
    destMat[1][0] = 0.0f - lightPos[0] * planeCoeff[1];
    destMat[2][0] = 0.0f - lightPos[0] * planeCoeff[2];
    destMat[3][0] = 0.0f - lightPos[0] * planeCoeff[3];
    destMat[0][1] = 0.0f - lightPos[1] * planeCoeff[0];
    destMat[1][1] = dot - lightPos[1] * planeCoeff[1];
    destMat[2][1] = 0.0f - lightPos[1] * planeCoeff[2];
    destMat[3][1] = 0.0f - lightPos[1] * planeCoeff[3];
    destMat[0][2] = 0.0f - lightPos[2] * planeCoeff[0];
    destMat[1][2] = 0.0f - lightPos[2] * planeCoeff[1];
    destMat[2][2] = dot - lightPos[2] * planeCoeff[2];
    destMat[3][2] = 0.0f - lightPos[2] * planeCoeff[3];
    destMat[0][3] = 0.0f - lightPos[3] * planeCoeff[0];
    destMat[1][3] = 0.0f - lightPos[3] * planeCoeff[1];
    destMat[2][3] = 0.0f - lightPos[3] * planeCoeff[2];
    destMat[3][3] = dot - lightPos[3] * planeCoeff[3];
}
