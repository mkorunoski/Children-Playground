#include <GL\glut.h>

void ReduceToUnit(float vector[3]);
void calcNormal(float v[3][3], float out[3]);
void MakeShadowMatrix(GLfloat points[3][3],
                      GLfloat lightPos[4],
                      GLfloat destMat[4][4]);
