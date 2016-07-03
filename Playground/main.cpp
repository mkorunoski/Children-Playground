#include <iostream>
#include <math.h>
#include "utils.h"
#include "utils_shadowing.h"
#include "Camera.h"

using namespace std;

enum Objects
{
    FLOOR = 0,
    WATER,
    POOL,
    BASE,
    DOOR,
    SWING_IYPE1_FIRST,
    SWING_IYPE2_FIRST,
    SWING_IYPE2_SECOND,
    SWING_IYPE1_SECOND,
    CAROUSEL,
    PLAY_HOUSE,
    BENCH,
    PILAR
};

bool SelectedObjects[13] = {false};

#define PI 3.1415926535898
#define Cos(th) cos(PI/180*(th))
#define Sin(th) sin(PI/180*(th))

int wndWidth    = 960;
int wndHeight   = 720;
char * wndName  = "Mladen Korunoski 141032 - Hit SPACEBAR to toggle FPS mode";

Camera g_camera;
bool g_key[256];
bool g_fps_mode = false;

// Movement settings
const float g_translation_speed = 0.1;
const float g_rotation_speed = PI/180*0.2;

float shadowDensity = 0.3;

static GLuint sandTexture;
static GLuint wallTexture;
static GLuint plasticTexture;
static GLuint woodTexture;
static GLuint concreteTexture;
static GLuint metalTexture;
static GLuint roofTexture;

static GLuint tilesFloor;
static GLuint tilesSSide;
static GLuint tilesLSide;

static GLuint water;

static void drawBox(GLfloat size, GLenum type, bool texure, GLuint texName);
void APIENTRY glutSolidCube(GLdouble size, bool texure, GLuint texName);

void DrawShadows();
void DrawReflections();
void MoveLight();

void MouseClick();
void MouseDown(int x, int y, int but);
void Select(int x, int y);
void ListHits(GLint hits, GLuint *names);

void SelectObjects();

float CalculateNormal(float x, float x1, float x2, float y1, float y2)
{
    float m = (y2 - y1)/(x2 - x1);
    return y1 - (x - x1)/m;
}

void LoadTexture(char* fname, GLuint* texName)
{
    Image *img;

    img = (Image *) malloc(sizeof(Image));
    if (img == NULL)
    {
        printf("Error allocating space for image");
        exit(0);
    }

    if (!ImageLoad(fname, img))
    {
        exit(1);
    }

    glGenTextures(1, texName);
    glBindTexture(GL_TEXTURE_2D, *texName);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, img->sizeX, img->sizeY,
                 0, GL_RGB, GL_UNSIGNED_BYTE, img->data);
}

void Init()
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glEnable(GL_DEPTH_TEST);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glFrontFace(GL_CCW);

    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_RGB);

    LoadTexture("sand_texture.bmp", &sandTexture);
    LoadTexture("wall_texture.bmp", &wallTexture);
    LoadTexture("plastic_texture.bmp", &plasticTexture);
    LoadTexture("wood_texture.bmp", &woodTexture);
    LoadTexture("concrete_texture.bmp", &concreteTexture);
    LoadTexture("metal_texture.bmp", &metalTexture);
    LoadTexture("roof_texture.bmp", &roofTexture);
    LoadTexture("tiles_texture_f.bmp", &tilesFloor);
    LoadTexture("tiles_texture_l.bmp", &tilesLSide);
    LoadTexture("tiles_texture_s.bmp", &tilesSSide);
    LoadTexture("water_texture.bmp", &water);

}

void SetPerspectiveProjection()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float fp = 0.1;
    glFrustum(-fp, fp, -fp, fp, fp, 100);
}

void SetLighting()
{
    glEnable(GL_LIGHTING);

    glEnable(GL_NORMALIZE);
    glFrontFace(GL_CCW);

    /* light model - ambient */
    float intensity = 0.3;
    float lm_ambient[] = {intensity, intensity, intensity, 1.0};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lm_ambient);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    glEnable(GL_LIGHT3);
    glEnable(GL_LIGHT4);
    glEnable(GL_LIGHT5);
    glEnable(GL_LIGHT6);
    glEnable(GL_LIGHT7);
}

/* Material Properties */
float
	red[]      		= {0.5, 0.0, 0.0, 1.0},
	green[]    		= {0.0, 0.5, 0.0, 1.0},
	blue[]     		= {0.1, 0.5, 0.8, 1.0},
	yellow[]   		= {1.0, 0.9, 0.5, 1.0},
	white[]    		= {1.0, 1.0, 1.0, 1.0},
	white_trans[] 	= {1.0, 1.0, 1.0, 0.9},
    gray[]      	= {0.2, 0.2, 0.2, 1.0};
float
	m_specular[] 	= {1.0, 1.0, 1.0, 1.0};


/* Vertex Normals */

typedef struct Point3
{
    float x, y, z;
} Point3;
typedef struct Vector
{
    float x, y, z;
} Vector;
Vector Point2Vector(const Point3 * p1, const Point3 * p2)
{
    Vector A;
    A.x = p2->x - p1->x;
    A.y = p2->y - p1->y;
    A.z = p2->z - p1->z;
    return A;
}
void Normalize(Vector * A)
{
    float vectorMagnitude = sqrt( pow(A->x, 2) + pow(A->y, 2) + pow(A->z, 2) );
    A->x /= vectorMagnitude;
    A->y /= vectorMagnitude;
    A->z /= vectorMagnitude;
}
Vector CrossProduct(const Vector * A, const Vector * B)
{
    Vector product;
    product.x = A->y * B->z - B->y * A->z;
    product.y = A->z * B->x - B->z * A->x;
    product.z = A->x * B->y - B->x * A->y;
    Normalize(&product);
    return product;
}

/* End Vertex Normals */

void DrawPath()
{
    float u = 0.25, d = u + 0.125;
    glBegin(GL_LINE_STRIP);
        float i;
        int p;
        for(i = -1, p = 0; i <= 6; i += 1, p++)
        {
            if(p % 2 == 0)
                glVertex3f(u, 0.0, i);
            else
                glVertex3f(d, 0.0, i);
        }
    glEnd();
    glBegin(GL_LINE_STRIP);
        for(i = -1, p = 0; i <= 6; i += 1, p++)
        {
            if(p % 2 == 0)
                glVertex3f(-u, 0.0, i);
            else
                glVertex3f(-d, 0.0, i);
        }
    glEnd();
}

void DrawChain(int noRings = 10)
{
    double outerR = 0.25, innerR = 0.0625;
    int sides = 3, rings = 5;
    for(int i = 0; i < noRings; i++)
    {
        glPushMatrix();
        glTranslated(0.0, -i*(innerR + 2*outerR), 0.0);
        if(i % 2 != 0)
        {
            glRotated(90, 0.0, 1.0, 0.0);
        }
        glutSolidTorus(innerR, outerR, sides, rings);
        glPopMatrix();
    }
}

float rotAngleSwingT1F = 0;
float rotAngleSwingT1S = 0;

void PositionChain(int noRings = 80,
                   double Tx = 0.0, double Ty = 0.0, double Tz = 0.0,
                   double Sx = 0.02, double Sy = 0.02, double Sz = 0.02)
{
    glPushMatrix();
        glTranslated(Tx, Ty, Tz);
        glScaled(Sx, Sy, Sz);
        DrawChain(noRings);
    glPopMatrix();
}

void DrawCylinder(float radius = 0.1, float height = 1, int slices = 12, int stacks = 12)
{
    float i, k;
    int j;
    float l = height/slices;
    k = 360/stacks;
    for(i = 0; i < height - l; i += l)
    {
        for(j = 0; j < 360; j += k)
        {
            glBegin(GL_QUADS);
                Point3 p1, p2, p4;
                p1.x = Sin(j)*radius;
                p1.y = i;
                p1.z = Cos(j)*radius;
                p2.x = Sin(j+k)*radius;
                p2.y = i;
                p2.z = Cos(j+k)*radius;
                p4.x = Sin(j)*radius;
                p4.y = i+1;
                p4.z = Cos(j)*radius;
                Vector A, B;
                A = Point2Vector(&p1, &p2);
                B = Point2Vector(&p1, &p4);
                Vector product = CrossProduct(&A, &B);
                glNormal3f(product.x, product.y, product.z);
                glVertex3f(Sin(j)*radius, i, Cos(j)*radius);
                glVertex3f(Sin(j+k)*radius, i, Cos(j+k)*radius);
                glVertex3f(Sin(j+k)*radius, i+l, Cos(j+k)*radius);
                glVertex3f(Sin(j)*radius, i+l, Cos(j)*radius);
            glEnd();
        }
    }
    //Bottom Cap
    for(i = 0; i < 360; i += k)
    {
        glBegin(GL_TRIANGLES);
            //Normals
            glNormal3f(0, -1, 0);

            glVertex3f(0.0, 0.0, 0.0);
            glVertex3f(Cos(i)*radius, 0.0, Sin(i)*radius);
            glVertex3f(Cos(i+k)*radius, 0.0, Sin(i+k)*radius);
        glEnd();
    }
    //Top Cap
    for(i = 0; i < 360; i += k)
    {
        glBegin(GL_TRIANGLES);
            //Normals
            glNormal3f(0, 1, 0);

            glVertex3f(0.0, height, 0.0);
            glVertex3f(Cos(i+k)*radius, height, Sin(i+k)*radius);
            glVertex3f(Cos(i)*radius, height, Sin(i)*radius);
        glEnd();
    }
}
void PositionCylinder(double Tx = 0.0, double Ty = 0.0, double Tz = 0.0,
                      int angle = 0.0, double Rx = 0.0, double Ry = 0.0, double Rz = 0.0,
                      double Sx = 1.0, double Sy = 1.0, double Sz = 1.0)
{
    glPushMatrix();
    glTranslated(Tx, Ty, Tz);
    glRotated(angle, Rx, Ry, Rz);
    glScaled(Sx, Sy, Sz);
    DrawCylinder();
    glPopMatrix();
}

void DrawBase(float height = 0.125)
{
    glMaterialfv(GL_FRONT, GL_AMBIENT, white);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
    int i;
    const int corners = 6;
    float baseL[corners][3] =
    {
        {-0.25,     0.0,    -1.0},
        {-2.00,     0.0,     0.0},
        {-2.00,     0.0,     6.0},
        { 2.00,     0.0,     6.0},
        { 2.00,     0.0,     0.0},
        {+0.25,     0.0,    -1.0}
    };
    float baseH[corners][3] =
    {
        {-0.25,     height,     -1.0},
        {-2.00,     height,      0.0},
        {-2.00,     height,      6.0},
        { 2.00,     height,      6.0},
        { 2.00,     height,      0.0},
        {+0.25,     height,     -1.0}
    };

    glBindTexture(GL_TEXTURE_2D, wallTexture);

    float z, x;
    glBegin(GL_QUADS);
    //Normals
    z = baseL[0][2] + 0.1;
    x = CalculateNormal(z, baseL[0][2], baseL[1][2], baseL[0][0], baseL[1][0]);
    glNormal3f(x, 0, z);

    glTexCoord2f(0.0f,0.0f);
    glVertex3fv(baseL[1]);
    glTexCoord2f(1.0f,0.0f);
    glVertex3fv(baseL[0]);
    glTexCoord2f(1.0f,1.0f);
    glVertex3fv(baseH[0]);
    glTexCoord2f(0.0f,1.0f);
    glVertex3fv(baseH[1]);

    glEnd();
    glBegin(GL_QUADS);
    //Normals
    glNormal3f(1, 0, 0);

    glTexCoord2f(0.0f,0.0f);
    glVertex3fv(baseL[2]);
    glTexCoord2f(1.0f,0.0f);
    glVertex3fv(baseL[1]);
    glTexCoord2f(1.0f,1.0f);
    glVertex3fv(baseH[1]);
    glTexCoord2f(0.0f,1.0f);
    glVertex3fv(baseH[2]);

    glEnd();
    glBegin(GL_QUADS);
    //Normals
    glNormal3f(0, 0, -1);

    glTexCoord2f(0.0f,0.0f);
    glVertex3fv(baseL[3]);
    glTexCoord2f(1.0f,0.0f);
    glVertex3fv(baseL[2]);
    glTexCoord2f(1.0f,1.0f);
    glVertex3fv(baseH[2]);
    glTexCoord2f(0.0f,1.0f);
    glVertex3fv(baseH[3]);

    glEnd();
    glBegin(GL_QUADS);
    //Normals
    glNormal3f(-1, 0, 0);

    glTexCoord2f(0.0f,0.0f);
    glVertex3fv(baseL[4]);
    glTexCoord2f(1.0f,0.0f);
    glVertex3fv(baseL[3]);
    glTexCoord2f(1.0f,1.0f);
    glVertex3fv(baseH[3]);
    glTexCoord2f(0.0f,1.0f);
    glVertex3fv(baseH[4]);
    glEnd();

    glBegin(GL_QUADS);
    //Normals
    glNormal3f(-x, 0, z);

    glTexCoord2f(0.0f,0.0f);
    glVertex3fv(baseL[5]);
    glTexCoord2f(1.0f,0.0f);
    glVertex3fv(baseL[4]);
    glTexCoord2f(1.0f,1.0f);
    glVertex3fv(baseH[4]);
    glTexCoord2f(0.0f,1.0f);
    glVertex3fv(baseH[5]);

    glEnd();

    glBindTexture(GL_TEXTURE_2D, 0);
}

typedef struct Point
{
    float x, z;
} Point;
float LineEquation2P(float z, Point P1, Point P2)
{
    float m = (P2.x - P1.x)/(P2.z - P1.z);
    return (P1.x + m*(z - P1.z));
}

int k = 16;

void DrawFloor1()
{
    glMaterialfv(GL_FRONT, GL_AMBIENT, white);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
    const int corners = 4;
    float base[corners][3] =
    {
        //A
        {-0.25,     0.0,    -1.0},
        //B
        {-2.00,     0.0,     0.0},
        //C
        { 2.00,     0.0,     0.0},
        //D
        { 0.25,     0.0,    -1.0}
    };

    float L1 = base[3][0] - base[0][0];
    float L2 = base[2][0] - base[1][0];
    float L3 = base[1][2] - base[0][2];

    float L1Segments[k+1], L2Segments[k+1];
    float horLSegments[k+1];

    int i, j;
    float u, v;
    for(i = 0; i < k+1; i++)
    {
        L1Segments[i]   = base[0][0] + i * L1/k;
        L2Segments[i]   = base[1][0] + i * L2/k;
        horLSegments[i] = base[0][2] + i * L3/k;
    }
    float POINTS[k+1][k+1][3];
    for(i = 0; i < k+1; i++)
    {
        Point P1, P2;
        P1.x = L1Segments[i], P1.z = horLSegments[0];
        P2.x = L2Segments[i], P2.z = horLSegments[k];
        glPointSize(5);
        for(j = 0; j < k+1; j++)
        {
            float x = LineEquation2P(horLSegments[j], P1, P2);
            POINTS[i][j][0] = x;
            POINTS[i][j][1] = 0.0 ;
            POINTS[i][j][2] = horLSegments[j];
        }
    }

    float inc = 1.0/k;

    glBindTexture(GL_TEXTURE_2D, concreteTexture);
    for(i = 0, u = 0; i < k; ++i, u+=inc)
    {
        for(j = 0, v = 0; j < k; ++j, v+=inc)
        {
            glBegin(GL_QUADS);
            glNormal3f(0.0, 1.0, 0.0);

            glTexCoord2f(u, v);
            glVertex3f(POINTS[i][j][0],     POINTS[i][j][1],        POINTS[i][j][2]);

            glTexCoord2f(u, v+inc);
            glVertex3f(POINTS[i][j+1][0],   POINTS[i][j+1][1],      POINTS[i][j+1][2]);

            glTexCoord2f(u+inc, v+inc);
            glVertex3f(POINTS[i+1][j+1][0], POINTS[i+1][j+1][1],    POINTS[i+1][j+1][2]);

            glTexCoord2f(u+inc, v);
            glVertex3f(POINTS[i+1][j][0],   POINTS[i+1][j][1],      POINTS[i+1][j][2]);

            glEnd();
        }
    }
    glBindTexture(GL_TEXTURE_2D, 0);
}

void DrawFloor2()
{
    glMaterialfv(GL_FRONT, GL_AMBIENT, yellow);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, yellow);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
    const int corners = 4;
    float base[corners][3] =
    {
        //A base[0]
        {-2.00,     0.0,     0.0},
        //B base[1]
        {-2.00,     0.0,     6.0},
        //C base[2]
        { 2.00,     0.0,     6.0},
        //D base[3]
        { 2.00,     0.0,     0.0}
    };

    glBindTexture(GL_TEXTURE_2D, sandTexture);

    float L1 = base[3][0] - base[0][0];
    float L2 = base[1][2] - base[0][2];
    float i, j;
    for(i = base[0][0]; i < base[3][0]; i += L1/k)
    {
        for(j = base[0][2]; j < base[1][2]; j += L2/k)
        {
            if((i >= -1.5 && i < -0.5) &&
                    (j >= 3.375 && j < 5.625))
                continue;
            glBegin(GL_QUADS);
            glNormal3f(0.0, 1.0, 0.0);
            glTexCoord2f(translate(i, base[0][0],base[3][0], 0,1),
                         translate(j, base[0][2],base[1][2], 0,1));
            glVertex3f(i , 0.0 , j);
            glTexCoord2f(translate(i, base[0][0],base[3][0], 0,1),
                         translate(j + L2/k, base[0][2],base[1][2], 0,1));
            glVertex3f(i , 0.0 , j + L2/k);
            glTexCoord2f(translate(i + L1/k, base[0][0],base[3][0], 0,1),
                         translate(j + L2/k, base[0][2],base[1][2], 0,1));
            glVertex3f(i + L1/k , 0.0 , j + L2/k);
            glTexCoord2f(translate(i + L1/k, base[0][0],base[3][0], 0,1),
                         translate(j, base[0][2],base[1][2], 0,1));
            glVertex3f(i + L1/k , 0.0 , j);
            glEnd();
        }
    }
    glBindTexture(GL_TEXTURE_2D, 0);
}

void DrawPool()
{
    glMaterialfv(GL_FRONT, GL_AMBIENT, white_trans);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, white_trans);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
    glBindTexture(GL_TEXTURE_2D, tilesFloor);
    //  length = 2.25
    //  width = 1
    float corners[4][2] =
    {
        {-1.5, 5.625},
        {-0.5, 5.625},
        {-0.5, 3.375},
        {-1.5, 3.375}
    };

    float height = 0.25;


    glBegin(GL_QUADS);
    glNormal3f(0.0, 1.0, 0.0);

    glTexCoord2f(0.0, 1.0);
    glVertex3f(corners[0][0], -height, corners[0][1]);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(corners[1][0], -height, corners[1][1]);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(corners[2][0], -height, corners[2][1]);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(corners[3][0], -height, corners[3][1]);
    glEnd();
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindTexture(GL_TEXTURE_2D, tilesLSide);
    glBegin(GL_QUADS);
    glNormal3f(1.0, 0.0, 0.0);

    glTexCoord2f(1.0, 1.0);
    glVertex3f(corners[0][0], -height, corners[0][1]);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(corners[3][0], -height, corners[3][1]);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(corners[3][0], 0.0, corners[3][1]);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(corners[0][0], 0.0, corners[0][1]);
    glEnd();
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindTexture(GL_TEXTURE_2D, tilesLSide);
    glBegin(GL_QUADS);
    glNormal3f(-1.0, 0.0, 0.0);

    glTexCoord2f(1.0, 1.0);
    glVertex3f(corners[1][0], -height, corners[1][1]);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(corners[2][0], -height, corners[2][1]);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(corners[2][0], 0.0, corners[2][1]);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(corners[1][0], 0.0, corners[1][1]);
    glEnd();
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindTexture(GL_TEXTURE_2D, tilesSSide);
    glBegin(GL_QUADS);
    glNormal3f(0.0, 0.0, -1.0);

    glTexCoord2f(0.0, 1.0);
    glVertex3f(corners[0][0], -height, corners[0][1]);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(corners[1][0], -height, corners[1][1]);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(corners[1][0], 0.0, corners[1][1]);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(corners[0][0], 0.0, corners[0][1]);
    glEnd();
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindTexture(GL_TEXTURE_2D, tilesSSide);
    glBegin(GL_QUADS);
    glNormal3f(0.0, 0.0, 1.0);

    glTexCoord2f(0.0, 1.0);
    glVertex3f(corners[2][0], -height, corners[2][1]);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(corners[3][0], -height, corners[3][1]);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(corners[3][0], 0.0, corners[3][1]);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(corners[2][0], 0.0, corners[2][1]);
    glEnd();
    glBindTexture(GL_TEXTURE_2D, 0);
}

float offset    = 0;
bool waves      = false;

float CalculateHeight(float x, float y)
{
    float height = ( (cos(PI*(x-offset)) * sin(2*PI*(y+offset)) - sin(2*PI*(x-offset)) * cos(PI*(y+offset))) / 10.0 );
    if(height > 0.0)
        return 0.0;
    return height;
}

void DrawWater()
{
	white_trans[3] = 0.6;
    glMaterialfv(GL_FRONT, GL_AMBIENT, white_trans);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, white_trans);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
    glBindTexture(GL_TEXTURE_2D, water);
    float height = 0.05;

    /*
    float k = 0.5;
    float p = 0.25;

    for(float j = 5.625; j > 3.375; j -= p)
    {
        for(float i = -1.5; i < -0.5; i += k)
        {
            glBegin(GL_QUADS);
                glNormal3f(0.0, 1.0, 0.0);

                glTexCoord2f(0.0, 1.0); glVertex3f(i,       CalculateHeight(i, j) - height,              j);
                glTexCoord2f(1.0, 1.0); glVertex3f(i + k,   CalculateHeight(i+k, j) - height,            j);
                glTexCoord2f(1.0, 0.0); glVertex3f(i + k,   CalculateHeight(i+k, j-p) - height,      j - p);
                glTexCoord2f(0.0, 0.0); glVertex3f(i,       CalculateHeight(i, j-p) - height,        j - p);
            glEnd();
        }
    }
    */

    /* DRAW WATER USING BEZIER CURVES */

    float eqdH[4];
    float eqdV[4];
    int index = 0;
    for(float i = -1.5, j = 3.375; index < 4 ; i += 1/3.0, j += 2.25/3.0, ++index)
    {
        eqdH[index] = i;
        eqdV[index] = j;
    };

    float ctrlPoints[4][4][3] =
    {
        {
            {eqdH[0], CalculateHeight(eqdH[0], eqdV[0])-height, eqdV[0]},
            {eqdH[1], CalculateHeight(eqdH[1], eqdV[0])-height, eqdV[0]},
            {eqdH[2], CalculateHeight(eqdH[2], eqdV[0])-height, eqdV[0]},
            {eqdH[3], CalculateHeight(eqdH[3], eqdV[0])-height, eqdV[0]},
        },
        {
            {eqdH[0], CalculateHeight(eqdH[0], eqdV[1])-height, eqdV[1]},
            {eqdH[1], CalculateHeight(eqdH[1], eqdV[1])-height, eqdV[1]},
            {eqdH[2], CalculateHeight(eqdH[2], eqdV[1])-height, eqdV[1]},
            {eqdH[3], CalculateHeight(eqdH[2], eqdV[1])-height, eqdV[1]},
        },
        {
            {eqdH[0], CalculateHeight(eqdH[0], eqdV[2])-height, eqdV[2]},
            {eqdH[1], CalculateHeight(eqdH[1], eqdV[2])-height, eqdV[2]},
            {eqdH[2], CalculateHeight(eqdH[2], eqdV[2])-height, eqdV[2]},
            {eqdH[3], CalculateHeight(eqdH[3], eqdV[2])-height, eqdV[2]},
        },
        {
            {eqdH[0], CalculateHeight(eqdH[0], eqdV[3])-height, eqdV[3]},
            {eqdH[1], CalculateHeight(eqdH[1], eqdV[3])-height, eqdV[3]},
            {eqdH[2], CalculateHeight(eqdH[2], eqdV[3])-height, eqdV[3]},
            {eqdH[3], CalculateHeight(eqdH[3], eqdV[3])-height, eqdV[3]},
        },
    };

    GLfloat texPts[2][2][2] =
    {
        {{0.0, 0.0}, {0.0, 1.0}},
        {{1.0, 0.0}, {1.0, 1.0}}
    };

    //  INIT
    glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, 4,
           0, 1, 12, 4, &ctrlPoints[0][0][0]);
    glMap2f(GL_MAP2_TEXTURE_COORD_2, 0, 1, 2, 2,
           0, 1, 4, 2, &texPts[0][0][0]);
    glEnable(GL_MAP2_TEXTURE_COORD_2);
    glEnable(GL_MAP2_VERTEX_3);
    glMapGrid2f(20, 0.0, 1.0, 20, 0.0, 1.0);
    //  DISPLAY
    glPushMatrix();
    glEvalMesh2(GL_FILL, 0, 20, 0, 20);
    /*
    for (int j = 0; j <= 8; j++)
    {
        glBegin(GL_LINE_STRIP);
            for (int i = 0; i <= 30; i++)
                glEvalCoord2f((GLfloat)i/30.0, (GLfloat)j/8.0);
        glEnd();
        glBegin(GL_LINE_STRIP);
            for (int i = 0; i <= 30; i++)
                glEvalCoord2f((GLfloat)j/8.0, (GLfloat)i/30.0);
        glEnd();
   }
   */
   glPopMatrix();
   glBindTexture(GL_TEXTURE_2D, 0);
}

void DrawDoor()
{
    glMaterialfv(GL_FRONT, GL_AMBIENT, green);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
    glMaterialf(GL_FRONT, GL_SHININESS, 50);
    PositionCylinder(-0.25,0.0,-1, 0.0,0.0,0.0,0.0, 0.2,0.8,0.2);
    PositionCylinder( 0.25,0.0,-1, 0.0,0.0,0.0,0.0, 0.2,0.8,0.2);
    PositionCylinder( 0.25,0.8,-1,  90,0.0,0.0,1.0, 0.2,0.5,0.2);
}

void DrawSeat(double length = 1)
{
    glPushMatrix();
        glScaled(length, 0.1, 0.2);
        glutSolidCube(1, true, woodTexture);
        glTranslated(0.0, 0.0, 1.5);
        glutSolidCube(1, true, woodTexture);
        glTranslated(0.0, 0.0, -3);
        glutSolidCube(1, true, woodTexture);
    glPopMatrix();
    glPushMatrix();
        glTranslated(length/2.5, -0.1, 0.0);
        glRotated(90, 0.0, 1.0, 0.0);
        glScaled(0.8, 0.1, 0.2);
        glutSolidCube(1, true, woodTexture);
    glPopMatrix();
    glPushMatrix();
        glTranslated(-length/2.5, -0.1, 0.0);
        glRotated(90, 0.0, 1.0, 0.0);
        glScaled(0.8, 0.1, 0.2);
        glutSolidCube(1, true, woodTexture);
    glPopMatrix();

}
void PositionSeat(double length,
                  double Tx = 0.0, double Ty = 0.0, double Tz = 0.0,
                  double Sx = 1.0, double Sy = 1.0, double Sz = 1.0)
{
    glPushMatrix();
        glTranslated(Tx, Ty, Tz);
        glScaled(Sx, Sy, Sz);
        DrawSeat(length);
    glPopMatrix();
}

bool swing1 = false;
bool swing2 = false;
bool swing3 = false;
bool swing4 = false;

void DrawSwingFirst(bool boja)
{
    int angle = 8;
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        if(SelectedObjects[SWING_IYPE1_FIRST] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, green);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);

    PositionCylinder( 0.5,0.0,-0.125,  angle,1.0,0.0,0.0, 0.1,0.8,0.1);
    PositionCylinder( 0.5,0.0, 0.125, -angle,1.0,0.0,0.0, 0.1,0.8,0.1);
    PositionCylinder(-0.5,0.0,-0.125,  angle,1.0,0.0,0.0, 0.1,0.8,0.1);
    PositionCylinder(-0.5,0.0, 0.125, -angle,1.0,0.0,0.0, 0.1,0.8,0.1);
    PositionCylinder( 0.5,0.8, 0.000,     90,0.0,0.0,1.0, 0.1,1.0,0.1);

    int noRings = 55;
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        glMaterialf(GL_FRONT, GL_SHININESS, 50);
        if(SelectedObjects[SWING_IYPE1_FIRST] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, gray);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, gray);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);

    glPushMatrix();
        glTranslated(0.0, 0.8, 0.0);
        glRotated(rotAngleSwingT1F, 1.0, 0.0, 0.0);
        PositionChain(noRings,  0.33);
        PositionChain(noRings,  0.17);
    glPopMatrix();

    glPushMatrix();
        glTranslated(0.0, 0.8, 0.0);
        glRotated(-rotAngleSwingT1F, 1.0, 0.0, 0.0);
        PositionChain(noRings, -0.33);
        PositionChain(noRings, -0.17);
    glPopMatrix();



    int length = 1;
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        if(SelectedObjects[SWING_IYPE1_FIRST] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, white);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);

    glPushMatrix();
        glTranslated(0.0, 0.8, 0.0);
        glRotated(rotAngleSwingT1F, 1.0, 0.0, 0.0);
        PositionSeat(length,  0.25,-0.6,0.0, 0.2,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
        glTranslated(0.0, 0.8, 0.0);
        glRotated(-rotAngleSwingT1F, 1.0, 0.0, 0.0);
        PositionSeat(length, -0.25,-0.6,0.0, 0.2,0.2,0.2);
    glPopMatrix();
}
void DrawSwingSecond(bool boja)
{
    int angle = 8;
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        if(SelectedObjects[SWING_IYPE1_SECOND] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, green);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);

    PositionCylinder( 0.5,0.0,-0.125,  angle,1.0,0.0,0.0, 0.1,0.8,0.1);
    PositionCylinder( 0.5,0.0, 0.125, -angle,1.0,0.0,0.0, 0.1,0.8,0.1);
    PositionCylinder(-0.5,0.0,-0.125,  angle,1.0,0.0,0.0, 0.1,0.8,0.1);
    PositionCylinder(-0.5,0.0, 0.125, -angle,1.0,0.0,0.0, 0.1,0.8,0.1);
    PositionCylinder( 0.5,0.8, 0.000,     90,0.0,0.0,1.0, 0.1,1.0,0.1);

    int noRings = 55;
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        glMaterialf(GL_FRONT, GL_SHININESS, 50);
        if(SelectedObjects[SWING_IYPE1_SECOND] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, gray);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, gray);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);

    glPushMatrix();
        glTranslated(0.0, 0.8, 0.0);
        glRotated(rotAngleSwingT1S, 1.0, 0.0, 0.0);
        PositionChain(noRings,  0.33);
        PositionChain(noRings,  0.17);
    glPopMatrix();

    glPushMatrix();
        glTranslated(0.0, 0.8, 0.0);
        glRotated(-rotAngleSwingT1S, 1.0, 0.0, 0.0);
        PositionChain(noRings, -0.33);
        PositionChain(noRings, -0.17);
    glPopMatrix();

    int length = 1;
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        if(SelectedObjects[SWING_IYPE1_SECOND] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, white);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);

    glPushMatrix();
        glTranslated(0.0, 0.8, 0.0);
        glRotated(rotAngleSwingT1S, 1.0, 0.0, 0.0);
        PositionSeat(length,  0.25,-0.6,0.0, 0.2,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
        glTranslated(0.0, 0.8, 0.0);
        glRotated(-rotAngleSwingT1S, 1.0, 0.0, 0.0);
        PositionSeat(length, -0.25,-0.6,0.0, 0.2,0.2,0.2);
    glPopMatrix();
}

float
    a1 = 20.0,
    a2 = -30.0,
    a3 = 15.0,
    a4 = 20.0,
    a5 = -30.0,
    a6 = 15.0;

void DrawSwing2First(bool boja)
{
    double outerR = 0.25, innerR = 0.03125;
    int sides = 5, rings = 15;

    glPushMatrix();
        if(boja)
        {
            if(SelectedObjects[SWING_IYPE2_FIRST] == false)
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, green);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
            }
            else
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
            }
        }
        else
            glColor4f(0.0, 0.0, 0.0, shadowDensity);
        glTranslated(0.5, 0.0, 0.0);
        PositionCylinder(0.25,0.6,0.0,  90,0.0,0.0,1.0, 0.2,0.5,0.2);
        PositionCylinder(0.25,0.0,0.0, 0.0,0.0,0.0,0.0, 0.2,0.6,0.2);
        glPushMatrix();
            if(boja)
            {
                if(SelectedObjects[SWING_IYPE2_FIRST] == false)
                {
                    glMaterialfv(GL_FRONT, GL_AMBIENT, red);
                    glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
                }
                else
                {
                    glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                    glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                }
            }
            else
                glColor4f(0.0, 0.0, 0.0, shadowDensity);
            glTranslated(0.0, 0.6, 0.0);
            glRotated(a1, 1.0, 0.0, 0.0);
            PositionCylinder(0.0,0.0,-1.0,  90,1.0,0.0,0.0, 0.2,2.0,0.2);
                glPushMatrix();
                    glPushMatrix();
                        glScaled(0.5, 0.5, 0.5);
                        glPushMatrix();
                            glTranslated(0.0, 0.2, -1.5);
                            glRotated(-15, 1.0, 0.0, 0.0);
                            glutSolidTorus(innerR, outerR, sides, rings);
                        glPopMatrix();
                        glPushMatrix();
                            glTranslated(0.0, 0.2, 1.5);
                            glRotated(15, 1.0, 0.0, 0.0);
                            glutSolidTorus(innerR, outerR, sides, rings);
                        glPopMatrix();
                    glPopMatrix();
                glPopMatrix();
                glPushMatrix();
                    if(boja)
                    {
                        if(SelectedObjects[SWING_IYPE2_FIRST] == false)
                        {
                            glMaterialfv(GL_FRONT, GL_AMBIENT, white);
                            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
                        }
                        else
                        {
                            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                        }
                    }
                    else
                        glColor4f(0.0, 0.0, 0.0, shadowDensity);
                    glScaled(0.15, 0.05, 0.15);
                    glPushMatrix();
                        glTranslated(0.0, 0, 6.5);
                        glutSolidCube(1, true, woodTexture);
                    glPopMatrix();
                    glPushMatrix();
                        glTranslated(0.0, 0, -6.5);
                        glutSolidCube(1, true, woodTexture);
                    glPopMatrix();
                glPopMatrix();
        glPopMatrix();
    glPopMatrix();

    if(boja)
    {
        if(SelectedObjects[SWING_IYPE2_FIRST] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, green);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);
    PositionCylinder(-0.25,0.0,0.0, 0.0,0.0,0.0,0.0, 0.2,0.8,0.2);
    PositionCylinder( 0.25,0.0,0.0, 0.0,0.0,0.0,0.0, 0.2,0.8,0.2);
    PositionCylinder( 0.25,0.8,0.0,  90,0.0,0.0,1.0, 0.2,0.5,0.2);
    glPushMatrix();
        if(boja)
        {
            if(SelectedObjects[SWING_IYPE2_FIRST] == false)
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, red);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
            }
            else
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
            }
        }
        else
            glColor4f(0.0, 0.0, 0.0, shadowDensity);
        glTranslated(0.0, 0.8, 0.0);
        glRotated(a2, 1.0, 0.0, 0.0);
        PositionCylinder(0.0,0.0,-1.0,  90,1.0,0.0,0.0, 0.2,2.0,0.2);
            glPushMatrix();
                glPushMatrix();
                    glScaled(0.5, 0.5, 0.5);
                    glPushMatrix();
                        glTranslated(0.0, 0.2, -1.5);
                        glRotated(-15, 1.0, 0.0, 0.0);
                        glutSolidTorus(innerR, outerR, sides, rings);
                    glPopMatrix();
                    glPushMatrix();
                        glTranslated(0.0, 0.2, 1.5);
                        glRotated(15, 1.0, 0.0, 0.0);
                        glutSolidTorus(innerR, outerR, sides, rings);
                    glPopMatrix();
                glPopMatrix();
                glPushMatrix();
                    if(boja)
                    {
                        if(SelectedObjects[SWING_IYPE2_FIRST] == false)
                        {
                            glMaterialfv(GL_FRONT, GL_AMBIENT, white);
                            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
                        }
                        else
                        {
                            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                        }
                    }
                    else
                        glColor4f(0.0, 0.0, 0.0, shadowDensity);
                    glScaled(0.15, 0.05, 0.15);
                    glPushMatrix();
                        glTranslated(0.0, 0, 6.5);
                        glutSolidCube(1, true, woodTexture);
                    glPopMatrix();
                    glPushMatrix();
                        glTranslated(0.0, 0, -6.5);
                        glutSolidCube(1, true, woodTexture);
                    glPopMatrix();
                glPopMatrix();
            glPopMatrix();
        glPopMatrix();


    glPushMatrix();
        if(boja)
        {
            if(SelectedObjects[SWING_IYPE2_FIRST] == false)
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, green);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
            }
            else
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
            }
        }
        else
            glColor4f(0.0, 0.0, 0.0, shadowDensity);
        glTranslated(-0.5, 0.0, 0.0);
        PositionCylinder(-0.25,0.0,0.0, 0.0,0.0,0.0,0.0, 0.2,0.4,0.2);
        PositionCylinder( 0.25,0.4,0.0,  90,0.0,0.0,1.0, 0.2,0.5,0.2);
        glPushMatrix();
            if(boja)
            {
                if(SelectedObjects[SWING_IYPE2_FIRST] == false)
                {
                    glMaterialfv(GL_FRONT, GL_AMBIENT, red);
                    glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
                }
                else
                {
                    glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                    glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                }
            }
            else
                glColor4f(0.0, 0.0, 0.0, shadowDensity);
            glTranslated(0.0, 0.4, 0.0);
            glRotated(a3, 1.0, 0.0, 0.0);
            PositionCylinder(0.0,0.0,-1.0,  90,1.0,0.0,0.0, 0.2,2.0,0.2);
                glPushMatrix();
                    glPushMatrix();
                        glScaled(0.5, 0.5, 0.5);
                        glPushMatrix();
                            glTranslated(0.0, 0.2, -1.5);
                            glRotated(-15, 1.0, 0.0, 0.0);
                            glutSolidTorus(innerR, outerR, sides, rings);
                        glPopMatrix();
                        glPushMatrix();
                            glTranslated(0.0, 0.2, 1.5);
                            glRotated(15, 1.0, 0.0, 0.0);
                            glutSolidTorus(innerR, outerR, sides, rings);
                        glPopMatrix();
                    glPopMatrix();
                    glPushMatrix();
                        if(boja)
                        {
                            if(SelectedObjects[SWING_IYPE2_FIRST] == false)
                            {
                                glMaterialfv(GL_FRONT, GL_AMBIENT, white);
                                glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
                            }
                            else
                            {
                                glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                                glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                            }
                        }
                        else
                            glColor4f(0.0, 0.0, 0.0, shadowDensity);
                        glScaled(0.15, 0.05, 0.15);
                        glPushMatrix();
                            glTranslated(0.0, 0, 6.5);
                            glutSolidCube(1, true, woodTexture);
                        glPopMatrix();
                        glPushMatrix();
                            glTranslated(0.0, 0, -6.5);
                            glutSolidCube(1, true, woodTexture);
                        glPopMatrix();
                    glPopMatrix();
                glPopMatrix();
        glPopMatrix();

    glPopMatrix();
}

void DrawSwing2Second(bool boja)
{
    double outerR = 0.25, innerR = 0.03125;
    int sides = 5, rings = 15;

    glPushMatrix();
        if(boja)
        {
            if(SelectedObjects[SWING_IYPE2_SECOND] == false)
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, green);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
            }
            else
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
            }
        }
        else
            glColor4f(0.0, 0.0, 0.0, shadowDensity);
        glTranslated(0.5, 0.0, 0.0);
        PositionCylinder(0.25,0.6,0.0,  90,0.0,0.0,1.0, 0.2,0.5,0.2);
        PositionCylinder(0.25,0.0,0.0, 0.0,0.0,0.0,0.0, 0.2,0.6,0.2);
        glPushMatrix();
            if(boja)
            {
                if(SelectedObjects[SWING_IYPE2_SECOND] == false)
                {
                    glMaterialfv(GL_FRONT, GL_AMBIENT, red);
                    glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
                }
                else
                {
                    glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                    glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                }
            }
            else
                glColor4f(0.0, 0.0, 0.0, shadowDensity);
            glTranslated(0.0, 0.6, 0.0);
            glRotated(a4, 1.0, 0.0, 0.0);
            PositionCylinder(0.0,0.0,-1.0,  90,1.0,0.0,0.0, 0.2,2.0,0.2);
                glPushMatrix();
                    glPushMatrix();
                        glScaled(0.5, 0.5, 0.5);
                        glPushMatrix();
                            glTranslated(0.0, 0.2, -1.5);
                            glRotated(-15, 1.0, 0.0, 0.0);
                            glutSolidTorus(innerR, outerR, sides, rings);
                        glPopMatrix();
                        glPushMatrix();
                            glTranslated(0.0, 0.2, 1.5);
                            glRotated(15, 1.0, 0.0, 0.0);
                            glutSolidTorus(innerR, outerR, sides, rings);
                        glPopMatrix();
                    glPopMatrix();
                glPopMatrix();
                glPushMatrix();
                    if(boja)
                    {
                        if(SelectedObjects[SWING_IYPE2_SECOND] == false)
                        {
                            glMaterialfv(GL_FRONT, GL_AMBIENT, white);
                            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
                        }
                        else
                        {
                            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                        }
                    }
                    else
                        glColor4f(0.0, 0.0, 0.0, shadowDensity);
                    glScaled(0.15, 0.05, 0.15);
                    glPushMatrix();
                        glTranslated(0.0, 0, 6.5);
                        glutSolidCube(1, true, woodTexture);
                    glPopMatrix();
                    glPushMatrix();
                        glTranslated(0.0, 0, -6.5);
                        glutSolidCube(1, true, woodTexture);
                    glPopMatrix();
                glPopMatrix();
        glPopMatrix();
    glPopMatrix();

    if(boja)
    {
        if(SelectedObjects[SWING_IYPE2_SECOND] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, green);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);
    PositionCylinder(-0.25,0.0,0.0, 0.0,0.0,0.0,0.0, 0.2,0.8,0.2);
    PositionCylinder( 0.25,0.0,0.0, 0.0,0.0,0.0,0.0, 0.2,0.8,0.2);
    PositionCylinder( 0.25,0.8,0.0,  90,0.0,0.0,1.0, 0.2,0.5,0.2);
    glPushMatrix();
        if(boja)
        {
            if(SelectedObjects[SWING_IYPE2_SECOND] == false)
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, red);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
            }
            else
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
            }
        }
        else
            glColor4f(0.0, 0.0, 0.0, shadowDensity);
        glTranslated(0.0, 0.8, 0.0);
        glRotated(a5, 1.0, 0.0, 0.0);
        PositionCylinder(0.0,0.0,-1.0,  90,1.0,0.0,0.0, 0.2,2.0,0.2);
            glPushMatrix();
                glPushMatrix();
                    glScaled(0.5, 0.5, 0.5);
                    glPushMatrix();
                        glTranslated(0.0, 0.2, -1.5);
                        glRotated(-15, 1.0, 0.0, 0.0);
                        glutSolidTorus(innerR, outerR, sides, rings);
                    glPopMatrix();
                    glPushMatrix();
                        glTranslated(0.0, 0.2, 1.5);
                        glRotated(15, 1.0, 0.0, 0.0);
                        glutSolidTorus(innerR, outerR, sides, rings);
                    glPopMatrix();
                glPopMatrix();
                glPushMatrix();
                    if(boja)
                    {
                        if(SelectedObjects[SWING_IYPE2_SECOND] == false)
                        {
                            glMaterialfv(GL_FRONT, GL_AMBIENT, white);
                            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
                        }
                        else
                        {
                            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                        }
                    }
                    else
                        glColor4f(0.0, 0.0, 0.0, shadowDensity);
                    glScaled(0.15, 0.05, 0.15);
                    glPushMatrix();
                        glTranslated(0.0, 0, 6.5);
                        glutSolidCube(1, true, woodTexture);
                    glPopMatrix();
                    glPushMatrix();
                        glTranslated(0.0, 0, -6.5);
                        glutSolidCube(1, true, woodTexture);
                    glPopMatrix();
                glPopMatrix();
            glPopMatrix();
        glPopMatrix();


    glPushMatrix();
        if(boja)
        {
            if(SelectedObjects[SWING_IYPE2_SECOND] == false)
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, green);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
            }
            else
            {
                glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
            }
        }
        else
            glColor4f(0.0, 0.0, 0.0, shadowDensity);
        glTranslated(-0.5, 0.0, 0.0);
        PositionCylinder(-0.25,0.0,0.0, 0.0,0.0,0.0,0.0, 0.2,0.4,0.2);
        PositionCylinder( 0.25,0.4,0.0,  90,0.0,0.0,1.0, 0.2,0.5,0.2);
        glPushMatrix();
            if(boja)
            {
                if(SelectedObjects[SWING_IYPE2_SECOND] == false)
                {
                    glMaterialfv(GL_FRONT, GL_AMBIENT, red);
                    glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
                }
                else
                {
                    glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                    glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                }
            }
            else
                glColor4f(0.0, 0.0, 0.0, shadowDensity);
            glTranslated(0.0, 0.4, 0.0);
            glRotated(a6, 1.0, 0.0, 0.0);
            PositionCylinder(0.0,0.0,-1.0,  90,1.0,0.0,0.0, 0.2,2.0,0.2);
                glPushMatrix();
                    glPushMatrix();
                        glScaled(0.5, 0.5, 0.5);
                        glPushMatrix();
                            glTranslated(0.0, 0.2, -1.5);
                            glRotated(-15, 1.0, 0.0, 0.0);
                            glutSolidTorus(innerR, outerR, sides, rings);
                        glPopMatrix();
                        glPushMatrix();
                            glTranslated(0.0, 0.2, 1.5);
                            glRotated(15, 1.0, 0.0, 0.0);
                            glutSolidTorus(innerR, outerR, sides, rings);
                        glPopMatrix();
                    glPopMatrix();
                    glPushMatrix();
                        if(boja)
                        {
                            if(SelectedObjects[SWING_IYPE2_SECOND] == false)
                            {
                                glMaterialfv(GL_FRONT, GL_AMBIENT, white);
                                glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
                            }
                            else
                            {
                                glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
                                glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
                            }
                        }
                        else
                            glColor4f(0.0, 0.0, 0.0, shadowDensity);
                        glScaled(0.15, 0.05, 0.15);
                        glPushMatrix();
                            glTranslated(0.0, 0, 6.5);
                            glutSolidCube(1, true, woodTexture);
                        glPopMatrix();
                        glPushMatrix();
                            glTranslated(0.0, 0, -6.5);
                            glutSolidCube(1, true, woodTexture);
                        glPopMatrix();
                    glPopMatrix();
                glPopMatrix();
        glPopMatrix();

    glPopMatrix();
}

void DrawCarousel(bool boja)
{
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);
    PositionCylinder(0.0,0.05,0.0, 0.0,0.0,0.0,0.0, 3.0,0.05,3.0);
    PositionCylinder(0.0,0.30,0.0, 0.0,0.0,0.0,0.0, 2.0,0.05,2.0);
    PositionCylinder(0.0,0.05,0.0, 0.0,0.0,0.0,0.0, 3.0,0.05,3.0);
    PositionCylinder(0.0,0.05,0.0, 0.0,0.0,0.0,0.0, 0.2,0.50,0.2);
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        glMaterialf(GL_FRONT, GL_SHININESS, 50);
        if(SelectedObjects[CAROUSEL] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, gray);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, gray);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);
    int step = 90;
    for(int i = 0; i < 360; i += step)
    {
        glPushMatrix();
        glTranslated(Cos(i)/4, 0.0, Sin(i)/4);
        PositionCylinder(0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 0.1,0.5,0.1);
        glPopMatrix();
    }
    PositionCylinder(0.25,0.5,0.0,   90,0.0,0.0,1.0, 0.1,0.5,0.1);
    PositionCylinder(0.0,0.5,-0.25, 180,0.0,1.0,1.0, 0.1,0.5,0.1);
}

void DrawBench(bool boja)
{
    glPushMatrix();
    glTranslated(0.0, 1.0, 0.0);
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        if(SelectedObjects[BENCH] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, white);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);
    DrawSeat(2);
    if(boja)
    {

        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        glMaterialf(GL_FRONT, GL_SHININESS, 50);
        if(SelectedObjects[BENCH] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, yellow);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, yellow);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, shadowDensity);
    PositionCylinder(-(float)12/16,-1,-0.25, 0.0,0.0,0.0,0.0, 1.0,0.9,1.0);
    PositionCylinder(-(float)12/16,-1, 0.25, 0.0,0.0,0.0,0.0, 1.0,0.9,1.0);
    PositionCylinder( (float)12/16,-1, 0.25, 0.0,0.0,0.0,0.0, 1.0,0.9,1.0);
    PositionCylinder( (float)12/16,-1,-0.25, 0.0,0.0,0.0,0.0, 1.0,0.9,1.0);
    glPopMatrix();
}

void DrawPlayHouse(bool boja)
{
    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        glMaterialf(GL_FRONT, GL_SHININESS, 50);
        if(SelectedObjects[PLAY_HOUSE] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, green);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, 0.0);
    PositionCylinder(-0.25,0.0,-0.25, 0.0,0.0,0.0,0.0, 0.1,1.0,0.1);
    PositionCylinder(-0.25,0.0, 0.25, 0.0,0.0,0.0,0.0, 0.1,1.0,0.1);
    PositionCylinder( 0.25,0.0, 0.25, 0.0,0.0,0.0,0.0, 0.1,1.0,0.1);
    PositionCylinder( 0.25,0.0,-0.25, 0.0,0.0,0.0,0.0, 0.1,1.0,0.1);

    float thiknes = 0.05, length = 0.5;
    PositionCylinder(-0.25,1.0,-0.25,  90,1,0.0,0.0, thiknes,length,thiknes);
    PositionCylinder(-0.25,1.0, 0.25, -90,0.0,0.0,1, thiknes,length,thiknes);
    PositionCylinder( 0.25,1.0, 0.25, -90,1,0.0,0.0, thiknes,length,thiknes);
    PositionCylinder( 0.25,1.0,-0.25,  90,0.0,0.0,1, thiknes,length,thiknes);

    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        glMaterialf(GL_FRONT, GL_SHININESS, 50);
        if(SelectedObjects[PLAY_HOUSE] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, red);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, 0.0);
    PositionCylinder(-(float)3/8,1.0,0.0, -90,0.0,0.0,1.0, 0.1,0.125,0.1);
    PositionCylinder(-(float)3/8,0.0,0.0, 0.0,0.0,0.0,0.0, 0.1,1.000,0.1);
    PositionCylinder(0.25,0.8,0.25, 90,0.0,0.0,1, thiknes,length,thiknes);
    PositionCylinder(0.25,0.7,0.25, 90,0.0,0.0,1, thiknes,length,thiknes);

    PositionCylinder(-0.25,0.0,-0.6, 30,1,0.0,0.0, thiknes,0.69,thiknes);
    PositionCylinder(  0.0,0.0,-0.6, 30,1,0.0,0.0, thiknes,0.69,thiknes);
    PositionCylinder( 0.25,0.0,-0.6, 30,1,0.0,0.0, thiknes,0.69,thiknes);

    for(float i = -0.6; i < -0.25; i += 0.05)
    {
        float y = 0.6/0.35*(i + 0.6);
        PositionCylinder(0.25,y,i, 90,0.0,0.0,1.0, thiknes,length,thiknes);
    }

    int len = 7;
    float surface[2*len][2];
    int i;
    float value;
    for(i = 2*len + 1, value = len; i >= 0; --i, value -= 0.5)
    {
        surface[i][1] = (float)value/10;
        surface[i][0] = 1/(8*pow(surface[i][1], 2) + 1);
    }

    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        glMaterialf(GL_FRONT, GL_SHININESS, 50);
        if(SelectedObjects[PLAY_HOUSE] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, white);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, 0.0);

    glBindTexture(GL_TEXTURE_2D, plasticTexture);

    float step = 1.0/(2*len - 1);

    for(i = 0; i <= 2*len - 2; ++i)
    {
        glBegin(GL_QUADS);
        //Normals
        Point3 p1, p2, p4;
        p1.x = surface[i][0];
        p1.y = surface[i][1];
        p1.z = -0.14;
        p2.x = surface[i+1][0];
        p2.y = surface[i+1][1];
        p2.z = -0.14;
        p4.x = surface[i][0];
        p4.y = surface[i][1];
        p4.z = 0.14;
        Vector A, B;
        A = Point2Vector(&p1, &p2);
        B = Point2Vector(&p1, &p4);
        Vector product = CrossProduct(&A, &B);
        glNormal3f(product.x, product.y, product.z);

        glTexCoord2f(0.0f, i*step);
        glVertex3f(surface[i][0],       surface[i][1],      -0.14);

        glTexCoord2f(0.0f, (i+1)*step);
        glVertex3f(surface[i+1][0],     surface[i+1][1],    -0.14);

        glTexCoord2f(1.0f, (i+1)*step);
        glVertex3f(surface[i+1][0],     surface[i+1][1],     0.14);

        glTexCoord2f(1.0f, i*step);
        glVertex3f(surface[i][0],       surface[i][1],       0.14);
        glEnd();
    }

    glBindTexture(GL_TEXTURE_2D, 0);

    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        glMaterialf(GL_FRONT, GL_SHININESS, 50);
        if(SelectedObjects[PLAY_HOUSE] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, red);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, 0.0);
    glPushMatrix();
        glTranslated(0.0, 0.6, 0.0);
        glScaled(0.5, 0.01, 0.5);
        glutSolidCube(1, true, metalTexture);
    glPopMatrix();

    if(boja)
    {
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
        glMaterialf(GL_FRONT, GL_SHININESS, 50);
        if(SelectedObjects[PLAY_HOUSE] == false)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, white);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
        }
        else
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
        }
    }
    else
        glColor4f(0.0, 0.0, 0.0, 0.0);
    /* ROOF */
    glPushMatrix();
        glTranslated(0.0, 1.05, 0.2);
        glRotated(37.5, 1, 0.0, 0.0);
        glScaled(0.5, 0.005, 0.5);
        glutSolidCube(1, true, roofTexture);
    glPopMatrix();
    glPushMatrix();
        glTranslated(0.0, 1.05, -0.2);
        glRotated(-37.5, 1, 0.0, 0.0);
        glScaled(0.5, 0.005, 0.5);
        glutSolidCube(1, true, roofTexture);
    glPopMatrix();
}

void DrawPilar()
{
    const int len = 13;
    float contour[len][2] =
    {
        {0.10, 0.00},
        {0.09, 0.20},
        {0.08, 0.40},
        {0.07, 0.60},
        {0.06, 0.80},
        {0.05, 1.00},
        {0.11, 1.20},
        {0.12, 1.30},
        {0.13, 1.40},
        {0.14, 1.50},
        {0.07, 1.55},
        {0.07, 1.60},
        {0.15, 1.70}
    };
    int i, j;
    int step = 10;
    glMaterialfv(GL_FRONT, GL_AMBIENT, gray);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, gray);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
    glMaterialf(GL_FRONT, GL_SHININESS, 50);
    for(int i = 0; i < len - 1; i++)
    {
        for(j = 0; j < 360; j += step)
        {
            glBegin(GL_QUADS);
            //Normals
            glNormal3f(Cos(j), 0.0, Sin(j));

            glVertex3f(Cos(j+step)*contour[i][0],   contour[i][1],      Sin(j+step)*contour[i][0]);
            glVertex3f(Cos(j)*contour[i][0],        contour[i][1],      Sin(j)*contour[i][0]);
            glVertex3f(Cos(j)*contour[i+1][0],      contour[i+1][1],    Sin(j)*contour[i+1][0]);
            glVertex3f(Cos(j+step)*contour[i+1][0], contour[i+1][1],    Sin(j+step)*contour[i+1][0]);
            glEnd();
        }
    }
}

void BlackSquare(void)
{
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(0.0, 0.0, 0.0, shadowDensity);
    glBegin(GL_QUADS);
        glVertex3f(-1000,   0.01,   -1000);
        glVertex3f( 1000,   0.01,   -1000);
        glVertex3f( 1000,   0.01,    1000);
        glVertex3f(-1000,   0.01,    1000);
    glEnd();
    glDisable(GL_BLEND);

    glEnable(GL_TEXTURE_2D);
    glEnable(GL_LIGHTING);
    glPopAttrib();
}

/* Light Properties */
float
    ml_x = 0.0, ml_z = 0.0;
float
    amb_inten   = 0.2,
    diff_inten  = 0.8,
    spec_inten  = 1.0;
float
    l_ambient[]     = {amb_inten, amb_inten, amb_inten, 1.0},
    l_diffuse[]     = {diff_inten, diff_inten, diff_inten, 1.0},
    l_specular[]    = {spec_inten, spec_inten, spec_inten, 1.0};
float
    spot_pos[][4]  =
    {
        {1, 3, 0.5, 1.0},
        {1, 2.5, 2, 1.0},
        {1, 2.5, 3, 1.0},
        {1, 3, 4.5, 1.0},
        {-1, 3.0, 0.5, 1.0},
        {-1, 4.5, 2.5, 1.0},
        {-1, 3.0, 4.5, 1.0},
        {ml_x, 3.0, ml_z, 1.0}
    },
    spot_dir[] = {0.0, -1.0, 0.0};
/* End Light Properties */

void SetLights()
{
    for(int i = 0; i < 7; ++i)
    {
        glLightfv(GL_LIGHT0 + i, GL_AMBIENT, l_ambient);
        glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, l_diffuse);
        glLightfv(GL_LIGHT0 + i, GL_SPECULAR, l_specular);
        glLightf(GL_LIGHT0 + i, GL_SPOT_CUTOFF, 15);
        glLightf(GL_LIGHT0 + i, GL_SPOT_EXPONENT, 270);
        glLightfv(GL_LIGHT0 + i, GL_POSITION, spot_pos[i]);
        glLightfv(GL_LIGHT0 + i, GL_SPOT_DIRECTION, spot_dir);
    }
}

void Draw(GLenum mode)
{
    if (mode == GL_RENDER)
        SetPerspectiveProjection();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    g_camera.Refresh();

        SetLights();
        /* LIGHT7 - Rotating */
        glLightfv(GL_LIGHT7, GL_AMBIENT, l_ambient);
        glLightfv(GL_LIGHT7, GL_DIFFUSE, l_diffuse);
        glLightfv(GL_LIGHT7, GL_SPECULAR, l_specular);
        glLightf(GL_LIGHT7, GL_SPOT_CUTOFF, 45);
        glLightf(GL_LIGHT7, GL_SPOT_EXPONENT, 270);
        glLightfv(GL_LIGHT7, GL_POSITION, spot_pos[7]);
        glLightfv(GL_LIGHT7, GL_SPOT_DIRECTION, spot_dir);

        /* DRAWING */

        glEnable(GL_TEXTURE_2D);

        if (mode == GL_SELECT) glLoadName(FLOOR);
        DrawFloor1();
        DrawFloor2();

        /*  Reflections  */
        if(mode == GL_RENDER)
        {
            glDisable(GL_DEPTH_TEST);
                glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

                glEnable(GL_STENCIL_TEST);

                    glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE);
                    glStencilFunc(GL_ALWAYS, 1, 0xffffffff);

                    DrawWater();

                    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
                    glEnable(GL_DEPTH_TEST);

                    glStencilFunc(GL_EQUAL, 1, 0xffffffff);
                    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

                    glPushMatrix();
                        glScalef(1.0, -1.0, 1.0);

                        glCullFace(GL_FRONT);
                        DrawReflections();
                        glCullFace(GL_BACK);

                    glPopMatrix();

                glDisable(GL_STENCIL_TEST);
            glEnable(GL_DEPTH_TEST);
        }
        /*  End Reflections  */

        /*  Shadows  */
        if(mode == GL_RENDER)
        {
            glClear(GL_STENCIL_BUFFER_BIT);
            glDisable(GL_DEPTH_TEST);
            glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

            glEnable(GL_STENCIL_TEST);

            glStencilFunc(GL_ALWAYS, 3, 0xffffffff);
            glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE);

            DrawShadows();

            glStencilFunc(GL_ALWAYS, 1, 0xffffffff);
            glStencilOp(GL_DECR, GL_DECR, GL_DECR);
            glEnable(GL_DEPTH_TEST);

            DrawFloor1();
            DrawFloor2();

            glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

            glStencilFunc(GL_EQUAL, 2, 0xffffffff);
            glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
            BlackSquare();
            glDisable(GL_STENCIL_TEST);
        }
        /*  End Shadows  */

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_DST_COLOR);
            if (mode == GL_SELECT) glLoadName(POOL);
            DrawPool();

            if (mode == GL_SELECT) glLoadName(WATER);
            DrawWater();
        glDisable(GL_BLEND);

        if (mode == GL_SELECT) glLoadName(BASE);
        DrawBase();
        //DrawPath();

        if (mode == GL_SELECT) glLoadName(DOOR);
        DrawDoor();

        glPushMatrix();
            glTranslated(1, 0.0, 0.5);
            glRotated(90, 0.0, 1, 0.0);
            if (mode == GL_SELECT) glLoadName(SWING_IYPE1_FIRST);
            DrawSwingFirst(true);
        glPopMatrix();
        glPushMatrix();
            glTranslated(1, 0.0, 2);
            glRotated(90, 0.0, 1, 0.0);
            glScaled(0.5, 0.5, 0.5);
            if (mode == GL_SELECT) glLoadName(SWING_IYPE2_FIRST);
            DrawSwing2First(true);
        glPopMatrix();
        glPushMatrix();
            glTranslated(1, 0.0, 3);
            glRotated(90, 0.0, 1, 0.0);
            glScaled(0.5, 0.5, 0.5);
            if (mode == GL_SELECT) glLoadName(SWING_IYPE2_SECOND);
            DrawSwing2Second(true);
        glPopMatrix();
        glPushMatrix();
            glTranslated(1, 0.0, 4.5);
            glRotated(90, 0.0, 1, 0.0);
            if (mode == GL_SELECT) glLoadName(SWING_IYPE1_SECOND);
            DrawSwingSecond(true);
        glPopMatrix();
        glPushMatrix();
            glTranslated(-1, 0.0, 0.5);
            glRotated(90, 0.0, 1, 0.0);
            glScaled(0.8, 0.8, 0.8);
            if (mode == GL_SELECT) glLoadName(CAROUSEL);
            DrawCarousel(true);
        glPopMatrix();

        glPushMatrix();
            glTranslated(-1, 0.0, 2.5);
            if (mode == GL_SELECT) glLoadName(PLAY_HOUSE);
            DrawPlayHouse(true);
        glPopMatrix();
        /*
        glPushMatrix();
            glTranslated(-1, 0.0, 4.5);
            glRotated(90, 0.0, 1, 0.0);
            glScaled(0.8, 0.8, 0.8);
            DrawCarousel(true);
        glPopMatrix();
        */
        if (mode == GL_SELECT) glLoadName(BENCH);
        float benchDistance[2] = {0.5, 4.5};
        for(int i = 0; i < 2; i++)
        {
            glPushMatrix();
                glTranslated((float)5/3, 0.0, benchDistance[i]);
                glRotated(90, 0.0, 1, 0.0);
                glScaled(0.2, 0.2, 0.2);
                DrawBench(true);
            glPopMatrix();
            glPushMatrix();
                glTranslated(-(float)5/3, 0.0, benchDistance[i]);
                glRotated(90, 0.0, 1, 0.0);
                glScaled(0.2, 0.2, 0.2);
                DrawBench(true);
            glPopMatrix();
        }
        float pilarPos[4][3] =
        {
            {-2, 0.0, 0.0},
            {-2, 0.0, 6.0},
            { 2, 0.0, 6.0},
            { 2, 0.0, 0.0},
        };
        if (mode == GL_SELECT) glLoadName(PILAR);
        for(int i = 0; i < 4; i++)
        {
            glPushMatrix();
                glTranslated(pilarPos[i][0], pilarPos[i][1], pilarPos[i][2]);
                glScaled(0.5, 0.9, 0.5);
                DrawPilar();
            glPopMatrix();
        }

        /* Select Objects */
        SelectObjects();
        /* End Select Objects */
}

void Display()
{
    Draw(GL_RENDER);
    glFlush();
    glutSwapBuffers();
    glDisable(GL_TEXTURE_2D);
}

void DrawShadows()
{
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);

    glPushMatrix();
        GLfloat shadowMat[4][4];
        GLfloat points[4][3] =
        {
            {-2.00,     0.001,     0.0},
            {-2.00,     0.001,     6.0},
            { 2.00,     0.001,     6.0},
            { 2.00,     0.001,     0.0}
        };

        MakeShadowMatrix(points, spot_pos[7], shadowMat);
        glMultMatrixf((GLfloat *)shadowMat);
        // Swing1 1 shadow
        glPushMatrix();
            glTranslated(1, 0.0, 0.5);
            glRotated(90, 0.0, 1, 0.0);
            DrawSwingFirst(false);
        glPopMatrix();
        // Swing2 1 shadow
        glPushMatrix();
            glTranslated(1, 0.0, 2);
            glRotated(90, 0.0, 1, 0.0);
            glScaled(0.5, 0.5, 0.5);
            DrawSwing2First(false);
        glPopMatrix();
        // Swing2 2 shadow
        glPushMatrix();
            glTranslated(1, 0.0, 3);
            glRotated(90, 0.0, 1, 0.0);
            glScaled(0.5, 0.5, 0.5);
            DrawSwing2Second(false);
        glPopMatrix();
        // Swing1 2 shadow
        glPushMatrix();
            glTranslated(1, 0.0, 4.5);
            glRotated(90, 0.0, 1, 0.0);
            DrawSwingSecond(false);
        glPopMatrix();
        // Carousel1
        glPushMatrix();
            glTranslated(-1, 0.0, 0.5);
            glRotated(90, 0.0, 1, 0.0);
            glScaled(0.8, 0.8, 0.8);
            DrawCarousel(false);
        glPopMatrix();
        // Play House
        glPushMatrix();
            glTranslated(-1, 0.0, 2.5);
            DrawPlayHouse(false);
        glPopMatrix();
        // Carousel2
        /*
            glPushMatrix();
                glTranslated(-1, 0.0, 4.5);
                glRotated(90, 0.0, 1, 0.0);
                glScaled(0.8, 0.8, 0.8);
                DrawCarousel(false);
            glPopMatrix();
        */
        // Bench
        float benchDistance[2] = {0.5, 4.5};
        for(int i = 0; i < 2; i++)
        {
            glPushMatrix();
                glTranslated((float)5/3, 0.0, benchDistance[i]);
                glRotated(90, 0.0, 1, 0.0);
                glScaled(0.2, 0.2, 0.2);
                DrawBench(false);
            glPopMatrix();
            glPushMatrix();
                glTranslated(-(float)5/3, 0.0, benchDistance[i]);
                glRotated(90, 0.0, 1, 0.0);
                glScaled(0.2, 0.2, 0.2);
                DrawBench(false);
            glPopMatrix();
        }
    glPopMatrix();

    glPopAttrib();
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
}

void DrawReflections()
{
    glPushMatrix();
        glTranslated(-1, 0.0, 2.5);
        DrawPlayHouse(true);
    glPopMatrix();
    glPushMatrix();
        glTranslated(-(float)5/3, 0.0, 4.5);
        glRotated(90, 0.0, 1, 0.0);
        glScaled(0.2, 0.2, 0.2);
        DrawBench(true);
    glPopMatrix();
    glPushMatrix();
        glTranslated(1, 0.0, 4.5);
        glRotated(90, 0.0, 1, 0.0);
        DrawSwingSecond(true);
    glPopMatrix();
    glPushMatrix();
        glTranslated(1, 0.0, 3);
        glRotated(90, 0.0, 1, 0.0);
        glScaled(0.5, 0.5, 0.5);
        DrawSwing2Second(true);
    glPopMatrix();
}

void MoveLight()
{
    spot_pos[7][0] = ml_x;
    spot_pos[7][2] = ml_z * 2 + 2.5;
    glLightfv(GL_LIGHT7, GL_POSITION, spot_pos[7]);
}

int factor  = 0;
int dAlphaF  = 0;
int dAlphaS  = 0;

int
    dA1 = 90,
    dA2 = 90,
    dA3 = 90,
    dA4 = 90,
    dA5 = 90,
    dA6 = 90;

void Idle(void)
{
    /*  Animate light */
    factor += 1;
    ml_x = Cos(factor);
    ml_z = Sin(factor);
    MoveLight();
    /* Animate water */
    if(waves)
        offset += 0.01;
    /* Animate Swing */
    if(swing1)
    {
        rotAngleSwingT1F = 30 * Sin(dAlphaF);
        dAlphaF += 10;
        if(dAlphaF >= 360)
            dAlphaF = 0;
    }
    if(swing2)
    {
        a1 = 20 * Sin(dA1);
        dA1 += 10;
        if(dA1 >= 360)
            dA1 = 0;

        a2 = -30 * Sin(dA2);
        dA2 += 10;
        if(dA2 >= 360)
            dA2 = 0;

        a3 = 15 * Sin(dA3);
        dA3 += 10;
        if(dA3 >= 360)
            dA3 = 0;
    }
    if(swing3)
    {
        a4 = 20 * Sin(dA4);
        dA4 += 10;
        if(dA4 >= 360)
            dA4 = 0;

        a5 = -30 * Sin(dA5);
        dA5 += 10;
        if(dA5 >= 360)
            dA5 = 0;

        a6 = 15 * Sin(dA6);
        dA6 += 10;
        if(dA6 >= 360)
            dA6 = 0;
    }
    if(swing4)
    {
        rotAngleSwingT1S = 30 * Sin(dAlphaS);
        dAlphaS += 10;
        if(dAlphaS >= 360)
            dAlphaS = 0;
    }

    glutPostRedisplay();
}

void Reshape(int width, int height)
{
    glViewport(0, 0, width, height);
    glutReshapeWindow(wndWidth, wndHeight);
}

void WindowKeyboard(unsigned char key, int x, int y)
{
    if (key == ' ')
    {
        g_fps_mode = !g_fps_mode;

        if(g_fps_mode) {
            glutSetCursor(GLUT_CURSOR_NONE);
            glutWarpPointer(wndWidth/2, wndHeight/2);
        }
        else {
            glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
        }
    }

    g_key[key] = true;

    glutPostRedisplay();
}

void WindowKeyboardUp(unsigned char key, int x, int y)
{
    g_key[key] = false;
}

void MouseMotion(int x, int y)
{
    static bool just_warped = false;

    if(just_warped) {
        just_warped = false;
        return;
    }

    if(g_fps_mode)
    {
        int dx = (x - wndWidth/2);
        int dy = -(y - wndHeight/2);

        if(dx)
            g_camera.RotateYaw(g_rotation_speed*dx);

        if(dy)
            g_camera.RotatePitch(g_rotation_speed*dy);

        glutWarpPointer(wndWidth/2, wndHeight/2);

        just_warped = true;
    }
}

void MouseClick(int button, int state, int x, int y)
{
    if(state == GLUT_DOWN)
    {
        if(button == GLUT_LEFT_BUTTON)
        {
            MouseDown(x, y, button);
        }
    }
}

void MouseDown(int x, int y, int button)
{
    cout << "Mouse button " << button << " pressed at " << x << " " << y << endl;
    Select(x, y);
}

void Select(int x, int y)
{
    GLuint buff[128] = {0};
    GLint hits, view[4];

    glSelectBuffer(128, buff);

    glGetIntegerv(GL_VIEWPORT, view);

    glRenderMode(GL_SELECT);

    glInitNames();
    glPushName(0);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
        glLoadIdentity();

        gluPickMatrix(x, view[3] - y, 0.5, 0.5, view);
        float fp = 0.1;
        glFrustum(-fp, fp, -fp, fp, fp, 100);

        glMatrixMode(GL_MODELVIEW);

        glutSwapBuffers();
        Draw(GL_SELECT);

        glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    hits = glRenderMode(GL_RENDER);

    ListHits(hits, buff);

    glMatrixMode(GL_MODELVIEW);
}

void ListHits(GLint hits, GLuint *names)
{
    for(int i = 0; i < hits; ++i)
    {
        cout << "Selected: " << (Objects)names[i * 4 + 3] << endl;
        switch(names[i * 4 + 3])
        {
            case SWING_IYPE1_FIRST:
            {
                swing1 = !swing1;
                SelectedObjects[SWING_IYPE1_FIRST] = !SelectedObjects[SWING_IYPE1_FIRST];
                break;
            }
            case SWING_IYPE2_FIRST:
            {
                swing2 = !swing2;
                SelectedObjects[SWING_IYPE2_FIRST] = !SelectedObjects[SWING_IYPE2_FIRST];
                break;
            }
            case SWING_IYPE2_SECOND:
            {
                swing3 = !swing3;
                SelectedObjects[SWING_IYPE2_SECOND] = !SelectedObjects[SWING_IYPE2_SECOND];
                break;
            }
            case SWING_IYPE1_SECOND:
            {
                swing4 = !swing4;
                SelectedObjects[SWING_IYPE1_SECOND] = !SelectedObjects[SWING_IYPE1_SECOND];
                break;
            }
            case CAROUSEL:
            {
                SelectedObjects[CAROUSEL] = !SelectedObjects[CAROUSEL];
                break;
            }
            case PLAY_HOUSE:
            {
                SelectedObjects[PLAY_HOUSE] = !SelectedObjects[PLAY_HOUSE];
                break;
            }
            case POOL:
            {
                waves = !waves;
                SelectedObjects[POOL] = !SelectedObjects[POOL];
                break;
            }
            case BENCH:
            {
                SelectedObjects[BENCH] = !SelectedObjects[BENCH];
                break;
            }
        }
    }
    cout << endl;
}

void SelectObjects()
{
    glMaterialfv(GL_FRONT, GL_AMBIENT, white);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular);
    if(SelectedObjects[SWING_IYPE1_FIRST] == true)
    {
        glPushMatrix();
            glScaled(0.6, 1.25, 1.25);
            glTranslated(1.65, 0.3, 0.4);
            glutWireCube(1);
        glPopMatrix();
    }
    if(SelectedObjects[SWING_IYPE2_FIRST] == true)
    {
        glPushMatrix();
            glScaled(1.25, 0.9, 0.8);
            glTranslated(0.75, 0.25, 2.5);
            glutWireCube(1);
        glPopMatrix();
    }
    if(SelectedObjects[SWING_IYPE2_SECOND] == true)
    {
        glPushMatrix();
            glScaled(1.25, 0.9, 0.8);
            glTranslated(0.75, 0.25, 3.7);
            glutWireCube(1);
        glPopMatrix();
    }
    if(SelectedObjects[SWING_IYPE1_SECOND] == true)
    {
        glPushMatrix();
            glScaled(0.6, 1.25, 1.25);
            glTranslated(1.65, 0.3, 3.6);
            glutWireCube(1);
        glPopMatrix();
    }
    if(SelectedObjects[CAROUSEL] == true)
    {
        glPushMatrix();
            glScaled(0.6, 0.75, 0.6);
            glTranslated(-1.65, 0.2, 0.85);
            glutWireCube(1);
        glPopMatrix();
    }
    if(SelectedObjects[PLAY_HOUSE] == true)
    {
        glPushMatrix();
            glScaled(1.0, 1.5, 1.0);
            glTranslated(-1.0, 0.35, 2.5);
            glutWireCube(1);
        glPopMatrix();
    }
    if(SelectedObjects[POOL] == true)
    {
        glPushMatrix();
            glScaled(1.0, 0.75, 2.5);
            glTranslated(-1.0, 0.2, 1.8);
            glutWireCube(1);
        glPopMatrix();
    }
    if(SelectedObjects[BENCH] == true)
    {
        float pos[4][2] =
        {
            {-4.8, 0.85},
            {-4.8, 7.5},
            {4.8, 7.5},
            {4.8, 0.85}
        };
        for(int i = 0; i < 4; ++i)
        {
            glPushMatrix();
                glScaled(0.35, 0.6, 0.6);
                glTranslated(pos[i][0], 0.15, pos[i][1]);
                glutWireCube(1);
            glPopMatrix();
        }
    }
}

void Timer(int value)
{
    if(g_fps_mode) {
        if(g_key['w'] || g_key['W']) {
            g_camera.Move(g_translation_speed);
        }
        else if(g_key['s'] || g_key['S']) {
            g_camera.Move(-g_translation_speed);
        }
        else if(g_key['a'] || g_key['A']) {
            g_camera.Strafe(g_translation_speed);
        }
        else if(g_key['d'] || g_key['D']) {
            g_camera.Strafe(-g_translation_speed);
        }
    }

    glutTimerFunc(1, Timer, 0);
}

int main(int argc, char ** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA | GLUT_STENCIL);
    glutInitWindowSize(wndWidth, wndHeight);
    glutInitWindowPosition(100, 100);
    glutCreateWindow(wndName);

    Init();
    SetLighting();

    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(WindowKeyboard);
    glutKeyboardUpFunc(WindowKeyboardUp);

    glutMouseFunc(MouseClick);
    glutMotionFunc(MouseMotion);
    glutPassiveMotionFunc(MouseMotion);

    glutIdleFunc(Idle);

    glutTimerFunc(1, Timer, 0);

    glutMainLoop();
    return 0;
}

static void
drawBox(GLfloat size, GLenum type, bool texture, GLuint texName)
{
    static GLfloat n[6][3] =
    {
        {-1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, -1.0, 0.0},
        {0.0, 0.0, 1.0},
        {0.0, 0.0, -1.0}
    };
    static GLint faces[6][4] =
    {
        {0, 1, 2, 3},
        {3, 2, 6, 7},
        {7, 6, 5, 4},
        {4, 5, 1, 0},
        {5, 6, 2, 1},
        {7, 4, 0, 3}
    };
    GLfloat v[8][3];
    GLint i;

    v[0][0] = v[1][0] = v[2][0] = v[3][0] = -size / 2;
    v[4][0] = v[5][0] = v[6][0] = v[7][0] =  size / 2;
    v[0][1] = v[1][1] = v[4][1] = v[5][1] = -size / 2;
    v[2][1] = v[3][1] = v[6][1] = v[7][1] =  size / 2;
    v[0][2] = v[3][2] = v[4][2] = v[7][2] = -size / 2;
    v[1][2] = v[2][2] = v[5][2] = v[6][2] =  size / 2;

    if(texture)
    {
        glBindTexture(GL_TEXTURE_2D, texName);
        for (i = 5; i >= 0; i--)
        {
            glBegin(type);
            glNormal3fv(&n[i][0]);

            glTexCoord2f(0.0, 0.0);
            glVertex3fv(&v[faces[i][0]][0]);
            glTexCoord2f(0.0, 1.0);
            glVertex3fv(&v[faces[i][1]][0]);
            glTexCoord2f(1.0, 1.0);
            glVertex3fv(&v[faces[i][2]][0]);
            glTexCoord2f(1.0, 0.0);
            glVertex3fv(&v[faces[i][3]][0]);
            glEnd();
        }
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    else
    {
        for (i = 5; i >= 0; i--)
        {
            glBegin(type);
            glNormal3fv(&n[i][0]);
            glVertex3fv(&v[faces[i][0]][0]);
            glVertex3fv(&v[faces[i][1]][0]);
            glVertex3fv(&v[faces[i][2]][0]);
            glVertex3fv(&v[faces[i][3]][0]);
            glEnd();
        }
    }
}

void APIENTRY
glutSolidCube(GLdouble size, bool texture, GLuint texName)
{
    drawBox(size, GL_QUADS, texture, texName);
}
