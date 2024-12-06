#include <iostream>
#include <cmath>
#include "Helpers.h"
#include "Translation.h"
#include "Scaling.h"
#include "Rotation.h"
/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v)
{
    double d = magnitudeOfVec3(v);
    return Vec3(v.x / d, v.y / d, v.z / d);
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v)
{
    return Vec3(-v.x, -v.y, -v.z);
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c)
{
    return Vec3(v.x * c, v.y * c, v.z * c);
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v)
{
    std::cout << "(" << v.x << "," << v.y << "," << v.z << ")" << std::endl;
}

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b)
{

    /* if x difference, y difference and z difference is smaller than threshold, then they are equal */
    if ((ABS((a.x - b.x)) < EPSILON) && (ABS((a.y - b.y)) < EPSILON) && (ABS((a.z - b.z)) < EPSILON))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
 */
Matrix4 getIdentityMatrix()
{
    Matrix4 result;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                result.values[i][j] = 1.0;
            }
            else
            {
                result.values[i][j] = 0.0;
            }
        }
    }

    return result;
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2)
{
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            total = 0;
            for (int k = 0; k < 4; k++)
            {
                total += m1.values[i][k] * m2.values[k][j];
            }

            result.values[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = 0;
        for (int j = 0; j < 4; j++)
        {
            total += m.values[i][j] * v.getNthComponent(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}



Matrix4 translate(Translation *translation, Matrix4 matrix){
    Matrix4 translationMatrix = getIdentityMatrix();

    translationMatrix.values[0][3] = translation->tx;
    translationMatrix.values[1][3] = translation->ty;
    translationMatrix.values[2][3] = translation->tz;

    return translationMatrix*matrix;
}

Matrix4 scale(Scaling *scaling, Matrix4 matrix){
    Matrix4 scalingMatrix = getIdentityMatrix();
    scalingMatrix.values[0][0] = scaling->sx;
    scalingMatrix.values[1][1] = scaling->sy;
    scalingMatrix.values[2][2] = scaling->sz;

    return scalingMatrix*matrix;
}


Matrix4 rotate(Rotation *rotation, Matrix4 matrix){
    Vec3 u(rotation->ux, rotation->uy, rotation->uz), v, w;
    u = normalizeVec3(u);
    double angle_radians = (rotation->angle * M_PI)/180.0;
    double smallest = std::min({abs(u.x), abs(u.y), abs(u.z)});

    if (smallest == abs(u.x)){
        v.x = 0.0;
        v.y = -u.z;
        v.z = u.y;
    }else if(smallest == abs(u.y)){
        v.x = -u.z;
        v.y = 0.0;
        v.z = u.x;
    } else{
        v.x = u.y;
        v.y = -u.x;
        v.z = 0.0;
    }

    w = crossProductVec3(u, v);
    v = normalizeVec3(v);
    w = normalizeVec3(w);

    double m_inverse_values[4][4] = {{u.x, v.x, w.x, 0}, {u.y, v.y, w.y, 0}, {u.z, v.z, w.z, 0}, {0, 0, 0, 1}};
    Matrix4 m_inverse(m_inverse_values);

    double m_values[4][4] = {{u.x, u.y, u.z, 0}, {v.x, v.y, v.z, 0}, {w.x, w.y, w.z, 0}, {0, 0, 0, 1}};
    Matrix4 m(m_values);

    Matrix4 rotationMatrix = getIdentityMatrix();
    rotationMatrix.values[1][1] = cos(angle_radians);
    rotationMatrix.values[1][2] = sin(angle_radians)*(-1.0);
    rotationMatrix.values[2][1] = sin(angle_radians);
    rotationMatrix.values[2][2] = cos(angle_radians);

    return m_inverse*(rotationMatrix*(m*matrix));

}
