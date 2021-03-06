#ifndef SETTINGS_H
#define SETTINGS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef double REAL;

typedef Eigen::Matrix<REAL, 2,  2>  MATRIX2;
typedef Eigen::Matrix<REAL, 3,  3>  MATRIX3;
typedef Eigen::Matrix<REAL, 4,  4>  MATRIX4;
typedef Eigen::Matrix<REAL, 6,  6>  MATRIX6;
typedef Eigen::Matrix<REAL, 9,  9>  MATRIX9;
typedef Eigen::Matrix<REAL, 4,  6>  MATRIX4x6;
typedef Eigen::Matrix<REAL, 3,  12> MATRIX3x12;
typedef Eigen::Matrix<REAL, 9,  12> MATRIX9x12;
typedef Eigen::Matrix<REAL, 12, 12> MATRIX12;
typedef Eigen::Matrix<REAL, 2,  1>  VECTOR2;
typedef Eigen::Matrix<REAL, 3,  1>  VECTOR3;
typedef Eigen::Matrix<REAL, 4,  1>  VECTOR4;
typedef Eigen::Matrix<REAL, 6,  1>  VECTOR6;
typedef Eigen::Matrix<REAL, 9,  1>  VECTOR9;
typedef Eigen::Matrix<REAL, 12, 1>  VECTOR12;

typedef Eigen::Matrix<int, 2, 1> VECTOR2I;
typedef Eigen::Matrix<int, 3, 1> VECTOR3I;
typedef Eigen::Matrix<int, 4, 1> VECTOR4I;

typedef Eigen::Matrix<REAL, Eigen::Dynamic, 1> VECTOR;
typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> MATRIX;
typedef Eigen::SparseMatrix<REAL> SPARSE_MATRIX;

#endif
