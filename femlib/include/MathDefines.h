// MathDefines.h
//
// Breannan Smith
// Last updated: 09/03/2015

#ifndef MATH_DEFINES_H
#define MATH_DEFINES_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

#include <climits>

namespace MathDefines
{
  template<typename T> inline T PI()
  {
    return T(3.14159265358979323846);
  }
}

// Special scalar values (infinity, not-a-number, etc)
// TODO: Actually check that infinity and not-a-number are supported by selected type
#define SCALAR_INFINITY std::numeric_limits<scalar>::infinity()
#define SCALAR_NAN std::numeric_limits<scalar>::signaling_NaN()

// TODO: Some of these are incosistent with others (e.g. MatrixXs vs Matrix3Xic), fix that

typedef double scalar;

// Define integer-valued version of Eigen types
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Matrix<unsigned,1,1> Vector1u;
typedef Eigen::Matrix<unsigned,3,1> Vector3u;
typedef Eigen::Matrix<unsigned,4,1> Vector4u;
typedef Eigen::VectorXi VectorXi;
typedef Eigen::Matrix<unsigned,Eigen::Dynamic,1> VectorXu;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> MatrixXXic;
typedef Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> MatrixX3ir;
typedef Eigen::Matrix<int,2,Eigen::Dynamic,Eigen::ColMajor> Matrix2Xic;
typedef Eigen::Matrix<int,3,Eigen::Dynamic,Eigen::ColMajor> Matrix3Xic;
typedef Eigen::Matrix<unsigned,3,Eigen::Dynamic,Eigen::ColMajor> Matrix3Xuc;

typedef Eigen::Array<unsigned, 2, 1> Array2u;
typedef Eigen::Array<unsigned, 3, 1> Array3u;
typedef Eigen::Array<int, 3, 1> Array3i;

typedef Eigen::Matrix<scalar,2,2> Matrix2s;
typedef Eigen::Matrix<scalar,3,3> Matrix3s;
typedef Eigen::Matrix<scalar,4,4> Matrix4s;
typedef Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> MatrixXs;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> MatrixXXfc;

// Define scalar-valued versions of the Eigen vector types
typedef Eigen::Matrix<scalar, 2, 1> Vector2s;
typedef Eigen::Matrix<scalar, 3, 1> Vector3s;
typedef Eigen::Matrix<scalar, 4, 1> Vector4s;
typedef Eigen::Matrix<scalar, 6, 1> Vector6s;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> VectorXs;

typedef Eigen::Array<scalar, 2, 1> Array2s;
typedef Eigen::Array<scalar, 3, 1> Array3s;
typedef Eigen::Array<scalar, 4, 1> Array4s;
typedef Eigen::Array<scalar, Eigen::Dynamic, 1> ArrayXs;


typedef Eigen::Quaternion<scalar> Quaternions;

// TODO: Clean this up :)
// Define scalar-valued versions of the Eigen matrix types. r postfix denotes row-major, c postfix denotes column-major.
typedef Eigen::Matrix<scalar,2,2,Eigen::RowMajor> Matrix22sr;
typedef Eigen::Matrix<scalar,2,2,Eigen::ColMajor> Matrix22sc;
typedef Eigen::Matrix<scalar,3,3,Eigen::RowMajor> Matrix33sr;
typedef Eigen::Matrix<scalar,3,3,Eigen::ColMajor> Matrix33sc;
typedef Eigen::Matrix<scalar,6,6,Eigen::ColMajor> Matrix66sc;
typedef Eigen::Matrix<scalar,2,Eigen::Dynamic,Eigen::ColMajor> Matrix2Xsc;
typedef Eigen::Matrix<scalar,3,Eigen::Dynamic,Eigen::ColMajor> Matrix3Xsc;
typedef Eigen::Matrix<scalar,Eigen::Dynamic,2,Eigen::RowMajor> MatrixX2sr;
typedef Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatrixXXsr;
typedef Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> MatrixXXsc;

// Define scalar-valued versions of the Eigen sparse matrix types
typedef Eigen::SparseMatrix<scalar,Eigen::ColMajor> SparseMatrixsc;
typedef Eigen::SparseMatrix<scalar,Eigen::RowMajor> SparseMatrixsr;

// Define some transformation geometry
typedef Eigen::AngleAxis<scalar> AnglesAxis3s;

// Permutation matrix
typedef Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,int> PermutationMatrix;

#endif
