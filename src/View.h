#pragma once

#include "QuSim.h"
#include "eigen/Eigen/Dense"
#include "Linear.h"

template<class T>
VectorView<T> View(std::vector<T> const &v)
{
    return VectorView<T>(v.data(), v.size());
}

template<class T>
VectorView<T> View(T *data, size_t n)
{
    return VectorView<T>(data, n);
}

template<class T>
VectorView<T> View(T const * data, size_t n)
{
    return VectorView<T>(data, n);
}

template<class T>
inline std::vector<T> ToVector(VectorView<T> const &v)
{
    return std::vector<T>(v.begin(), v.end());
}


inline MatrixView<Complex> View(Eigen::MatrixXcd const &v)
{
    return MatrixView<Complex>(v.data(), v.rows(), v.cols());
}

inline MatrixView<Real> View(Eigen::MatrixXd const &v)
{
    return MatrixView<Real>(v.data(), v.rows(), v.cols());
}

inline MatrixView<Real> View(Real const * d, size_t sz1, size_t sz2)
{
    return MatrixView<Real>(d, sz1, sz2);
}

inline MatrixView<Complex> View(Complex const* d, size_t sz1, size_t sz2)
{
    return MatrixView<Complex>(d, sz1, sz2);
}

inline Tensor3View<Complex> View(Complex const * d, size_t sz1, size_t sz2, size_t sz3)
{
    return Tensor3View<Complex>(d, sz1, sz2, sz3);
}

inline Tensor3View<Real> View(Real const * d, size_t sz1, size_t sz2, size_t sz3)
{
    return Tensor3View<Real>(d, sz1, sz2, sz3);
}

inline Eigen::MatrixXd ToEigen(MatrixView<Real> const &v)
{
    using Map = Eigen::Map< Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > >;
    using Index = Map::Index;
    return Map((Real*)v.data(), (Index)v.rows(), (Index)v.cols());
}

inline Eigen::MatrixXcd ToEigen(MatrixView<Complex> const& v)
{
    using Map = Eigen::Map< Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > >;
    using Index = Map::Index;
    return Map((Complex*)v.data(), (Index)v.rows(), (Index)v.cols());
}

