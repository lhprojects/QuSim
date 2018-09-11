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

inline Eigen::MatrixXd ToEigen(MatrixView<Real> const &v)
{
	return Eigen::Map< Eigen::MatrixXd >((double*)v.data(), (Eigen::Map< Eigen::MatrixXd >::Index)v.rows(),
		(Eigen::Map< Eigen::MatrixXd >::Index)v.cols());
}

inline Eigen::MatrixXcd ToEigen(MatrixView<Complex> const &v)
{
	return Eigen::Map<Eigen::MatrixXcd>((Complex*)v.data(), v.rows(), v.cols());
}

