#ifndef MATRIX2_H
#define MATRIX2_H

#include <stdint.h>

template<class Scalar>
struct ColVec2 {

	ColVec2(double a, double b) : x(a), y(b) {}
	Scalar &operator()(size_t i)
	{
		assert(i < 2);
		if (i == 0) return x;
		else return y;
	}
	Scalar const &operator()(size_t i) const
	{
		assert(i < 2);
		if (i == 0) return x;
		else return y;
	}
private:
	Scalar x;
	Scalar y;
};

template<class Scalar>
struct Mat2Identity {
};

template<class Scalar>
struct Mat2IdentityK {
	Mat2IdentityK(double k_) : k(k_) {}
	Scalar k;
};


template<class Scalar>
struct Mat2 {

	size_t rows() const { return 2; }
	size_t cols() const { return 2; }
	Mat2() {}
	Mat2(Scalar a11, Scalar a12, Scalar a21, Scalar a22) : m11(a11), m12(a12), m21(a21), m22(a22)
	{
	}

	ColVec2<Scalar> operator*(ColVec2<Scalar> const &cv) const
	{
		return ColVec2<Scalar>(m11*cv.x + m12 * cv.y, m21*cv.x + m22 * cv.y);
	}

	Scalar &operator()(size_t i, size_t j)
	{
		assert(i < 2 && j < 2);
		if (i == 0 && j == 0) return m11;
		if (i == 0 && j == 1) return m12;
		if (i == 1 && j == 0) return m21;
		return m22;
	}

	Scalar const &operator()(size_t i, size_t j) const
	{
		assert(i < 2 && j < 2);
		if (i == 0 && j == 0) return m11;
		if (i == 0 && j == 1) return m12;
		if (i == 1 && j == 0) return m21;
		return m22;
	}

	inline Mat2<Scalar> inverse() const
	{
		Scalar det = Scalar(1) / (m11*m22 - m12 * m21);
		return Mat2<Scalar>(m22*det, -m12 * det, -m21 * det, m11*det);
	}

	static Mat2Identity<Scalar> Identity()
	{
		return Mat2Identity<Scalar>();
	}

	Scalar m11;
	Scalar m12;
	Scalar m21;
	Scalar m22;
};

template<class Scalar>
struct AntiDiagonalMatrix2
{
	AntiDiagonalMatrix2()
	{
	}

	AntiDiagonalMatrix2(Scalar v12, Scalar v21)
	{
		m12 = v12;
		m21 = v21;
	}

	Scalar &operator()(size_t i, size_t j)
	{
		assert((i == 1 && j == 0) || (i == 0 && j == 1));
		if (i == 0 && j == 1) return m12;
		else return m21;
		
	}

	Scalar const &operator()(size_t i, size_t j) const
	{
		assert((i == 1 && j == 0) || (i == 0 && j == 1));
		if (i == 0 && j == 1) return m12;
		else return m21;
	}

	inline AntiDiagonalMatrix2<Scalar> inverse() const
	{
		Scalar inv_m_det = Scalar(1)/(m12 * m21);
		return AntiDiagonalMatrix2<Scalar>(m12*inv_m_det, m21*inv_m_det);
	}
	Scalar m12;
	Scalar m21;
};

template<class Scalar>
struct DiagonalMatrix2 {
	DiagonalMatrix2()
	{
	}

	DiagonalMatrix2(Scalar v11, Scalar v22)
	{
		m11 = v11;
		m22 = v22;
	}

	Scalar &operator()(size_t i, size_t j)
	{
		assert((i == 0 && j == 0) || (i == 1 && j == 1));
		if (i == 0 && j == 0) return m11;
		else return m22;

	}

	Scalar const &operator()(size_t i, size_t j) const
	{
		assert((i == 0 && j == 0) || (i == 1 && j == 1));
		if (i == 0 && j == 0) return m11;
		else return m22;
	}

	inline DiagonalMatrix2<Scalar> inverse() const
	{
		Scalar inv_det = Scalar(1) /(m11 * m22);
		return AntiDiagonalMatrix2<Scalar>(m22*inv_det, m11*inv_det);
	}

	Scalar m11;
	Scalar m22;
};

template<class Scalar>
inline AntiDiagonalMatrix2<Scalar> operator+(AntiDiagonalMatrix2<Scalar> const &l,
	AntiDiagonalMatrix2<Scalar> const &r)
{
	return AntiDiagonalMatrix2<Scalar>(l.m12 + r.m12, l.m21 + r.m21);
}

template<class Scalar>
inline Mat2<Scalar> operator+(Mat2Identity<Scalar> const &l,
	AntiDiagonalMatrix2<Scalar> const &r)
{
	return Mat2<Scalar>(Scalar(1), r.m12, r.m21, Scalar(1));
}

template<class Scalar>
inline Mat2<Scalar> operator-(Mat2Identity<Scalar> const &l,
	AntiDiagonalMatrix2<Scalar> const &r)
{
	return Mat2<Scalar>(Scalar(1), -r.m12, -r.m21, Scalar(1));
}

template<class Scalar>
inline Mat2<Scalar> operator+(AntiDiagonalMatrix2<Scalar> const &l,
	Mat2<Scalar> const &r)
{
	return Mat2<Scalar>(r.m11, l.m12 + r.m12, l.m21 + r.m21, r.m22);
}

template<class Scalar>
inline Mat2<Scalar> operator+(Mat2<Scalar> const &l,
	AntiDiagonalMatrix2<Scalar> const &r)
{
	return Mat2<Scalar>(l.m11, l.m12 + r.m12, l.m21 + r.m21, l.m22);
}

template<class Scalar>
inline Mat2<Scalar> operator-(Mat2<Scalar> const &l,
	AntiDiagonalMatrix2<Scalar> const &r)
{
	return Mat2<Scalar>(l.m11, l.m12 - r.m12, l.m21 - r.m21, l.m22);
}

template<class Scalar>
inline Mat2<Scalar> operator+(Mat2<Scalar> const &l,
	DiagonalMatrix2<Scalar> const &r)
{
	return Mat2<Scalar>(l.m11 + r.m11, l.m12, l.m21, l.m22 + r.m22);
}

template<class Scalar>
inline AntiDiagonalMatrix2<Scalar> operator*(Scalar l,
	AntiDiagonalMatrix2<Scalar> const &r)
{
	return AntiDiagonalMatrix2<Scalar>(l * r.m12, l * r.m21);
}

template<class Scalar>
inline DiagonalMatrix2<Scalar> operator*(AntiDiagonalMatrix2<Scalar> const &l,
	AntiDiagonalMatrix2<Scalar> const &r)
{
	return DiagonalMatrix2<Scalar>(l.m12 * r.m21, l.m21 * r.m12);
}

template<class Scalar>
inline Mat2<Scalar> operator*(AntiDiagonalMatrix2<Scalar> const &l,
	Mat2<Scalar> const &r)
{
	return Mat2<Scalar>(l.m12 * r.m21,
		l.m12 * r.m22,
		l.m21 * r.m11,
		l.m21 * r.m12);
}


template<class Scalar>
inline Mat2<Scalar> operator*(Mat2<Scalar> const &l, Mat2<Scalar> const &r)
{
	return Mat2<Scalar>(l.m11*r.m11 + l.m12 * r.m21,
		l.m11*r.m12 + l.m12 * r.m22,
		l.m21*r.m11 + l.m22 * r.m21,
		l.m21*r.m12 + l.m22 * r.m22);
}

template<class Scalar>
Mat2<Scalar> operator*(Mat2<Scalar> const &m, Mat2IdentityK<Scalar> k)
{
	return Mat2<Scalar>(m.m11 *k.k, m.m12 *k.k, m.m21 *k.k, m.m22 *k.k);
}

template<class Scalar>
Mat2<Scalar> operator*(Mat2IdentityK<Scalar> k, Mat2<Scalar> const &m)
{
	return Mat2<Scalar>(k.k*m.m11, k.k*m.m12, k.k*m.m21, k.k*m.m22);
}

template<class Scalar>
Mat2<Scalar> operator*(Scalar k, Mat2<Scalar> const &m)
{
	return Mat2<Scalar>(k*m.m11, k*m.m12, k*m.m21, k*m.m22);
}

template<class Scalar>
Mat2<Scalar> operator*(Mat2<Scalar> const &m, Scalar k)
{
	return Mat2<Scalar>(m.m11 *k, m.m12 *k, m.m21 *k, m.m22 *k);
}

template<class Scalar>
Mat2<Scalar> operator*(Mat2<Scalar> const &m, Mat2Identity<Scalar> &)
{
	return m;
}

template<class Scalar>
Mat2<Scalar> operator*(Mat2Identity<Scalar> &, Mat2<Scalar> const &m)
{
	return m;
}

template<class Scalar>
Mat2IdentityK<Scalar> operator*(Mat2IdentityK<Scalar> const &k1, Mat2IdentityK<Scalar> const &k2)
{
	return Mat2IdentityK<Scalar>(k1.k * k2.k);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator*(Mat2IdentityK<Scalar> const &k, Scalar k2)
{
	return Mat2IdentityK<Scalar>(k.k*k2);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator*(Scalar k2, Mat2IdentityK<Scalar> const &k)
{
	return Mat2IdentityK<Scalar>(k2 * k.k);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator*(Scalar k, Mat2Identity<Scalar> const &)
{
	return Mat2IdentityK<Scalar>(k);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator*(Mat2Identity<Scalar> const &, Scalar k)
{
	return Mat2IdentityK<Scalar>(k);
}







template<class Scalar>
Mat2<Scalar> operator+(Mat2<Scalar> const &m, Mat2<Scalar> const &r)
{
	return Mat2<Scalar>(m.m11 + r.m11, m.m12 + r.m12, m.m21 + r.m21, m.m22 + r.m22);
}

template<class Scalar>
Mat2<Scalar> operator+(Mat2<Scalar> const &m, Mat2IdentityK<Scalar> const &k)
{
	return Mat2<Scalar>(m.m11 + Scalar(k.k), m.m12, m.m21, m.m22 + Scalar(k.k));
}

template<class Scalar>
Mat2<Scalar> operator+(Mat2IdentityK<Scalar> const &k, Mat2<Scalar> const &m)
{
	return Mat2<Scalar>(Scalar(k.k) + m.m11, m.m12, m.m21, Scalar(k.k) + m.m22);
}

template<class Scalar>
Mat2<Scalar> operator+(Mat2<Scalar> const &m, Mat2Identity<Scalar> const &)
{
	return Mat2<Scalar>(m.m11 + Scalar(1), m.m12, m.m21, m.m22 + Scalar(1));
}

template<class Scalar>
Mat2<Scalar> operator+(Mat2Identity<Scalar> const &, Mat2<Scalar> const &m)
{
	return Mat2<Scalar>(Scalar(1) + m.m11, m.m12, m.m21, Scalar(1) + m.m22);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator+(Mat2IdentityK<Scalar> const &k1, Mat2IdentityK<Scalar> const &k2)
{
	return Mat2Identity<Scalar>(k1 + k2);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator+(Mat2Identity<Scalar> const &, Mat2IdentityK<Scalar> const &k2)
{
	return Mat2Identity<Scalar>(Scalar(1) + k2);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator+(Mat2IdentityK<Scalar> const &k1, Mat2Identity<Scalar> const &k2)
{
	return Mat2Identity<Scalar>(k2 + Scalar(1));
}






template<class Scalar>
Mat2<Scalar> operator-(Mat2<Scalar> const &m, Mat2<Scalar> const &r)
{
	return Mat2<Scalar>(m.m11 - r.m11, m.m12 - r.m12, m.m21 - r.m21, m.m22 - r.m22);
}

template<class Scalar>
Mat2<Scalar> operator-(Mat2<Scalar> const &m, Mat2IdentityK<Scalar> const &k)
{
	return Mat2<Scalar>(m.m11 - Scalar(k.k), m.m12, m.m21, m.m22 - Scalar(k.k));
}

template<class Scalar>
Mat2<Scalar> operator-(Mat2IdentityK<Scalar> const &k, Mat2<Scalar> const &m)
{
	return Mat2<Scalar>(Scalar(k.k) - m.m11, -m.m12, -m.m21, Scalar(k.k) - m.m22);
}

template<class Scalar>
Mat2<Scalar> operator-(Mat2<Scalar> const &m, Mat2Identity<Scalar> const &)
{
	return Mat2<Scalar>(m.m11 - Scalar(1), m.m12, m.m21, m.m22 - Scalar(1));
}

template<class Scalar>
Mat2<Scalar> operator-(Mat2Identity<Scalar> const &, Mat2<Scalar> const &m)
{
	return Mat2<Scalar>(Scalar(1) - m.m11, -m.m12, -m.m21, Scalar(1) - m.m22);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator-(Mat2IdentityK<Scalar> const &k1, Mat2IdentityK<Scalar> const &k2)
{
	return Mat2Identity<Scalar>(k1 - k2);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator-(Mat2Identity<Scalar> const &, Mat2IdentityK<Scalar> const &k2)
{
	return Mat2Identity<Scalar>(Scalar(1) - k2);
}

template<class Scalar>
Mat2IdentityK<Scalar> operator-(Mat2IdentityK<Scalar> const &k1, Mat2Identity<Scalar> const &k2)
{
	return Mat2Identity<Scalar>(k2 - Scalar(1));
}


#endif
