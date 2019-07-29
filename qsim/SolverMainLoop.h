

void SolverImpl1D:: LOOP_FUNC_NAME ()
{
	Complex psi = fPsi[0];
	Complex psiPrime = fPsiPrime[0];

	Real e = -2 * fMass / (fHbar*fHbar) * fE;
	Real vk = -2 * fMass / (fHbar*fHbar);

#ifdef SMALL_ROUND_ERROR
	Matrix big;
	Matrix little;
	big(0, 0) = fTMat(0, 0);
	big(0, 1) = fTMat(0, 1);
	big(1, 0) = fTMat(1, 0);
	big(1, 1) = fTMat(1, 1);

	little(0, 0) = 0;
	little(0, 1) = 0;
	little(1, 0) = 0;
	little(1, 1) = 0;
#else
	Matrix mat;
	mat(0, 0) = fTMat(0, 0);
	mat(0, 1) = fTMat(0, 1);
	mat(1, 0) = fTMat(1, 0);
	mat(1, 1) = fTMat(1, 1);
#endif

#ifdef Tail
#undef Tail
#endif

#if !defined(SMALL_ROUND_ERROR)

#define Tail()											   \
	mat = tr * mat;                                        \
	Complex psi_ = psi;									   \
	Complex psiPrime_ = psiPrime;						   \
	psi = psi_ * tr.Get11() + psiPrime_ * tr.Get12();	   \
	psiPrime = psi_ * tr.Get21() + psiPrime_ * tr.Get22(); \
														   \
	fPsi[i + 1] = psi;									   \
	fPsiPrime[i + 1] = psiPrime;

#else

#define Tail()		                       \
	little = tr * little + (tr_mI*big - ((tr*big) - big)); \
	big = tr * big;                                        \
	Complex psi_ = psi;									   \
	Complex psiPrime_ = psiPrime;						   \
	psi = psi_ * tr.Get11() + psiPrime_ * tr.Get12();	   \
	psiPrime = psi_ * tr.Get21() + psiPrime_ * tr.Get22(); \
														   \
	fPsi[i + 1] = psi;									   \
	fPsiPrime[i + 1] = psiPrime;

#endif


	//                                    | 0        1 |
	//       d{psi, psiPrime}^T /dx =     |            |  {psi, psiPrime}^T
	//                                    | e-vk*v   0 |
	if (fMethod == SolverMethod::ImplicitMidpointMethod) {
		for (size_t i = 0; i < fNBins; ++i) {
			// K = A (1 + 1/2 K)
			// tr = 1 + 1.0 K

			Real a = e - vk * fV[i];
			AntiDiagonalMatrix2<Real> A(fDx, a*fDx);
			Matrix tr_mI = A * (Matrix::Identity() - 0.5*A).inverse();
			Matrix tr = Matrix::Identity() + tr_mI;

			Tail();
		}

	} else if (fMethod == SolverMethod::ExplicitRungeKuttaO4Classical) {
		for (size_t i = 0; i < fNBins; ++i) {

			Real a1 = e - vk * fV[3 * i];
			Real a2 = e - vk * fV[3 * i + 1];
			//Real a3 = a2;
			Real a4 = e - vk * fV[3 * i + 2];

			AntiDiagonalMatrix2<Real> const A1(fDx, a1*fDx);
			AntiDiagonalMatrix2<Real> const A2(fDx, a2*fDx);
			AntiDiagonalMatrix2<Real> const &A3 = A2;
			AntiDiagonalMatrix2<Real> const A4(fDx, a4*fDx);

			AntiDiagonalMatrix2<Real> const &K1 = A1;
			Matrix const K2 = A2 * (Matrix::Identity() + 1. / 2 * K1);
			Matrix const K3 = A3 * (Matrix::Identity() + 1. / 2 * K2);
			Matrix const K4 = A4 * (Matrix::Identity() + K3);

			Matrix tr_mI = 1. / 6 * (K1 + 2. * K2 + 2. * K3 + K4);
			Matrix tr = Matrix::Identity() + tr_mI;

			Tail();
		}

	} else if (fMethod == SolverMethod::GaussLegendreO4) {

		// K1 = A1(1 + a11*K1 + a12*K2)
		// K2 = A2(1 + a21*K1 + a22*K2)
		// K1 = (1 - A2 a22 - A1 a11 + (a22 a11-a12 a21)A2A1)^-1 (A1 - A2A1 a22 + A2A1 a12)
		// K2 = (1 - A2 a22 - A1 a11 + (a22 a11-a12 a21)A1A2)^-1 (A2 - A1A2 a11 + A1A2 a21)
		Real const a12 = 1. / 4 - sqrt(3) / 6;
		Real const a21 = 1. / 4 + sqrt(3) / 6;
		Real const a11 = 1. / 4;
		Real const a22 = 1. / 4;
		Real const C2 = a22*a11 - a12*a21;

		for (size_t i = 0; i < fNBins; ++i) {
			Real a1 = e - vk * fV[2 * i];
			Real a2 = e - vk * fV[2 * i + 1];

			AntiDiagonalMatrix2<Real> A1(fDx, fDx*a1);
			AntiDiagonalMatrix2<Real> A2(fDx, fDx*a2);
			double dx2a1 = fDx*fDx*a1;
			double dx2a2 = fDx*fDx*a2;
			DiagonalMatrix2<Real> const A1A2(dx2a1, dx2a2);
			DiagonalMatrix2<Real> const A2A1(dx2a2, dx2a1);

			AntiDiagonalMatrix2<Real> T = a11*A1 + a22*A2;
			Matrix K1 = A1 * (Matrix::Identity() + C2 * A2A1 - T).inverse() * (Matrix::Identity() + (a12 - a22)*A2);
			Matrix K2 = A2 * (Matrix::Identity() + C2 * A1A2 - T).inverse() * (Matrix::Identity() + (a21 - a11)*A1);

			Matrix tr_mI = 0.5*(K1 + K2);
			Matrix tr = Matrix::Identity() + tr_mI;

			Tail();
		}
	} else if (fMethod == SolverMethod::ExplicitRungeKuttaO6Luther1967) {

		Real const c31 = 3. / 8;
		Real const c32 = 1. / 8;

		Real const c41 = 8. / 27.;
		Real const c42 = 2. / 27;
		Real const c43 = 8. / 27;

		Real const c51 = 3 * (3 * sqrt(21) - 7) / 392;
		Real const c52 = -8 * (7 - sqrt(21)) / 392;
		Real const c53 = 48 * (7 - sqrt(21)) / 392;
		Real const c54 = -3 * (21 - sqrt(21)) / 392;
		//Real const one1 = c51 + c52 + c53 + c54;

		Real const c61 = -5 * (231 + 51 * sqrt(21)) / 1960;
		Real const c62 = -40 * (7 + sqrt(21)) / 1960;
		Real const c63 = -320 * sqrt(21) / 1960;
		Real const c64 = 3 * (21 + 121 * sqrt(21)) / 1960;
		Real const c65 = 392 * (6 + sqrt(21)) / 1960;
		//Real const one2 = c61 + c62 + c63 + c64 + c65;

		Real const c71 = 15 * (22 + 7 * sqrt(21)) / 180;
		Real const c72 = 120. / 180;
		Real const c73 = 40 * (7 * sqrt(21) - 5) / 180;
		Real const c74 = -63 * (3 * sqrt(21) - 2) / 180;
		Real const c75 = -14 * (49 + 9 * sqrt(21)) / 180;
		Real const c76 = 70 * (7 - sqrt(21)) / 180;
		//Real const one3 = c71 + c72 + c73 + c74 + c75 + c76;

		Real const f1 = 0;
		Real const f2 = 1;
		Real const f3 = 1. / 2;
		Real const f4 = 2. / 3;
		Real const f5 = (7 - sqrt(21)) / 14;
		Real const f6 = (7 + sqrt(21)) / 14;
		Real const f7 = 1;

		Real const b1 = 9. / 180;
		Real const b3 = 64. / 180;
		Real const b5 = 49. / 180;
		Real const b6 = 49. / 180;
		Real const b7 = 9. / 180;;


		for (size_t i = 0; i < fNBins; ++i) {

			Real a1 = e - vk * fV[7 * i];
			Real a2 = e - vk * fV[7 * i + 1];
			Real a3 = e - vk * fV[7 * i + 2];
			Real a4 = e - vk * fV[7 * i + 3];
			Real a5 = e - vk * fV[7 * i + 4];
			Real a6 = e - vk * fV[7 * i + 5];
			Real a7 = e - vk * fV[7 * i + 6];

			AntiDiagonalMatrix2<Real> A1(fDx, a1*fDx);
			AntiDiagonalMatrix2<Real> A2(fDx, a2*fDx);
			AntiDiagonalMatrix2<Real> A3(fDx, a3*fDx);
			AntiDiagonalMatrix2<Real> A4(fDx, a4*fDx);
			AntiDiagonalMatrix2<Real> A5(fDx, a5*fDx);
			AntiDiagonalMatrix2<Real> A6(fDx, a6*fDx);
			AntiDiagonalMatrix2<Real> A7(fDx, a7*fDx);

			AntiDiagonalMatrix2<Real> const &K1 = A1;
			Matrix const K2 = A2 * (Matrix::Identity() + K1);
			Matrix const K3 = A3 * (Matrix::Identity() + c31 * K1 + c32 * K2);
			Matrix const K4 = A4 * (Matrix::Identity() + c41 * K1 + c42 * K2 + c43 * K3);
			Matrix const K5 = A5 * (Matrix::Identity() + c51 * K1 + c52 * K2 + c53 * K3 + c54 * K4);
			Matrix const K6 = A6 * (Matrix::Identity() + c61 * K1 + c62 * K2 + c63 * K3 + c64 * K4 + c65 * K5);
			Matrix const K7 = A7 * (Matrix::Identity() + c71 * K1 + c72 * K2 + c73 * K3 + c74 * K4 + c75 * K5 + c76 * K6);

			Matrix tr_mI = b1 * K1 + b3 * K3 + b5 * K5 + b6 * K6 + b7 * K7;
			Matrix tr = Matrix::Identity() + tr_mI;

			Tail();
		}

	}

#ifdef SMALL_ROUND_ERROR
	Matrix mat = big + little;
#endif

	fTMat(0, 0) = mat(0, 0);
	fTMat(0, 1) = mat(0, 1);
	fTMat(1, 0) = mat(1, 0);
	fTMat(1, 1) = mat(1, 1);

}
