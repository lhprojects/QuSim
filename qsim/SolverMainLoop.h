

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
#define Tail()											   \
	 mat = tr * mat;                                       \
	Complex psi_ = psi;									   \
	Complex psiPrime_ = psiPrime;						   \
	psi = psi_ * tr(0, 0) + psiPrime_ * tr(0, 1);		   \
	psiPrime = psi_ * tr(1, 0) + psiPrime_ * tr(1, 1);	   \
														   \
	fPsi[i + 1] = psi;									   \
	fPsiPrime[i + 1] = psiPrime;

#define Tail_SMALLROUNDERR()		                       \
	little = tr * little + (tr_mI*big - ((tr*big) - big)); \
	big = tr * big;                                        \
	Complex psi_ = psi;									   \
	Complex psiPrime_ = psiPrime;						   \
	psi = psi_ * tr(0, 0) + psiPrime_ * tr(0, 1);		   \
	psiPrime = psi_ * tr(1, 0) + psiPrime_ * tr(1, 1);	   \
														   \
	fPsi[i + 1] = psi;									   \
	fPsiPrime[i + 1] = psiPrime;



	if (fMethod == SolverMethod::ImplicitMidpointMethod) {
		for (size_t i = 0; i < fNBins; ++i) {
			Matrix tr;
			Real a = e - vk * fV[i];

			AntiDiagonalMatrix2<Real> A;
			A(0, 1) = 1; A(1, 0) = a;

			// K = dx A (1 + 1/2 K)
			AntiDiagonalMatrix2<Real> dxA = fDx * A;

#ifdef SMALL_ROUND_ERROR
				Matrix tr_mI = dxA * (Matrix::Identity() - 0.5*dxA).inverse();
				tr = Matrix::Identity() + tr_mI;
				Tail_SMALLROUNDERR();
#else
				tr = Matrix::Identity() + dxA * (Matrix::Identity() - 0.5*dxA).inverse();
				Tail();
#endif
		}

	} else if (fMethod == SolverMethod::ExplicitRungeKuttaO4Classical) {
		for (size_t i = 0; i < fNBins; ++i) {
			Matrix tr;
			Real a1 = e - vk * fV[3 * i];
			Real a2 = e - vk * fV[3 * i + 1];
			//Real a3 = a2;
			Real a4 = e - vk * fV[3 * i + 2];

			AntiDiagonalMatrix2<Real> A1;
			A1(0, 1) = 1; A1(1, 0) = a1;
			AntiDiagonalMatrix2<Real> A2;
			A2(0, 1) = 1; A2(1, 0) = a2;

			AntiDiagonalMatrix2<Real> const &A3 = A2;
			AntiDiagonalMatrix2<Real> A4;
			A4(0, 1) = 1; A4(1, 0) = a4;

			AntiDiagonalMatrix2<Real> const &K1 = A1;
			Matrix K2 = A2 * (Matrix::Identity() + 1. / 2 * fDx*K1);
			Matrix K3 = A3 * (Matrix::Identity() + 1. / 2 * fDx*K2);
			Matrix K4 = A4 * (Matrix::Identity() + fDx * K3);

#ifdef SMALL_ROUND_ERROR
				Matrix tr_mI = 1. / 6 * fDx * (
					K1 + 2. * K2 + 2. * K3 + K4
					);
				tr = Matrix::Identity() + tr_mI;
				Tail_SMALLROUNDERR();
#else
				tr = Matrix::Identity() + 1. / 6 * fDx * (
					K1 + 2. * K2 + 2. * K3 + K4
					);
				Tail();
#endif
		}

	} else if (fMethod == SolverMethod::GaussLegendreO4) {

		Real const a12 = 1. / 4 - sqrt(3) / 6;
		Real const a21 = 1. / 4 + sqrt(3) / 6;
		Real const BB = -a12 * fDx;
		Real const CC = -a21 * fDx;
		Real const FF = 1. / 16 * fDx*fDx - CC * BB;

		for (size_t i = 0; i < fNBins; ++i) {
			Matrix tr;
			Real a1 = e - vk * fV[2 * i];
			Real a2 = e - vk * fV[2 * i + 1];


			AntiDiagonalMatrix2<Real> A1;
			A1(0, 1) = 1; A1(1, 0) = a1;
			AntiDiagonalMatrix2<Real> A2;
			A2(0, 1) = 1; A2(1, 0) = a2;

			AntiDiagonalMatrix2<Real> tmp = 0.25*fDx*(A1 + A2);
			Matrix kk1 = A1 * (Matrix::Identity() - tmp + FF * A2 * A1).inverse() * (Matrix::Identity() - (1. / 4 * fDx + BB)*A2);
			Matrix kk2 = A2 * (Matrix::Identity() - tmp + FF * A1 * A2).inverse() * (Matrix::Identity() - (1. / 4 * fDx + CC)*A1);

#ifdef SMALL_ROUND_ERROR
			Matrix tr_mI = fDx * 0.5*(kk1 + kk2);
			tr = Matrix::Identity() + tr_mI;
			Tail_SMALLROUNDERR();
#else
			Matrix tr_mI = fDx * 0.5*(kk1 + kk2);
				tr = Matrix::Identity() + tr_mI;
				//tr = Matrix::Identity() + fDx * 0.5*(kk1 + kk2);
				Tail();
#endif
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
			Matrix tr;
			Real a1 = e - vk * fV[7 * i];
			Real a2 = e - vk * fV[7 * i + 1];
			Real a3 = e - vk * fV[7 * i + 2];
			Real a4 = e - vk * fV[7 * i + 3];
			Real a5 = e - vk * fV[7 * i + 4];
			Real a6 = e - vk * fV[7 * i + 5];
			Real a7 = e - vk * fV[7 * i + 6];

			AntiDiagonalMatrix2<Real> A1;
			AntiDiagonalMatrix2<Real> A2;
			AntiDiagonalMatrix2<Real> A3;
			AntiDiagonalMatrix2<Real> A4;
			AntiDiagonalMatrix2<Real> A5;
			AntiDiagonalMatrix2<Real> A6;
			AntiDiagonalMatrix2<Real> A7;
			A1(0, 1) = 1; A1(1, 0) = a1;
			A2(0, 1) = 1; A2(1, 0) = a2;
			A3(0, 1) = 1; A3(1, 0) = a3;
			A4(0, 1) = 1; A4(1, 0) = a4;
			A5(0, 1) = 1; A5(1, 0) = a5;
			A6(0, 1) = 1; A6(1, 0) = a6;
			A7(0, 1) = 1; A7(1, 0) = a7;


			AntiDiagonalMatrix2<Real> K1 = fDx * A1;
			Matrix K2 = fDx * A2 * (Matrix::Identity() + K1);
			Matrix K3 = fDx * A3 * (Matrix::Identity() + c31 * K1 + c32 * K2);
			Matrix K4 = fDx * A4 * (Matrix::Identity() + c41 * K1 + c42 * K2 + c43 * K3);
			Matrix K5 = fDx * A5 * (Matrix::Identity() + c51 * K1 + c52 * K2 + c53 * K3 + c54 * K4);
			Matrix K6 = fDx * A6 * (Matrix::Identity() + c61 * K1 + c62 * K2 + c63 * K3 + c64 * K4 + c65 * K5);
			Matrix K7 = fDx * A7 * (Matrix::Identity() + c71 * K1 + c72 * K2 + c73 * K3 + c74 * K4 + c75 * K5 + c76 * K6);

#ifdef SMALL_ROUND_ERROR
			Matrix tr_mI = b1 * K1 + b3 * K3 + b5 * K5 + b6 * K6 + b7 * K7;
			tr = Matrix::Identity() + tr_mI;
			Tail_SMALLROUNDERR();
#else
			tr = Matrix::Identity() + b1 * K1 + b3 * K3 + b5 * K5 + b6 * K6 + b7 * K7;
			Tail();
#endif
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
