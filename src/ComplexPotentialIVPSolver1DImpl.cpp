#include "ComplexPotentialIVPSolver1DImpl.h"
#include "Matrix2.h"
#define SMALL_ROUND_ERROR

typedef Mat2<Complex> Matrix;
void ComplexPotentialIVPSolver1DImpl::Compute()
{

    auto met = fMethod;
    auto x0 = fX0;

    if (met == SolverMethod::ImplicitMidpointMethod) {
        fV.resize(fNBins);
        for (size_t i = 0; i < fNBins; ++i) {
            Real x = x0 + (i + 0.5)*fDx;
            fV[i] = fVFunc(x);
        }
    } else if (met == SolverMethod::ExplicitRungeKuttaO4Classical) {
        fV.resize(3 * fNBins);
        for (size_t i = 0; i < fNBins; ++i) {
            Real x1 = x0 + (i + 0.0)*fDx;
            Real x2 = x0 + (i + 0.5)*fDx;
            Real x3 = x0 + (i + 1.0)*fDx;
            fV[3 * i] = fVFunc(x1);
            fV[3 * i + 1] = fVFunc(x2);
            fV[3 * i + 2] = fVFunc(x3);
        }
    } else if (met == SolverMethod::GaussLegendreO4) {
        fV.resize(2 * fNBins);
        Real f1 = 1. / 2 - sqrt(3) / 6;
        Real f2 = 1. / 2 + sqrt(3) / 6;
        for (size_t i = 0; i < fNBins; ++i) {
            Real x1 = x0 + (i + f1)*fDx;
            Real x2 = x0 + (i + f2)*fDx;
            fV[2 * i] = fVFunc(x1);
            fV[2 * i + 1] = fVFunc(x2);
        }
    } else if (met == SolverMethod::ExplicitRungeKuttaO6Luther1967) {
        fV.resize(7 * fNBins);
        Real f1 = 0;
        Real f2 = 1;
        Real f3 = 1. / 2;
        Real f4 = 2. / 3;
        Real f5 = (7 - sqrt(21)) / 14;
        Real f6 = (7 + sqrt(21)) / 14;
        Real f7 = 1;

        for (size_t i = 0; i < fNBins; ++i) {
            Real x1 = x0 + (i + f1)*fDx;
            Real x2 = x0 + (i + f2)*fDx;
            Real x3 = x0 + (i + f3)*fDx;
            Real x4 = x0 + (i + f4)*fDx;
            Real x5 = x0 + (i + f5)*fDx;
            Real x6 = x0 + (i + f6)*fDx;
            Real x7 = x0 + (i + f7)*fDx;
            fV[7 * i] = fVFunc(x1);
            fV[7 * i + 1] = fVFunc(x2);
            fV[7 * i + 2] = fVFunc(x3);
            fV[7 * i + 3] = fVFunc(x4);
            fV[7 * i + 4] = fVFunc(x5);
            fV[7 * i + 5] = fVFunc(x6);
            fV[7 * i + 6] = fVFunc(x7);
        }
    } else {
        throw std::runtime_error("Unsupported method");
    }


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
        //Real const OneFourthDxDx = 0.25 * fDx*fDx;
        for (size_t i = 0; i < fNBins; ++i) {
            Matrix tr;
            Complex a = e - vk * fV[i];

            AntiDiagonalMatrix2<Complex> A;
            A(0, 1) = 1; A(1, 0) = a;

            // K = dx A (1 + 1/2 K)
            AntiDiagonalMatrix2<Complex> dxA = Complex(fDx) * A;

#ifdef SMALL_ROUND_ERROR
            Matrix tr_mI = dxA * (Matrix::Identity() - Complex(0.5)*dxA).inverse();
            tr = Matrix::Identity() + tr_mI;
            Tail_SMALLROUNDERR();
#else
            tr = Matrix::Identity() + dxA * (Matrix::Identity() - Complex(0.5)*dxA).inverse();
            Tail();
#endif
        }

    } else if (fMethod == SolverMethod::ExplicitRungeKuttaO4Classical) {
        for (size_t i = 0; i < fNBins; ++i) {
            Matrix tr;
            Complex a1 = e - vk * fV[3 * i];
            Complex a2 = e - vk * fV[3 * i + 1];
            //Real a3 = a2;
            Complex a4 = e - vk * fV[3 * i + 2];

            AntiDiagonalMatrix2<Complex> A1;
            A1(0, 1) = 1; A1(1, 0) = a1;
            AntiDiagonalMatrix2<Complex> A2;
            A2(0, 1) = 1; A2(1, 0) = a2;

            AntiDiagonalMatrix2<Complex> const &A3 = A2;
            AntiDiagonalMatrix2<Complex> A4;
            A4(0, 1) = 1; A4(1, 0) = a4;

            AntiDiagonalMatrix2<Complex> const &K1 = A1;
            Matrix K2 = A2 * (Matrix::Identity() + Complex(1. / 2 * fDx)*K1);
            Matrix K3 = A3 * (Matrix::Identity() + Complex(1. / 2 * fDx)*K2);
            Matrix K4 = A4 * (Matrix::Identity() + Complex(fDx) * K3);

#ifdef SMALL_ROUND_ERROR
            Matrix tr_mI = Complex(1. / 6) * fDx * (
                K1 + Complex(2.) * K2 + Complex(2.) * K3 + K4
                );
            tr = Matrix::Identity() + tr_mI;
            Tail_SMALLROUNDERR();
#else
            tr = Matrix::Identity() + Complex(1. / 6) * fDx * (
                K1 + Complex(2.) * K2 + Complex(2.) * K3 + K4
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
            Complex a1 = e - vk * fV[2 * i];
            Complex a2 = e - vk * fV[2 * i + 1];


            AntiDiagonalMatrix2<Complex> A1;
            A1(0, 1) = 1; A1(1, 0) = a1;
            AntiDiagonalMatrix2<Complex> A2;
            A2(0, 1) = 1; A2(1, 0) = a2;

            AntiDiagonalMatrix2<Complex> tmp = Complex(0.25*fDx)*(A1 + A2);
            Matrix kk1 = A1 * (Matrix::Identity() - tmp + Complex(FF) * A2 * A1).inverse() * (Matrix::Identity() - Complex(1. / 4 * fDx + BB)*A2);
            Matrix kk2 = A2 * (Matrix::Identity() - tmp + Complex(FF) * A1 * A2).inverse() * (Matrix::Identity() - Complex(1. / 4 * fDx + CC)*A1);

#ifdef SMALL_ROUND_ERROR
            Matrix tr_mI = fDx * Complex(0.5)*(kk1 + kk2);
            tr = Matrix::Identity() + tr_mI;
            Tail_SMALLROUNDERR();
#else
            Matrix tr_mI = Complex(fDx * 0.5)*(kk1 + kk2);
            tr = Matrix::Identity() + tr_mI;
            //tr = Matrix::Identity() + fDx * 0.5*(kk1 + kk2);
            Tail();
#endif
        }
    } else if (fMethod == SolverMethod::ExplicitRungeKuttaO6Luther1967) {

        Complex const c31 = 3. / 8;
        Complex const c32 = 1. / 8;

        Complex const c41 = 8. / 27.;
        Complex const c42 = 2. / 27;
        Complex const c43 = 8. / 27;

        Complex const c51 = 3 * (3 * sqrt(21) - 7) / 392;
        Complex const c52 = -8 * (7 - sqrt(21)) / 392;
        Complex const c53 = 48 * (7 - sqrt(21)) / 392;
        Complex const c54 = -3 * (21 - sqrt(21)) / 392;
        //Real const one1 = c51 + c52 + c53 + c54;

        Complex const c61 = -5 * (231 + 51 * sqrt(21)) / 1960;
        Complex const c62 = -40 * (7 + sqrt(21)) / 1960;
        Complex const c63 = -320 * sqrt(21) / 1960;
        Complex const c64 = 3 * (21 + 121 * sqrt(21)) / 1960;
        Complex const c65 = 392 * (6 + sqrt(21)) / 1960;
        //Real const one2 = c61 + c62 + c63 + c64 + c65;

        Complex const c71 = 15 * (22 + 7 * sqrt(21)) / 180;
        Complex const c72 = 120. / 180;
        Complex const c73 = 40 * (7 * sqrt(21) - 5) / 180;
        Complex const c74 = -63 * (3 * sqrt(21) - 2) / 180;
        Complex const c75 = -14 * (49 + 9 * sqrt(21)) / 180;
        Complex const c76 = 70 * (7 - sqrt(21)) / 180;
        //Real const one3 = c71 + c72 + c73 + c74 + c75 + c76;

        //Complex const f1 = 0;
        //Complex const f2 = 1;
        //Complex const f3 = 1. / 2;
        //Complex const f4 = 2. / 3;
        //Complex const f5 = (7 - sqrt(21)) / 14;
        //Complex const f6 = (7 + sqrt(21)) / 14;
        //Complex const f7 = 1;

        Complex const b1 = 9. / 180;
        Complex const b3 = 64. / 180;
        Complex const b5 = 49. / 180;
        Complex const b6 = 49. / 180;
        Complex const b7 = 9. / 180;;


        for (size_t i = 0; i < fNBins; ++i) {
            Matrix tr;
            Complex a1 = e - vk * fV[7 * i];
            Complex a2 = e - vk * fV[7 * i + 1];
            Complex a3 = e - vk * fV[7 * i + 2];
            Complex a4 = e - vk * fV[7 * i + 3];
            Complex a5 = e - vk * fV[7 * i + 4];
            Complex a6 = e - vk * fV[7 * i + 5];
            Complex a7 = e - vk * fV[7 * i + 6];

            AntiDiagonalMatrix2<Complex> A1;
            AntiDiagonalMatrix2<Complex> A2;
            AntiDiagonalMatrix2<Complex> A3;
            AntiDiagonalMatrix2<Complex> A4;
            AntiDiagonalMatrix2<Complex> A5;
            AntiDiagonalMatrix2<Complex> A6;
            AntiDiagonalMatrix2<Complex> A7;
            A1(0, 1) = 1; A1(1, 0) = a1;
            A2(0, 1) = 1; A2(1, 0) = a2;
            A3(0, 1) = 1; A3(1, 0) = a3;
            A4(0, 1) = 1; A4(1, 0) = a4;
            A5(0, 1) = 1; A5(1, 0) = a5;
            A6(0, 1) = 1; A6(1, 0) = a6;
            A7(0, 1) = 1; A7(1, 0) = a7;


            AntiDiagonalMatrix2<Complex> K1 = Complex(fDx) * A1;
            Matrix K2 = Complex(fDx) * A2 * (Matrix::Identity() + K1);
            Matrix K3 = Complex(fDx) * A3 * (Matrix::Identity() + c31 * K1 + c32 * K2);
            Matrix K4 = Complex(fDx) * A4 * (Matrix::Identity() + c41 * K1 + c42 * K2 + c43 * K3);
            Matrix K5 = Complex(fDx) * A5 * (Matrix::Identity() + c51 * K1 + c52 * K2 + c53 * K3 + c54 * K4);
            Matrix K6 = Complex(fDx) * A6 * (Matrix::Identity() + c61 * K1 + c62 * K2 + c63 * K3 + c64 * K4 + c65 * K5);
            Matrix K7 = Complex(fDx) * A7 * (Matrix::Identity() + c71 * K1 + c72 * K2 + c73 * K3 + c74 * K4 + c75 * K5 + c76 * K6);

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
    CalculateFinalJFromPsi();

    {
        Real k0 = sqrt(2 * fMass*(fE - fV0)) / fHbar;
        Real k1 = sqrt(2 * fMass*(fE - fV1)) / fHbar;
        ColVec2<Complex> inc(Complex(1.), I * k0);
        ColVec2<Complex> ref(Complex(1.), -I * k0);
        ColVec2<Complex> trn(Complex(1.), I * k1);

        // mat * (inc + ref*r) == trn*t
        inc = mat * inc;
        ref = mat * ref;
        // inc + ref*r == trn*t
        // inc = (-ref | trn) . (r |t )^T
        Matrix RefTrn(-ref(0), trn(0), -ref(1), trn(1));
        ColVec2<Complex> rt = RefTrn.inverse() * inc;

        fR = abs2(rt(0));
        fT = abs2(rt(1));
    }


}
