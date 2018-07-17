#include "SolverImpl.h"

SolverImpl::SolverImpl() : fMass(0), fHbar(0), fE(0)
{
}

void SolverImpl::initSystem(
	Real en,
	Real mass,
	Real hbar)
{
	const_cast<Real&>(fE) = en;
	const_cast<Real&>(fMass) = mass;
	const_cast<Real&>(fHbar) = hbar;
}


SolverImpl1D::SolverImpl1D()
{
}

void SolverImpl1D::initSystem1D(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n,
	Real en, Complex initPsi, Complex initPsiPrime,
	SolverMethod met,
	Real mass, Real hbar)
{
	initSystem(en, mass, hbar);

	fMethod = met;
	fNPoints = n;
	fNBins = n - 1;
	fX0 = x0;
	fDx = (x1 - x0) / fNBins;
	fVFunc = v;

	if (met == SolverMethod::ImplicitMidpointMethod) {
		fV.resize(fNBins);
		for (size_t i = 0; i < fNBins; ++i) {
			Real x = x0 + (i + 0.5)*fDx;
			fV[i] = fVFunc(x).real();
		}
	} else if(met == SolverMethod::ExplicitRungeKuttaO4Classical) {
		fV.resize(3 * fNBins);
		for (size_t i = 0; i < fNBins; ++i) {
			Real x1 = x0 + (i + 0.0)*fDx;
			Real x2 = x0 + (i + 0.5)*fDx;
			Real x3 = x0 + (i + 1.0)*fDx;
			fV[3 * i] = fVFunc(x1).real();
			fV[3 * i + 1] = fVFunc(x2).real();
			fV[3 * i + 2] = fVFunc(x3).real();
		}
	} else if (met == SolverMethod::GaussLegendreO4) {
		fV.resize(2 * fNBins);
		Real f1 = 1. / 2 - sqrt(3) / 6;
		Real f2 = 1. / 2 + sqrt(3) / 6;
		for (size_t i = 0; i < fNBins; ++i) {
			Real x1 = x0 + (i + f1)*fDx;
			Real x2 = x0 + (i + f2)*fDx;
			fV[2 * i] = fVFunc(x1).real();
			fV[2 * i + 1] = fVFunc(x2).real();
		}
	} else {
		throw std::runtime_error("Unsupported method");
	}
	
	fPsi.resize(fNPoints);
	fPsiPrime.resize(fNPoints);

	fPsi[0] = initPsi;
	fPsiPrime[0] = initPsiPrime;

	fMass = mass;
	fHbar = hbar;
	fE = en;

	T11 = 1;
	T12 = 0;
	T21 = 0;
	T22 = 1;

	fV0 = fVFunc(x0).real();
	fV1 = fVFunc(x1).real();
}

void SolverImpl1D::Calculate()
{

	Complex psi = fPsi[0];
	Complex psiPrime = fPsiPrime[0];
	Real e = -2 * fMass / (fHbar*fHbar) * fE;
	Real vk = -2 * fMass / (fHbar*fHbar);

	fInitJ = (std::conj(psi)*psiPrime - psi * std::conj(psiPrime)).imag();
	//          |  psi       |
	//  Psi =   |            |
	//          |  psi^prime |
	//
	//   d Psi / dx = A Psi
	//  
	//               | 0    1 |  
	//   A   =       |        |  
	//               | a    0 |  
	//
	//   a = -2m/hbar^2 (E - V)
	//
	//   2rd Order
	//                       1 + 1/2 A(x + 1/2dx) dx
	//      Psi(x + dx) =  -------------------------- Psi(x)
	//                       1 - 1/2 A(x + 1/2dx) dx
	//   Classical Runge Kutta 4rd Order
	//                      
	//      Psi(x + dx) =   () Psi(x)
	//                       

	//Real amax = -10000;
	//Real amin = +10000;

	for (size_t i = 0; i < fNBins; ++i) {

		Real f_11;
		Real f_12;
		Real f_21;
		Real f_22;
		if (fMethod == SolverMethod::ImplicitMidpointMethod) {
			Real a = e - vk * fV[i];
			//amax = std::max(a, amax);
			//amin = std::min(a, amin);

			Real c1 = 1 + 0.25 * a * fDx*fDx;
			Real c2 = 1 - 0.25 * a * fDx*fDx;
			
			f_11 = c1 / c2;
			f_12 = fDx / (c2);
			f_21 = a * f_12;
			f_22 = f_11;


		} else if (fMethod == SolverMethod::ExplicitRungeKuttaO4Classical) {
			Real a1 = e - vk * fV[3 * i];
			Real a2 = e - vk * fV[3 * i + 1];
			Real a3 = a2;
			Real a4 = e - vk * fV[3 * i + 2];
			//amax = std::max(a, amax);
			//amin = std::min(a, amin);

			Eigen::Matrix2d A1;
			A1(0, 0) = A1(1, 1) = 0; A1(0, 1) = 1; A1(1, 0) = a1;
			Eigen::Matrix2d A2;
			A2(0, 0) = A2(1, 1) = 0; A2(0, 1) = 1; A2(1, 0) = a2;
			Eigen::Matrix2d const &A3 = A2;
			Eigen::Matrix2d A4;
			A4(0, 0) = A4(1, 1) = 0; A4(0, 1) = 1; A4(1, 0) = a4;

			Eigen::Matrix2d const &K1 = A1;
			Eigen::Matrix2d K2 = A2 * (Eigen::Matrix2d::Identity() + 1. / 2 * fDx*K1);
			Eigen::Matrix2d K3 = A3 * (Eigen::Matrix2d::Identity() + 1. / 2 * fDx*K2);
			Eigen::Matrix2d K4 = A4 * (Eigen::Matrix2d::Identity() + fDx*K3);

			Eigen::Matrix2d ARK4 = Eigen::Matrix2d::Identity() + 1. / 6 * fDx * (
				K1 + 2 * K2 + 2 * K3 + K4
				);

			f_11 = ARK4(0, 0);
			f_12 = ARK4(0, 1);
			f_21 = ARK4(1, 0);
			f_22 = ARK4(1, 1);


		} else if(fMethod == SolverMethod::GaussLegendreO4){
			Real a1 = e - vk * fV[2 * i];
			Real a2 = e - vk * fV[2 * i + 1];

			Eigen::Matrix2d A1;
			A1(0, 0) = A1(1, 1) = 0; A1(0, 1) = 1; A1(1, 0) = a1;
			Eigen::Matrix2d A2;
			A2(0, 0) = A2(1, 1) = 0; A2(0, 1) = 1; A2(1, 0) = a2;

			Real a12 = 1. / 4 - sqrt(3) / 6;
			Real a21 = 1. / 4 + sqrt(3) / 6;
			Eigen::Matrix2d AA = A1.inverse() - 1. / 4 * fDx * Eigen::Matrix2d::Identity();
			//Eigen::Matrix2d AA_inv = AA.inverse();
			Real BB = -a12 * fDx;
			Real CC = -a21 * fDx;
			Eigen::Matrix2d DD = A2.inverse() - 1. / 4 * fDx * Eigen::Matrix2d::Identity();
			//Eigen::Matrix2d DD_inv = DD.inverse();

			//Eigen::Matrix2d kk1 = (AA / BB - DD_inv * CC).inverse() * (Eigen::Matrix2d::Identity() / BB - DD_inv);
			//Eigen::Matrix2d kk2 = (AA_inv * BB - DD / CC).inverse() * (AA_inv - Eigen::Matrix2d::Identity()/CC);
			Eigen::Matrix2d kk1 = (DD * AA - CC * BB * Eigen::Matrix2d::Identity()).inverse() * (DD - Eigen::Matrix2d::Identity() * BB);
			Eigen::Matrix2d kk2 = (BB * CC * Eigen::Matrix2d::Identity() - AA * DD).inverse() * (CC * Eigen::Matrix2d::Identity() - AA);

			Eigen::Matrix2d glo4 = Eigen::Matrix2d::Identity() + fDx * 0.5*(kk1 + kk2);

			f_11 = glo4(0, 0);
			f_12 = glo4(0, 1);
			f_21 = glo4(1, 0);
			f_22 = glo4(1, 1);

		}
		
		if (0) {
			Eigen::Matrix2cd U;
			Eigen::Matrix2cd O;
			U(0, 0) = f_11;
			U(0, 1) = f_12;
			U(1, 0) = f_21;
			U(1, 1) = f_22;
			O(0, 0) = O(1, 1) = 0;
			O(0, 1) = 1;
			O(1, 0) = -1;
			
			Eigen::Matrix2cd err = U.transpose() * O * U - O;
			if (err.norm() > 1E-10) {
				err(0, 0);
			}
		}

		Real t11 = T11;
		Real t12 = T12;
		Real t21 = T21;
		Real t22 = T22;
		
		T11 = f_11 * t11 + f_12 * t21;
		T12 = f_11 * t12 + f_12 * t22;
		T21 = f_21 * t11 + f_22 * t21;
		T22 = f_21 * t12 + f_22 * t22;

		Complex psi_ = psi;
		Complex psiPrime_  = psiPrime;
		psi = psi_ * f_11 + psiPrime_ * f_12;
		psiPrime = psi_ * f_21 + psiPrime_ * f_22;

		fPsi[i + 1] = psi;
		fPsiPrime[i + 1] = psiPrime;
	}
	fFinalJ = (std::conj(psi)*psiPrime - psi * std::conj(psiPrime)).imag();

	fTMat(0, 0) = T11;
	fTMat(0, 1) = T12;
	fTMat(1, 0) = T21;
	fTMat(1, 1) = T22;

	do {
		fT = 0;
		fR = 0;
		Real kin1 = fE - fV0;
		if (kin1 <= 0) break;

		Real kin2 = fE - fV1;
		if (kin2 <= 0) break;

		Real k1 = sqrt(2 * fMass*kin1) / fHbar;
		Real k2 = sqrt(2 * fMass*kin2) / fHbar;

		//t (1 , I k2)^T =  TMat (1 + r, i k1 - i k1 r)^T
		// t = T11 + T12 i k1 + r (T11 - T12 i k1)
		// i k2 t = i k2 (T11 + T12 i k1) + r i k2 (T11 - T12 i k1)
		// i k2 t =  T21 + T22 i k1 + r (T21 - T22 i k1)
		Complex r = -(I * k2*(T11 + T12 * I*k1) - (T21 + T22 * I*k1)) / (I*k2*(T11 - T12*I*k1) - (T21 - T22*I*k1));
		
		Complex One = (-T12 * T21 + T11 * T22);
		One = 1;
		Complex t = 2 * k1 * One / (k2*T11 - I * k1*k2*T12 + I * T21 + k1 * T22);

		fT = abs2(t) * k2 / k1;
		fR = abs2(r);
	} while (false);

}
