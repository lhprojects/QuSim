
#include "SplittingMethod1D.h"
#include "SplittingUtils.h"

void SplittingMethod1D::InitSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization,
    Complex dt, bool force_normalization_each_step,
    std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
    BoundaryCondition b, SolverMethod solver,
    Real mass, Real hbar,
    OptionsImpl const &opts)
{
    QuEvolver1DImpl::InitSystem1D(psi, force_normalization,
        dt, force_normalization_each_step,
        vs, x0, x1, n, b, solver,
        mass, hbar, opts);

    fFourierTransformOptions.Init(opts, fDeviceType);

    if (fBoundaryCondition == BoundaryCondition::Period) {
        fFTPsi = fDevice->Alloc<Complex>(n);

        fft_N.reset(FourierTransform1D::Create(fN, false, fFourierTransformOptions.fLib));
        inv_fft_N.reset(FourierTransform1D::Create(fN, true, fFourierTransformOptions.fLib));

    } else if (fBoundaryCondition == BoundaryCondition::InfiniteWall) {
        mutable_cast(fIWPsi) = fDevice->Alloc<Complex>(2 * n);
        mutable_cast(fIWKPsi) = fDevice->Alloc<Complex>(2 * n);

        inv_fft_2N.reset(FourierTransform1D::Create(2 * fN, true, fFourierTransformOptions.fLib));
    }
}


void SplittingMethod1D::UpdatePsi()
{
    QuUpdatePsi(this);
}


void SplittingMethod1D::ExpT(Complex* psi, Real tt) const
{
    if (fBoundaryCondition == BoundaryCondition::Period) {

        double n2 = fDevice->Norm2(psi, fN);
        fft_N->Transform(psi, fFTPsi);

        // Dk = 2 Pi / (a * N)
        // k = i Dk
        // K = k^2/(2m)
        // exp(-I *K dt /hbar)
        // exp(f*i^2)
        Real Dk = 2 * Pi / (fDx * fN);
        Real DT = QuSqr(Dk * fHbar) / (2 * fMass);
        Complex f = -I * DT * fDt * tt / fHbar;

        // exp(-I K Dt tt/hbar)
        fDevice->MulExpK1D(fFTPsi, f, fN);

        inv_fft_N->Transform(fFTPsi, psi);
        fDevice->Scale(psi, 1. / fN, fN);
        double n3 = fDevice->Norm2(psi, fN);
        double n4 = fDevice->Norm2(psi, fN);

    } else if (fBoundaryCondition == BoundaryCondition::InfiniteWall) {

        fDevice->Copy(fIWPsi, psi, fN);
        fDevice->Set(fIWPsi, 0, 0.);
        fDevice->Set(fIWPsi, fN, 0.);
        fDevice->CopyReverseMinu(fIWPsi + 2 * fN - 1, fIWPsi + 1, fN - 1);

        inv_fft_2N->Transform(fIWPsi, fIWKPsi);


        RealType const Dp = fHbar * Pi / (fDx * fN);
        RealType const DT = Dp * Dp / (2 * fMass);
        fDevice->MulExpIdx2(fIWKPsi, -I * DT * fDt * tt / fHbar, fN);

        fDevice->Set(fIWKPsi,  0, 0.);
        fDevice->Set(fIWKPsi, fN, 0.);
        fDevice->CopyReverseMinu(fIWKPsi + 2 * fN - 1, fIWKPsi + 1, fN - 1);

        inv_fft_2N->Transform(fIWKPsi, fIWPsi);
        fDevice->Scale(psi, fIWPsi, -1. / (2. * fN), fN);
    }
}

// vpsi = psi exp(-i/hbar Dt V)
void SplittingMethod1D::ExpV(Complex *psi, double t) const
{
    Complex f = -I / fHbar * fDt * t;
    fDevice->MulExp(psi, f, fV, fN);
}

Real SplittingMethod1D::CalKinEn() const
{
    if (fBoundaryCondition == BoundaryCondition::Period) {
        //Zero(fFTPsi);

        fft_N->Transform(fPsi, fFTPsi);

        // k = Dk * i
        // T = i^2 Dk**2/(2m) 
        // T = i^2 DT^2
        Real Dk = 2 * Pi / (fDx * fN);
        Real DT = Dk * Dk / (2 * fMass);

        Real Kin = fDevice->Abs2K1D(fFTPsi, DT, fN);
        Real N2 = fDevice->Norm2(fFTPsi, fN);
        Kin = Kin / N2;
        return Kin;


    } else if (fBoundaryCondition == BoundaryCondition::InfiniteWall) {

        fDevice->Copy(fIWPsi, fPsi, fN);
        fDevice->Set(fIWPsi,  0, 0.);
        fDevice->Set(fIWPsi, fN, 0.);
        fDevice->CopyReverseMinu(fIWPsi + 2 * fN - 1, fIWPsi + 1, fN - 1);
        inv_fft_2N->Transform(fIWPsi, fIWKPsi);

        RealType const Dk = Pi / (fDx * fN);
        RealType const DT = QuSqr(Dk * fHbar) / (2 * fMass);
        RealType const kinen = fDevice->Abs2Idx2(fIWKPsi, fN) * DT;
        RealType const norm2 = fDevice->Norm2(fIWKPsi, fN);

        RealType j2 = fDevice->Abs2Idx2(fIWKPsi, fN);
        RealType j1 = fDevice->Abs2Idx(fIWKPsi, fN);
        RealType kv2 = sqrt(j2 / norm2)*Dk;
        RealType kv = Dk * j1 / norm2;

        return kinen / norm2;

    }
    return 0;
}
