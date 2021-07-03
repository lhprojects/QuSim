#include "FourierTransform.h"
#include "kissfft.hh"

#include "../third_party/fftw-3.3.5-dll64/fftw3.h"
#include <utility>

#ifdef QUSIM_USE_CUDA
#include <cufft.h>
#include "CudaUtility.h"

inline std::string my_itoa(int i)
{
    char b[10];
    sprintf(b, "%d", i);
    return b;
}
#define check_err(err) do {\
int err_code = err;\
if(err_code !=cudaSuccess)\
throw std::runtime_error("cuda error: " #err " : " + std::to_string(err_code));\
} while(0)

#endif

#include <algorithm>

struct KissFourierTransform1D : FourierTransform1D {
    kissfft<double> impl64;
    kissfft<float> impl32;
    bool inverse;
    size_t n;
    KissFourierTransform1D(size_t n, bool inverse) :
        impl64(n, inverse),
        impl32(n, inverse),
        inverse(inverse),
        n(n)
    {

    }
    void Transform(Complex64 const *from, Complex64*to) override
    {
        impl64.transform(from, to);
    }

    void Transform(Complex32 const* from, Complex32* to)
    {
        impl32.transform(from, to);
    }

};

struct FFTWFourierTransform1D : FourierTransform1D {
    fftw_plan plan;
    fftw_complex* in;
    fftw_complex* out;
    size_t fn;
    bool finverse;
    FFTWFourierTransform1D(size_t n, bool inverse)
    {
        in = fftw_alloc_complex(n);
        out = fftw_alloc_complex(n);
        plan = fftw_plan_dft_1d((int)n, in, out, (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);
        fn = n;
        finverse = inverse;
    }

    void Transform(std::complex<double> const* from, std::complex<double>* to) override
    {
        memcpy(in, from, sizeof(std::complex<double>) * fn);
        fftw_execute(plan);
        std::copy((std::complex<double>*)out, (std::complex<double>*)out + fn, to);
    }
        

    ~FFTWFourierTransform1D()
    {
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
    }
};

#if defined(QUSIM_USE_CUDA) && QUSIM_USE_CUDA
struct CudaFourierTransform1D : FourierTransform1D
{
    bool fInverse;
    cufftHandle fPlan32;
    cufftHandle fPlan64;

    CudaFourierTransform1D(size_t nx, bool inverse)
    {

        fInverse = inverse;
        cufftResult_t res32 = cufftPlan1d(&fPlan32, nx, CUFFT_C2C, 1);
        if (res32 != CUFFT_SUCCESS) {
            throw std::exception("can't create 1d plan");
        }
        cufftResult_t res64 = cufftPlan1d(&fPlan64, nx, CUFFT_Z2Z, 1);
        if (res64 != CUFFT_SUCCESS) {
            throw std::exception("can't create 1d plan");
        }
    }

    void Transform(Complex64 const* from, Complex64* to) override
    {
        auto res = cufftExecZ2Z(fPlan64, (cufftDoubleComplex*)from, (cufftDoubleComplex*)to,
            (fInverse ? CUFFT_INVERSE : CUFFT_FORWARD));
        if (res != CUFFT_SUCCESS) {
            throw std::runtime_error("can't fft! " + std::to_string(res));
        }
    }

    void Transform(Complex32 const* from, Complex32* to)
    {
        auto res = cufftExecC2C(fPlan32, (cufftComplex*)from, (cufftComplex*)to,
            (fInverse ? CUFFT_INVERSE : CUFFT_FORWARD));
        if (res != CUFFT_SUCCESS) {
            throw std::exception("can't fft!");
        }
    }

    ~CudaFourierTransform1D()
    {
        cufftDestroy(fPlan32);
        cufftDestroy(fPlan64);
    }
};
#endif

FourierTransform1D* FourierTransform1D::Create(size_t n, bool inverse, FourierTransformLibrary lib)
{
    //return CreateCUDAFourierTransform(n, inverse);
    if (lib == FourierTransformLibrary::KISS) {
        return new KissFourierTransform1D(n, inverse);
    } else if (lib == FourierTransformLibrary::FFTW) {
        return new FFTWFourierTransform1D(n, inverse);
    } else if (lib == FourierTransformLibrary::CUDA) {
#if defined(QUSIM_USE_CUDA) && QUSIM_USE_CUDA
        return new CudaFourierTransform1D(n, inverse);
#else
        throw std::runtime_error("cuda not compiled");
#endif
    } else {
        throw std::runtime_error("not supported fft library");
    }
}

struct FFT2DTransY {

    template<size_t block_size_>
    void trans(size_t nx, size_t ny,
        size_t j,
        kissfft<double> &impl1,
        Complex* to,
        Complex *t2, Complex *t3)
    {
        for (size_t i = 0; i < nx; ++i) {
            size_t const idx = CalGlobalIdx(i, j, nx, ny);
            for (int b = 0; b < block_size_; ++b) {
                t2[i + b * nx] = to[idx + b];
            }
        }
        for (int b = 0; b < block_size_; ++b) {
            impl1.transform(t2 + nx * b, t3 + nx * b);
        }
        for (size_t i = 0; i < nx; ++i) {
            size_t const idx = CalGlobalIdx(i, j, nx, ny);
            for (int b = 0; b < block_size_; ++b) {
                to[idx + b] = t3[i + b * nx];
            }
        }

    }
};

struct KissFourierTransform2D : FourierTransform2D {
    kissfft<double> impl1;
    kissfft<double> impl2;
    std::vector<std::complex<double> > tmp2;
    std::vector<std::complex<double> > tmp3;
    size_t rows;
    size_t cols;
    bool inverse;

    static size_t const fNfold = 8; // <= 3 please

    KissFourierTransform2D(size_t rows, size_t cols, bool inverse) :
        impl1(rows, inverse), impl2(cols, inverse), inverse(inverse)
    {
        this->rows = rows;
        this->cols = cols;
        tmp2.resize(fNfold * rows);
        tmp3.resize(fNfold * rows);
    }

    void Transform(std::complex<double> const *from, std::complex<double> *to) override
    {

        std::complex<double> *t2 = tmp2.data();
        std::complex<double> *t3 = tmp3.data();

#if 0
        for (size_t i = 0; i < rows; ++i) {
            impl2.transform(from + i * cols, to + i * cols);
        }



        for (size_t j = 0; j < cols; j += fNfold) {
                    
            size_t nfold = std::min(fNfold, cols - j);

            if (nfold == 1) {

                for (size_t i = 0; i < rows; ++i) {
                    t2[i] = to[i * cols + j];
                }
                impl1.transform(t2, t3);
                for (size_t i = 0; i < rows; ++i) {
                    to[i * cols + j] = t3[i];
                }

            } else if(nfold == 2) {

                auto t2_0 = t2;
                auto t2_1 = t2 + rows;
                auto t3_0 = t3;
                auto t3_1 = t3 + rows;
                size_t j1 = j + 0;
                size_t j2 = j + 1;

                for (size_t i = 0; i < rows; ++i) {
                    t2_0[i] = to[i * cols + j1];
                    t2_1[i] = to[i * cols + j2];
                }

                impl1.transform(t2_0, t3_0);
                impl1.transform(t2_1, t3_1);

                for (size_t i = 0; i < rows; ++i) {
                    to[i * cols + j1] = t3_0[i];
                    to[i * cols + j2] = t3_1[i];
                }
            }
        }

#else
        size_t const nx = rows;
        size_t const ny = cols;

        // transform for x
        for (size_t i = 0; i < nx; ++i) {
            impl2.transform(from + CalGlobalIdx(i, 0, nx, ny),
                to + CalGlobalIdx(i, 0, nx, ny));
        }

        // transform for y
        for (size_t j = 0; j < ny;) {

            if (ny - j >= fNfold) {
                FFT2DTransY().trans<fNfold>(nx, ny, j, impl1, to, t2, t3);
                j += fNfold;
            } else {
                size_t block_size = std::min(ny - j, fNfold);
                if (block_size >= 4) {
                    FFT2DTransY().trans<4>(nx, ny, j, impl1, to, t2, t3);
                    j += 4;
                    block_size -= 4;
                }
                if (block_size >= 2) {
                    FFT2DTransY().trans<2>(nx, ny, j, impl1, to, t2, t3);
                    j += 2;
                    block_size -= 2;
                }
                if (block_size >= 1) {
                    FFT2DTransY().trans<1>(nx, ny, j, impl1, to, t2, t3);
                    j += 1;
                    block_size -= 1;
                }
            }
        }


#endif
    }

    ~KissFourierTransform2D()
    {
    }

};


struct FFTWFourierTransform2D : FourierTransform2D {

    fftw_plan plan;
    fftw_complex* in;
    fftw_complex* out;
    size_t fn1;
    size_t fn2;
    bool finverse;

    FFTWFourierTransform2D(size_t n1, size_t n2, bool inverse)
    {
        in = (fftw_complex*)fftw_malloc(n1 * n2 * sizeof(fftw_complex));
        out = (fftw_complex*)fftw_malloc(n1 * n2 * sizeof(fftw_complex));
        // we are using row major, fftw too
        plan = fftw_plan_dft_2d((int)n1, (int)n2, in, out,
            (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);
        fn1 = n1;
        fn2 = n2;
        finverse = inverse;
    }

    void Transform(std::complex<double> const* from, std::complex<double>* to) override
    {
        
        //memcpy(in, from, sizeof(std::complex<double>) * fn1 * fn2);

        std::copy((std::complex<double>*)from, (std::complex<double>*)from + fn1 * fn2, (std::complex<double>*)in);
        fftw_execute(plan);
        std::copy((std::complex<double>*)out, (std::complex<double>*)out + fn1 * fn2, to);
        
    }

    ~FFTWFourierTransform2D()
    {
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
    }

};

#if QUSIM_USE_CUDA
struct CudaFourierTransform2D : FourierTransform2D
{
    bool fInverse;
    cufftHandle fPlan;
    CudaFourierTransform2D(size_t n1, size_t n2, bool inverse)
    {
        fInverse = inverse;
        cufftResult_t res = cufftPlan2d(&fPlan, n1, n2, CUFFT_Z2Z);
        if (res != CUFFT_SUCCESS) {
            throw std::exception("can't create 2d plan");
        }
    }

    void Transform(std::complex<double> const* from, std::complex<double>* to) override
    {
        check_err(cufftExecZ2Z(fPlan, (cufftDoubleComplex*)from, (cufftDoubleComplex*)to,
            (fInverse ? CUFFT_INVERSE : CUFFT_FORWARD)));
    }

    ~CudaFourierTransform2D()
    {
        cufftDestroy(fPlan);
    }
};
#endif

FourierTransform2D * FourierTransform2D::Create(size_t n1, size_t n2, bool inverse, FourierTransformLibrary lib)
{
    if (lib == FourierTransformLibrary::KISS) {
        return new KissFourierTransform2D(n1, n2, inverse);
    } else if (lib == FourierTransformLibrary::FFTW) {
        return new FFTWFourierTransform2D(n1, n2, inverse);
    } else if (lib == FourierTransformLibrary::CUDA) {
#if QUSIM_USE_CUDA
        return new CudaFourierTransform2D(n1, n2, inverse);
#else
        throw std::runtime_error("cuda not compiled");
#endif
    } else {
        throw std::runtime_error("not supported fft library");
    }
}

struct KissFourierTransform3D : FourierTransform3D {
    kissfft<double> impl1;
    kissfft<double> impl2;
    kissfft<double> impl3;
    std::vector<std::complex<double> > tmp2;
    std::vector<std::complex<double> > tmp2_out;
    std::vector<std::complex<double> > tmp1;
    std::vector<std::complex<double> > tmp1_out;
    size_t fn1;
    size_t fn2;
    size_t fn3;

    KissFourierTransform3D(size_t n1, size_t n2, size_t n3, bool inverse) :
        impl1(n1, inverse),
        impl2(n2, inverse),
        impl3(n3, inverse)
    {
        this->fn1 = n1;
        this->fn2 = n2;
        this->fn3 = n3;
        tmp1.resize(fn1);
        tmp1_out.resize(fn1);
        tmp2.resize(fn2);
        tmp2_out.resize(fn2);
    }
    void Transform(std::complex<double> const *from, std::complex<double> *to) override
    {
        size_t const nx = fn1;
        size_t const ny = fn2;
        size_t const nz = fn3;
        auto* t1 = tmp1.data();
        auto* t1_out = tmp1_out.data();
        auto* t2 = tmp2.data();
        auto* t2_out = tmp2_out.data();

        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                size_t const idx = CalGlobalIdx(i, j, 0, nx, ny, nz);
                impl3.transform(from + idx, to + idx);
            }
        }

        for (size_t i = 0; i < nx; ++i) {
            for (size_t k = 0; k < nz; ++k) {

                for (size_t j = 0; j < ny; ++j) {
                    size_t idx = CalGlobalIdx(i, j, k, nx, ny, nz);
                    t2[j] = to[idx];
                }
                impl2.transform(t2, t2_out);
                for (size_t j = 0; j < ny; ++j) {
                    size_t idx = CalGlobalIdx(i, j, k, nx, ny, nz);
                    to[idx] = t2_out[j];
                }
            }
        }

        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                for (size_t i = 0; i < nx; ++i) {
                    t1[i] = to[CalGlobalIdx(i, j, k, nx, ny, nz)];
                }
                impl1.transform(t1, t1_out);
                for (size_t i = 0; i < nx; ++i) {
                    to[CalGlobalIdx(i, j, k, nx, ny, nz)] = t1_out[i];
                }
            }
        }
    }

    ~KissFourierTransform3D()
    {
    }

};

struct FFTWFourierTransform3D : FourierTransform3D {

    fftw_plan plan;
    fftw_complex* in;
    fftw_complex* out;
    size_t fn1;
    size_t fn2;
    size_t fn3;
    bool finverse;

    FFTWFourierTransform3D(size_t n1, size_t n2, size_t n3, bool inverse)
    {
        in = (fftw_complex*)fftw_malloc(n1 * n2 * n3 * sizeof(fftw_complex));
        out = (fftw_complex*)fftw_malloc(n1 * n2 * n3 * sizeof(fftw_complex));

        plan = fftw_plan_dft_3d((int)n1, (int)n2, (int)n3, in, out,
            (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);
        fn1 = n1;
        fn2 = n2;
        fn3 = n3;
        finverse = inverse;
    }

    void Transform(std::complex<double> const* from, std::complex<double>* to) override
    {
        std::copy((std::complex<double>*)from, (std::complex<double>*)from + fn1 * fn2 * fn3, (std::complex<double>*)in);
        fftw_execute(plan);
        std::copy((std::complex<double>*)out, (std::complex<double>*)out + fn1 * fn2 * fn3, to);
    }

    ~FFTWFourierTransform3D()
    {
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
    }

};

#if QUSIM_USE_CUDA
struct CudaFourierTransform3D : FourierTransform3D
{
    bool fInverse;
    cufftHandle fPlan;
    CudaFourierTransform3D(size_t n1, size_t n2, size_t n3, bool inverse)
    {
        fInverse = inverse;
        cufftResult_t res = cufftPlan3d(&fPlan, n1, n2, n3, CUFFT_Z2Z);
        if (res != CUFFT_SUCCESS) {
            throw std::exception("can't create 3d plan");
        }
    }

    void Transform(std::complex<double> const* from, std::complex<double>* to) override
    {
        check_err(cufftExec(fPlan, (cufftDoubleComplex*)from, (cufftDoubleComplex*)to,
            (fInverse ? CUFFT_INVERSE : CUFFT_FORWARD)));
    }

    ~CudaFourierTransform3D()
    {
        cufftDestroy(fPlan);
    }
};
#endif

FourierTransform3D * FourierTransform3D::Create(size_t n1, size_t n2, size_t n3, bool inverse, FourierTransformLibrary lib)
{
    if (lib == FourierTransformLibrary::KISS) {
        return new KissFourierTransform3D(n1, n2, n3, inverse);
    } else if (lib == FourierTransformLibrary::FFTW) {
        return new FFTWFourierTransform3D(n1, n2, n3, inverse);
    } else if (lib == FourierTransformLibrary::CUDA) {
#if defined(QUSIM_USE_CUDA) && QUSIM_USE_CUDA
        return new CudaFourierTransform3D(n1, n2, n3, inverse);
#else
        throw std::runtime_error("cuda not compiled");
#endif
    } else {
        throw std::runtime_error("not supported fft library");
    }
}
