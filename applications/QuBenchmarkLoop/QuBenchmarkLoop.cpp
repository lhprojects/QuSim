#include <chrono>
#include <algorithm>
#include <complex>
#include <execution>
#include "../QuBenchmark/Benchmark.h"

template<class Op>
void timeit(std::string const& title, Op op)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    op();
    auto t1 = std::chrono::high_resolution_clock::now();
    printf("%20s %9d us\n", title.c_str(),
        (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());

}

template<class R>
R Sqr(R r) { return r * r; }

double Abs2(double r) { return r * r; }
double Abs2(std::complex<double> r) { return r.real() * r.real() + r.imag() * r.imag(); }

template<class Integer>
struct int_iter {

    using difference_type = ptrdiff_t;
    using value_type = Integer;
    using pointer = Integer*;
    using reference = Integer&;
    using iterator_category = std::random_access_iterator_tag;

    Integer fInt;

    int_iter() : fInt() {}
    int_iter(Integer i) : fInt(i) {}

    int_iter operator+(difference_type i) const
    {
        return int_iter(fInt + i);
    }

    int_iter& operator+=(difference_type i)
    {
        fInt += i;
        return *this;
    }

    int_iter operator-(difference_type i) const
    {
        return int_iter(fInt - i);
    }

    int_iter& operator-=(difference_type i)
    {
        fInt -= i;
        return *this;
    }

    difference_type operator-(int_iter r) const
    {
        return difference_type(fInt - r.fInt);
    }

    int_iter& operator++()
    {
        fInt += 1;
        return *this;
    }
    int_iter& operator--()
    {
        fInt -= 1;
        return *this;
    }

    Integer operator*() const
    {
        return fInt;
    }

    bool operator!=(int_iter r) const
    {
        return fInt != r.fInt;
    }
    bool operator==(int_iter r) const
    {
        return fInt == r.fInt;
    }
};


template<class S>
struct Trait {
    static std::string title(std::string str)
    {
        return str + " ( double)";
    }

    static S OneORI()
    {
        return 1;
    }
};

template<>
struct Trait<std::complex<double> > {
    static std::string title(std::string str)
    {
        return str + " (complex)";
    }

    static std::complex<double> OneORI()
    {
        return std::complex<double>(0, 1);
    }
};

size_t nx = 10000;
size_t ny = 10000;

double alpha1 = 1.2;
double alpha2 = 1.3;


struct op_SetOne
{
    template<class U>
    auto operator()(U u)
    {
        return 1.;
    }
};
struct op_Exp {
    template<class U>
    auto operator()(U u)
    {
        return exp(u);
    }
};
struct op_Add {
    template<class U, class V>
    auto operator()(U u, V v)
    {
        return u + v;
    }
};
struct op_Mul {
    template<class U, class V>
    auto operator()(U u, V v)
    {
        return u * v;
    }
};
volatile double dont_optimize_out;

template<class S>
void benchmark_loop()
{

    begin_section("benchmark_loop");
    S* in = new S[nx * ny];
    std::unique_ptr<S> auto_free(in);

    std::fill(in, in + nx * ny, 1.);
    using Tr = Trait<S>;

    timeit(Tr::title("nestedloop   assign"), [in]() {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                size_t idx = j + i * ny;
                in[idx] = (1. * i + 1. * j);
            }
        }
        });

    timeit(Tr::title("nestedloop      mul"), [in, alpha1 = alpha1, alpha2 = alpha2]() {
        double r = 0;
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                size_t idx = j + i * ny;
                r += Abs2(in[idx]) * (alpha1 * Sqr(1. * i) + alpha2 * Sqr(1. * j));
            }
        }
        dont_optimize_out = r;
        });

    timeit(Tr::title("nestedloop      exp"), [in, alpha1 = alpha1, alpha2 = alpha2]() {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                size_t idx = j + i * ny;
                in[idx] *= exp(Tr::OneORI() * (alpha1 * Sqr(1. * i) + alpha2 * Sqr(1. * j)));
            }
        }
        });


    timeit(Tr::title("flatteloop   assign"), [in]() {
        for (size_t idx = 0; idx < nx * ny; ++idx) {
            size_t i = idx / ny;
            size_t j = idx % ny;
            in[idx] = (1. * i + 1. * j);
        }
        });

    timeit(Tr::title("flatteloop      mul"), [in, alpha1 = alpha1, alpha2 = alpha2]() {
        double r = 0;
        for (size_t idx = 0; idx < nx * ny; ++idx) {
            size_t i = idx / ny;
            size_t j = idx % ny;
            r += Abs2(in[idx]) * (alpha1 * Sqr(1. * i) + alpha2 * Sqr(1. * j));
        }
        dont_optimize_out = r;
        });

    timeit(Tr::title("flatteloop      exp"), [in, alpha1 = alpha1, alpha2 = alpha2]() {
        for (size_t idx = 0; idx < nx * ny; ++idx) {
            size_t i = idx / ny;
            size_t j = idx % ny;
            in[idx] *= exp(Tr::OneORI() * (alpha1 * Sqr(1. * i) + alpha2 * Sqr(1. * j)));
        }
        });
    timeit(Tr::title(" transform_red  mul"), [in, alpha1 = alpha1, alpha2 = alpha2]() {
        dont_optimize_out = std::transform_reduce(std::execution::par, in, in + nx * ny, int_iter<size_t>(0), double(0),
            op_Add(),
            [alpha1, alpha2](S s, size_t idx) {
                size_t i = idx / ny;
                size_t j = idx % ny;
                return Abs2(s) * (alpha1 * Sqr(1. * i) + alpha2 * Sqr(1. * j));
            });
        });
    timeit(Tr::title(" transform      exp"), [in, alpha1 = alpha1, alpha2 = alpha2]() {
        std::transform(std::execution::par, in, in + nx * ny, int_iter<size_t>(0), in,
            [alpha1, alpha2](S s, size_t idx) {
                size_t i = idx / ny;
                size_t j = idx % ny;
                return s * exp(Tr::OneORI() * (alpha1 * Sqr(1. * i) + alpha2 * Sqr(1. * j)));
            });
        });

    end_section();
}


int nbt = 10000000;

template<class R>
struct QuComp {
    QuComp() : fReal(), fImag() { }
    QuComp(R real, R imag) : fReal(real), fImag(imag) { }
    R fReal;
    R fImag;
};

template<class R>
inline QuComp<R> operator*(QuComp<R> l, QuComp<R> r)
{
    return QuComp<R>(l.fReal * r.fReal - l.fImag * r.fImag,
        l.fReal * r.fImag + l.fImag * r.fReal);
}

void benchmark_transform_fill()
{
    begin_section("benchmark_transform_fill");
    double* in = new double[nbt];
    std::unique_ptr<double> auto_free(in);
    timeit("      seq", [in]() {
        std::transform(std::execution::seq, in, in + nbt, in, op_SetOne());
        });
    timeit("      par", [in]() {
        std::transform(std::execution::par, in, in + nbt, in, op_SetOne());
        });
    timeit("par_unseq", [in]() {
        std::transform(std::execution::par_unseq, in, in + nbt, in, op_SetOne());
        });
    timeit("par_unseq", [in]() {
        std::transform(std::execution::par_unseq, in, in + nbt, in, op_SetOne());
        });
    timeit("par_unseq", [in]() {
        std::transform(std::execution::par_unseq, in, in + nbt, in, op_SetOne());
        });
    timeit("      par", [in]() {
        std::transform(std::execution::par, in, in + nbt, in, op_SetOne());
        });
    timeit("      seq", [in]() {
        std::transform(std::execution::seq, in, in + nbt, in, op_SetOne());
        });
    timeit("      par", [in]() {
        std::transform(std::execution::par, in, in + nbt, in, op_SetOne());
        });
    timeit("par_unseq", [in]() {
        std::transform(std::execution::par_unseq, in, in + nbt, in, op_SetOne());
        });
    end_section();
}

void benchmark_transform_std_complex_mul()
{
    begin_section("benchmark_transform_std_complex_mul");
    std::complex<double>* in1 = new std::complex<double>[nbt];
    QuComp<double>* in2 = new QuComp<double>[nbt];

    std::unique_ptr<std::complex<double> > auto_free1(in1);
    std::unique_ptr<QuComp<double> > auto_free2(in2);

    for (int i = 0; i < 3; ++i) {
        timeit("      seq & std::complex", [in1]() {
            std::transform(std::execution::seq, in1, in1 + nbt, in1, in1, op_Mul());
            });
        timeit("      par & std::complex", [in1]() {
            std::transform(std::execution::par, in1, in1 + nbt, in1, in1, op_Mul());
            });
        timeit("par_unseq & std::complex", [in1]() {
            std::transform(std::execution::par_unseq, in1, in1 + nbt, in1, in1, op_Mul());
            });
    }
    for (int i = 0; i < 3; ++i) {
        timeit("      seq &       QuComp", [in2]() {
            std::transform(std::execution::seq, in2, in2 + nbt, in2, in2, op_Mul());
            });
        timeit("      par &       QuComp", [in2]() {
            std::transform(std::execution::par, in2, in2 + nbt, in2, in2, op_Mul());
            });
        timeit("par_unseq &       QuComp", [in2]() {
            std::transform(std::execution::par_unseq, in2, in2 + nbt, in2, in2, op_Mul());
            });
    }
    end_section();
}

void benchmark_transform_Exp()
{
    begin_section("benchmark_transform_Exp");
    double* in = new double[nbt];
    std::unique_ptr<double> auto_free(in);
    timeit("      seq", [in]() {
        std::transform(std::execution::seq, in, in + nbt, in, op_Exp());
        });
    timeit("      par", [in]() {
        std::transform(std::execution::par, in, in + nbt, in, op_Exp());
        });
    timeit("par_unseq", [in]() {
        std::transform(std::execution::par_unseq, in, in + nbt, in, op_Exp());
        });
    timeit("par_unseq", [in]() {
        std::transform(std::execution::par_unseq, in, in + nbt, in, op_Exp());
        });
    timeit("par_unseq", [in]() {
        std::transform(std::execution::par_unseq, in, in + nbt, in, op_Exp());
        });
    timeit("      par", [in]() {
        std::transform(std::execution::par, in, in + nbt, in, op_Exp());
        });
    timeit("      seq", [in]() {
        std::transform(std::execution::seq, in, in + nbt, in, op_Exp());
        });
    timeit("      par", [in]() {
        std::transform(std::execution::par, in, in + nbt, in, op_Exp());
        });
    timeit("par_unseq", [in]() {
        std::transform(std::execution::par_unseq, in, in + nbt, in, op_Exp());
        });
    end_section();
}

void benchmark_transform_Dot()
{
    begin_section("benchmark_transform_Dot");
    double* in = new double[nbt];
    std::unique_ptr<double> auto_free(in);
    timeit("      seq", [in]() {
        dont_optimize_out = std::transform_reduce(std::execution::seq, in, in + nbt, in, double(), op_Add(), op_Mul());
        });
    timeit("      par", [in]() {
        dont_optimize_out = std::transform_reduce(std::execution::par, in, in + nbt, in, double(), op_Add(), op_Mul());
        });
    timeit("par_unseq", [in]() {
        dont_optimize_out = std::transform_reduce(std::execution::par_unseq, in, in + nbt, in, double(), op_Add(), op_Mul());
        });
    timeit("par_unseq", [in]() {
        dont_optimize_out = std::transform_reduce(std::execution::par_unseq, in, in + nbt, in, double(), op_Add(), op_Mul());
        });
    timeit("par_unseq", [in]() {
        dont_optimize_out = std::transform_reduce(std::execution::par_unseq, in, in + nbt, in, double(), op_Add(), op_Mul());
        });
    timeit("      par", [in]() {
        dont_optimize_out = std::transform_reduce(std::execution::par, in, in + nbt, in, double(), op_Add(), op_Mul());
        });
    timeit("      seq", [in]() {
        dont_optimize_out = std::transform_reduce(std::execution::seq, in, in + nbt, in, double(), op_Add(), op_Mul());
        });
    timeit("      par", [in]() {
        dont_optimize_out = std::transform_reduce(std::execution::par, in, in + nbt, in, double(), op_Add(), op_Mul());
        });
    timeit("par_unseq", [in]() {
        dont_optimize_out = std::transform_reduce(std::execution::par_unseq, in, in + nbt, in, double(), op_Add(), op_Mul());
        });
    end_section();
}

#include <Windows.h>
void benchmark_for_each_big()
{
    begin_section("benchmark_for_each sleep(1 second)");

    static int const sleep = 1;

    timeit("par", [] {
        std::for_each(std::execution::par, int_iter<size_t>(0), int_iter<size_t>(16), [](size_t) {
            std::this_thread::sleep_for(std::chrono::seconds(sleep));
            });
        });
    timeit("seq", [] {
        std::for_each(std::execution::seq, int_iter<size_t>(0), int_iter<size_t>(16), [](size_t) {
            std::this_thread::sleep_for(std::chrono::seconds(sleep));
            });
        });

    timeit("par", [] {
        std::for_each(std::execution::par, int_iter<size_t>(0), int_iter<size_t>(16), [](size_t) {
            std::this_thread::sleep_for(std::chrono::seconds(sleep));
            });
        });
    timeit("seq", [] {
        std::for_each(std::execution::seq, int_iter<size_t>(0), int_iter<size_t>(16), [](size_t) {
            std::this_thread::sleep_for(std::chrono::seconds(sleep));
            });
        });

    timeit("par", [] {
        std::for_each(std::execution::par, int_iter<size_t>(0), int_iter<size_t>(16), [](size_t) {
            std::this_thread::sleep_for(std::chrono::seconds(sleep));
            });
        });
    timeit("seq", [] {
        std::for_each(std::execution::seq, int_iter<size_t>(0), int_iter<size_t>(16), [](size_t) {
            std::this_thread::sleep_for(std::chrono::seconds(sleep));
            });
        });
    end_section();
}



#include "../src/kissfft.hh"
#include "../src/FourierTransform.h"

#include <memory>

void benchmark_fouier_2d()
{
    size_t nx = 10000;
    size_t ny = 10000;
    begin_section("benchmark_fouier fft + invfft 10000x10000");
    Complex* in = new Complex[nx * ny];
    Complex* out = new Complex[nx * ny];
    std::unique_ptr<Complex> autofree1(in);
    std::unique_ptr<Complex> autofree2(out);

    std::fill(in, in + nx * ny, 1.);

    std::unique_ptr<FourierTransform2D> kiss(FourierTransform2D::Create(nx, ny, false, FourierTransformLibrary::KISS));
    std::unique_ptr<FourierTransform2D> ikiss(FourierTransform2D::Create(nx, ny, true, FourierTransformLibrary::KISS));

    std::unique_ptr<FourierTransform2D> fftw(FourierTransform2D::Create(nx, ny, false, FourierTransformLibrary::FFTW));
    std::unique_ptr<FourierTransform2D> ifftw(FourierTransform2D::Create(nx, ny, true, FourierTransformLibrary::FFTW));

    timeit("kiss", [&]() {
        kiss->Transform(in, out);
        ikiss->Transform(out, in);
        });

    timeit("FFTW", [&]() {
        fftw->Transform(in, out);
        ifftw->Transform(out, in);
        });

    timeit("kiss", [&]() {
        kiss->Transform(in, out);
        ikiss->Transform(out, in);
        });

    timeit("FFTW", [&]() {
        fftw->Transform(in, out);
        ifftw->Transform(out, in);
        });

    timeit("kiss", [&]() {
        kiss->Transform(in, out);
        ikiss->Transform(out, in);
        });

    timeit("FFTW", [&]() {
        fftw->Transform(in, out);
        ifftw->Transform(out, in);
        });
    timeit("kiss", [&]() {
        kiss->Transform(in, out);
        ikiss->Transform(out, in);
        });

    timeit("FFTW", [&]() {
        fftw->Transform(in, out);
        ifftw->Transform(out, in);
        });

    end_section();
}

void benchmark_fouier_3d()
{
    begin_section("benchmark_fouier fft + invfft 500x500x500");
    size_t nx = 500;
    size_t ny = 500;
    size_t nz = 500;
    Complex* in = new Complex[nx * ny * nz];
    Complex* out = new Complex[nx * ny * nz];
    std::unique_ptr<Complex> autofree1(in);
    std::unique_ptr<Complex> autofree2(out);

    std::fill(in, in + nx * ny, 1.);

#ifdef QUSIM_USE_CUDA
    Device::Create();
    Complex* in_d = new Complexp[];
#endif

    std::unique_ptr<FourierTransform3D> kiss(FourierTransform3D::Create(nx, ny, nz, false, FourierTransformLibrary::KISS));
    std::unique_ptr<FourierTransform3D> ikiss(FourierTransform3D::Create(nx, ny, nz, true, FourierTransformLibrary::KISS));

    std::unique_ptr<FourierTransform3D> fftw(FourierTransform3D::Create(nx, ny, nz, false, FourierTransformLibrary::FFTW));
    std::unique_ptr<FourierTransform3D> ifftw(FourierTransform3D::Create(nx, ny, nz, true, FourierTransformLibrary::FFTW));

    timeit("kiss", [&]() {
        kiss->Transform(in, out);
        ikiss->Transform(out, in);
        });

    timeit("FFTW", [&]() {
        fftw->Transform(in, out);
        ifftw->Transform(out, in);
        });

    timeit("kiss", [&]() {
        kiss->Transform(in, out);
        ikiss->Transform(out, in);
        });

    timeit("FFTW", [&]() {
        fftw->Transform(in, out);
        ifftw->Transform(out, in);
        });

    timeit("kiss", [&]() {
        kiss->Transform(in, out);
        ikiss->Transform(out, in);
        });

    timeit("FFTW", [&]() {
        fftw->Transform(in, out);
        ifftw->Transform(out, in);
        });
    timeit("kiss", [&]() {
        kiss->Transform(in, out);
        ikiss->Transform(out, in);
        });

    timeit("FFTW", [&]() {
        fftw->Transform(in, out);
        ifftw->Transform(out, in);
        });

    end_section();
}

int main()
{

    benchmark_transform_std_complex_mul();
    benchmark_transform_Dot();
    benchmark_transform_Exp();
    benchmark_transform_fill();
    benchmark_loop<double>();
    benchmark_loop<std::complex<double> >();
    benchmark_for_each_big();

    benchmark_fouier_2d();
    benchmark_fouier_3d();
    return 0;
}
