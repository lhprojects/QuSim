#ifndef QUSIM_UTILS_H
#define QUSIM_UTILS_H

#include "QuSim.h"
#include <stdint.h>
#include <stddef.h>

template<class T>
T& mutable_cast(T const& r)
{
    return const_cast<T&>(r);
}

template<class T>
T* &mutable_ptr_cast(T const * &r)
{
    return const_cast<T*&>(r);
}

template<class T>
T*& mutable_ptr_cast(T *& r)
{
    return const_cast<T*&>(r);
}

template<class T>
T* const& mutable_ptr_cast(T const* const& r)
{
    return const_cast<T* const &>(r);
}

template<class T>
T* const& mutable_ptr_cast(T* const& r)
{
    return const_cast<T* const &>(r);
}

inline Real Abs2(Complex c)
{
    return real(c) * real(c) + imag(c) * imag(c);
}

inline Real Abs2(Real c)
{
    return c * c;
}

template<class R>
inline R QuGetReal(std::complex<R> c)
{
    return c.real();
}

inline Real32 QuGetReal(Real32 c)
{
    return c;
}

inline Real64 QuGetReal(Real64 c)
{
    return c;
}

template<class R>
inline auto make_complex(R r, R i)
{
    return std::complex<R>(r, i);
}

constexpr inline ptrdiff_t QuIdx(ptrdiff_t i, ptrdiff_t j, ptrdiff_t nx, ptrdiff_t ny)
{
    return j + ny * i;
}

constexpr inline ptrdiff_t QuIdxFold(ptrdiff_t i, ptrdiff_t j, ptrdiff_t nx, ptrdiff_t ny)
{
    if (i < 0) i += nx;
    else if (i >= nx) i -= nx;

    if (j < 0) j += ny;
    else if (j >= ny) j -= ny;

    return QuIdx(i, j, nx, ny);
}

template<class C>
struct ScalarTrait {
};

template<class S>
inline S QuSqr(S s)
{
    return s * s;
}

template<>
struct ScalarTrait<Complex32> {
    using RealType = Real32;
    using ComplexType = Complex32;
    static Complex32 One() { return 1.; }
    static RealType RealOne() { return 1.; }
};

template<>
struct ScalarTrait<Real32> {
    using RealType = Real32;
    using ComplexType = Complex32;
    static Real32 One() { return 1.; }
    static RealType RealOne() { return 1.; }
};

template<>
struct ScalarTrait<Complex64> {
    using RealType = Real64;
    using ComplexType = Complex64;
    static Complex64 One() { return 1.; }
    static RealType RealOne() { return 1.; }
};

template<>
struct ScalarTrait<Real64> {
    using RealType = Real64;
    using ComplexType = Complex64;
    static Real64 One() { return 1.; }
    static RealType RealOne() { return 1.; }
};

struct Delayed {

    Delayed(std::function<void()> f) : f(f)
    {
    }
    ~Delayed()
    {
        f();
    }
    std::function<void()> f;
};

// https://godbolt.org/#g:!((g:!((g:!((h:codeEditor,i:(fontScale:14,fontUsePx:'0',j:1,lang:c%2B%2B,selection:(endColumn:1,endLineNumber:2,positionColumn:1,positionLineNumber:2,selectionStartColumn:1,selectionStartLineNumber:2,startColumn:1,startLineNumber:2),source:'%0Avoid+foo(double+%26a)+%7B%0A++++a+*%3D+1%3B%0A%7D%0A%0A%0Avoid+bar(double+%26a)+%7B%0A++++a+%2B%3D+0.%3B%0A%7D%0A'),l:'5',n:'0',o:'C%2B%2B+source+%231',t:'0')),k:53.22175732217573,l:'4',n:'0',o:'',s:0,t:'0'),(g:!((h:compiler,i:(compiler:g111,filters:(b:'0',binary:'1',commentOnly:'0',demangle:'0',directives:'0',execute:'1',intel:'0',libraryCode:'0',trim:'1'),fontScale:14,fontUsePx:'0',j:1,lang:c%2B%2B,libs:!(),options:'',selection:(endColumn:1,endLineNumber:1,positionColumn:1,positionLineNumber:1,selectionStartColumn:1,selectionStartLineNumber:1,startColumn:1,startLineNumber:1),source:1),l:'5',n:'0',o:'x86-64+gcc+11.1+(Editor+%231,+Compiler+%231)+C%2B%2B',t:'0')),k:46.77824267782427,l:'4',n:'0',o:'',s:0,t:'0')),l:'2',n:'0',o:'',t:'0')),version:4
struct Zero {

    template<class U>
    U operator+(U u) const
    {
        return +u;
    }

    template<class U>
    U operator-(U u) const
    {
        return -u;
    }

    template<class U>
    U operator*(U u) const
    {
        return U();
    }

    template<class U>
    U operator/(U u) const
    {
        return U();
    }
};

struct One;
struct Minus1 {

    template<class U>
    U operator*(U u) const
    {
        return -u;
    }

    One operator-() const;
};

struct One {

    template<class U>
    U operator*(U u) const
    {
        return u;
    }
};

inline One Minus1::operator-() const
{
    return One();
}

template<class T>
inline void SafeFree(T *&p)
{
    if (p) {
        free(p);
        p = nullptr;
    }
}
#endif
