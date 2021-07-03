#pragma once

#include "Device.h"
#include <memory>
#include <complex>
#include <algorithm>
#include <execution>
#include <utility>
#include <tuple>
#include <cstddef>
#include <cstring> // for memmove

#ifndef QUSIM_HOST_DEVICE
#define QUSIM_HOST_DEVICE
#endif



template<class Integer>
struct int_iter {

    using difference_type = ptrdiff_t;
    using value_type = Integer;
    using pointer = Integer*;
    using reference = const Integer; // uh..., you can't change the value the iterator points to
    using iterator_category = std::random_access_iterator_tag;

    Integer fInt;

    int_iter() : fInt() {}
    int_iter(Integer i) : fInt(i) {}

    Integer operator[](size_t n) const
    {
        return fInt + n;
    }

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
        return difference_type(fInt - r);
    }

    int_iter& operator++()
    {
        fInt += 1;
        return *this;
    }
    int_iter& operator--()
    {
        fInt += 1;
        return *this;
    }

    Integer operator*() const
    {
        return fInt;
    }

};

template<class Iter1, class Iter2>
struct zip_iter
{
    Iter1 iter1;
    Iter2 iter2;
    using value_type1 = typename std::iterator_traits<Iter1>::value_type;
    using value_type2 = typename std::iterator_traits<Iter2>::value_type;
    using reference1 = typename std::iterator_traits<Iter1>::reference;
    using reference2 = typename std::iterator_traits<Iter2>::reference;

    using zip_value = std::pair<value_type1, value_type2>;
    using reference = std::pair<reference1, reference2>;

    using iterator_category = std::random_access_iterator_tag;
    using difference_type = ptrdiff_t;
    using value_type = zip_value;
    using pointer = value_type*;

    zip_iter() : iter1(), iter2()
    {
    }

    zip_iter(Iter1 iter1, Iter2 iter2) : iter1(iter1), iter2(iter2)
    {
    }

    zip_iter& operator+=(difference_type n)
    {
        iter1 += n;
        iter2 += n;
    }

    zip_iter operator+(difference_type n) const
    {
        return zip_iter(iter1 + n, iter2 + n);
    }

    zip_iter& operator-=(difference_type n)
    {
        iter1 -= n;
        iter2 -= n;
    }

    zip_iter operator-(difference_type n) const
    {
        return zip_iter(iter1 - n, iter2 - n);
    }

    difference_type operator-(zip_iter r) const
    {
        return r.iter1 - iter1;
    }

    zip_iter& operator++()
    {
        iter1 += 1;
        iter2 += 1;
        return *this;
    }
    zip_iter& operator--()
    {
        iter1 -= 1;
        iter2 -= 1;
        return *this;
    }

    reference operator*() const
    {
        return reference(*iter1, *iter2);
    }

    reference operator[](size_t n) const
    {
        return reference(iter1[n], iter2[n]);
    }
};

template<class Iter1, class Iter2, class Iter3>
struct zip_iter3
{
    Iter1 iter1;
    Iter2 iter2;
    Iter3 iter3;
    using value_type1 = typename std::iterator_traits<Iter1>::value_type;
    using value_type2 = typename std::iterator_traits<Iter2>::value_type;
    using value_type3 = typename std::iterator_traits<Iter3>::value_type;
    using reference1 = typename std::iterator_traits<Iter1>::reference;
    using reference2 = typename std::iterator_traits<Iter2>::reference;
    using reference3 = typename std::iterator_traits<Iter3>::reference;

    using zip_value = std::tuple<value_type1, value_type2, value_type3>;
    using reference = std::tuple<reference1, reference2, reference3>;

    using iterator_category = std::random_access_iterator_tag;
    using difference_type = ptrdiff_t;
    using value_type = zip_value;
    using pointer = value_type*;

    zip_iter3() : iter1(), iter2(), iter3()
    {
    }

    zip_iter3(Iter1 iter1, Iter2 iter2, Iter3 iter3) : iter1(iter1), iter2(iter2), iter3(iter3)
    {
    }

    zip_iter3& operator+=(difference_type n)
    {
        iter1 += n;
        iter2 += n;
        iter3 += n;
    }

    zip_iter3 operator+(difference_type n) const
    {
        return zip_iter3(iter1 + n, iter2 + n, iter3 + n);
    }

    zip_iter3& operator-=(difference_type n)
    {
        iter1 -= n;
        iter2 -= n;
        iter3 -= n;
    }

    zip_iter3 operator-(difference_type n) const
    {
        return zip_iter(iter1 - n, iter2 - n, iter3);
    }

    difference_type operator-(zip_iter3 r) const
    {
        return r.iter3 - iter3, r.iter2 - iter2, r.iter1 - iter1;
    }

    zip_iter3& operator++()
    {
        iter1 += 1;
        iter2 += 1;
        iter3 += 1;
        return *this;
    }
    zip_iter3& operator--()
    {
        iter1 -= 1;
        iter2 -= 1;
        iter3 -= 1;
        return *this;
    }

    reference operator*() const
    {
        return reference(*iter1, *iter2, *iter3);
    }

    reference operator[](size_t n) const
    {
        return reference(iter1[n], iter2[n], iter3[n]);
    }
};


template<class E>
struct ExecTrait {

    using Exec = E;
    static const bool cOnMainMemory = true;

    template<class T>
    using IntIter = int_iter<T>;

    template<class T>
    using ReverseIter = std::reverse_iterator<T>;

    template<class T>
    using DevicePtr = T*;

    using ComplexType = std::complex<double>;
    using RealType = double;

    template<class R>
    struct ScalarTrait {
        using DeviceType = R;
        using ComplexType = std::complex<R>;
        using RealType = R;
        static RealType real_max()
        {
            return std::numeric_limits<R>::max();
        }
    };

    template<class R>
    struct ScalarTrait<const R> {
        using DeviceType = const R;
        using ComplexType = const std::complex<R>;
        using RealType = const R;
        static RealType real_max()
        {
            return std::numeric_limits<R>::max();
        }
    };

    template<class R>
    struct ScalarTrait<std::complex<R>> {
        using DeviceType = std::complex<R>;
        using ComplexType = std::complex<R>;
        using RealType = R;
        static RealType real_max()
        {
            return std::numeric_limits<R>::max();
        }
    };

    template<class R>
    struct ScalarTrait<const std::complex<R>> {
        using DeviceType = const  std::complex<R>;
        using ComplexType = const std::complex<R>;
        using RealType = const R;
        static RealType real_max()
        {
            return std::numeric_limits<R>::max();
        }
    };

    static void MemToDevice(Exec& exec, void* out, void const* in, size_t bytes)
    {
        std::copy(exec, (std::byte*)in, (std::byte*)in + bytes, (std::byte*)out);
    }

    static void MemToHost(Exec& exec, void* out, void const* in, size_t bytes)
    {
        std::copy(exec, (std::byte*)in, (std::byte*)in + bytes, (std::byte*)out);
    }

    static void MoveToDevice(void* out, void const* in, size_t bytes)
    {
        memmove(out, in, bytes);
    }

    static void MoveToHost(void* out, void const* in, size_t bytes)
    {
        memmove(out, in, bytes);
    }

    static void* Alloc(size_t bytes)
    {
        void *p = malloc(bytes);
        if (!p) {
            throw std::runtime_error("can't alloc " + std::to_string(bytes) + " bytes");
        }
        return p;
    }

    static void Free(void* p)
    {
        return free(p);
    }

    template<class Iter, class T>
    static void Fill(Exec &exec, Iter b1, Iter e1, T const &v)
    {
        std::fill(exec, b1, e1, v);
    }

    template<class Iter, class OutIter>
    static void Copy(Exec& exec, Iter b1, Iter e1, OutIter b2)
    {
        std::copy(exec, b1, e1, b2);
    }

    template<class Iter, class OutIter, class UnitaryOp>
    static void Transform(Exec& exec, Iter b1, Iter e1, OutIter outIter, UnitaryOp op)
    {
        std::transform(exec, b1, e1, outIter, op);
    }

    template<class Iter,  class Iter2, class OutIter, class BinaryOp>
    static void Transform(Exec& exec, Iter b1, Iter e1, Iter2 b2, OutIter outIter, BinaryOp op)
    {
        std::transform(exec, b1, e1, b2, outIter, op);
    }

    template<class TransOp>
    struct op_Expand {

        op_Expand(TransOp op) : transOp(op)
        {
        }

        TransOp transOp;
        template<class U, class V>
        auto operator()(U u, V v)
        {
            return transOp(u, std::get<0>(v), std::get<1>(v));
        }
    };

    template<class TransOp>
    struct op_Expand3 {

        TransOp transOp;
        op_Expand3(TransOp op) : transOp(op)
        {
        }

        template<class U, class V>
        auto operator()(U u, V v)
        {
            return transOp(u, std::get<0>(v), std::get<1>(v), std::get<2>(v));
        }
    };

    template<class Iter, class Iter2, class Iter3, class OutIter, class TransOp>
    static void Transform(Exec& exec, Iter b1, Iter e1, Iter2 b2, Iter3 b3, OutIter outIter, TransOp op)
    {
        using ZipIter = zip_iter<Iter2, Iter3>;
        using OpExpand = op_Expand<TransOp>;
        std::transform(exec, b1, e1, ZipIter(b2, b3), outIter, OpExpand(op));
    }

    template<class Iter, class ReduceOp, class T, class TransOp>
    static T TransformReduce(Exec& exec, Iter b1, Iter e1, T init, ReduceOp reduce_op, TransOp transOp)
    {
        return std::transform_reduce(exec, b1, e1, init, reduce_op, transOp);
    }

    template<class Iter, class Iter2, class T, class ReduceOp, class TransOp>
    static T TransformReduce(Exec& exec, Iter b1, Iter e1, Iter2 b2, T init, ReduceOp reduceOp, TransOp transOp)
    {
        return std::transform_reduce(exec, b1, e1, b2, init, reduceOp, transOp);
    }

    template<class Iter, class Iter2, class Iter3, class T, class ReduceOp, class TransOp>
    static T TransformReduce(Exec& exec, Iter b1, Iter e1, Iter2 b2, Iter3 b3,
        T init, ReduceOp reduceOp, TransOp transOp)
    {
        using ZipIter = zip_iter<Iter2, Iter3>;
        using OpExpand2 = op_Expand<TransOp>;
        return std::transform_reduce(exec, b1, e1, ZipIter(b2, b3),
            init, reduceOp, OpExpand2(transOp));
    }

    template<class Iter, class Iter2, class Iter3, class Iter4, class T, class ReduceOp, class TransOp>
    static T TransformReduce(Exec& exec, Iter b1, Iter e1, Iter2 b2, Iter3 b3, Iter4 b4,
        T init, ReduceOp reduceOp, TransOp transOp)
    {
        using ZipIter = zip_iter3<Iter2, Iter3, Iter4>;
        using OpExpand = op_Expand3<TransOp>;
        return std::transform_reduce(exec, b1, e1, ZipIter(b2, b3, b4),
            init, reduceOp, OpExpand(transOp));
    }
    static inline QUSIM_HOST_DEVICE double Abs2(double s)
    {
        return s * s;
    }
    static inline QUSIM_HOST_DEVICE float Abs2(float s)
    {
        return s * s;
    }

    template<class R>
    static inline QUSIM_HOST_DEVICE R Abs2(std::complex<R> s)
    {
        return s.imag() * s.imag() + s.real() * s.real();
    }

    static inline QUSIM_HOST_DEVICE std::complex<double> IMul(double c)
    {
        return std::complex<double>(double(0), c);
    }

    static inline QUSIM_HOST_DEVICE std::complex<double> IMul(std::complex<double> c)
    {
        return std::complex<double>(-c.imag(), c.real());
    }

    static inline QUSIM_HOST_DEVICE double GetReal(std::complex<double> c)
    {
        return c.real();
    }

    static inline QUSIM_HOST_DEVICE double GetImag(std::complex<double> c)
    {
        return c.imag();
    }

};

template<class ExecPar>
struct DeviceTemplate : Device 
{

    using Trait = ExecTrait<ExecPar>;
    using ExecType = typename Trait::Exec;

    template<class T>
    using IntIter = typename Trait::template IntIter<T>;

    template<class T>
    using ReverseIter = typename Trait::template ReverseIter<T>;

    template<class T>
    using DevicePtr = typename Trait::template DevicePtr<T>;

    using ComplexType = typename Trait::ComplexType;
    using RealType = typename Trait::RealType;

    template<class T>
    using ScalarTrait = typename Trait::ScalarTrait<T>;

    DevicePtr<RealType> DevicePtrCast(DReal* p)
    {
        return DevicePtr<RealType>((RealType*)p);
    }

    DevicePtr<const RealType> DevicePtrCast(DReal const * p)
    {
        return DevicePtr<const RealType>((RealType*)p);
    }

    DevicePtr<ComplexType> DevicePtrCast(DComplex* p)
    {
        return DevicePtr<ComplexType>((ComplexType*)p);
    }

    DevicePtr<const ComplexType> DevicePtrCast(DComplex const* p)
    {
        return DevicePtr<const ComplexType>((ComplexType*)p);
    }

    template<class Iter>
    ReverseIter<Iter> MakeReverseIter(Iter iter)
    {
        return ReverseIter<Iter>(iter);
    }

    ExecType fExecPar;
    QUSIM_HOST_DEVICE ExecType& Exec()
    {
        return fExecPar;
    }


    template<class S>
    static inline QUSIM_HOST_DEVICE S QuSqr(S s)
    {
        return s * s;
    }
    template<class R>
    struct op_LinearUpdate {

        QUSIM_HOST_DEVICE op_LinearUpdate(R gamma) : gamma(gamma) {}
        R const gamma;

        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V v) const
        {
            return u * (1 - gamma) + v * gamma;
        }
    };

    QUSIM_HOST_DEVICE static inline size_t fold_half(size_t i, size_t N)
    {
        return (i > N / 2) ? (N - i) : i;
    }

    QUSIM_HOST_DEVICE static inline void cal_index_2d(size_t &ix, size_t &iy, size_t idx, size_t nx, size_t ny)
    {
        auto const ix_ = idx / ny;
        auto const iy_ = idx % ny;
        ix = ix_;
        iy = iy_;
    }

    QUSIM_HOST_DEVICE static inline void cal_index_3d(size_t &ix, size_t &iy, size_t &iz,
        size_t idx, size_t nx, size_t ny, size_t nz)
    {
        size_t const k = idx % nz;
        size_t const idx2 = idx / nz;
        size_t const j = idx2 % ny;
        size_t const i = idx2 / ny;

        ix = i;
        iy = j;
        iz = k;
    }

    template<class P>
    struct op_ExpI {

        QUSIM_HOST_DEVICE op_ExpI(P a) : alpha(a) {}

        P const alpha;
        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U o1) const
        {
            return exp(Trait::IMul(o1 * alpha));
        }
    };

    struct op_Abs2 {
        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U o1) const
        {
            return Trait::Abs2(o1);
        }
    };

    template<class P>
    struct op_MulExpI {

        P const alpha;
        QUSIM_HOST_DEVICE op_MulExpI(P a) : alpha(a) {}

        template<class V, class U>
        QUSIM_HOST_DEVICE auto operator()(V o1, U o2) const
        {
            return o1 * Trait::IMul(alpha*o2);
        }
    };

    template<class P>
    struct op_ExpMul {

        QUSIM_HOST_DEVICE op_ExpMul(P a) : alpha(a) {}

        P const alpha;
        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U o1) const
        {
            return exp(alpha * o1);
        }
    };

    template<class P>
    struct op_MulExpMul {

        QUSIM_HOST_DEVICE op_MulExpMul(P a) : alpha(a) {}

        P const alpha;
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U o1, V o2) const
        {
            return o1 * exp(alpha * o2);
        }
    };

    template<class P>
    struct op_Scale {
        QUSIM_HOST_DEVICE op_Scale(P a) : alpha(a) {}
        P const alpha;

        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U o1) const
        {
            return alpha * o1;
        }

    };

    struct op_Min {
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V v) const
        {
            return u < v ? u : v;
        }
    };

    struct op_Max {
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V v) const
        {
            return u >= v ? u : v;
        }
    };

    struct op_Abs2Mul {
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V v) const
        {
            return Trait::Abs2(u) * v;
        }
    };

    struct op_GetReal {

        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U u) const
        {
            return Trait::GetReal(u);
        }
    };

    struct op_GetImag {

        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U u) const
        {
            return Trait::GetImag(u);
        }
    };

    struct op_Add {
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V v) const
        {
            return u + v;
        }
    };

    struct op_Sub {
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V v) const
        {
            return u - v;
        }
    };

    struct op_Mul {
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U o1, V o2) const
        {
            return o1 * o2;
        }
    };

    struct op_MulR {
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U o1, V o2) const
        {
            return o1 * Trait::GetReal(o2);
        }
    };

    template<class P>
    struct op_MulMul {

        QUSIM_HOST_DEVICE op_MulMul(P a) : alpha(a) {}
        P const alpha;
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U o1, V o2) const
        {
            return o1 * alpha * o2;
        }
    };

    template<class R>
    struct op_Abs2K1D {

        QUSIM_HOST_DEVICE op_Abs2K1D(size_t sz, R ax) : Nx(sz), ax(ax) {}

        R const ax;
        size_t const Nx;

        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U u, size_t idx) const
        {
            size_t i = fold_half(idx, Nx);
            using RealType = typename ScalarTrait<U>::RealType;
            return Trait::Abs2(u) * QuSqr(RealType(i)) * ax;
        }
    };

    template<class R, class S>
    struct op_MulExpK1D {

        QUSIM_HOST_DEVICE op_MulExpK1D(size_t sz, S ax) : Nx(sz), ax(ax) {}

        S const ax;
        size_t const Nx;

        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V idx) const
        {
            size_t const i = fold_half(idx, Nx);
            using RealType = typename ScalarTrait<U>::RealType;
            return u * exp(ax * QuSqr(RealType(i)));
        }
    };

    template<class R>
    struct op_Abs2K2D {

        QUSIM_HOST_DEVICE op_Abs2K2D(size_t nx, size_t ny, R ax, R ay) : Nx(nx), Ny(ny), ax(ax), ay(ay) {}

        R const ax;
        R const ay;
        size_t const Nx;
        size_t const Ny;

        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V idx) const
        {
            size_t i, j;
            cal_index_2d(i, j, idx, Nx, Ny);

            i = fold_half(i, Nx);
            j = fold_half(j, Ny);
            using RealType = typename ScalarTrait<U>::RealType;
            return Trait::Abs2(u) * (QuSqr(RealType(i)) * ax + QuSqr(RealType(j)) * ay);
        }
    };

    template<class R, class S>
    struct op_MulExpK2D {

        QUSIM_HOST_DEVICE op_MulExpK2D(size_t nx, size_t ny, R  ax, R ay, S alpha)
            : Nx(nx), Ny(ny), ax(ax), ay(ay), alpha(alpha)
        {
        }

        S const alpha;
        R const ax;
        R const ay;
        size_t const Nx;
        size_t const Ny;

        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V idx) const
        {
            size_t i, j;
            cal_index_2d(i, j, idx, Nx, Ny);

            i = fold_half(i, Nx);
            j = fold_half(j, Ny);
            using RealType = typename ScalarTrait<U>::RealType;
            return u * exp(alpha * (QuSqr(RealType(i)) * ax + QuSqr(RealType(j)) * ay));
        }
    };

    template<class R>
    struct op_Abs2K3D {

        QUSIM_HOST_DEVICE op_Abs2K3D(size_t nx, size_t ny, size_t nz, R ax, R ay, R az)
            : Nx(nx), Ny(ny), Nz(nz), ax(ax), ay(ay), az(az)
        {
        }

        R const ax;
        R const ay;
        R const az;
        size_t const Nx;
        size_t const Ny;
        size_t const Nz;

        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V idx) const
        {
            size_t i, j, k;
            cal_index_3d(i, j, k, idx, Nx, Ny, Nz);
            i = fold_half(i, Nx);
            j = fold_half(j, Ny);
            k = fold_half(k, Nz);

            using RealType = typename ScalarTrait<U>::RealType;
            return Trait::Abs2(u) * (QuSqr(RealType(i)) * ax + QuSqr(RealType(j)) * ay + QuSqr(RealType(k)) * az);
        }
    };

    template<class R, class S>
    struct op_MulExpK3D {

        QUSIM_HOST_DEVICE op_MulExpK3D(size_t nx, size_t ny, size_t nz, R ax, R ay, R az, S alpha)
            : Nx(nx), Ny(ny), Nz(nz), ax(ax), ay(ay), az(az), alpha(alpha)
        {
        }

        S const alpha;
        R const ax;
        R const ay;
        R const az;
        size_t const Nx;
        size_t const Ny;
        size_t const Nz;

        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V idx) const
        {
            size_t i, j, k;
            cal_index_3d(i, j, k, idx, Nx, Ny, Nz);
            i = fold_half(i, Nx);
            j = fold_half(j, Ny);
            k = fold_half(k, Nz);

            using RealType = typename ScalarTrait<U>::RealType;
            return u * exp(QuSqr(RealType(i)) * ax + QuSqr(RealType(j)) * ay + QuSqr(RealType(k)) * az);
        }
    };

    template<class R>
    struct op_G01D {

        QUSIM_HOST_DEVICE op_G01D(size_t n, R en, R eps, R dk2div2m) :
            N(n), E(en), Epsilon(eps), Dpx2Div2M(dk2div2m)
        {
        }

        size_t const N;
        R const E;
        R const Epsilon;
        R const Dpx2Div2M;

        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U u, size_t v) const
        {
            using ComplexType = U;
            using RealType = RealType;

            v = fold_half(v, N);
            RealType t = QuSqr(RealType(v)) * Dpx2Div2M;
            return u / ComplexType(E - t, Epsilon);
        }
    };

    template<class R>
    struct op_G02D {

        QUSIM_HOST_DEVICE op_G02D(size_t nx, size_t ny,
            R en, R eps,
            R Dpx2Div2M, R Dpy2Div2M) :
            Nx(nx), Ny(ny),
            E(en), Epsilon(eps),
            Dpx2Div2M(Dpx2Div2M),
            Dpy2Div2M(Dpy2Div2M)
        {
        }

        size_t const Nx;
        size_t const Ny;
        R const E;
        R const Epsilon;
        R const Dpx2Div2M;
        R const Dpy2Div2M;

        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U psi, size_t idx) const
        {
            size_t i, j;
            cal_index_2d(i, j, idx, Nx, Ny);
            i = fold_half(i, Nx);
            j = fold_half(j, Ny);
            RealType t = QuSqr(RealType(i)) * Dpx2Div2M + QuSqr(RealType(j)) * Dpy2Div2M;
            return psi / ComplexType(E - t, Epsilon);
        }
    };

    template<class R>
    struct op_G03D {

        QUSIM_HOST_DEVICE op_G03D(size_t nx, size_t ny, size_t nz,
            R en, R eps,
            R dkx2div2m, R dky2div2m, R dkz2div2m) :
            Nx(nx), Ny(ny), Nz(nz),
            E(en), Epsilon(eps),
            Dkx2HbarDiv2M(dkx2div2m),
            Dky2HbarDiv2M(dky2div2m),
            Dkz2HbarDiv2M(dkz2div2m)
        {
        }

        size_t const Nx;
        size_t const Ny;
        size_t const Nz;
        R const E;
        R const Epsilon;
        R const Dkx2HbarDiv2M;
        R const Dky2HbarDiv2M;
        R const Dkz2HbarDiv2M;

        template<class C>
        QUSIM_HOST_DEVICE auto operator()(C psi, size_t idx) const
        {
            using ComplexType = C;
            size_t k = idx % Nz;
            size_t idx2 = idx / Nz;
            size_t j = idx2 % Ny;
            size_t i = idx2 / Ny;

            i = fold_half(i, Nx);
            j = fold_half(j, Ny);
            k = fold_half(k, Nz);
            RealType t = QuSqr(RealType(i)) * Dkx2HbarDiv2M + QuSqr(RealType(j)) * Dky2HbarDiv2M
                + QuSqr(RealType(k)) * Dkz2HbarDiv2M;
            return psi / ComplexType(E - t, Epsilon);
        }
    };


    template<class P>
    struct op_MulExpSquare {
        QUSIM_HOST_DEVICE op_MulExpSquare(P a) : alpha(a) {}

        P const alpha;

        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U u, size_t idx) const
        {
            using RealType = typename ScalarTrait<P>::RealType;
            return u * exp(QuSqr(RealType(idx)) * alpha);
        }
    };

    struct op_Abs2MulSquare {

        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V v) const
        {
            using RealType = typename ScalarTrait<U>::RealType;
            return Trait::Abs2(u) * QuSqr(RealType(v));
        }
    };


    struct op_AddMul {
        template<class U, class V, class P>
        QUSIM_HOST_DEVICE auto operator()(U delpsi, V psi0, P v)
        {
            return (delpsi + psi0) * v;
        }
    };

    struct op_Minus {
        template<class U>
        QUSIM_HOST_DEVICE U operator()(U u) const
        {
            return -u;
        }
    };

    template<class R>
    struct op_Vellekoop {

        QUSIM_HOST_DEVICE op_Vellekoop(R slow, R epsilon) : slow(slow), epsilon(epsilon)
        {
        }

        R const slow;
        R const epsilon;
        template<class U, class V, class P>
        QUSIM_HOST_DEVICE U operator()(U oldpsi, V newpsi, P v_) const
        {
            using ComplexType = typename ScalarTrait<R>::ComplexType;
            ComplexType gamma = 1. - Trait::IMul(v_) / epsilon;
            gamma *= slow;
            return (1. - gamma) * oldpsi + gamma * newpsi;
        }
    };

    template<class R>
    struct op_Hao1 {

        QUSIM_HOST_DEVICE op_Hao1(R slow, R epsilon) : slow(slow), epsilon(epsilon)
        {
        }

        R const slow;
        R const epsilon;
        template<class U, class V, class P>
        QUSIM_HOST_DEVICE U operator()(U oldpsi, V newpsi, P v_) const
        {
            using ComplexType = typename ScalarTrait<R>::ComplexType;
            ComplexType gamma = 1. - Trait::IMul(ComplexType(v_.real(), -v_.imag())) / epsilon;
            gamma *= slow;
            return (1. - gamma) * oldpsi + gamma * newpsi;
        }
    };

    template<class R>
    struct op_Hao2 {

        QUSIM_HOST_DEVICE op_Hao2(R slow, R epsilon) : slow(slow), epsilon(epsilon)
        {
        }

        R const slow;
        R const epsilon;

        template<class U, class V, class P>
        QUSIM_HOST_DEVICE U operator()(U oldpsi, V newpsi, P v_) const
        {

            using ComplexType = typename ScalarTrait<R>::ComplexType;
            ComplexType gamma = 2. / (1. + Trait::IMul(v_) / epsilon);
            gamma *= slow;
            return (1. - gamma) * oldpsi + gamma * newpsi;
        }
    };

    template<class R>
    struct op_Add3 {

        QUSIM_HOST_DEVICE op_Add3(R epsilon) : epsilon(epsilon) {}

        R const epsilon;

        template<class U, class V, class P>
        QUSIM_HOST_DEVICE U operator()(U u, V psi0, P v_) const
        {
            using ComplexType = typename ScalarTrait<R>::ComplexType;
            return ComplexType(v_.real(), v_.imag() + epsilon) * u + v_.real() * psi0;
        }

    };


    struct op_MinusAbs2 {
        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V v) const
        {
            return Trait::Abs2(u - v);
        }
    };

    template<class R>
    struct op_Xsec3D {
        R kx, ky, kz;
        R dx, dy, dz;
        R x0, y0, z0;
        size_t nx, ny, nz;

        op_Xsec3D(R kx, R ky, R kz,
            R dx, R dy, R dz,
            R x0, R y0, R z0,
            size_t nx, size_t ny, size_t nz) :
            kx(kx), ky(ky), kz(kz),
            dx(dx), dy(dy), dz(dz),
            x0(x0), y0(y0), z0(z0),
            nx(nx), ny(ny), nz(nz)
        {
        }

        template<class U>
        QUSIM_HOST_DEVICE auto operator()(U psi0, U psi, U v, size_t idx) const
        {
            size_t ix, iy, iz;
            cal_index_3d(ix, iy, iz, idx, nx, ny, nz);
            auto phase = (x0 + RealType(ix) * dx) * kx +
                (y0 + RealType(iy) * dy) * ky +
                (z0 + RealType(iz) * dz) * kz;

            return (psi0 + psi) * Trait::GetReal(v)
                * exp(-Trait::IMul(phase));
        }
    };

    bool OnMainMem() override
    { 
        return Trait::cOnMainMemory;
    }
       
    void MemToDevice(void* out, void const* in, size_t bytes) override
    {
        Trait::MemToDevice(Exec(), out, in, bytes);
    }

    void MemToHost(void* out, void const* in, size_t bytes) override
    {
        Trait::MemToHost(Exec(), out, in, bytes);
    }
    void* AllocMem(size_t bytes) override
    {
        return Trait::Alloc(bytes);
    }

    DComplex At(DComplex const* in, size_t idx, size_t n) override
    {
        DComplex ans;
        Trait::MoveToHost(&ans, in + idx, sizeof(DComplex));
        return ans;
    }

    DReal At(DReal const* in, size_t idx, size_t n) override
    {
        DReal ans;
        Trait::MoveToHost(&ans, in + idx, sizeof(DReal));
        return ans;
    }

    void Set(DComplex* in, size_t idx, DComplex alpha)
    {
        Trait::MoveToDevice(in + idx, &alpha, sizeof(DComplex));
    }

    void Set(DReal* in, size_t idx, DReal alpha)
    {
        Trait::MoveToDevice(in + idx, &alpha, sizeof(DReal));
    }

    void Free(void* ptr) override
    {
        return Trait::Free(ptr);
    }

    void SetZero(DComplex* in, size_t n) override
    {
        Trait::Fill(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + n),
            DComplex(0));
    }

    void SetOne(DComplex* in, size_t n) override
    {
        Trait::Fill(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + n),
            DComplex(1));
    }

    DReal Abs2Idx(DComplex const* in, size_t n) override
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + n),
            IntIter<size_t>(0),
            RealType(0),
            op_Add(),
            op_Abs2Mul());
    }

    DReal Abs2Idx(DReal const* in, size_t n) override
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + n),
            IntIter<size_t>(0),
            RealType(0),
            op_Add(),
            op_Abs2Mul());
    }

    DReal Norm2(DComplex const* in, size_t nx)
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + nx),
            RealType(0),
            op_Add(), op_Abs2());
    }

    DReal Min(DComplex const* in, size_t n) override
    {
        return sqrt(Trait::TransformReduce(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + n),
            ScalarTrait<DReal>::real_max(),
            op_Min(),
            op_Abs2()));
    }

    DReal  Max(DComplex const* in, size_t n) override
    {
        return sqrt(Trait::TransformReduce(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + n),
            DReal(0),
            op_Max(),
            op_Abs2()));
    }

    DReal SumReal(DComplex const* in, size_t n)
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + n),
            DReal(0),
            op_Add(),
            op_GetReal()
            );
    }

    DReal SumImag(DComplex const* in, size_t n)
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + n),
            DReal(0),
            op_Add(),
            op_GetImag()
            );
    }

    DReal MinusNorm2(DComplex const* in1, DComplex const* in2, size_t n) override
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(in1),
            DevicePtrCast(in1 + n),
            DevicePtrCast(in2),
            RealType(0),
            op_Add(),
            op_MinusAbs2());
    }

    DReal Norm2(DReal const* in, size_t nx)
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(in), 
            DevicePtrCast(in + nx),
            RealType(0),
            op_Add(),
            op_Abs2());
    }

    void Exp(DComplex* to, DReal const* from, DComplex alpha, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from),
            DevicePtrCast(from + nx),
            DevicePtrCast(to),
            op_ExpMul<ComplexType>(alpha));
    }

    void ExpI(DComplex* to, DReal const* from, DReal alpha, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from),
            DevicePtrCast(from + nx),
            DevicePtrCast(to),
            op_ExpI<RealType>(alpha));
    }

    void Exp(DComplex* to, DComplex const* from, DComplex alpha, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from), 
            DevicePtrCast(from + nx), 
            DevicePtrCast(to),
            op_ExpMul<ComplexType>(alpha));
    }

    void CopyReverseMinu(DComplex* out, DComplex const* in, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + nx),
            MakeReverseIter(DevicePtrCast(out + 1)),
            op_Minus());
    }

    void Copy(DComplex* out, DComplex const* in, size_t nx) override
    {
        Trait::Copy(Exec(),
            DevicePtrCast(in),
            DevicePtrCast(in + nx),
            DevicePtrCast(out));
    }

    void Copy(DReal* out, DReal const* in, size_t nx) override
    {
        Trait::Copy(Exec(), 
            DevicePtrCast(in), 
            DevicePtrCast(in + nx),
            DevicePtrCast(out));
    }

    void Scale(DReal* fromto, DReal scale, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(fromto),
            DevicePtrCast(fromto + nx),
            DevicePtrCast(fromto),
            op_Scale<RealType>(scale));
    }
    void Scale(DComplex* fromto, DReal scale, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(fromto),
            DevicePtrCast(fromto + nx),
            DevicePtrCast(fromto),
            op_Scale<RealType>(scale));
    }
    void Scale(DComplex* fromto, DComplex scale, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(fromto),
            DevicePtrCast(fromto + nx),
            DevicePtrCast(fromto),
            op_Scale<ComplexType>(scale));
    }

    void Scale(DReal* to, DReal const* from, DReal scale, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from), 
            DevicePtrCast(from + nx), 
            DevicePtrCast(to),
            op_Scale<RealType>(scale));
    }

    void Scale(DComplex* to, DComplex const* from, DReal scale, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from),
            DevicePtrCast(from + nx),
            DevicePtrCast(to),
            op_Scale<RealType>(scale));
    }

    void Scale(DComplex* to, std::complex<double> const* from, DComplex scale, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from),
            DevicePtrCast(from + nx),
            DevicePtrCast(to),
            op_Scale<ComplexType>(scale));
    }

    void Sub(DComplex* to, DComplex const* in1, DComplex const* in2, size_t n)
    {
        Trait::Transform(Exec(),
            DevicePtrCast(in1),
            DevicePtrCast(in1 + n),
            DevicePtrCast(in2),
            DevicePtrCast(to),
            op_Sub());
    }

    void MulR(DComplex* to, DComplex const* from1, DComplex const* from2, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from1),
            DevicePtrCast(from1 + nx),
            DevicePtrCast(from2),
            DevicePtrCast(to),
            op_MulR());
    }

    void Mul(DComplex* to, DComplex const* from1, DComplex const* from2, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from1),
            DevicePtrCast(from1 + nx),
            DevicePtrCast(from2),
            DevicePtrCast(to),
            op_Mul());
    }

    void Mul(DComplex* to, DReal const* from1, DComplex const* from2, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from1),
            DevicePtrCast(from1 + nx),
            DevicePtrCast(from2),
            DevicePtrCast(to),
            op_Mul());
    }

    void Mul(DComplex* to, DReal const* from1, DComplex alpha, DComplex const* from2, size_t nx) override
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(from1), 
            DevicePtrCast(from1 + nx),
            DevicePtrCast(from2), 
            DevicePtrCast(to),
            op_MulMul<ComplexType>(alpha));
    }

    void Mul(DComplex* to, DComplex const* from2, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from2),
            DevicePtrCast(from2 + nx),
            DevicePtrCast(to),
            DevicePtrCast(to),
            op_Mul());
    }

    void MulExp(DComplex* to, DComplex const* from1, DComplex alpha, DComplex const* from2, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(from1),
            DevicePtrCast(from1 + nx),
            DevicePtrCast(from2),
            DevicePtrCast(to),
            op_MulExpMul<ComplexType>(alpha));
    }

    void MulExp(DComplex* fromto, DComplex alpha, DComplex const* from2, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(fromto),
            DevicePtrCast(fromto + nx),
            DevicePtrCast(from2),
            DevicePtrCast(fromto),
            op_MulExpMul<ComplexType>(alpha));
    }

    void MulExp(DComplex* fromto, DComplex alpha, DReal const* from2, size_t nx) override
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(fromto),
            DevicePtrCast(fromto + nx),
            DevicePtrCast(from2),
            DevicePtrCast(fromto),
            op_MulExpMul<ComplexType>(alpha));
    }

    void MulExpI(DComplex* fromto, DReal alpha, DReal const* from2, size_t nx) override
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(fromto),
            DevicePtrCast(fromto + nx),
            DevicePtrCast(from2),
            DevicePtrCast(fromto),
            op_MulExpI<RealType>(alpha));
    }

    void MulExpIdx2(DComplex* psi, DComplex alpha, size_t nx) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(psi),
            DevicePtrCast(psi + nx),
            IntIter<size_t>(0),
            DevicePtrCast(psi),
            op_MulExpSquare<ComplexType>(alpha));
    }

    DReal Abs2Idx2(DComplex const* psi, size_t nx) override
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(psi),
            DevicePtrCast(psi + nx),
            IntIter<size_t>(0),
            RealType(0), op_Add(), op_Abs2MulSquare());
    }

    DReal Abs2K1D(DComplex const* psi, DReal ax, size_t nx) override
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(psi),
            DevicePtrCast(psi + nx),
            IntIter<size_t>(0),
            RealType(0),
            op_Add(),
            op_Abs2K1D<RealType>(nx, ax));
    }


    DReal Abs2K2D(DComplex const* psi, DReal ax, size_t nx, DReal ay, size_t ny) override
    {
        return Trait::TransformReduce(Exec(), 
            DevicePtrCast(psi),
            DevicePtrCast(psi + nx * ny),
            IntIter<size_t>(0),
            RealType(0),
            op_Add(),
            op_Abs2K2D<RealType>(nx, ny, ax, ay));
    }

    DReal Abs2K3D(DComplex const* psi, DReal ax, size_t nx, DReal ay, size_t ny, DReal az, size_t nz) override
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(psi),
            DevicePtrCast(psi + nx * ny * nz),
            IntIter<size_t>(0),
            RealType(0),
            op_Add(),
            op_Abs2K3D<RealType>(nx, ny, nz, ax, ay, az));
    }

    void MulExpK1D(DComplex* psi, DComplex ax, size_t nx) override
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(psi),
            DevicePtrCast(psi + nx),
            IntIter<size_t>(0),
            DevicePtrCast(psi),
            op_MulExpK1D<RealType, ComplexType>(nx, ax));
    }

    void MulExpK2D(DComplex* psi, DComplex alpha, DReal ax, size_t nx, DReal ay, size_t ny) override
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(psi),
            DevicePtrCast(psi + nx * ny),
            IntIter<size_t>(0),
            DevicePtrCast(psi),
            op_MulExpK2D<RealType, ComplexType>(nx, ny, ax, ay, alpha));
    }

    void MulExpK3D(DComplex* psi, DComplex alpha, DReal ax, size_t nx, DReal ay, size_t ny, DReal az, size_t nz) override
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(psi),
            DevicePtrCast(psi + nx * ny * nz),
            IntIter<size_t>(0),
            DevicePtrCast(psi),
            op_MulExpK3D<RealType, ComplexType>(nx, ny, nz, ax, ay, az, alpha));
    }

    DComplex Abs2Mul(DComplex const* in1, DComplex const* in2, size_t nx) override
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(in1),
            DevicePtrCast(in1 + nx),
            DevicePtrCast(in2),
            ComplexType(0),
            op_Add(), op_Abs2Mul());
    }

    DReal Abs2Mul(DComplex const* in1, DReal const* in2, size_t nx) override
    {
        return Trait::TransformReduce(Exec(), 
            DevicePtrCast(in1), 
            DevicePtrCast(in1 + nx),
            DevicePtrCast(in2),
            RealType(0),
            op_Add(), op_Abs2Mul());
    }


    void G01D(DComplex* psik, DReal E, DReal epsilon, DReal Dpx2Div2M, size_t n) override
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(psik),
            DevicePtrCast(psik + n),
            IntIter<size_t>(0),
            DevicePtrCast(psik),
            op_G01D<RealType>(n, E, epsilon, Dpx2Div2M));
    }

    void G02D(DComplex* psik, DReal E, DReal epsilon,
        DReal Dkx2HbarDiv2M,
        DReal Dky2HbarDiv2M,
        size_t nx, size_t ny) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(psik),
            DevicePtrCast(psik + nx * ny),
            IntIter<size_t>(0),
            DevicePtrCast(psik),
            op_G02D<RealType>(nx, ny, E, epsilon, Dkx2HbarDiv2M, Dky2HbarDiv2M));
    }

    void G03D(DComplex* psik, DReal E, DReal epsilon, DReal a1, DReal a2, DReal a3,
        size_t nx, size_t ny, size_t nz)
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(psik),
            DevicePtrCast(psik + nx * ny * nz),
            IntIter<size_t>(0),
            DevicePtrCast(psik),
            op_G03D<RealType>(nx, ny, nz, E, epsilon, a1, a2, a3));
    }

    virtual DComplex Xsection3D(DComplex const* psi0, DComplex const* psi, DComplex const* v,
        DReal kx, DReal ky, DReal kz,
        DReal x0, DReal y0, DReal z0,
        DReal dx, DReal dy, DReal dz,
        size_t nx, size_t ny, size_t nz)
    {
        return Trait::TransformReduce(Exec(),
            DevicePtrCast(psi0),
            DevicePtrCast(psi0 + nx * ny * nz),
            DevicePtrCast(psi),
            DevicePtrCast(v),
            IntIter<size_t>(0),
            ComplexType(0),
            op_Add(),
            op_Xsec3D<DReal>(kx, ky, kz,
                dx, dy, dz,
                x0, y0, z0,
                nx, ny, nz)
        );
    }

    void AddMul(DComplex* psi, DComplex const* psi0, DComplex const* v, size_t n) override
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(psi), 
            DevicePtrCast(psi + n),
            DevicePtrCast(psi0),
            DevicePtrCast(v),
            DevicePtrCast(psi),
            op_AddMul());
    }

    void Add3(DComplex* deltaPsix, DComplex const* psi0x, DComplex const* v, DReal epsilon, size_t n) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(deltaPsix),
            DevicePtrCast(deltaPsix + n),
            DevicePtrCast(psi0x),
            DevicePtrCast(v),
            DevicePtrCast(deltaPsix),
            op_Add3<RealType>(epsilon));
    }

    void Add3(DComplex* add3, DComplex const* deltaPsix, DComplex const* psi0x, DComplex const* v, DReal epsilon, size_t n) override
    {
        Trait::Transform(Exec(),
            DevicePtrCast(deltaPsix),
            DevicePtrCast(deltaPsix + n),
            DevicePtrCast(psi0x),
            DevicePtrCast(v),
            DevicePtrCast(add3),
            op_Add3<RealType>(epsilon));
    }

    void LinearUpdate(DComplex* psi, DComplex const* newPsi, DReal gamma, size_t n) override
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(psi), 
            DevicePtrCast(psi + n),
            DevicePtrCast(newPsi), 
            DevicePtrCast(psi),
            op_LinearUpdate<RealType>(gamma));
    }

    void Vellekoop(DComplex* psi, DComplex const* newPsi, DComplex const* V, DReal slow, DReal epsilon, size_t n)
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(psi), 
            DevicePtrCast(psi + n),
            DevicePtrCast(newPsi), 
            DevicePtrCast(V),
            DevicePtrCast(psi),
            op_Vellekoop<RealType>(slow, epsilon));
    }
    void Hao1(DComplex* psi, DComplex const* newPsi, DComplex const* V, DReal slow, DReal epsilon, size_t n)
    {
        Trait::Transform(Exec(), 
            DevicePtrCast(psi),
            DevicePtrCast(psi + n),
            DevicePtrCast(newPsi),
            DevicePtrCast(V), 
            DevicePtrCast(psi),
            op_Hao1<RealType>(slow, epsilon));
    }
    void Hao2(DComplex* psi, DComplex const* newPsi, DComplex const* V, DReal slow, DReal epsilon, size_t n)
    {
        Trait::Transform(Exec(),
            DevicePtrCast(psi),
            DevicePtrCast(psi + n),
            DevicePtrCast(newPsi),
            DevicePtrCast(V),
            DevicePtrCast(psi),
            op_Hao2<RealType>(slow, epsilon));
    }
};

