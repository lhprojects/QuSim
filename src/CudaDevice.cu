#if defined(QUSIM_USE_CUDA) && QUSIM_USE_CUDA

#include <cuda_runtime.h>
#include <thrust/transform.h>
#include <thrust/device_ptr.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/reverse.h>
#include <thrust/fill.h>
#include <thrust/transform_reduce.h>
#include <thrust/zip_function.h>
#include <thrust/complex.h>

#define QUSIM_HOST_DEVICE __host__ __device__

template<class R>
static inline QUSIM_HOST_DEVICE R Abs2(thrust::complex<R> s)
{
    return s.imag() * s.imag() + s.real() * s.real();
}

#include "DeviceTemplate.h"

struct CudaExec { };

inline std::string my_itoa(int i)
{
    char b[10];
    sprintf(b, "%d", i);
    return b;
}
#define check_err(err) do { if(err!=cudaSuccess) throw std::runtime_error("cuda error: " #err " : " + my_itoa((int)err)); } while(0)

template<>
struct ExecTrait<CudaExec>
{
    using Exec = std::execution::parallel_unsequenced_policy;
    static const bool cOnMainMemory = false;

    template<class T>
    using IntIter = thrust::counting_iterator<T>;

    template<class T>
    using ReverseIter = thrust::reverse_iterator<T>;

    template<class T>
    using DevicePtr = thrust::device_ptr<T>;

    using ComplexType = thrust::complex<double>;
    using RealType = double;

    template<class T>
    struct ScalarTrait {
    };

    template<>
    struct ScalarTrait<double> {
        using DeviceType = double;
        using ComplexType = thrust::complex<double>;
        using RealType = double;
        static RealType real_max()
        {
            return std::numeric_limits<double>::max();
        }
    };
    template<>
    struct ScalarTrait<const double> {
        using DeviceType = const double;
        using ComplexType = const thrust::complex<double>;
        using RealType = const double;
        static RealType real_max()
        {
            return std::numeric_limits<double>::max();
        }
    };

    template<>
    struct ScalarTrait<thrust::complex<double> > {
        using DeviceType = ComplexType;
        using ComplexType = thrust::complex<double>;
        using RealType = double;
        static RealType real_max()
        {
            return std::numeric_limits<double>::max();
        }
    };

    template<>
    struct ScalarTrait<const thrust::complex<double> > {
        using DeviceType = const ComplexType;
        using ComplexType = const thrust::complex<double>;
        using RealType = const double;
        static RealType real_max()
        {
            return std::numeric_limits<double>::max();
        }
    };

    static void Free(void* ptr)
    {
        check_err(cudaFree(ptr));
    }
    static void* Alloc(size_t bytes)
    {
        void* ptr;
        check_err(cudaMalloc(&ptr, bytes));
        return (void*)ptr;
    }

    static void MemToDevice(Exec& exec, void* out, void const* in, size_t bytes)
    {
        check_err(cudaMemcpy(out, in, bytes, cudaMemcpyDefault));
    }

    static void MemToHost(Exec& exec, void* out, void const* in, size_t bytes)
    {
        check_err(cudaMemcpy(out, in, bytes, cudaMemcpyDefault));
    }

    static void MoveToDevice(void* out, void const* in, size_t bytes)
    {
        check_err(cudaMemcpy(out, in, bytes, cudaMemcpyDefault));
    }

    static void MoveToHost(void* out, void const* in, size_t bytes)
    {
        check_err(cudaMemcpy(out, in, bytes, cudaMemcpyDefault));
    }

    template<class Exec, class Iter1, class OutIter, class UnitaryOp>
    static void Transform(Exec& exec, Iter1 b1, Iter1 e1, OutIter b2, UnitaryOp unitaryOp)
    {
        thrust::transform(b1, e1, b2, unitaryOp);
    }

    template<class Exec, class Iter1, class Iter2, class OutIter, class BinaryOp>
    static void Transform(Exec& exec, Iter1 b1, Iter1 e1, Iter2 b2, OutIter outIter, BinaryOp binaryOp)
    {
#if 1
        thrust::transform(b1, e1, b2, outIter, binaryOp);
#endif
    }


    template<class TransOp>
    struct op_Expand {

        TransOp transOp;

        QUSIM_HOST_DEVICE op_Expand(TransOp transOp) : transOp(transOp) {
        }

        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(U u, V v) {
            return transOp(u, thrust::get<0>(v), thrust::get<1>(v));
        }
        
    };

    template<class Exec, class Iter1, class Iter2, class Iter3, class OutIter, class TransOp>
    static void Transform(Exec& exec, Iter1 b1, Iter1 e1, Iter2 b2, Iter3 b3, OutIter outIter, TransOp transOp)
    {
#if 1
        thrust::transform(b1, e1, thrust::make_zip_iterator(thrust::make_tuple(b2, b3)),
        outIter, op_Expand<TransOp>(transOp));
#endif
    }

    template<class Iter1, class T>
    static void Fill(Exec& exec, Iter1 b1, Iter1 e1, T const& v)
    {
#if 1
        thrust::fill(b1, e1, v);
#endif
    }

    template<class Iter1, class OutIter>
    static void Copy(Exec& exec, Iter1 b1, Iter1 e1, OutIter outIter)
    {
#if 1
        thrust::copy(b1, e1, outIter);
#endif
    }

    template<class Iter, class T, class BinaryOp, class UnitaryOp>
    static T TransformReduce(Exec& exec, Iter b1, Iter e1, T init, BinaryOp binaryOp, UnitaryOp unitaryOp)
    {
#if 1
        return thrust::transform_reduce(b1, e1, unitaryOp, init, binaryOp);
#else
        return T();
#endif
    }

    template<class V>
    QUSIM_HOST_DEVICE static inline V value(V v) {
        return v;
    }

    template<class V>
    QUSIM_HOST_DEVICE static inline V value(thrust::device_reference<V> v) {
        return v;
    }

    template<class BinaryOp>
    struct op_Expand2 {

        BinaryOp binaryOp;

        QUSIM_HOST_DEVICE op_Expand2(BinaryOp binaryOp) : binaryOp(binaryOp) {
        }

        template<class U, class V>
        QUSIM_HOST_DEVICE auto operator()(thrust::tuple<U,V> v) {
            auto v1 = value(thrust::get<0>(v));
            auto v2 = value(thrust::get<1>(v));
            return binaryOp(v1, v2);
        }        
    };

    template<class TriOp>
    struct op_Expand3 {

        TriOp triOp;

        QUSIM_HOST_DEVICE op_Expand3(TriOp triOp) : triOp(triOp) {
        }

        template<class U, class V, class P>
        QUSIM_HOST_DEVICE auto operator()(thrust::tuple<U,V,P> v) {
            auto v1 = value(thrust::get<0>(v));
            auto v2 = value(thrust::get<1>(v));
            auto v3 = value(thrust::get<2>(v));
            return triOp(v1, v2, v3);
        }        
    };
    
    template<class QuadOp>
    struct op_Expand4 {

        QuadOp quadOp;

        QUSIM_HOST_DEVICE op_Expand4(QuadOp quadOp) : quadOp(quadOp) {
        }

        template<class U, class V, class P, class R>
        QUSIM_HOST_DEVICE auto operator()(thrust::tuple<U,V,P,R> v) {
            auto v1 = value(thrust::get<0>(v));
            auto v2 = value(thrust::get<1>(v));
            auto v3 = value(thrust::get<2>(v));
            auto v4 = value(thrust::get<3>(v));
            return quadOp(v1, v2, v3, v4);
        }
    };

    template<class Iter, class Iter2, class T, class BinaryOp, class TransOp>
    static T TransformReduce(Exec& exec, Iter b1, Iter e1, Iter2 b2, T init, BinaryOp binaryOp, TransOp trans)
    {
#if 1
        auto B1 = thrust::make_zip_iterator(thrust::make_tuple(b1, b2));
        auto E1 = thrust::make_zip_iterator(thrust::make_tuple(e1, b2 + (e1 - b1)));

        using V1 = typename thrust::iterator_traits<Iter>::value_type;
        using V2 = typename thrust::iterator_traits<Iter2>::value_type;

        return thrust::transform_reduce(B1, E1, op_Expand2<TransOp>(trans), init, binaryOp);
#else
        return T();
#endif
    }

    template<class Iter, class Iter2, class Iter3, class T, class ReduceOp, class TransOp>
    static T TransformReduce(Exec& exec, Iter b1, Iter e1, Iter2 b2, Iter3 b3,
        T init, ReduceOp reduceOp, TransOp trans)
    {
#if 1
        auto B1 = thrust::make_zip_iterator(thrust::make_tuple(b1, b2, b3));
        auto E1 = thrust::make_zip_iterator(thrust::make_tuple(e1, b2 + (e1 - b1), b3 + (e1 - b1)));

        using V1 = typename thrust::iterator_traits<Iter>::value_type;
        using V2 = typename thrust::iterator_traits<Iter2>::value_type;

        return thrust::transform_reduce(B1, E1, op_Expand3<TransOp>(trans), init, reduceOp);
#else
        return T();
#endif
    }
    template<class Iter, class Iter2, class Iter3, class Iter4, class T, class ReduceOp, class TransOp>
    static T TransformReduce(Exec& exec, Iter b1, Iter e1, Iter2 b2, Iter3 b3, Iter4 b4,
        T init, ReduceOp reduceOp, TransOp trans)
    {
#if 1
        auto B1 = thrust::make_zip_iterator(thrust::make_tuple(b1, b2, b3, b4));
        auto E1 = thrust::make_zip_iterator(thrust::make_tuple(e1, b2 + (e1 - b1), b3 + (e1 - b1), b4 + (e1 - b1)));

        using V1 = typename thrust::iterator_traits<Iter>::value_type;
        using V2 = typename thrust::iterator_traits<Iter2>::value_type;

        return thrust::transform_reduce(B1, E1, op_Expand4<TransOp>(trans), init, reduceOp);
#else
        return T();
#endif
    }

    static inline QUSIM_HOST_DEVICE thrust::complex<double> IMul(double c)
    {
        return thrust::complex<double>(double(0), -c);
    }

    static inline QUSIM_HOST_DEVICE thrust::complex<double> IMul(thrust::complex<double> c)
    {
        return thrust::complex<double>(-c.real(), c.real());
    }

    static inline QUSIM_HOST_DEVICE double Abs2(double c)
    {
        return c*c;
    }

    static inline QUSIM_HOST_DEVICE double Abs2(thrust::complex<double> c)
    {
        return c.real() * c.real() + c.imag() * c.imag();
    }

    static inline QUSIM_HOST_DEVICE double GetReal(thrust::complex<double> c)
    {
        return c.real();
    }

    static inline QUSIM_HOST_DEVICE double GetImag(thrust::complex<double> c)
    {
        return c.imag();
    }

};

template class DeviceTemplate<CudaExec>;

Device* CreateCudaDevice()
{
    return new DeviceTemplate<CudaExec>();
}


#endif
