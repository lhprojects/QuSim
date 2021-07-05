#ifndef QUSIM_DEVICE
#define QUSIM_DEVICE

#include <complex>
#include <memory>

enum class DeviceType {
    CPU_SEQ,
    CPU_PAR,
    CPU_VEC,
    CPU_PAR_VEC,
    GPU_CUDA,
};

struct Device;
struct DeviceAutoFree {

    DeviceAutoFree(Device &dev, void *mem) : fDev(dev), fMem(mem) { }
    ~DeviceAutoFree();
    
private:
    Device& fDev;
    void* const fMem;
};

template<class T>
struct DeviceAutoFreeT {

    DeviceAutoFreeT(Device& dev, T* mem) : fDev(dev), fMem(mem) {}
    ~DeviceAutoFreeT();
    T* Get() const { return fMem; }
private:
    Device& fDev;
    T* const fMem;
};

struct Device
{

    using DReal = double;
    using DComplex = std::complex<DReal>;
    
    using DReal64 = double;
    using DComplex64 = std::complex<DReal>;
    
    using DReal32 = float;
    using DComplex32 = std::complex<float>;


    static std::unique_ptr<Device> Create(DeviceType);
    virtual ~Device() {}

    virtual bool OnMainMem() = 0;
    virtual void* AllocMem(size_t bytes) = 0;
    virtual void Free(void*) = 0;

    template<class T>
    void SafeFree(T * &p);

    template<class T>
    T* Alloc(size_t bytes)
    {
        return (T*)AllocMem(sizeof(T) * bytes);
    }

    DeviceAutoFree AutoFree(void* f);

    template<class T>
    DeviceAutoFreeT<T> AllocAutoFree(size_t bytes)
    {
        return DeviceAutoFreeT<T>(*this, (T*)AllocMem(sizeof(T) * bytes));
    }

    void ToDevice(DComplex* out, DComplex const* in, size_t n);
    void ToDevice(DReal* out, DReal const* in, size_t n);
    virtual void MemToDevice(void* out, void const* in, size_t bytes) = 0;

    void ToHost(DComplex *out, DComplex const* in, size_t n);
    void ToHost(DReal* out, DReal const* in, size_t n);
    virtual void MemToHost(void* out, void const* in, size_t bytes) = 0;

    virtual void SetZero(DComplex* in, size_t n) = 0;
    virtual void SetOne(DComplex* in, size_t n) = 0;

    // sum_i |in_i|^2 i
    virtual DReal Abs2Idx(DComplex const* in, size_t n) = 0;
    virtual DReal Abs2Idx(DReal const* in, size_t n) = 0;

    virtual void Set(DComplex* in, size_t idx, DComplex alpha) = 0;
    virtual void Set(DReal* in, size_t idx, DReal alpha) = 0;

    virtual DComplex At(DComplex const* in, size_t idx, size_t n) = 0;
    virtual DReal At(DReal const* in, size_t idx, size_t n) = 0;

    virtual DReal Min(DComplex const* in, size_t n) = 0;
    virtual DReal Max(DComplex const* in, size_t n) = 0;

    virtual DReal SumReal(DComplex const* in, size_t n) = 0;
    virtual DReal SumImag(DComplex const* in, size_t n) = 0;

    virtual DReal MinusNorm2(DComplex const* in1, DComplex const* in2, size_t n) = 0;
    virtual DReal Norm2(DReal const* in, size_t nx) = 0;
    virtual DReal Norm2(DComplex const* in, size_t nx) = 0;
    virtual void Exp(DComplex* to, DComplex const* from, DComplex alpha, size_t nx) = 0;
    virtual void Exp(DComplex* to, DReal const* from, DComplex alpha, size_t nx) = 0;
    virtual void ExpI(DComplex* to, DReal const* from, DReal alpha, size_t nx) = 0;

    virtual void Copy(DComplex *to, DComplex const* from, size_t n) = 0;
    virtual void Copy(DReal* to, DReal const* from, size_t n) = 0;

    virtual void Scale(double* fromto, double scale, size_t nx) = 0;
    virtual void Scale(DComplex* fromto, double scale, size_t nx) = 0;
    virtual void Scale(DComplex* fromto, DComplex scale, size_t nx) = 0;

    virtual void Scale(DReal* to, DReal const* from, DReal scale, size_t nx) = 0;
    virtual void Scale(DComplex* to, DComplex const* from, double scale, size_t nx) = 0;
    virtual void Scale(DComplex* to, DComplex const* from, DComplex scale, size_t nx) = 0;


    virtual void Sub(DComplex* to, DComplex const* in1, DComplex const* in2, size_t n) = 0;
    virtual void Mul(DComplex* to, DComplex const* from1, DComplex const* from2, size_t nx) = 0;
    virtual void MulR(DComplex* to, DComplex const* from1, DComplex const* from2, size_t nx) = 0;
    virtual void Mul(DComplex* to, DReal const* from1, DComplex const* from2, size_t nx) = 0;
    virtual void Mul(DComplex* to, DReal const* from1, DComplex alpha, DComplex const* from2, size_t nx) = 0;
    virtual void Mul(DComplex* to, DComplex const* from2, size_t nx) = 0;

    virtual void MulExp(DComplex* to, DComplex const* from1, DComplex alpha, DComplex const* from2, size_t nx) = 0;
    virtual void MulExp(DComplex* fromto, DComplex alpha, DComplex const* from2, size_t nx) = 0;
    virtual void MulExp(DComplex* fromto, DComplex alpha, DReal const* from2, size_t nx) = 0;
    virtual void MulExpI(DComplex* fromto, DReal alpha, DReal const* from2, size_t nx) = 0;

    virtual void CopyReverseMinu(DComplex* to, DComplex const *from, size_t nx) = 0;

    // sum_i |psi_i|^2 * i^2 alpha

    virtual DReal Abs2Idx2(DComplex const* psi, size_t nx) = 0;
    virtual void MulExpIdx2(DComplex * psi, DComplex alpha, size_t nx) = 0;

    virtual DReal Abs2K1D(DComplex const* psi, DReal ax, size_t nx) = 0;
    virtual DReal Abs2K2D(DComplex const* psi, DReal ax, size_t nx, DReal ay, size_t ny) = 0;
    virtual DReal Abs2K3D(DComplex const* psi, DReal ax, size_t nx, DReal ay, size_t ny, DReal az, size_t nz) = 0;

    // psi *= psi * exp(I * ax * idx^2)
    virtual void MulExpK1D(DComplex* psi, DComplex ax, size_t nx) = 0;
    virtual void MulExpK2D(DComplex* psi, DComplex alpha, DReal ax, size_t nx, DReal ay, size_t ny) = 0;
    virtual void MulExpK3D(DComplex* psi, DComplex alpha, DReal ax, size_t nx, DReal ay, size_t ny, DReal az, size_t nz) = 0;

    // sum_i in1_i^2 *in2_i
    virtual DComplex Abs2Mul(DComplex const* in1, DComplex const* in2, size_t nx) = 0;
    virtual DReal Abs2Mul(DComplex const* in1, DReal const* in2, size_t nx) = 0;

    virtual void G01D(DComplex *psik, DReal E, DReal epsilon, DReal alpha, size_t n) = 0;
    virtual void G02D(DComplex* psik, DReal E, DReal epsilon, DReal a1, DReal a2, size_t nx, size_t ny) = 0;
    virtual void G03D(DComplex* psik, DReal E, DReal epsilon, DReal a1, DReal a2, DReal a3, size_t nx, size_t ny, size_t nz) = 0;

    virtual DComplex Xsection3D(DComplex const* psi0,
        DComplex const* psi,
        DComplex const* v,
        DReal kx, DReal ky, DReal kz,
        DReal x0, DReal y0, DReal z0,
        DReal dx, DReal dy, DReal dz,
        size_t nx, size_t ny, size_t nz) = 0;

    virtual void AddMul(DComplex* psi, DComplex const* psi0, DComplex const* v, size_t n) = 0;

    virtual void Add3(DComplex* deltaPsix, DComplex const* psi0x, DComplex const* V, DReal epsilon, size_t n) = 0;
    virtual void Add3(DComplex* add3, DComplex const* deltaPsix, DComplex const* psi0x, DComplex const* V, DReal epsilon, size_t n) = 0;

    virtual void LinearUpdate(DComplex* psi, DComplex const* newPsi, DReal gamma, size_t n) = 0;
    virtual void Vellekoop(DComplex* psi, DComplex const* newPsi, DComplex const* V, DReal slow, DReal epsilon, size_t n) = 0;
    virtual void Hao1(DComplex* psi, DComplex const* newPsi, DComplex const* V, DReal slow, DReal epsilon, size_t n) = 0;
    virtual void Hao2(DComplex* psi, DComplex const* newPsi, DComplex const* V, DReal slow, DReal epsilon, size_t n) = 0;
};

inline DeviceAutoFree Device::AutoFree(void* f)
{
    return DeviceAutoFree(*this, f);
}

inline void Device::ToDevice(Device::DComplex* out, Device::DComplex const* in, size_t n)
{
    MemToDevice(out, in, n * sizeof(Device::DComplex));
}

inline void Device::ToDevice(Device::DReal* out, Device::DReal const* in, size_t n)
{
    MemToDevice(out, in, n * sizeof(Device::DReal));
}

inline void Device::ToHost(DComplex* out, DComplex const* in, size_t n)
{
    MemToHost(out, in, sizeof(DComplex) * n);
}

inline void Device::ToHost(DReal* out, DReal const* in, size_t n)
{
    MemToHost(out, in, sizeof(DComplex) * n);
}


inline DeviceAutoFree::~DeviceAutoFree()
{
    if(fMem) fDev.Free(fMem);
}

template<class T>
inline DeviceAutoFreeT<T>::~DeviceAutoFreeT() { 
    if (fMem) fDev.Free(fMem);
}

template<class T>
void Device::SafeFree(T*& p)
{
    if (p) {
        Free(p);
        p = nullptr;
    }
}

#endif
