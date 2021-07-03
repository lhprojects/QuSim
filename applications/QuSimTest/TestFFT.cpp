#include "../../src/FourierTransform.h"
#include "../../src/Device.h"
#include "../QuTest/Test.h"


template<class CreateFFT>
void testFFT(CreateFFT createFFT, int n)
{
    auto cpu = Device::Create(DeviceType::CPU_SEQ);
    auto cuda = Device::Create(DeviceType::GPU_CUDA);

    auto cpu_o = cpu->Alloc<Complex>(n);
    for (int i = 0; i < n; ++i) {
        cpu_o[i] = i;
    }

    auto cpu_k1 = cpu->Alloc<Complex>(n);
    auto cpu_k2 = cpu->Alloc<Complex>(n);
    auto cpu_k3 = cpu->Alloc<Complex>(n);
    auto cuda_k3 = cuda->Alloc<Complex>(n);

    auto cpu_x1 = cpu->Alloc<Complex>(n);
    auto cpu_x2 = cpu->Alloc<Complex>(n);
    auto cpu_x3 = cpu->Alloc<Complex>(n);
    auto cuda_x3 = cuda->Alloc<Complex>(n);

    auto fft1 = createFFT(n, false, FourierTransformLibrary::KISS);
    auto fft2 = createFFT(n, false, FourierTransformLibrary::FFTW);
    auto fft3 = createFFT(n, false, FourierTransformLibrary::CUDA);

    auto invfft1 = createFFT(n, true, FourierTransformLibrary::KISS);
    auto invfft2 = createFFT(n, true, FourierTransformLibrary::FFTW);
    auto invfft3 = createFFT(n, true, FourierTransformLibrary::CUDA);

    fft1->Transform(cpu_o, cpu_k1);
    fft2->Transform(cpu_o, cpu_k2);

    cuda->ToDevice(cuda_x3, cpu_o, n);
    fft3->Transform(cuda_x3, cuda_k3);
    cuda->ToHost(cpu_k3, cuda_k3, n);

    auto e12 = cpu->MinusNorm2(cpu_k1, cpu_k2, n);
    auto e13 = cpu->MinusNorm2(cpu_k1, cpu_k3, n);
    auto e23 = cpu->MinusNorm2(cpu_k2, cpu_k3, n);

    invfft1->Transform(cpu_k1, cpu_x1);
    invfft2->Transform(cpu_k2, cpu_x2);

    invfft3->Transform(cuda_k3, cuda_x3);
    cuda->ToHost(cpu_x3, cuda_x3, n);

    cpu->Scale(cpu_x1, 1. / n, n);
    cpu->Scale(cpu_x2, 1. / n, n);
    cpu->Scale(cpu_x3, 1. / n, n);

    auto e12_x = cpu->MinusNorm2(cpu_x1, cpu_x2, n);
    auto e13_x = cpu->MinusNorm2(cpu_x1, cpu_x3, n);
    auto e23_x = cpu->MinusNorm2(cpu_x2, cpu_x3, n);

    auto e01_x = cpu->MinusNorm2(cpu_o, cpu_x1, n);
    auto e02_x = cpu->MinusNorm2(cpu_o, cpu_x2, n);
    auto e03_x = cpu->MinusNorm2(cpu_o, cpu_x3, n);

    QU_TEST_LESS(abs(e12), 1E-8);
    QU_TEST_LESS(abs(e23), 1E-8);
    QU_TEST_LESS(abs(e13), 1E-8);

    QU_TEST_LESS(abs(e12_x), 1E-8);
    QU_TEST_LESS(abs(e23_x), 1E-8);
    QU_TEST_LESS(abs(e13_x), 1E-8);

    QU_TEST_LESS(abs(e01_x), 1E-8);
    QU_TEST_LESS(abs(e02_x), 1E-8);
    QU_TEST_LESS(abs(e03_x), 1E-8);


    //printf("%8.1E %8.1E %8.1E\n", e12, e12, e23);
    //printf("%8.1E %8.1E %8.1E\n", e12_x, e12_x, e23_x);
}


void testFFT()
{
    testFFT([](int n, bool inv, FourierTransformLibrary lib) {
        return FourierTransform1D::Create(1234, inv, lib);
        }, 1234);
    testFFT([](int n, bool inv, FourierTransformLibrary lib) {
        return FourierTransform2D::Create(123, 234, inv, lib);
        }, 123*234);
    testFFT([](int n, bool inv, FourierTransformLibrary lib) {
        return FourierTransform3D::Create(12, 23, 34, inv, lib);
        }, 12 * 23 * 34);
}
