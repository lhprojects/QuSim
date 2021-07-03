#include "Device.h"
#include "DeviceTemplate.h"

extern Device* CreateCudaDevice();

std::unique_ptr<Device> Device::Create(DeviceType type)
{
    Device* d = nullptr;

    if (type == DeviceType::CPU_SEQ) {
        d = new DeviceTemplate<std::execution::sequenced_policy>();
    } else if (type == DeviceType::CPU_PAR) {
        d = new DeviceTemplate<std::execution::parallel_policy>();
    } else if (type == DeviceType::CPU_PAR_VEC) {
        d = new DeviceTemplate<std::execution::parallel_unsequenced_policy>();
    } else if (type == DeviceType::CPU_VEC) {

#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201902L
        d = new DeviceTemplate<std::execution::unsequenced_policy>();
#else
        d = new DeviceTemplate<std::execution::parallel_unsequenced_policy>();
#endif

    } else if (type == DeviceType::GPU_CUDA) {

#if defined(QUSIM_USE_CUDA) && QUSIM_USE_CUDA

        d = CreateCudaDevice();
#else
        throw std::invalid_argument("Cuda not compiled");
#endif

    } else {
        throw std::invalid_argument("Unkown device type");
    }

    return std::unique_ptr<Device>(d);
}
