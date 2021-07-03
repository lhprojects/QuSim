#pragma once

#include "Device.h"
#include "OptionsImpl.h"

inline DeviceType GetDeviceType(OptionsImpl const &opts)
{
    DeviceType devType;
    std::string device = opts.GetDefaut("device", std::string("CPU_SEQ"));
    if (device == "GPU_CUDA") {
        devType = DeviceType::GPU_CUDA;
    } else if (device == "CPU_SEQ") {
        devType = DeviceType::CPU_SEQ;
    } else if (device == "CPU_PAR") {
        devType = DeviceType::CPU_PAR;
    } else if (device == "CPU_PAR_VEC") {
        devType = DeviceType::CPU_PAR_VEC;
    } else if (device == "CPU_VEC") {
        devType = DeviceType::CPU_VEC;
    } else {
        throw std::invalid_argument("invalid device: " + device);
    }
    return devType;
}
