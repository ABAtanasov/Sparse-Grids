// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <fstream>
#include <memory>
#include <string>

#include "sgpp/base/exception/operation_exception.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/base/opencl/OCLDevice.hpp"
#include "sgpp/base/opencl/KernelSourceBuilderBase.hpp"

namespace SGPP {
namespace datadriven {
namespace StreamingModOCLMaskMultiPlatform {

template <typename T>
class KernelSourceBuilderMult : public base::KernelSourceBuilderBase<T> {
 private:
  std::shared_ptr<base::OCLDevice> device;

  json::Node &kernelConfiguration;

  size_t dims;

  size_t localWorkgroupSize;
  bool useLocalMemory;

 public:
  KernelSourceBuilderMult(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration,
                          size_t dims)
      : device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
    localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
  }

  std::string generateSource() {
    if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
      return this->reuseSource("streamingModOCLMaskMP_mult.cl");
    }

    std::stringstream sourceStream;

    sourceStream << "// platform: " << device->platformName << " device: " << device->deviceName
                 << std::endl
                 << std::endl;

    if (std::is_same<T, double>::value) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl
                   << std::endl;
    }

    sourceStream << "__kernel" << std::endl;
    sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))"
                 << std::endl;
    sourceStream << "void multOCLMask(__global const " << this->floatType() << "* ptrLevel,"
                 << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrIndex," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrMask,"
                 << std::endl;  // not needed for this kernel, but there for uniformity
    sourceStream << "           __global const " << this->floatType() << "* ptrOffset,"
                 << std::endl;  // not needed for this kernel, but there for uniformity
    sourceStream << "           __global const " << this->floatType() << "* ptrData," << std::endl;
    sourceStream << "           __global const " << this->floatType() << "* ptrAlpha," << std::endl;
    sourceStream << "           __global       " << this->floatType() << "* ptrResult,"
                 << std::endl;
    sourceStream << "           uint resultSize," << std::endl;
    sourceStream << "           uint start_grid," << std::endl;
    sourceStream << "           uint end_grid) " << std::endl;
    sourceStream << "{" << std::endl;
    sourceStream << "   int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << "   int localIdx = get_local_id(0);" << std::endl;
    sourceStream << std::endl;

    if (useLocalMemory) {
      sourceStream << "   __local " << this->floatType() << " locLevel["
                   << dims * localWorkgroupSize << "];" << std::endl;
      sourceStream << "   __local " << this->floatType() << " locIndex["
                   << dims * localWorkgroupSize << "];" << std::endl;
      sourceStream << "   __local " << this->floatType() << " locMask[" << dims * localWorkgroupSize
                   << "];" << std::endl;
      sourceStream << "   __local " << this->floatType() << " locOffset["
                   << dims * localWorkgroupSize << "];" << std::endl;
      sourceStream << "   __local " << this->floatType() << " locAlpha[" << localWorkgroupSize
                   << "];" << std::endl;
      sourceStream << std::endl;
    }

    sourceStream << "   " << this->floatType()
                 << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl
                 << std::endl;
    sourceStream << "   " << this->floatType() << " myResult = ptrResult[globalIdx];" << std::endl
                 << std::endl;
    sourceStream << "   // Create registers for the data" << std::endl;

    for (size_t d = 0; d < dims; d++) {
      sourceStream << " " << this->floatType() << " data_" << d
                   << " = ptrData[globalIdx+(resultSize*" << d << ")];" << std::endl;
    }

    sourceStream << std::endl;
    if (useLocalMemory) {
      sourceStream << "   // Iterate over all grid points (fast ones, with cache)" << std::endl;
      sourceStream << " uint chunkSizeGrid = end_grid - start_grid;" << std::endl;
      sourceStream << " uint fastChunkSizeGrid = (chunkSizeGrid / " << localWorkgroupSize << ") * "
                   << localWorkgroupSize << ";" << std::endl;
      sourceStream << " for(int j = start_grid; j < start_grid + fastChunkSizeGrid; j+="
                   << localWorkgroupSize << ")" << std::endl;
      sourceStream << "   {" << std::endl;

      for (size_t d = 0; d < dims; d++) {
        sourceStream << "     locLevel[(localIdx*" << dims << ")+" << d
                     << "] = ptrLevel[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
        sourceStream << "     locIndex[(localIdx*" << dims << ")+" << d
                     << "] = ptrIndex[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
        sourceStream << "     locMask[(localIdx*" << dims << ")+" << d
                     << "] = ptrMask[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
        sourceStream << "     locOffset[(localIdx*" << dims << ")+" << d
                     << "] = ptrOffset[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
      }

      sourceStream << "       locAlpha[localIdx] = ptrAlpha[j+localIdx];" << std::endl;
      sourceStream << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      sourceStream << std::endl;
      sourceStream << "       for(int k = 0; k < " << localWorkgroupSize << "; k++)" << std::endl;
      sourceStream << "       {" << std::endl;
      sourceStream << "           curSupport = locAlpha[k];" << std::endl
                   << std::endl;

      for (size_t d = 0; d < dims; d++) {
        sourceStream << "         eval = ((locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d
                     << "));" << std::endl;
        sourceStream << "         index_calc = eval - (locIndex[(k*" << dims << ")+" << d << "]);"
                     << std::endl;
        sourceStream << "         abs = as_" << this->floatType() << "(as_" << this->intType()
                     << "(index_calc) | as_" << this->intType() << "(locMask[(k*" << dims << ")+"
                     << d << "]));" << std::endl;
        sourceStream << "         last = locOffset[(k*" << dims << ")+" << d << "] + abs;"
                     << std::endl;
        sourceStream << "         localSupport = fmax(last, 0.0" << this->constSuffix() << ");"
                     << std::endl;
        sourceStream << "         curSupport *= localSupport;" << std::endl
                     << std::endl;
      }

      sourceStream << "           myResult += curSupport;" << std::endl;
      sourceStream << "       }" << std::endl;
      sourceStream << std::endl;
      sourceStream << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      sourceStream << "   }" << std::endl;
      sourceStream << std::endl;
      sourceStream << "   // Iterate over all grid points (slow ones, without cache)" << std::endl;
      sourceStream << " for(int m = start_grid + fastChunkSizeGrid; m < end_grid; m++)"
                   << std::endl;
      sourceStream << "   {" << std::endl;
      sourceStream << "       curSupport = ptrAlpha[m];" << std::endl
                   << std::endl;

      for (size_t d = 0; d < dims; d++) {
        sourceStream << "     eval = ((ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d
                     << "));" << std::endl;
        sourceStream << "     index_calc = eval - (ptrIndex[(m*" << dims << ")+" << d << "]);"
                     << std::endl;
        sourceStream << "     abs = as_" << this->floatType() << "(as_" << this->intType()
                     << "(index_calc) | as_" << this->intType() << "(ptrMask[(m*" << dims << ")+"
                     << d << "]));" << std::endl;
        sourceStream << "     last = ptrOffset[(m*" << dims << ")+" << d << "] + abs;" << std::endl;
        sourceStream << "     localSupport = fmax(last, 0.0" << this->constSuffix() << ");"
                     << std::endl;
        sourceStream << "     curSupport *= localSupport;" << std::endl
                     << std::endl;
      }

      sourceStream << "       myResult += curSupport;" << std::endl;
      sourceStream << "   }" << std::endl;
    } else {
      sourceStream << "   // Iterate over all grid points (without cache)" << std::endl;
      sourceStream << " for(int m = start_grid; m < end_grid; m++)" << std::endl;
      sourceStream << "   {" << std::endl;
      sourceStream << "       curSupport = ptrAlpha[m];" << std::endl
                   << std::endl;

      for (size_t d = 0; d < dims; d++) {
        sourceStream << "     eval = ((ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d
                     << "));" << std::endl;
        sourceStream << "     index_calc = eval - (ptrIndex[(m*" << dims << ")+" << d << "]);"
                     << std::endl;
        sourceStream << "     abs = as_" << this->floatType() << "(as_" << this->intType()
                     << "(index_calc) | as_" << this->intType() << "(ptrMask[(m*" << dims << ")+"
                     << d << "]));" << std::endl;
        sourceStream << "     last = ptrOffset[(m*" << dims << ")+" << d << "] + abs;" << std::endl;
        sourceStream << "     localSupport = fmax(last, 0.0" << this->constSuffix() << ");"
                     << std::endl;
        sourceStream << "     curSupport *= localSupport;" << std::endl
                     << std::endl;
      }

      sourceStream << "       myResult += curSupport;" << std::endl;
      sourceStream << "   }" << std::endl;
    }
    sourceStream << std::endl;
    sourceStream << "   ptrResult[globalIdx] = myResult;" << std::endl;
    sourceStream << "}" << std::endl;

    if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
      this->writeSource("streamingModOCLMaskMP_mult.cl", sourceStream.str());
    }

    return sourceStream.str();
  }
};

}  // namespace StreamingModOCLMaskMultiPlatform
}  // namespace datadriven
}  // namespace SGPP
