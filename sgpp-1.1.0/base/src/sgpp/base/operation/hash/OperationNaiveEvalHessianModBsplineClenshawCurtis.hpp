// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALHESSIANMODBSPLINECLENSHAWCURTIS_HPP
#define OPERATIONEVALHESSIANMODBSPLINECLENSHAWCURTIS_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessian.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace SGPP {
namespace base {

/**
 * Operation for evaluating modified Clenshaw-Curtis B-spline
 * linear combinations on Noboundary grids, their gradients,
 * and their Hessians.
 */
class OperationNaiveEvalHessianModBsplineClenshawCurtis :
  public OperationNaiveEvalHessian {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationNaiveEvalHessianModBsplineClenshawCurtis(
    GridStorage* storage, size_t degree) :
    storage(storage), base(degree) {
  }

  /**
   * Destructor.
   */
  ~OperationNaiveEvalHessianModBsplineClenshawCurtis() override {
  }

  /**
   * @param       alpha       coefficient vector
   * @param       point       evaluation point
   * @param[out]  gradient    gradient vector of linear combination
   * @param[out]  hessian     Hessian matrix of linear combination
   * @return                  value of linear combination
   */
  float_t evalHessian(const DataVector& alpha,
                      const DataVector& point,
                      DataVector& gradient,
                      DataMatrix& hessian) override;

 protected:
  /// storage of the sparse grid
  GridStorage* storage;
  /// 1D B-spline basis
  SBsplineModifiedClenshawCurtisBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONEVALHESSIANMODBSPLINECLENSHAWCURTIS_HPP */
