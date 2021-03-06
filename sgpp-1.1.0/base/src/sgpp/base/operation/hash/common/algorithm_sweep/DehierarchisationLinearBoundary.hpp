// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DEHIERARCHISATIONLINEARBOUNDARY_HPP
#define DEHIERARCHISATIONLINEARBOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {



/**
 * Class that implements the dehierarchisation on a linear sparse grid with boundaries. Therefore
 * the ()operator has to be implement in order to use the sweep algorithm for
 * the grid traversal
 */
class DehierarchisationLinearBoundary : public DehierarchisationLinear {
 public:
  /**
   * Constructor, must be bind to a grid
   *
   * @param storage the grid storage object of the the grid, on which the dehierarchisation should be executed
   */
  explicit DehierarchisationLinearBoundary(GridStorage* storage);

  /**
   * Destructor
   */
  ~DehierarchisationLinearBoundary() override;

  /**
   * Implements operator() needed by the sweep class during the grid traversal. This function
   * is applied to the whole grid.
   *
   * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
   * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
   * result)
   * So please assure that both functions do exist!
   *
   * @param source this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
   * @param result this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
   * @param index a iterator object of the grid
   * @param dim current fixed dimension of the 'execution direction'
   */
  virtual void operator()(DataVector& source, DataVector& result,
                          grid_iterator& index, size_t dim) override;
};

}  // namespace base
}  // namespace SGPP

#endif /* DEHIERARCHISATIONLINEARBOUNDARY_HPP */
