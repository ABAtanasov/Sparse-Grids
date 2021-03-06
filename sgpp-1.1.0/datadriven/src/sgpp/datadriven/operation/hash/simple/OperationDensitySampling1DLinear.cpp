// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationDensitySampling1DLinear.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>
#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>

namespace SGPP {
namespace datadriven {
/**
 * WARNING: the grid must be a 1D grid!
 */
OperationDensitySampling1DLinear::OperationDensitySampling1DLinear(base::Grid* grid) {
  this->grid = grid;
}

OperationDensitySampling1DLinear::~OperationDensitySampling1DLinear() {}

void OperationDensitySampling1DLinear::doSampling1D(base::DataVector* alpha, size_t num_samples,
                                                    base::DataVector*& samples,
                                                    unsigned int* seedp) {
  /***************** STEP 1. Compute CDF  ********************/

  // compute PDF, sort by coordinates
  std::multimap<float_t, float_t> coord_pdf, coord_cdf;
  std::multimap<float_t, float_t>::iterator it1, it2;

  base::GridStorage* gs = this->grid->getStorage();
  base::OperationEval* opEval = op_factory::createOperationEval(*(this->grid));
  base::DataVector coord(1);

  for (unsigned int i = 0; i < gs->size(); i++) {
    coord[0] = gs->get(i)->getCoord(0);
    coord_pdf.insert(std::pair<float_t, float_t>(coord[0], opEval->eval(*alpha, coord)));
    coord_cdf.insert(std::pair<float_t, float_t>(coord[0], i));
  }

  delete opEval;
  opEval = NULL;
  // include values at the boundary [0,1]
  coord_pdf.insert(std::pair<float_t, float_t>(0.0, 0.0));
  coord_pdf.insert(std::pair<float_t, float_t>(1.0, 0.0));
  coord_cdf.insert(std::pair<float_t, float_t>(0.0, 0.0));
  coord_cdf.insert(std::pair<float_t, float_t>(1.0, 1.0));

  // Composite rule: trapezoidal (b-a)/2 * (f(a)+f(b))
  it1 = coord_pdf.begin();
  it2 = coord_pdf.begin();
  std::vector<float_t> tmp;
  tmp.push_back(0.0);
  float_t sum = 0.0, area;

  for (++it2; it2 != coord_pdf.end(); ++it2) {
    // (*it).first : the coordinate
    // (*it).second : the function value
    area = ((*it2).first - (*it1).first) / 2 * ((*it1).second + (*it2).second);
    tmp.push_back(area);
    sum += area;
    ++it1;
  }

  // compute CDF
  float_t tmp_sum;
  unsigned int i = 0;

  for (it1 = coord_cdf.begin(); it1 != coord_cdf.end(); ++it1) {
    tmp_sum = 0.0;

    for (unsigned int j = 0; j <= i; ++j) tmp_sum += tmp[j];

    ++i;
    (*it1).second = tmp_sum / sum;
  }

  tmp.clear();
  coord_pdf.clear();
  /***************** STEP 1. Done  ********************/

  /***************** STEP 2. Sampling  ********************/
  samples = new base::DataVector(num_samples);

  float_t y, x1, x2, y1, y2;

  for (size_t i = 0; i < num_samples; i++) {
// drand48_r() is thread_safe
#ifdef _WIN32
    y = static_cast<float_t>(rand()) / RAND_MAX;
#else
    y = static_cast<float_t>(rand_r(seedp)) / RAND_MAX;
#endif

    // find cdf interval
    for (it1 = coord_cdf.begin(); it1 != coord_cdf.end(); ++it1)
      if ((*it1).second >= y) break;

    x2 = (*it1).first;
    y2 = (*it1).second;
    --it1;
    x1 = (*it1).first;
    y1 = (*it1).second;
    // find x (linear interpolation): (y-y1)/(x-x1) = (y2-y1)/(x2-x1)
    (*samples)[i] = (x2 - x1) / (y2 - y1) * (y - y1) + x1;
  }

  /***************** STEP 2. Done  ********************/

  coord_cdf.clear();
  return;
}  // end of compute_1D_cdf()
}  // namespace datadriven
}  // namespace SGPP
