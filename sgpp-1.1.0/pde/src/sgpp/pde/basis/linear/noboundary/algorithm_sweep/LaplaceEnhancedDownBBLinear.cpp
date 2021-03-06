// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedDownBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

LaplaceEnhancedDownBBLinear::LaplaceEnhancedDownBBLinear(SGPP::base::GridStorage* storage)
    : storage(storage),
      boundingBox(storage->getBoundingBox()),
      algoDims(storage->getAlgorithmicDimensions()),
      numAlgoDims_(storage->getAlgorithmicDimensions().size()),
      ptr_source_(NULL),
      ptr_result_(NULL),
      cur_algo_dim_(0),
      q_(0.0),
      t_(0.0)
#if 1
#if defined(__SSE3__) && USE_DOUBLE_PRECISION == 1
      ,
      half_in_(_mm_set1_pd(0.5)),
      twothird_(_mm_set1_pd(2.0 / 3.0))
#endif
#else
#ifdef __SSE3__
      ,
      half_in_(_mm_set1_pd(0.5)),
      twothird_(_mm_set1_pd(2.0 / 3.0))
#endif
#endif
//      ,h_table_(NULL),grad_table_(NULL)
{
}

LaplaceEnhancedDownBBLinear::~LaplaceEnhancedDownBBLinear() {}

void LaplaceEnhancedDownBBLinear::operator()(SGPP::base::DataMatrix& source,
                                             SGPP::base::DataMatrix& result, grid_iterator& index,
                                             size_t dim) {
  q_ = this->boundingBox->getIntervalWidth(this->algoDims[dim]);
  t_ = this->boundingBox->getIntervalOffset(this->algoDims[dim]);

  //    h_table_ = new float_t[MAX_TABLE_DEPTH+1];
  //    grad_table_ = new float_t[MAX_TABLE_DEPTH+1];
  //
  //    for (int i = 0; i <= MAX_TABLE_DEPTH; i++)
  //    {
  //        h_table_[i] = 1.0/static_cast<float_t>(1<<i);
  //        grad_table_[i] = (static_cast<float_t>(1<<(i+1)))/q;
  //    }

  ptr_source_ = source.getPointer();
  ptr_result_ = result.getPointer();
  cur_algo_dim_ = this->algoDims[dim];

  if (q_ != 1.0 || t_ != 0.0) {
    size_t i = 0;

    for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
      float_t fl = 0.0;
      float_t fr = 0.0;

      if (dim == i) {
        recBB_GL(fl, fr, i, index);
      } else if (dim == i + 1) {
        recBB_LG(fl, fr, i, index);
      } else {
        recBB_LL(fl, fr, fl, fr, i, index);
      }
    }

    for (; i < this->numAlgoDims_; i++) {
      float_t fl = 0.0;
      float_t fr = 0.0;

      if (dim == i) {
        recBB_grad(i, index);
      } else {
        recBB(fl, fr, i, index);
      }
    }
  } else {
    size_t i = 0;

    for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
      float_t fl = 0.0;
      float_t fr = 0.0;
#if 1
#if defined(__SSE3__) && USE_DOUBLE_PRECISION == 1
      __m128d fl_xmm = _mm_set1_pd(0.0);
      __m128d fr_xmm = _mm_set1_pd(0.0);
#endif
#else
#ifdef __SSE3__
      __m128d fl_xmm = _mm_set1_pd(0.0);
      __m128d fr_xmm = _mm_set1_pd(0.0);
#endif
#endif

      if (dim == i) {
        rec_GL(fl, fr, i, index);
      } else if (dim == i + 1) {
        rec_LG(fl, fr, i, index);
      } else {
#if 1
#if defined(__SSE3__) && USE_DOUBLE_PRECISION == 1
        rec_LL(fl_xmm, fr_xmm, i, index);
#else
        rec_LL(fl, fr, fl, fr, i, index);
#endif
#else
#ifdef __SSE3__
        rec_LL(fl_xmm, fr_xmm, i, index);
#else
        rec_LL(fl, fr, fl, fr, i, index);
#endif
#endif
      }
    }

    for (; i < this->numAlgoDims_; i++) {
      float_t fl = 0.0;
      float_t fr = 0.0;

      if (dim == i) {
        rec_grad(i, index);
      } else {
        rec(fl, fr, i, index);
      }
    }
  }

  //    delete[] h_table_;
  //    delete[] grad_table_;
}

void LaplaceEnhancedDownBBLinear::rec(float_t fl, float_t fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // mesh-width
  float_t h = 1.0 / static_cast<float_t>(1 << l);

  // L2 scalar product
  float_t tmp_m = ((fl + fr) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)));
  float_t fm = tmp_m + alpha_value;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      rec(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      rec(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::rec_grad(size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // Gradient just in selected dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (static_cast<float_t>(1 << (l + 1)) * alpha_value);

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      rec_grad(dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      rec_grad(dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

#if 1
#if defined(__SSE3__) && USE_DOUBLE_PRECISION == 1
void LaplaceEnhancedDownBBLinear::rec_LL(__m128d fl, __m128d fr, size_t dim, grid_iterator& index)
#else
void LaplaceEnhancedDownBBLinear::rec_LL(float_t fl, float_t fr, float_t fl2, float_t fr2,
                                         size_t dim, grid_iterator& index)
#endif
#else
#ifdef __SSE3__
void LaplaceEnhancedDownBBLinear::rec_LL(__m128d fl, __m128d fr, size_t dim, grid_iterator& index)
#else
void LaplaceEnhancedDownBBLinear::rec_LL(float_t fl, float_t fr, float_t fl2, float_t fr2,
                                         size_t dim, grid_iterator& index)
#endif
#endif
{
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  // mesh-width
  float_t h = 1.0 / static_cast<float_t>(1 << l);

// L2 scalar product
#if 1
#if defined(__SSE3__) && USE_DOUBLE_PRECISION == 1
  // with intrinsics
  __m128d h_in = _mm_loaddup_pd(&h);
  __m128d fl_in = fl;
  __m128d fr_in = fr;
  __m128d alpha = _mm_loadu_pd(&ptr_source_[(seq * this->numAlgoDims_) + dim]);
  __m128d tmp = _mm_mul_pd(_mm_add_pd(fl_in, fr_in), half_in_);
  __m128d res = _mm_add_pd(_mm_mul_pd(h_in, tmp), _mm_mul_pd(alpha, _mm_mul_pd(h_in, twothird_)));
  __m128d new_fm = _mm_add_pd(alpha, tmp);
  _mm_storeu_pd(&ptr_result_[(seq * this->numAlgoDims_) + dim], res);
#else
  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  float_t alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];
  float_t tmp_m = ((fl + fr) / 2.0);
  float_t tmp_m2 = ((fl2 + fr2) / 2.0);
  float_t res = ((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value));
  float_t res2 = ((h * tmp_m2) + (((2.0 / 3.0) * h) * alpha_value2));
  ptr_result_[(seq * this->numAlgoDims_) + dim] = res;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = res2;
  float_t fm = tmp_m + alpha_value;
  float_t fm2 = tmp_m2 + alpha_value2;
#endif
#else
#ifdef __SSE3__
  // with intrinsics
  __m128d h_in = _mm_loaddup_pd(&h);
  __m128d fl_in = fl;
  __m128d fr_in = fr;
  __m128d alpha = _mm_loadu_pd(&ptr_source_[(seq * this->numAlgoDims_) + dim]);
  __m128d tmp = _mm_mul_pd(_mm_add_pd(fl_in, fr_in), half_in_);
  __m128d res = _mm_add_pd(_mm_mul_pd(h_in, tmp), _mm_mul_pd(alpha, _mm_mul_pd(h_in, twothird_)));
  __m128d new_fm = _mm_add_pd(alpha, tmp);
  _mm_storeu_pd(&ptr_result_[(seq * this->numAlgoDims_) + dim], res);
#else
  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  float_t alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];
  float_t tmp_m = ((fl + fr) / 2.0);
  float_t tmp_m2 = ((fl2 + fr2) / 2.0);
  float_t res = ((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value));
  float_t res2 = ((h * tmp_m2) + (((2.0 / 3.0) * h) * alpha_value2));
  ptr_result_[(seq * this->numAlgoDims_) + dim] = res;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = res2;
  float_t fm = tmp_m + alpha_value;
  float_t fm2 = tmp_m2 + alpha_value2;
#endif
#endif

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
#if 1
#if defined(__SSE3__) && USE_DOUBLE_PRECISION == 1
      rec_LL(fl, new_fm, dim, index);
#else
      rec_LL(fl, fm, fl2, fm2, dim, index);
#endif
#else
#ifdef __SSE3__
      rec_LL(fl, new_fm, dim, index);
#else
      rec_LL(fl, fm, fl2, fm2, dim, index);
#endif
#endif
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
#if 1
#if defined(__SSE3__) && USE_DOUBLE_PRECISION == 1
      rec_LL(new_fm, fr, dim, index);
#else
      rec_LL(fm, fr, fm2, fr2, dim, index);
#endif
#else
#ifdef __SSE3__
      rec_LL(new_fm, fr, dim, index);
#else
      rec_LL(fm, fr, fm2, fr2, dim, index);
#endif
#endif
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::rec_LG(float_t fl, float_t fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  float_t alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  float_t h = 1.0 / static_cast<float_t>(1 << l);

  // L2 scalar product
  float_t tmp_m = ((fl + fr) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)));
  // Gradient in second dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      (static_cast<float_t>(1 << (l + 1)) * alpha_value2);

  float_t fm = tmp_m + alpha_value;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      rec_LG(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      rec_LG(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::rec_GL(float_t fl, float_t fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  float_t alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  float_t h = 1.0 / static_cast<float_t>(1 << l);

  float_t tmp_m = ((fl + fr) * 0.5);

  // Gradient in second dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      ((static_cast<float_t>(1 << (l + 1))) * alpha_value);
  // L2 scalar product
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value2)));

  float_t fm = tmp_m + alpha_value2;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      rec_GL(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      rec_GL(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB(float_t fl, float_t fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // mesh-width
  float_t h = 1.0 / static_cast<float_t>(1 << l);

  // L2 scalar product
  float_t tmp_m = ((fl + fr) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)) * q_);
  float_t fm = tmp_m + alpha_value;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB_LL(float_t fl, float_t fr, float_t fl2, float_t fr2,
                                           size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  float_t alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  float_t h = 1.0 / static_cast<float_t>(1 << l);

  // L2 scalar product
  float_t tmp_m = ((fl + fr) * 0.5);
  float_t tmp_m2 = ((fl2 + fr2) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)) * q_);
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      (((h * tmp_m2) + (((2.0 / 3.0) * h) * alpha_value2)) * q_);
  float_t fm = tmp_m + alpha_value;
  float_t fm2 = tmp_m2 + alpha_value2;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB_LL(fl, fm, fl2, fm2, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB_LL(fm, fr, fm2, fr2, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB_LG(float_t fl, float_t fr, size_t dim,
                                           grid_iterator& index) {
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  float_t alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  float_t h = 1.0 / static_cast<float_t>(1 << l);

  // L2 scalar product
  float_t tmp_m = ((fl + fr) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)) * q_);
  // Gradient in second dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      ((static_cast<float_t>(1 << (l + 1)) / q_) * alpha_value2);

  float_t fm = tmp_m + alpha_value;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB_LG(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB_LG(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB_GL(float_t fl, float_t fr, size_t dim,
                                           grid_iterator& index) {
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  float_t alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  float_t h = 1.0 / static_cast<float_t>(1 << l);

  float_t tmp_m = ((fl + fr) * 0.5);
  // Gradient in second dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      ((static_cast<float_t>(1 << (l + 1)) / q_) * alpha_value);
  // L2 scalar product
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value2)) * q_);
  float_t fm = tmp_m + alpha_value2;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB_GL(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB_GL(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB_grad(size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;
  index.get(cur_algo_dim_, l, i);

  float_t alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // Gradient just in selected dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      ((static_cast<float_t>(1 << (l + 1)) / q_) * alpha_value);

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB_grad(dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->end(index.seq())) {
      recBB_grad(dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

}  // namespace pde
}  // namespace SGPP
