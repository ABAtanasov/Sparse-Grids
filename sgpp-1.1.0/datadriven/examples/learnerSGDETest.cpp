// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <string>

#include "sgpp/globaldef.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/datadriven/application/LearnerSGDE.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/datadriven/application/RegularizationConfiguration.hpp"
#include "sgpp/datadriven/application/GaussianKDE.hpp"

using std::cout;
using std::endl;

int main(int argc, char** argv) {
  std::string filename = "../tests/data/friedman_4d_2000.arff";

  cout << "# loading file: " << filename << endl;
  SGPP::datadriven::Dataset dataset =
      SGPP::datadriven::ARFFTools::readARFF(filename);
  SGPP::base::DataMatrix& samples = dataset.getData();

  // configure grid
  cout << "# create grid config" << endl;
  SGPP::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = dataset.getDimension();
  gridConfig.level_ = 4;
  gridConfig.type_ = SGPP::base::GridType::Linear;

  // configure adaptive refinement
  cout << "# create adaptive refinement config" << endl;
  SGPP::base::AdpativityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.noPoints_ = 10;

  // configure solver
  cout << "# create solver config" << endl;
  SGPP::solver::SLESolverConfiguration solverConfig;
  solverConfig.type_ = SGPP::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-10;
  solverConfig.threshold_ = 1e-10;

  // configure regularization
  cout << "# create regularization config" << endl;
  SGPP::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.regType_ = SGPP::datadriven::RegularizationType::Laplace;

  // configure learner
  cout << "# create learner config" << endl;
  SGPP::datadriven::LearnerSGDEConfiguration learnerConfig;
  learnerConfig.doCrossValidation_ = true;
  learnerConfig.kfold_ = 3;
  learnerConfig.lambdaStart_ = 1e-1;
  learnerConfig.lambdaEnd_ = 1e-10;
  learnerConfig.lambdaSteps_ = 3;
  learnerConfig.logScale_ = true;
  learnerConfig.shuffle_ = true;
  learnerConfig.seed_ = 1234567;
  learnerConfig.silent_ = false;

  cout << "# creating the learner" << endl;
  SGPP::datadriven::LearnerSGDE learner(gridConfig, adaptConfig, solverConfig,
                                        regularizationConfig, learnerConfig);
  learner.initialize(samples);

  SGPP::datadriven::GaussianKDE kde(samples);
  SGPP::base::DataVector x(learner.getDim());

  for (size_t i = 0; i < x.getSize(); i++) {
    x[i] = 0.5;
  }

  cout << "--------------------------------------------------------" << endl;
  cout << (learner.getGrid())->getSize() << " -> "
       << (learner.getAlpha())->sum() << endl;
  cout << "pdf_SGDE(x) = " << learner.pdf(x) << " ~ " << kde.pdf(x)
       << " = pdf_KDE(x)" << endl;
  cout << "mean_SGDE(x) = " << learner.mean() << " ~ " << kde.mean()
       << " = mean_KDE(x)" << endl;
  cout << "var_SGDE(x) = " << learner.variance() << " ~ " << kde.variance()
       << " = var_KDE(x)" << endl;
}
