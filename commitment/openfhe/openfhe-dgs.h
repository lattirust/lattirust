#pragma once
#include "rust/cxx.h"
#include "openfhe/core/math/discretegaussiangeneratorgeneric.h"
#include "openfhe/core/math/discretegaussiangenerator-impl.h"
#include "openfhe/core/math/discretegaussiangenerator.h"

struct GaussianParameters;

// std::unique_ptr<DiscreteGaussianGenerator> new_gaussian_sampler(const GaussianParameters& params);