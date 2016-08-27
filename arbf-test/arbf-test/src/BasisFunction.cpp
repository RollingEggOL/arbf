//
// Created by Ke Liu on 8/5/16.
//

#include <cmath>
#include "../include/Common.h"
#include "../include/Config.h"
#include "../include/BasisFunction.h"
// BasisFunction implementations
BasisFunction::~BasisFunction() {}


// MQBasisFunction implementations
MQBasisFunction::~MQBasisFunction() {}

double MQBasisFunction::phi(double radial, double c) const {
    return std::sqrt(SQ(radial) + SQ(c));
}


// IMQBasisFunction implementations
IMQBasisFunction::~IMQBasisFunction() {}

double IMQBasisFunction::phi(double radial, double c) const {
    return 1.0 / (std::sqrt(SQ(radial) + SQ(c)));
}


// GaussianBasisFunction implementations
GaussianBasisFunction::~GaussianBasisFunction() {}

double GaussianBasisFunction::phi(double radial, double c) const {
    return std::exp(-SQ(c*radial));
}


// TPSBasisFunction implementations
TPSBasisFunction::~TPSBasisFunction() {}

double TPSBasisFunction::phi(double radial, double c) const {
    if (std::abs(radial) < Config::epsilon) {
        return 0;
    } else {
        return SQ(radial) * std::log(radial);
    }
}
