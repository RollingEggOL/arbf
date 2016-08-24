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

double MQBasisFunction::phi(double radial) const {
    return std::sqrt(SQ(radial) + SQ(0.01));
}


// IMQBasisFunction implementations
IMQBasisFunction::~IMQBasisFunction() {}

double IMQBasisFunction::phi(double radial) const {
    return 1.0 / (std::sqrt(SQ(radial) + SQ(1.8)));
}


// GaussianBasisFunction implementations
GaussianBasisFunction::~GaussianBasisFunction() {}

double GaussianBasisFunction::phi(double radial) const {
    double c = 0.01;
    return std::exp(-SQ(c*radial));
}


// TPSBasisFunction implementations
TPSBasisFunction::~TPSBasisFunction() {}

double TPSBasisFunction::phi(double radial) const {
    if (std::abs(radial) < Config::epsilon) {
        return 0;
    } else {
        return SQ(radial) * std::log(radial);
    }
}
