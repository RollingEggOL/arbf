//
// Created by Ke Liu on 8/5/16.
//

#include <memory>
#include "../include/BasisFunctionFactory.h"

// BasisFunctionFactory implementations
BasisFunctionFactory::~BasisFunctionFactory() {}


// MQBasisFunctionFactory implementations
MQBasisFunctionFactory::~MQBasisFunctionFactory() {}

std::unique_ptr<BasisFunction> MQBasisFunctionFactory::createBasisFunction() {
    std::unique_ptr<BasisFunction> ptr(new MQBasisFunction());
    return std::move(ptr);
}


// IMQBasisFunctionFactory implementations
IMQBasisFunctionFactory::~IMQBasisFunctionFactory() {}

std::unique_ptr<BasisFunction> IMQBasisFunctionFactory::createBasisFunction() {
    std::unique_ptr<BasisFunction> ptr(new IMQBasisFunction());
    return std::move(ptr);
}


// GaussianBasisFunctionFactory implementations
GaussianBasisFunctionFactory::~GaussianBasisFunctionFactory() {}

std::unique_ptr<BasisFunction> GaussianBasisFunctionFactory::createBasisFunction() {
    std::unique_ptr<BasisFunction> ptr(new GaussianBasisFunction());
    return std::move(ptr);
}


// TPSBasisFunctionFactory implementations
TPSBasisFunctionFactory::~TPSBasisFunctionFactory() {}

std::unique_ptr<BasisFunction> TPSBasisFunctionFactory::createBasisFunction() {
    std::unique_ptr<BasisFunction> ptr(new TPSBasisFunction());
    return std::move(ptr);
}
