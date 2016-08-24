//
// Created by Ke Liu on 8/5/16.
//

#ifndef ARBF_TEST_BASISFUNCTIONFACTORY_H
#define ARBF_TEST_BASISFUNCTIONFACTORY_H

#include <memory>
#include "BasisFunction.h"

/*
 * Abstract basis function factory
 */
class BasisFunctionFactory {
public:
    virtual ~BasisFunctionFactory();
    virtual std::unique_ptr<BasisFunction> createBasisFunction() = 0;
};

class MQBasisFunctionFactory: public BasisFunctionFactory {
public:
    ~MQBasisFunctionFactory();
    virtual std::unique_ptr<BasisFunction> createBasisFunction();
};

class IMQBasisFunctionFactory: public BasisFunctionFactory {
public:
    ~IMQBasisFunctionFactory();
    virtual std::unique_ptr<BasisFunction> createBasisFunction();
};

class GaussianBasisFunctionFactory: public BasisFunctionFactory {
public:
    ~GaussianBasisFunctionFactory();
    virtual std::unique_ptr<BasisFunction> createBasisFunction();
};

class TPSBasisFunctionFactory: public BasisFunctionFactory {
public:
    ~TPSBasisFunctionFactory();
    virtual std::unique_ptr<BasisFunction> createBasisFunction();
};

#endif //ARBF_TEST_BASISFUNCTIONFACTORY_H
