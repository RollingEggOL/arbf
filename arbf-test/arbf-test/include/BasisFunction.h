//
// Created by Ke Liu on 8/5/16.
//

#ifndef ARBF_TEST_BASISFUNCTION_H
#define ARBF_TEST_BASISFUNCTION_H

/*
 * Abstract basis function.
 */
class BasisFunction {
public:
    virtual ~BasisFunction();
    virtual double phi(double radial) const = 0;
};

class MQBasisFunction: public BasisFunction {
public:
    ~MQBasisFunction();
    virtual double phi(double radial) const;
};

class IMQBasisFunction: public BasisFunction {
public:
    ~IMQBasisFunction();
    virtual double phi(double radial) const;
};

class GaussianBasisFunction: public BasisFunction {
public:
    ~GaussianBasisFunction();
    virtual double phi(double radial) const;
};

class TPSBasisFunction: public BasisFunction {
public:
    ~TPSBasisFunction();
    virtual double phi(double radial) const;
};

#endif //ARBF_TEST_BASISFUNCTION_H
