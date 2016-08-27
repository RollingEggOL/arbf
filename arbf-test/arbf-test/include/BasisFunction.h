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
    virtual double phi(double radial, double c) const = 0;
};

class MQBasisFunction: public BasisFunction {
public:
    ~MQBasisFunction();
    virtual double phi(double radial, double c=0.01) const;
};

class IMQBasisFunction: public BasisFunction {
public:
    ~IMQBasisFunction();
    virtual double phi(double radial, double c=1.8) const;
};

class GaussianBasisFunction: public BasisFunction {
public:
    ~GaussianBasisFunction();
    virtual double phi(double radial, double c=0.01) const;
};

class TPSBasisFunction: public BasisFunction {
public:
    ~TPSBasisFunction();
    virtual double phi(double radial, double c=0) const;
};

#endif //ARBF_TEST_BASISFUNCTION_H
