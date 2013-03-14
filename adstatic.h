/*
 * adstatic.h
 *
 *  Created on: 14.03.2013
 *      Author: vladimir
 */

#ifndef ADSTATIC_H_
#define ADSTATIC_H_

#include <valarray>
#include <complex>

typedef std::complex<double> dcomplex;
typedef std::valarray<dcomplex> zvalarray;

const double W0 = 100.0 * M_PI;

struct ADParams
{
    int p;
    double Rs, Rr, Ls, Lr, Lm, U;
};

class ADStatic
{
public:
    ADStatic() : U(3), M(42) {init();}
    ADStatic(const ADParams &_params) : params(_params), U(3), M(42) {init();}
    virtual ~ADStatic() {}
    zvalarray Solve(double w);
private:
    ADParams params;
    zvalarray U;
    zvalarray M; //порядок расположения элементов матрицы - вертикальное
    void init();
    void FormMatrix(double w);
    zvalarray SolveMatrix();
};


#endif /* ADSTATIC_H_ */
