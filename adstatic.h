/*
 * adstatic.h
 *
 * Copyright (c) 2013 Vladimir Dabarov
 *
 * This file is part of adstatic.
 * 
 * Adstatic is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 * 
 * Adstatic is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Adstatic.  If not, see <http://www.gnu.org/licenses/>.
 *
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
