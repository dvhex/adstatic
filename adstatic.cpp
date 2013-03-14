/*
 * adstatic.cpp
 *
 *  Created on: 14.03.2013
 *      Author: vladimir
 */

#include "adstatic.h"
#include "atlas/clapack.h"
#include "cblas.h"

using std::slice;

zvalarray ADStatic::Solve(double w)
{
    FormMatrix(w);
    return SolveMatrix();
}

void ADStatic::init()
{
    U[0] = params.U;
    U[1] = std::polar(1.0, - 2. / 3. * M_PI) * params.U;
    U[2] = std::polar(1.0, 2. / 3. * M_PI) * params.U;
}

void ADStatic::FormMatrix(double w)
{
    /* Формирование матрицы коэффициентов для вычисления токов на
     * основе параметров АД.
     */
    //сначала везде, кроме последнего столбца -1/Lm
    M[slice(0, 36, 1)] = - params.Lm / 2.;
    //потом заменяем диагонали на (Lm + Ls), (Lm+Lr) и Lm
    M[slice(3, 3, 7)] = params.Lm;
    M[slice(3*6, 3, 7)] = params.Lm;
    M[slice(0, 3, 7)] = params.Lm + params.Ls;
    M[slice(3*7, 3, 7)] = params.Lm + params.Lr;
    //умножаем всё на jw
    M *= dcomplex(0, W0);
    //добавляем к диагоналям коэффициенты из уравнения двигателя
    M[slice(0, 3, 7)] -= params.Rs;
    M[slice(3*7, 3, 7)] -= params.Rr;
    double k1 = params.p * w * sqrt(3.) * params.Lm * .5;
    double k2 = params.p * w / sqrt(3.) * (1.5 * params.Lm + params.Ls);
    size_t i[4][3] = {{5, 9, 16}, {4, 11, 15}, {23, 27, 34}, {22, 29, 33}};
    std::valarray<size_t> index[4];
    for (size_t j=0; j<4; j++)
        index[j] = std::valarray<size_t>(i[j], 3);
    M[index[0]] -= k1;
    M[index[1]] += k1;
    M[index[2]] -= k2;
    M[index[3]] += k2;
    //добавляем напряжения
    M[slice(6*6, 3, 1)] = U;
}


zvalarray ADStatic::SolveMatrix()
{
    /* Решение СЛАУ */
    int ipiv[6];
    clapack_zgesv(CblasRowMajor, 6, 6, &M[0], 6, ipiv, &M[36], 6);
    return M[slice(36, 6, 1)];
}
