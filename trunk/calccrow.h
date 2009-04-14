#ifndef CALCCROW_H
#define CALCCROW_H

#include <QThread>

#include "common_math_tools.h"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class CalcCRow : public QThread
{
    Q_OBJECT

public:
    CalcCRow(CommonMathTools *cmtObj,
             const int & TSNum,
             const int & Nu,
             const int & P,
             const int & S,
             const int & M,
             const matrix<double> & Lags,
             const matrix<double> & Shifts,
             const int & row,
             const int & N,
             matrix<double> *C);

protected:
    void run();

private:
    double calcPhi(const int & index,
                   const int & num);

    CommonMathTools *cmtObj;
    int TSNum;
    int Nu;
    int P;
    int S;
    int M;
    matrix<double> Lags;
    matrix<double> Shifts;
    int row;
    int N;
    matrix<double> *C;

    int TSLen;
};

#endif // CALCCROW_H
