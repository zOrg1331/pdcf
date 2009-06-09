#ifndef CALCRBLOCK_H
#define CALCRBLOCK_H

#include <QThread>

#include "common_math_tools.h"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class CalcRBlock : public QThread
{
    Q_OBJECT

    public:
    CalcRBlock(const int & P,
               const int & N,
               const int & S,
               CommonMathTools *cmtObj,
               /* const matrix<double> & Lags, */
               /* const matrix<double> & Shifts, */
               const int & i,
               matrix<double> *R,
               matrix<double> *H);

protected:
    void run();

private:
    double calcCov(const int & j,
                   const int & k,
                   const int & l);

    int P;
    int N;
    int S;
    CommonMathTools *cmtObj;
    /* matrix<double> Lags; */
    /* matrix<double> Shifts; */
    int i;
    matrix<double> *R;
    matrix<double> *H;

    int TSLen;
};

#endif // CALCRBLOCK_H
