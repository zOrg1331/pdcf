#include "calcrblock.h"

CalcRBlock::CalcRBlock(const int & P_,
                       const int & N_,
                       const int & S_,
                       CommonMathTools *cmtObj_,
                       const matrix<double> & Lags_,
                       const matrix<double> & Shifts_,
                       const int & i_,
                       matrix<double> *R_,
                       matrix<double> *H_) {
    P = P_;
    N = N_;
    S = S_;
    cmtObj = cmtObj_;
    Lags = Lags_;
    Shifts = Shifts_;
    i = i_;
    R = R_;
    H = H_;
}

void CalcRBlock::run() {
    TSLen = cmtObj->getTSlen();
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < P; k++) {
            for (int l = 0; l < P; l++) {
                // считаем ковариацию между временными рядами i и j систем
                // при соответствующих размерностях k и l
                R->at_element(i*P+k, j*P+l) = calcCov(j, k+1, l+1);
                H->at_element(i*P+k, j*P+l) = 0;
            }
        }
    }
}

double CalcRBlock::calcCov(const int & j,
                           const int & k,
                           const int & l) {
    double eXi = 0;
    double eXj = 0;
    double eXji = 0;
    double x1 = 0, x2 = 0;
    for (int r = 0; r < TSLen; r++) {
        if (((r-k-S) >= 0) && ((r-l-S) >= 0)) {
//            x1 = cmtObj->getTSvalue(i, r-k-S/**Lags(i, j)-Shifts(i, j)*/);
            x1 = cmtObj->getTSvalueNorm(i, r-k-S/**Lags(i, j)-Shifts(i, j)*/);
            eXi += x1;
//            x2 = cmtObj->getTSvalue(j, r-l-S/**Lags(i, j)-Shifts(i, j)*/);
            x2 = cmtObj->getTSvalueNorm(j, r-l-S/**Lags(i, j)-Shifts(i, j)*/);
            eXj += x2;
            eXji += x1*x2;
        }
    }
    eXi /= TSLen;//-((k>l) ? k : l);
    eXj /= TSLen;//-((k>l) ? k : l);
    eXji /= TSLen;//-((k>l) ? k : l);
    return (eXji - eXi*eXj);
}
