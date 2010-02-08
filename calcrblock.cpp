#include "calcrblock.h"

CalcRBlock::CalcRBlock(const int & P_,
                       const int & N_,
                       const int & S_,
                       CommonMathTools *cmtObj_,
                       // const matrix<double> & Lags_,
                       // const matrix<double> & Shifts_,
                       const int & i_,
                       matrix<double> *R_,
                       matrix<double> *H_) {
    P = P_;
    N = N_;
    S = S_;
    cmtObj = cmtObj_;
    // Lags = Lags_;
    // Shifts = Shifts_;
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
    // расчет ковариации между i и j системами
    // с учетом возможного сдвига.
    double eXi = 0;
    double eXj = 0;
    double eXji = 0;
    double x1 = 0, x2 = 0;
    int r_start = (i != j) ? S : 0;
    int TSLen_real = (i != j) ? (TSLen-S) : TSLen;
    for (int r = r_start; r < TSLen; r++) {
        if (i != j) {
            x1 = cmtObj->getTSvalueNorm(i, r-k /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
            x2 = cmtObj->getTSvalueNorm(j, r-l -S /**Lags(i, j)-Shifts(i, j)*/);
        } else {
            x1 = cmtObj->getTSvalueNorm(i, r-k /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
            x2 = cmtObj->getTSvalueNorm(j, r-l /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
        }
        eXi += x1;
        eXj += x2;
    }
    eXi /= (double)TSLen_real;
    eXj /= (double)TSLen_real;
    for (int r = r_start; r < TSLen; r++) {
        if (i != j) {
            x1 = cmtObj->getTSvalueNorm(i, r-k /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
            x2 = cmtObj->getTSvalueNorm(j, r-l -S /**Lags(i, j)-Shifts(i, j)*/);
        } else {
            x1 = cmtObj->getTSvalueNorm(i, r-k /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
            x2 = cmtObj->getTSvalueNorm(j, r-l /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
        }
        eXji += (x1-eXi)*(x2-eXj);
    }
    eXji /= (double)TSLen_real-1;
    return eXji;
}
