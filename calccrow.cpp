#include "calccrow.h"

CalcCRow::CalcCRow(CommonMathTools *cmtObj_,
                   const int & TSNum_,
                   const int & Nu_,
                   const int & P_,
                   const int & S_,
                   const int & M_,
                   const matrix<double> & Lags_,
                   const matrix<double> & Shifts_,
                   const int & row_,
                   const int & N_,
                   matrix<double> *C_) {
    cmtObj = cmtObj_;
    TSNum = TSNum_;
    Nu = Nu_;
    P = P_;
    S = S_;
    M = M_;
    Lags = Lags_;
    Shifts = Shifts_;
    row = row_;
    N = N_;
    C = C_;

    TSLen = cmtObj->getTSlen();
}

double CalcCRow::calcPhi(const int & index,
                         const int & num) {

    if (Nu == 1) {

        double x1 = 0;
        if (index == -1) {
            //                x1 = cmtObj->getTSvalue(TSNum, num-S);
            x1 = cmtObj->getTSvalueNorm(TSNum, num/*-S*/);
            return x1;
        }

        int ip = index/P;
        //            x1 = cmtObj->getTSvalue(ip,
        //                                    num-
        //                                    //Lags(TSNum, ip)*
        //                                    (index-ip*P+1)/*-
        //                                    Shifts(TSNum, ip)*/-
        //                                                       S);
        if (ip == TSNum) {
            x1 = cmtObj->getTSvalueNorm(ip,
                                        num-
                                        //Lags(TSNum, ip)*
                                        (index-ip*P+1)/*-
                                    Shifts(TSNum, ip)*/);

            return x1;
        } else {
            x1 = cmtObj->getTSvalueNorm(ip,
                                        num-
                                        //Lags(TSNum, ip)*
                                        (index-ip*P+1)/*-
                                    Shifts(TSNum, ip)*/-S);

            return x1;
        }
    } else {
        return 0;
    }
}

void CalcCRow::run() {
    for (int col = 0; col < P*M; col++) {
        double summC = 0;
        double x1C, x2C;
        for (int n = 0; n < N; n++) {
            x1C = calcPhi(row, n);
            x2C = calcPhi(col, n);
            summC += x1C*x2C;
        }

        C->at_element(row, col) = summC;
    }
}
