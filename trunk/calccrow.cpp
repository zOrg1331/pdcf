#include "calccrow.h"

CalcCRow::CalcCRow(CommonMathTools *cmtObj_,
                   const int TSNum_,
                   const int Nu_,
                   const int P_,
                   const int S_,
                   const int M_,
                   // const matrix<double> & Lags_,
                   // const matrix<double> & Shifts_,
                   const int row_,
                   const int N_,
                   matrix<double> *C_) {
    cmtObj = cmtObj_;
    TSNum = TSNum_;
    Nu = Nu_;
    P = P_;
    S = S_;
    M = M_;
    // Lags = Lags_;
    // Shifts = Shifts_;
    row = row_;
    N = N_;
    C = C_;
}

double CalcCRow::calcPhi(const int index,
                         const int num) {

    if (Nu == 1) {

        if (index == -1) {
            //                x1 = cmtObj->getTSvalue(TSNum, num-S);
            return cmtObj->getTSvalueNorm(TSNum, num);
        }

        //            x1 = cmtObj->getTSvalue(ip,
        //                                    num-
        //                                    //Lags(TSNum, ip)*
        //                                    (index-ip*P+1)/*-
        //                                    Shifts(TSNum, ip)*/-
        //                                                       S);
        int ip = index/P;
        if (ip == TSNum) {
            return cmtObj->getTSvalueNorm(ip,
                                          num-
                                          //Lags(TSNum, ip)*
                                          (index-ip*P+1)/*-
                                                          Shifts(TSNum, ip)*/);
        } else {
            return cmtObj->getTSvalueNorm(ip,
                                          num-
                                          //Lags(TSNum, ip)*
                                          (index-ip*P+1)/*-
                                                          Shifts(TSNum, ip)*/-S);
        }
    } else {
        return 0;
    }
}

void CalcCRow::run() {
    double summC = 0;
    for (int col = 0; col < P*M; col++) {
        summC = 0;
        for (int n = 0; n < N; n++) {
            summC += calcPhi(row, n)*calcPhi(col, n);
        }

        C->at_element(row, col) = summC;
    }
}
