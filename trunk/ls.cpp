#include "ls.h"

#include "common_math_tools.h"

LS::LS()
{

}

void LS::setParams(const QString & baseDir_,
                   CommonMathTools *cmtObj_,
                   const int dimension_,
                   const int shift_,
                   VECTOR_M *arCoeffsVector_)
{
    baseDir = baseDir_;
    cmtObj = cmtObj_;
    dimension = dimension_;
    shift = shift_;
    arCoeffsVector = arCoeffsVector_;
}

void LS::startCalc()
{
    calcLS();
}

int LS::calcLS()
{
    int M = cmtObj->getTScount();

    int coeffsNum = dimension*M;

    VECTOR_M arCoeffsVectorTmp;
    arCoeffsVectorTmp.resize(dimension);
    for (unsigned int i = 0; i < arCoeffsVectorTmp.size(); i++) {
        MATRIX t(M, M);
        arCoeffsVectorTmp[i] = t;
    }

    for (int TSNum = 0; TSNum < M; TSNum++) {
        VECTOR_I tsIndexes;
        tsIndexes.resize(M-1);
        for (int i = 0, j = 0; i < M; i++) {
            if (i != TSNum) {
                tsIndexes[j] = i;
                j++;
            }
        }
        VECTOR_D arCoeffs;
        arCoeffs.resize(coeffsNum);

        cmtObj->lls_solve(TSNum, tsIndexes, dimension, &arCoeffs);

        emit infoMsg(QString("LS: P=%2, Sh=%3 solving matrix equation for TS = %1 finished")
                     .arg(TSNum).arg(dimension).arg(shift));
        for (int i = 0; i < coeffsNum; i++) {
            //                    if (i < P) {
            //                        Ar[i](TSNum, 0) = Res.at(i);
            //                    } else if (i >= P) {
            //                        Ar[i-(i/P)*P](TSNum, i/P) = Res.at(i);
            //                    }
//            if ((TSNum == 1)) {
//            }
//            int i_dim = i/dimension;
//            int i_ost = i - i_dim*dimension;
//            if (i_ost%M == 0) {
//                arCoeffsVectorTmp[i_dim](TSNum, TSNum) = arCoeffs[i];
//            } else {
//                int i1 = i_ost-M*M*(i_ost/(M*M));
//                if ((i1-(i1/M)*M) > TSNum) {
//                    arCoeffsVectorTmp[i_dim](TSNum, i1-(i1/M)*M) = arCoeffs[i];
//                } else {
//                    arCoeffsVectorTmp[i_dim](TSNum, (i1-(i1/M)*M)-1) = arCoeffs[i];
//                }
//            }
            //            arCoeffsVectorTmp[i%dimension](TSNum, i%M) = arCoeffs[i];

            int dim = i%dimension;
            int ts_1 = TSNum;
            int ts_2 = 0;
            if (i < dimension) {
                ts_2 = TSNum;
            } else {
                ts_2 = i/dimension;
                if (i/dimension > TSNum) {
                    ts_2 = i/dimension;
                } else {
                    ts_2 = i/dimension - 1;
                }
            }
//            qDebug() << "TSNum" << TSNum << "i" << i << "coeff" << arCoeffs[i];
//            qDebug() << "dim" << dim << "ts_1" << ts_1 << "ts_2" << ts_2;
            arCoeffsVectorTmp[dim](ts_1, ts_2) = arCoeffs[i];
            
//            if (i < dimension) {
//                arCoeffsVectorTmp[i](TSNum, TSNum) = arCoeffs[i];
//            } else if (i >= dimension) {
//                arCoeffsVectorTmp[i%dimension](TSNum, M-TSNum-1) = arCoeffs[i];
//            }
        }
    }
    VECTOR_I tsIndexes;
    tsIndexes.resize(M);
    for (int i = 0; i < M; i++) tsIndexes[i] = i;
    VECTOR_VD residuals;

    cmtObj->calcResiduals(arCoeffsVectorTmp, dimension, tsIndexes, shift, residuals);

    for (int ts = 0; ts < M; ts++) {
        QFile res_out(QString("./%1/residuals_p=%2_s=%3_ts=%4.txt")
                      .arg(baseDir)
                      .arg(dimension)
                      .arg(shift)
                      .arg(ts));
        if (!res_out.open(QIODevice::WriteOnly | QIODevice::Text))
            return -1;

        QTextStream out_res(&res_out);
        for (unsigned int i = 0; i < residuals[ts].size(); i++) {
            out_res << QString("%1\n").arg(residuals[ts][i], 13, 'E', 6, ' ');
        }
        res_out.close();
    }
    (*arCoeffsVector) = arCoeffsVectorTmp;
    return 0;
}
