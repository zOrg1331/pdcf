#include "ls.h"

#include "common_math_tools.h"

#include "calccrow.h"

LS::LS() {

}

void LS::setParams(const QString & baseDir_,
                   CommonMathTools *cmtObj_,
                   const int & P_from_,
                   const int & P_to_,
                   const int & P_inc_,
                   const int & S_from_,
                   const int & S_to_,
                   const int & S_inc_,
                   const matrix<double> & Lags_,
                   const matrix<double> & Shifts_,
                   QVector<QList<matrix<double> > > *ls_coeffs_list_) {
    baseDir = baseDir_;
    cmtObj = cmtObj_;
    P_from = P_from_;
    P_to = P_to_;
    P_inc = P_inc_;
    S_from = S_from_;
    S_to = S_to_;
    S_inc = S_inc_;
    Lags = Lags_;
    Shifts = Shifts_;
    ls_coeffs_list = ls_coeffs_list_;
}

void LS::run() {
    calcLS();
}

double LS::calcPhi(const int & TSNum,
                   const int & Nu,
                   const int & P,
                   const int & S,
                   const int & index,
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
        //                                    Shifts(TSNum, ip)*/-S);
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

int LS::calcLS() {

    int Pi = 0;
    for (int P = P_from; P <= P_to; P += P_inc) {
        for (int Sh = S_from; Sh != (S_to+S_inc); Sh += S_inc) {
            int M = cmtObj->getTScount();
            //    qDebug() << "M: " << M;
            int N = cmtObj->getTSlen();
            TSLen = N;
            //    qDebug() << "TSLen: " << N;

            int coeffsNum = P*M;
            QList<matrix<double> > Ar;
            for (int i = 0; i < P; i++) {
                matrix<double> Ari(M, M);
                for (int j = 0; j < M; j++) {
                    for (int k = 0; k < M; k++) {
                        Ari(j, k) = 0.0;
                    }
                }
                Ar.append(Ari);
            }


            for (int TSNum = 0; TSNum < M; TSNum++) {

                matrix<double> C (coeffsNum, coeffsNum);
                for (int j = 0; j < coeffsNum; j++) {
                    for (int k = 0; k < coeffsNum; k++) {
                        C(j, k) = 0.0;
                    }
                }
                vector<double> B (coeffsNum);
                for (int j = 0; j < coeffsNum; j++) {
                    B(j) = 0.0;
                }

                // генерируем матрицу C
                emit infoMsg(QString("LS: P=%1, Sh=%2 prepairing matrix").arg(P).arg(Sh));
                //    for (int row = 0; row < coeffsNum; row++) {
                //        for (int col = 0; col < coeffsNum; col++) {
                //            double summC = 0;
                //            double x1C, x2C;
                //            for (int n = 0; n < N; n++) {
                //                x1C = calcPhi(cmtObj, 0, 1, P, Lags, Shifts, row, n);
                //                x2C = calcPhi(cmtObj, 0, 1, P, Lags, Shifts, col, n);
                //                summC += x1C*x2C;
                //            }
                //
                //            C(row, col) = summC;
                //        }
                //        printf("\testimated: %d of %d coeffs\n", row+1, coeffsNum);
                //    }
                int threadCount = QThread::idealThreadCount();
                if (threadCount < 1) threadCount = 1;
                for (int row = 0; row < coeffsNum; row += threadCount) {
                    QVector<CalcCRow *> threads;
                    for (int thr_i = 0; thr_i < threadCount; thr_i++) {
                        if ((row+thr_i) < coeffsNum) {
                            CalcCRow *thr = new CalcCRow(cmtObj, TSNum, 1, P, Sh, M, Lags, Shifts,
                                                         row+thr_i,
                                                         N, &C);
                            threads.append(thr);
                            thr->start();
                        }
                    }
                    for (int thr_i = 0; thr_i < threads.count(); thr_i++) {
                        threads.at(thr_i)->wait();
                    }
//                    CalcCRow thr1(cmtObj, TSNum, 1, P, Sh, M, Lags, Shifts, row, N, &C);
//                    thr1.start();
//                    if ((row+1) < coeffsNum) {
//                        CalcCRow thr2(cmtObj, TSNum, 1, P, Sh, M, Lags, Shifts, row+1, N, &C);
//                        thr2.start();
//                        thr2.wait();
//                    }
//                    thr1.wait();
                    QString str = QString("LS: P=%3, Sh=%4, prepearing matrix for TS=%5: estimated %1 of %2 coeffs")
//                                  .arg(row+1)
                                  .arg(row+threadCount)
                                  .arg(coeffsNum)
                                  .arg(P)
                                  .arg(Sh)
                                  .arg(TSNum);
                    emit infoMsg(str);
                }


                QString str = QString("LS: P=%2, Sh=%3 working with TS: %1").arg(TSNum+1).arg(P).arg(Sh);
                emit infoMsg(str);


                for (int row = 0; row < coeffsNum; row++) {
                    double summB = 0;
                    double x1B, x2B;

                    for (int n = 0; n < N; n++) {
                        x1B = calcPhi(TSNum, 1, P, Sh, row, n);
                        x2B = calcPhi(TSNum, 1, P, Sh,  -1, n);
                        summB += x1B*x2B;
                    }
                    B(row) = summB;
                    QString str1 = QString("LS: P=%3, Sh=%4 estimated: %1 of %2 coeffs").arg(row+1).arg(coeffsNum).arg(P).arg(Sh);
                    emit infoMsg(str1);
                }

                vector<double> Res(coeffsNum);
                for (int i = 0; i < coeffsNum; i++) Res(i) = 0.0;

                str = QString("LS: P=%2, Sh=%3 solving matrix equation for TS = %1").arg(TSNum).arg(P).arg(Sh);
                emit infoMsg(str);
                if (cmtObj->gaussSolve(C, B, Res) == -1) {
                    emit infoMsg("cannot calc coeffs");
                    return -1;
                }
                str = QString("LS: P=%2, Sh=%3 solving matrix equation for TS = %1 finished").arg(TSNum).arg(P).arg(Sh);
                emit infoMsg(str);
                for (int i = 0; i < coeffsNum; i++) {
                    if (i < P) {
                        Ar[i](TSNum, 0) = Res(i);
                    } else if (i >= P) {
                        Ar[i-(i/P)*P](TSNum, i/P) = Res(i);
                    }
                }

            }
            QVector<int> Ps;
            for (int i = 0; i < P; i++) Ps << (i+1);
            QVector<int> Ss;
            for (int i = 0; i < M; i++) Ss << i;
            QVector<QVector<double> > Residuals;

            cmtObj->calcResiduals(Ar, Ps, Ss, Sh, Residuals);

            for (int ts = 0; ts < M; ts++) {
                QFile res_out(QString("./%1/residuals_p=%2_s=%3_ts=%4.txt")
                              .arg(baseDir)
                              .arg(P)
                              .arg(Sh)
                              .arg(ts));
                if (!res_out.open(QIODevice::WriteOnly | QIODevice::Text))
                    return -1;

                QTextStream out_res(&res_out);
                for (int i = 0; i < Residuals.at(ts).count(); i++) {
                    out_res << QString("%1\n").arg(Residuals.at(ts).at(i), 13, 'E', 6, ' ');
                }
                res_out.close();
            }
            (*ls_coeffs_list)[Pi] = Ar;
            Pi++;
        }
    }
    return 0;
}
