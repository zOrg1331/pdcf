#include "pdcf.h"

#include <boost/math/distributions/chi_squared.hpp>

#include <boost/numeric/ublas/lu.hpp>

#include "common_math_tools.h"

#include "calcrblock.h"

#define PI   3.14159265358979323846
#define ZERO 0.0000001

PDCF::PDCF() {

}

void PDCF::setParams(CommonMathTools *cmtObj_,
                     const matrix<double> & Lags_,
                     const matrix<double> & Shifts_,
                     const int & P_from_,
                     const int & P_to_,
                     const int & P_inc_,
                     const int & S_from_,
                     const int & S_to_,
                     const int & S_inc_,
                     const double & freqFrom_,
                     const double & freqTo_,
                     const QVector<QList<matrix<double> > >& Ar_,
                     QVector<QVector<QVector<QVector<double> > > > *PDCFRes_) {

    cmtObj = cmtObj_;
    Lags = Lags_;
    Shifts = Shifts_;
    P_from = P_from_;
    P_to = P_to_;
    P_inc = P_inc_;
    S_from = S_from_;
    S_to = S_to_;
    S_inc = S_inc_;
    freqFrom = freqFrom_;
    freqTo = freqTo_;
    Ar = Ar_;
    PDCFRes = PDCFRes_;
}

void PDCF::run() {
    calcPDCF();
}

void PDCF::calcPDCF() {

    int Ari = 0;
    for (int P = P_from; P <= P_to; P += P_inc) {
        for (int Sh = S_from; Sh != (S_to+S_inc); Sh += S_inc) {

            // количество систем
            int N = Ar.at(Ari).at(0).size1();

            int TSLen = cmtObj->getTSlen();

            // создаем совокупность данных, представляющую собой результат
            QVector<QVector<QVector<double> > > result;
            {
                // заполняем вектора нулями
                QVector<double> tmp1;
                tmp1 << 0 << 0;
                QVector<QVector<double> > tmp2;
                for (int j = 0; j < N*N; j++) {
                    tmp2 << tmp1;
                }

                for (int i = 0; i < FREQ_RES; i++) {
                    result << tmp2;
                }
            }

            // заранее расчитаем квантиль распределения Хи-квадрат
            double alpha = 0.95;
            boost::math::chi_squared dist(1);
            double q = boost::math::quantile(dist, alpha);

            emit infoMsg(QString("PDC: starting estimating covariance matrix of %1 system").arg(Ari));
            matrix<double> H(N*P, N*P);
            matrix<double> R(N*P, N*P);
            int threadCount = QThread::idealThreadCount();
            if (threadCount < 1) threadCount = 1;
            for (int i = 0; i < N; i += threadCount) {
                QVector<CalcRBlock *> threads;
                for (int thr_i = 0; thr_i < threadCount; thr_i++) {
                    if ((i+thr_i) < N) {
                        CalcRBlock *thr = new CalcRBlock(P, N, Sh, cmtObj, Lags, Shifts,
                                                       i+thr_i,
                                                       &R, &H);
                        threads.append(thr);
                        thr->start();
                    }
                }
                for (int thr_i = 0; thr_i < threads.count(); thr_i++) {
                    threads.at(thr_i)->wait();
                }

//                CalcRBlock thr1(P, N, Sh, cmtObj, Lags, Shifts, i, &R, &H);
//                thr1.start();
//                if ((i+1) < N) {
//                    CalcRBlock thr2(P, N, Sh, cmtObj, Lags, Shifts, i+1, &R, &H);
//                    thr2.start();
//                    thr1.wait();
//                    thr2.wait();
//                } else {
//                    thr1.wait();
//                }
                QString str = QString("PDC: estimating R: estimated %1 of %2 blocks of %3 system")
//                              .arg(i+1)
                              .arg(i+threadCount)
                              .arg(N)
                              .arg(Ari);
                emit infoMsg(str);
                //for (int j = 0; j < N; j += 2) {
                //                for (int k = 0; k < P; k++) {
                //                    for (int l = 0; l < P; l++) {
                //                        // считаем ковариацию между временными рядами i и j систем
                //                        // при соответствующих размерностях k и l
                //                        R(i*P+k, j*P+l) = calcCov(cmtObj, Lags, Shifts, i, j, k+1, l+1);
                //                        H(i*P+k, j*P+l) = 0;
                //                    }
                //                }
                //  }
                //}
            }

            // считаем обратную матрицу
            emit infoMsg(QString("PDC: inverting covariance matrix of %1 system").arg(Ari));
            cmtObj->calcInverseMatrix(R, H);
            emit infoMsg(QString("PDC: inverting covariance matrix finished of %1 system").arg(Ari));

            // считаем матрицу ковариаций шума
            emit infoMsg(QString("PDC: estimating noise covariance matrix of %1 system").arg(Ari));
            matrix<double> Z(N, N);
            {
                QVector<int> Ps;
                for (int i = 0; i < P; i++) Ps << (i+1);
                QVector<int> Ss;
                for (int i = 0; i < N; i++) Ss << i;
                QVector<QVector<double> > Res;

                cmtObj->calcResiduals(Ar.at(Ari), Ps, Ss, Sh, Res);

                for (int ts = 0; ts < N; ts++) {
                    double eXX = 0;
                    for (int i = 0; i < Res.at(ts).count(); i++) {
                        eXX += Res.at(ts).at(i)*Res.at(ts).at(i);
                    }
                    eXX /= (double)Res.at(ts).count();
                    Z(ts, ts) = eXX;
                    //                qDebug() << "Z" << ts << Z(ts, ts);
                }
            }
            emit infoMsg(QString("PDC: extimating noise covariance matrix finished of %1 system").arg(Ari));

            // пробегаемся по всему интервалу частот
            int freqI = 0;
            emit infoMsg(QString("PDC: estimating PDC started of %1 system").arg(Ari));
            for (double freq = freqFrom; freq < freqTo; freq += 0.5/FREQ_RES, freqI++) {
                // переводим заданную матрицу коэффициентов в частотную область
                matrix<std::complex<double> > Af;
                Af = calcAf(Ar.at(Ari), freq);

                // вычитаем полученную матрицу из единичной
                identity_matrix<std::complex<double> > I (N);
                matrix<std::complex<double> > AfI = I - Af;

                // пробегаемся по всем сочетаниям компонент и считаем pdcf
                for (int source = 0; source < N; source++) {
                    for (int dest = 0; dest < N; dest++) {
                        // считаем знаменатель формулы для расчета PDCF
                        double summ = 0;
                        for (int i = 0; i < N; i++) {
                            summ += ( AfI(i, source).real() * AfI(i, source).real() +
                                      AfI(i, source).imag() * AfI(i, source).imag() );
                        }
                        // считаем PDCF
                        double pdcf = sqrt( AfI(dest, source).real() * AfI(dest, source).real() +
                                            AfI(dest, source).imag() * AfI(dest, source).imag() ) /
                                      sqrt( summ );
                        result[freqI][dest+source*N][0] = pdcf;

                        double summ1 = 0;
                        for (int ki = 0; ki < P; ki++) {
                            for (int li = 0; li < P; li++) {
                                summ1 += H(source*P+ki, source*P+li)
                                         * (cos(ki*2*PI*freq)*cos(li*2*PI*freq)
                                            + sin(ki*2*PI*freq)*sin(li*2*PI*freq));
                            }
                        }
                        double C = Z(dest, dest)*summ1;
                        double p_level = sqrt((C * q)/(TSLen * summ));
                        result[freqI][dest+source*N][1] = p_level;
                    }
                }
                QString str = QString("PDC: estimating PDC: estimated %1 of %2 freqs of %3 system").arg(freqI).arg(FREQ_RES).arg(Ari);
                emit infoMsg(str);
            }

            (*PDCFRes)[Ari] = result;
            Ari++;
        }
    }
}

matrix<std::complex<double> > PDCF::calcAf(const QList<matrix<double> >& Ar,
                                           const double & freq) {
    int N = Ar.at(0).size1();
    // создаем матрицу комплексных чисел, тем же размером, что и матрица коэффициентов
    matrix<std::complex<double> > Af (N, N);
    // обязательной зануляем, инча вылазиют численные ошибки
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::complex<double> z(0, 0);
            Af(i, j) = z;
        }
    }
    for (int i = 0, r = 1; i < Ar.count(); i++, r++) {
        // создаем промежуточную матрицу комплексных чисел, идентичную матрице вещественных коэффициентов
        // фактически, делаем преобразование типов
        matrix<std::complex<double> > ArC = Ar.at(i);
        // показатель степени экспоненты
        double z0 = -2*PI*freq*r;
        // формула Эйлера
        std::complex<double> z(cos(z0), sin(z0));
        // считаем матрицу, которая равна сумме матриц на входе,
        // помноженных на экспоненту в комплексной степени
        Af += ArC * z;
        for (int k1 = 0; k1 < N; k1++) {
            for (int k2 = 0; k2 < N; k2++) {
                // обязательно сравниваем на оч. маленькое число, иначе вылазиют численные ошибки
                if (fabs(Af(k1, k2).real()) < ZERO) Af(k1, k2).real() = 0;
                if (fabs(Af(k1, k2).imag()) < ZERO) Af(k1, k2).imag() = 0;
            }
        }
    }
    return Af;
}

double PDCF::calcCov(CommonMathTools *cmtObj,
                     const matrix<double> & Lags,
                     const matrix<double> & Shifts,
                     const int & i, // система i
                     const int & j, // система j
                     const int & k,
                     const int & l) {
    qDebug() << "ACHTUNG!";
//    int TSLen = cmtObj->getTSlen();
//    double eXi = 0;
//    double eXj = 0;
//    double eXji = 0;
//    double x1 = 0, x2 = 0;
//    for (int r = ((k>l) ? k : l); r < TSLen; r++) {
//        x1 = cmtObj->getTSvalue(i, r-k*Lags(i, j)-Shifts(i, j));
//        eXi += x1;
//        x2 = cmtObj->getTSvalue(j, r-l*Lags(i, j)-Shifts(i, j));
//        eXj += x2;
//        eXji += x1*x2;
//    }
//    eXi /= TSLen;//-((k>l) ? k : l);
//    eXj /= TSLen;//-((k>l) ? k : l);
//    eXji /= TSLen;//-((k>l) ? k : l);
//    return (eXji - eXi*eXj);
}
