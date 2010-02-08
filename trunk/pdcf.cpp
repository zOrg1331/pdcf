#include "pdcf.h"

#include <boost/math/distributions/chi_squared.hpp>

#include <boost/numeric/ublas/lu.hpp>

#include "common_math_tools.h"

#define PI   M_PIl
#define ZERO 0.0000001

PDCF::PDCF()
{
}

void PDCF::setParams(CommonMathTools *cmtObj_,
                     const int dimension_,
                     const int shift_,
                     const VECTOR_M & arCoeffs_,
                     const int windowSize_,
                     VECTOR_VVD *PDCFRes_)
{
    cmtObj = cmtObj_;
    dimension = dimension_;
    shift = shift_;
    arCoeffs = arCoeffs_;
    windowSize = windowSize_;
    PDCFRes = PDCFRes_;
}

void PDCF::startCalc()
{
    calcPDCF();
}

void PDCF::calcPDCF()
{
    /*
      Далее идет расчет характеристики связанности набора систем
      частной направленной когерентности.
      Исходная информация полностью содержится в коэффициентах
      АР-моделей этих систем. Коэффициенты этих АР-моделей
      хранятся в Ar.
      Ar - вектор из N элементов, где N количество систем,
      для которых необходимо провести расчет (например при разных
      размерностях и сдвигах).
      Каждый элемент - список из P матриц, где P - размерность
      АР-моделей.
      Каждая матрица - коэффициенты совместных и индивидуальных
      АР-моделей. Элементы матрицы задаются в следующем виде:
          1    2  .. n
      1 1<-1 1<-2 .. 1<-n
      2 2<-1 2<-2 .. 2<-n
      .   .    .  ..  .
      n n<-1 n<-2 .. n<-n
      где n - количество исходных временных рядов.
    */

    // количество систем
    int M = arCoeffs[0].size1();
    // длины временных рядов (все одинаковые)
    TSLen = cmtObj->getTSlen();

    // создаем совокупность данных, представляющую собой результат
    /*
     Для системы Ari результат представляется в следующем виде:
     вектор из FREQ_RES векторов, каждый из которых сам по себе
     вектор из N*N элементов, где N - количество временных
     рядов на входе, каждый из которых вектор из двух
     элементов - характеристика связанности и ее значимость.
    */
    // весь этот набор векторов необходимо инициализировать
    VECTOR_VVD result;
    {
        result.resize((int)FREQ_RES);
        for (int i = 0; i < (int)FREQ_RES; i++) {
            result[i].resize(M*M);
            for (int j = 0; j < M*M; j++) {
                result[i][j].resize(2);
            }
        }
    }

    // заранее расчитаем 95% квантиль распределения Хи-квадрат
    // потребуется при расчете значимости
    double alpha = 0.95;
    if (windowSize > 0) {
        alpha = 1.00 - 0.05/(double)((double)cmtObj->getTSlenAbs()/(double)windowSize);
    }
    boost::math::chi_squared dist(1);
    double q = boost::math::quantile(dist, alpha);

    /*
     Для расчета значимости необходимо вычислить матрицу
     ковариаций остатков. Это квадратная N*P размерная
     матрица. Каждый элемент матрицы - ковариация
     между временными рядами i-й и j-й систем при
     размерности p.
     i: 1..N; j: 1..N; p: 1..P
    */
    emit infoMsg(QString("PDC: starting estimating covariance matrix"));
    MATRIX H(M*dimension, M*dimension);
    MATRIX R(M*dimension, M*dimension);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < dimension; k++) {
                for (int l = 0; l < dimension; l++) {
                    // считаем ковариацию между временными рядами i и j систем
                    // при соответствующих размерностях k и l
                    R.at_element(i*dimension+k, j*dimension+l) = calcCov(i, j, k+1, l+1);
                    H.at_element(i*dimension+k, j*dimension+l) = 0;
                }
            }
        }
    }

    // считаем обратную матрицу матрице ковариаций
    cmtObj->calcInverseMatrix(R, H);

    /*
     Считаем матрицу ковариаций шума. Т.к. исходного ряда шумового
     процесса нет, то используются остатки АР-моделей.
     В результате, в матрице Z будет храниться NxN элементов,
     каждый из которых - ковариация остатков АР модели
     восстановленного временного ряда i-й системы при учете j-й
    */
    emit infoMsg(QString("PDC: estimating noise covariance matrix"));
    MATRIX Z(M, M);
    {
        // считаем ковариации только для диагональных элементов матрицы,
        // остальные нам не понадобятся
        for (int ts = 0; ts < M; ts++) {
            // вектор номером систем, для которых надо считать остатки
            VECTOR_I tsIndexes;
            tsIndexes.resize(1);
            tsIndexes[0] = ts;
            // вектор из N элементов, каждый из которых - ряд остатков N-й модели
            VECTOR_VD residuals;

            cmtObj->calcResiduals(arCoeffs, dimension, tsIndexes, shift, residuals);

            double eXX = 0;
            double eX = 0;
            for (unsigned int i = 0; i < residuals[0].size(); i++) {
                eX += residuals[0][i];
            }
            eX /= (double)residuals[0].size();
            for (unsigned int i = 0; i < residuals[0].size(); i++) {
                eXX += (residuals[0][i]-eX)*(residuals[0][i]-eX);
            }
            eXX /= (double)(residuals[0].size()-1);
            Z(ts, ts) = eXX;
        }
    }
    emit infoMsg(QString("PDC: extimating noise covariance matrix finished"));

    // пробегаемся по всему интервалу частот
    int freqI = 0;
    emit infoMsg(QString("PDC: estimating PDC started"));
    for (double freq = 0; freq < 0.5; freq += 0.5/FREQ_RES, freqI++) {
        // переводим заданную матрицу коэффициентов в частотную область
        MATRIXcmplx Af;
        Af = calcAf(arCoeffs, freq);

        // вычитаем полученную матрицу из единичной
        boost::numeric::ublas::identity_matrix<std::complex<double> > I (M);
        MATRIXcmplx AfI = I - Af;

        // пробегаемся по всем сочетаниям компонент и считаем pdcf
        for (int source = 0; source < M; source++) {
            for (int dest = 0; dest < M; dest++) {
                // считаем знаменатель формулы для расчета PDCF
                double summ = 0;
                for (int i = 0; i < M; i++) {
                    summ += fabs( AfI(i, source).real() * AfI(i, source).real() +
                                  AfI(i, source).imag() * AfI(i, source).imag() );
                }
                // считаем PDCF
                double pdcf = sqrt( AfI(dest, source).real() * AfI(dest, source).real() +
                                    AfI(dest, source).imag() * AfI(dest, source).imag() ) /
                              sqrt( summ );
                result[freqI][dest+source*M][0] = pdcf;

                // считаем значимость
                double summ1 = 0;
                for (int ki = 1; ki <= dimension; ki++) {
                    for (int li = 1; li <= dimension; li++) {
                        summ1 += H(source*dimension+ki-1, source*dimension+li-1)
                                 * (cos(ki*2*PI*freq)*cos(li*2*PI*freq)
                                    + sin(ki*2*PI*freq)*sin(li*2*PI*freq));
                    }
                }
                double C = Z(dest, dest)*fabs(summ1);
                double p_level = sqrt((C * q)/(TSLen * summ));
                result[freqI][dest+source*M][1] = p_level;
            }
        }
        //        QString str = QString("PDC: estimating PDC: estimated %1 of %2 freqs of %3 system").arg(freqI).arg(FREQ_RES).arg(Ari);
        //        emit infoMsg(str);
    }

    (*PDCFRes) = result;
}

MATRIXcmplx PDCF::calcAf(const VECTOR_M & arCoeffs_,
                         const double freq)
{
    int N = arCoeffs_[0].size1();
    // создаем матрицу комплексных чисел, тем же размером, что и матрица коэффициентов
    MATRIXcmplx Af (N, N);
    // обязательной зануляем, инча вылазиют численные ошибки
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::complex<double> z(0, 0);
            Af(i, j) = z;
        }
    }
    for (unsigned int i = 0, r = 1; i < arCoeffs_.size(); i++, r++) {
        // создаем промежуточную матрицу комплексных чисел, идентичную матрице вещественных коэффициентов
        // фактически, делаем преобразование типов
        MATRIXcmplx ArC = arCoeffs_[i];
        // показатель степени экспоненты
        double z0 = -2*PI*freq*r;
        // формула Эйлера
        std::complex<double> z(cos(z0), sin(z0));
        // считаем матрицу, которая равна сумме матриц на входе,
        // помноженных на экспоненту в комплексной степени
        Af += ArC * z;
//        for (int k1 = 0; k1 < N; k1++) {
//            for (int k2 = 0; k2 < N; k2++) {
//                // обязательно сравниваем на оч. маленькое число, иначе вылазиют численные ошибки
//                if (fabs(Af(k1, k2).real()) < ZERO) Af(k1, k2).real() = 0;
//                if (fabs(Af(k1, k2).imag()) < ZERO) Af(k1, k2).imag() = 0;
//            }
//        }
    }
    return Af;
}

double PDCF::calcCov(const int i,
                     const int j,
                     const int k,
                     const int l)
{
    // расчет ковариации между i и j системами
    // с учетом возможного сдвига.
    double eXi = 0;
    double eXj = 0;
    double eXji = 0;
    double x1 = 0, x2 = 0;
    int r_start = (i != j) ? (shift + (k > l ? k : l)) : (k > l ? k : l);
    int TSMax = (i != j) ? (TSLen-shift) : TSLen;
    int TSLen_real = TSMax-r_start;
    for (int r = r_start; r < TSMax; r++) {
        if (i != j) {
            x1 = cmtObj->getTSvalueNorm(i, r-k /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
            x2 = cmtObj->getTSvalueNorm(j, r-l -shift /**Lags(i, j)-Shifts(i, j)*/);
        } else {
            x1 = cmtObj->getTSvalueNorm(i, r-k /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
            x2 = cmtObj->getTSvalueNorm(j, r-l /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
        }
        eXi += x1;
        eXj += x2;
    }
    eXi /= (double)TSLen_real;
    eXj /= (double)TSLen_real;
    for (int r = r_start; r < TSMax; r++) {
        if (i != j) {
            x1 = cmtObj->getTSvalueNorm(i, r-k /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
            x2 = cmtObj->getTSvalueNorm(j, r-l -shift /**Lags(i, j)-Shifts(i, j)*/);
        } else {
            x1 = cmtObj->getTSvalueNorm(i, r-k /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
            x2 = cmtObj->getTSvalueNorm(j, r-l /*-S*/ /**Lags(i, j)-Shifts(i, j)*/);
        }
        eXji += (x1-eXi)*(x2-eXj);
    }
    eXji /= (double)TSLen_real;
    return eXji;
}
