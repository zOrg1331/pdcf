#ifndef PDCF_H
#define PDCF_H

#include <QtCore>

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include "common_math_tools.h"

#define FREQ_RES 1000.0

class PDCF : public QThread
{
    Q_OBJECT

public:
    PDCF();

    void setParams(CommonMathTools *cmtObj,
                   const matrix<double> & Lags,
                   const matrix<double> & Shifts,
                   const int & P_from,
                   const int & P_to,
                   const int & P_inc,
                   const int & S_from,
                   const int & S_to,
                   const int & S_inc,
                   const double & freqFrom,
                   const double & freqTo,
                   const QVector<QList<matrix<double> > >& Ar,
                   QVector<QVector<QVector<QVector<double> > > > *PDCFRes,
                   const int & cpuCount);

protected:
    void run();

signals:
    void infoMsg(QString);

private:
    // функция возвращает указатель вектор векторов
    // 0й элемент внешнего вектора - соответствует 0й частоте
    // FREQ_RES элемент внешнего вектора - соответствует частоте 0.5
    // 0й элемент второго внешнего вектора - влияние 1->1
    // 1й элемент второго внешнего вектора - влияние 2->1 и т.д.
    // внутренний вектор состоит всегда из двух элементов:
    //   0й - значение pdcf
    //   1й - 95%й уровень значимости
    // на вход функция ожидает матрицу коэффициентов и длину временного ряда,
    // по которому она получена для расчета значимости
    void calcPDCF();

    // вспомогательная функция, переводит матрицу коэффициентов в частотную область
    matrix<std::complex<double> > calcAf(const QList<matrix<double> >& Ar,
                                         const double & freq);

    double calcCov(CommonMathTools *cmtObj,
                   const matrix<double> & Lags,
                   const matrix<double> & Shifts,
                   const int & i,
                   const int & j,
                   const int & k,
                   const int & l);

    CommonMathTools *cmtObj;
    matrix<double> Lags;
    matrix<double> Shifts;
    int P_from;
    int P_to;
    int P_inc;
    int S_from;
    int S_to;
    int S_inc;
    double freqFrom;
    double freqTo;
    QVector<QList<matrix<double> > > Ar;
    QVector<QVector<QVector<QVector<double> > > > *PDCFRes;
    int cpuCount;

};

#endif // PDCF_H
