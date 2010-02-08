#ifndef PDCF_H
#define PDCF_H

#include <QtCore>

#include "common_math_tools.h"

#define FREQ_RES 1000.0

class PDCF : public QObject
{
    Q_OBJECT

public:
    PDCF();

    void setParams(CommonMathTools *cmtObj,
                   const int dimension,
                   const int shift,
                   const VECTOR_M & arCoeffs,
                   const int windowSize,
                   VECTOR_VVD *PDCFRes);

    void startCalc();

signals:
    void infoMsg(QString);

private:
    double calcCov(const int i,
                   const int j,
                   const int k,
                   const int l);

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
    MATRIXcmplx calcAf(const VECTOR_M & arCoeffs,
                       const double freq);

    CommonMathTools *cmtObj;
    int dimension;
    int shift;
    VECTOR_M arCoeffs;
    VECTOR_VVD *PDCFRes;
    int windowSize;
    int TSLen;

};

#endif // PDCF_H
