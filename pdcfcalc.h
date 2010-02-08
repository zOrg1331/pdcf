#ifndef PDCFCALC_H
#define PDCFCALC_H

#include <QtCore>

#include "common_math_tools.h"

class PDCFcalc : public QThread
{
    Q_OBJECT
public:
    PDCFcalc(const QStringList & filesWithData,
             const QString dirName,
             const int currDataFrom, const int currDataTo,
             CommonMathTools *cmtObj,
             const int dim,
             const int shift,
             const int window);
    ~PDCFcalc();

protected:
    void run();

signals:
    void infoMsg(QString);

private:
    void estLS();
    void estLSFinished();
    void estPDCF();
    void estPDCFinished();
    void printResult();

    QStringList filesWithData;
    QString dirName;
    int currDataFrom;
    int currDataTo;
    CommonMathTools *cmtObj;
    int dim;
    int shift;
    int window;

    VECTOR_M Ar;
    VECTOR_VVD pdcfResult;

};

#endif // PDCFCALC_H
