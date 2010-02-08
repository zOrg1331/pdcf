#ifndef LS_H
#define LS_H

#include <QtCore>

#include "common_math_tools.h"

class LS : public QObject
{
    Q_OBJECT

    public:
    LS();

    void setParams(const QString & baseDir,
                   CommonMathTools *cmtObj,
                   const int dimension,
                   const int shift,
                   VECTOR_M *ls_coeffs_list);

    void startCalc();

signals:
    void infoMsg(QString);

private:
    int calcLS();

    QString baseDir;
    CommonMathTools *cmtObj;
    int dimension;
    int shift;
    VECTOR_M *arCoeffsVector;

};

#endif // LS_H
