#ifndef LS_H
#define LS_H

#include <QtCore>

#include "common_math_tools.h"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class LS : public QThread
{
    Q_OBJECT

public:
    LS();

    void setParams(const QString & baseDir,
                   CommonMathTools *cmtObj,
                   const int & P_from,
                   const int & P_to,
                   const int & P_inc,
                   const int & S_from,
                   const int & S_to,
                   const int & S_inc,
                   const matrix<double> & Lags,
                   const matrix<double> & Shifts,
                   QVector<QList<matrix<double> > > *ls_coeffs_list);

protected:
    void run();

signals:
    void infoMsg(QString);

private:
    int calcLS();

    double calcPhi(const int & TSNum,
                   const int & Nu,
                   const int & P,
                   const int & S,
                   const int & index,
                   const int & num);

    QString baseDir;
    CommonMathTools *cmtObj;
    int P_from;
    int P_to;
    int P_inc;
    int S_from;
    int S_to;
    int S_inc;
    matrix<double> Lags;
    matrix<double> Shifts;
    QVector<QList<matrix<double> > > *ls_coeffs_list;

    int TSLen;

};

#endif // LS_H
