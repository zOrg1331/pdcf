#ifndef PDCF_SHELL_H
#define PDCF_SHELL_H

#include <QtCore>

#include "common_math_tools.h"

#include "ls.h"
#include "pdcf.h"

class PdcfShell : public QObject {
    Q_OBJECT

public:
    PdcfShell(QStringList filesWithData,
              bool calcOnlyAr = false,
              int dimFrom = 1,
              int dimTo = 0,
              int dimStep = 0,
              int shiftFrom = 0,
              int shiftTo = 0,
              int shiftStep = 0,
              int window = 0,
              int dataFrom = 0,
              int dataTo = 0,
              int dataStep = 0);
    ~PdcfShell();

    void startCalc();

private:
    void printResult();
    void printReport();
    void estLS();
    void estPDCF();
    void incDataInterval();

    QString dirName;

    CommonMathTools *cmtObj;
    LS *lsObj;
    PDCF *pdcfObj;
    QVector<QList<boost::numeric::ublas::matrix<double> > > Ar;
    QVector<QVector<QVector<QVector<double> > > > pdcfResult;

    QTime t;

    QStringList filesWithData;
    bool calcOnlyAr;
    int dimFrom;
    int dimTo;
    int dimStep;
    int shiftFrom;
    int shiftTo;
    int shiftStep;
    int window;
    int dataFrom;
    int dataTo;
    int dataStep;

    int currDataFrom;
    int currDataTo;

    int time_elapsed;

    QFile log_out;
    QTextStream out_log;

private slots:
    void setStatusMsg(QString);
    void estLSFinished();
    void estPDCFinished();

};

#endif