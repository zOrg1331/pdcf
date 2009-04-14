#ifndef COMMONMATHTOOLS_H
#define COMMONMATHTOOLS_H

#include <QtCore>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

class CommonMathTools : public QObject
{
    Q_OBJECT

public:
    CommonMathTools();

    int loadTS(const QStringList & fileNames);

    int getTScount();

    int getTSlen();

    int getTSlenAbs();

    double getTSvalue(const int & TSNum,
                      const int & num);

    double getTSvalueNorm(const int & TSNum,
                          const int & num);

    int gaussSolve(const boost::numeric::ublas::matrix<double>& A,
                   const boost::numeric::ublas::vector<double>& B,
                   boost::numeric::ublas::vector<double>& X);

    int calcInverseMatrix(const boost::numeric::ublas::matrix<double>& A,
                          boost::numeric::ublas::matrix<double>& B);

    int calcDeterminant(const boost::numeric::ublas::matrix<double>& A,
                        double & Det);

    int calcResiduals(const QList<boost::numeric::ublas::matrix<double> >& Ar,
                      const QVector<int> & Ps,
                      const QVector<int> & Ss,
                      const int & Sh,
                      QVector<QVector<double> > & Residuals);

    double calcStdDeviation(const QVector<double> & ts);

    void setDataWindow(int dataFrom, int dataTo);

private:
    void normalizeTS();

    QStringList timeseries;
    QVector<QVector<double> > timeseriesValues;
    QVector<QVector<double> > timeseriesValuesNorm;
    QVector<double> stdDeviation;
    QVector<double> dispersion;
    QVector<double> mean;
    QVector<int> tsLen;
    int timeseriesCount;

    int currDataFrom;
    int currDataTo;

signals:
    void infoMsg(QString);

};

#endif // COMMONMATHTOOLS_H
