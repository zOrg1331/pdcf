#include "common_math_tools.h"

#include <boost/numeric/ublas/lu.hpp>

using namespace boost::numeric::ublas;

#define ZERO 0.0000001

CommonMathTools::CommonMathTools() {

}

int CommonMathTools::loadTS(const QStringList & fileNames) {
    timeseries = fileNames;
    timeseriesCount = fileNames.count();
    QVector<QFile *> files;
    files.resize(timeseriesCount);
    for (int i = 0; i < timeseriesCount; i++) {
        QFile *file = new QFile(fileNames.at(i));
        if (!file->open(QIODevice::ReadOnly)) {
            return -1;
        }
        files[i] = file;
    }
    stdDeviation.resize(timeseriesCount);
    dispersion.resize(timeseriesCount);
    mean.resize(timeseriesCount);
    tsLen.resize(timeseriesCount);
    timeseriesValues.resize(timeseriesCount);
    timeseriesValuesNorm.resize(timeseriesCount);
    for (int i = 0; i < timeseriesCount; i++) {
        int line = 0;
        QString str;
        bool ok;
        double currNum = 0;
        double summ1 = 0;
        double summ2 = 0;
        double dev = 0;
        double disp = 0;

        while (!files.at(i)->atEnd()) {
            str = files[i]->readLine();
            str = str.split(' ', QString::SkipEmptyParts).last();
            currNum = str.toDouble(&ok);
            if (!ok) {
                QString info = QString("garbage %1 at line %2 in file: %3").arg(str).arg(line).arg(fileNames.at(i));
                emit infoMsg(info);
                return -1;
            }
            timeseriesValues[i].append(currNum);
            summ1 += currNum;
            summ2 += currNum*currNum;
            line++;
        }
        summ1 /= line;
        summ2 /= line;
        disp = summ2 - summ1*summ1;
        dev = sqrt(disp);

        dispersion[i] = disp;
        stdDeviation[i] = dev;
        mean[i] = summ1;
        tsLen[i] = line;
    }
    // проверяем, одинаковой ли длины файлы
    for (int i = 1; i < timeseriesCount; i++) {
        if (tsLen.at(i) != tsLen.at(i-1)) {
            emit infoMsg("timeseries has different lengths");
            return -1;
        }
    }
    // отнормируем ряды
    normalizeTS();

    for (int i = 0; i < timeseriesCount; i++) {
        files[i]->close();
    }

    return 0;
}

void CommonMathTools::normalizeTS() {
    for (int i = 0; i < timeseriesCount; i++) {
        timeseriesValuesNorm[i].resize(tsLen.at(0));
    }
    // нормируем к нулевому среднему
    // нормируем к единичной дисперсии
    for (int j = 0; j < tsLen.at(0); j++) {
        for (int i = 0; i < timeseriesCount; i++) {
            timeseriesValuesNorm[i][j] = timeseriesValues.at(i).at(j) - mean.at(i);
            timeseriesValuesNorm[i][j] /= sqrt(dispersion.at(i));
        }
    }
}

int CommonMathTools::getTScount() {
    return timeseriesCount;
}

int CommonMathTools::getTSlen() {
    return tsLen.at(0);
}

double CommonMathTools::getTSvalue(const int & TSNum,
                                   const int & num) {
    return timeseriesValues.at(TSNum).at(num);
}

double CommonMathTools::getTSvalueNorm(const int & TSNum,
                                       const int & num) {
    return timeseriesValuesNorm.at(TSNum).at(num);
}

int CommonMathTools::gaussSolve(const matrix<double>& A,
                                const vector<double>& B,
                                vector<double>& X) {
    int N = A.size1();

    double det = 0;
    calcDeterminant(A, det);
    if (fabs(det) < ZERO) {
        return -1;
    }

    matrix<double> Ac = A;
    permutation_matrix<double> P(N);
    vector<double> x(N);

    lu_factorize(Ac, P);
    x = B;
    lu_substitute(Ac, P, x);
    X = x;
    return 0;
}

int CommonMathTools::calcInverseMatrix(const matrix<double>& A,
                                       matrix<double>& B) {
    // create a working copy of the input
    matrix<double> mLu(A);

    // perform LU-factorization
    lu_factorize(mLu);

    // create identity matrix of "inverse"
    B.assign(identity_matrix<double> (A.size1()));

    // backsubstitute to get the inverse
    lu_substitute<matrix<double> const, matrix<double> >(mLu, B);
    return 0;
}

int CommonMathTools::calcDeterminant(const matrix<double>& A,
                                     double & Det) {
    // create a working copy of the input
    matrix<double> aLu(A);
    permutation_matrix<std::size_t> pivots(A.size1());

    lu_factorize(aLu, pivots);

    double det = 1.0;

    for (std::size_t i = 0; i < pivots.size(); ++i) {
        if (pivots(i) != i)
            det *= -1.0;
        det *= aLu(i, i);
    }
    Det = det;
    return 0;
}

int CommonMathTools::calcResiduals(const QList<matrix<double> >& Ar,
                                   const QVector<int> & Ps,
                                   const QVector<int> & Ss,
                                   const int & Sh,
                                   QVector<QVector<double> > & Residuals) {
    int N = tsLen.at(0);

    // сортируем входной вектор размерностей, чтобы найти потом максимальную
    QVector<int> P(Ps);
    qSort(P.begin(), P.end());

    // готовим вектор векторов остатков моделей
    QVector<QVector<double> > Res;
    Res.resize(Ss.count());
    for (int i = 0; i < Ss.count(); i++) {
        Res[i].resize(N-P.last()-Sh);
    }

    for (int i = P.last()+Sh; i < N; i++) {
        for (int s = 0; s < Ss.count(); s++) {
            double x_mine = 0;
            for (int s1 = 0; s1 < Ss.count(); s1++) {
                for (int p = 0; p < Ps.count(); p++) {
                    x_mine += getTSvalue(Ss.at(s1), i-Ps.at(p))*Ar.at(p)(Ss.at(s), Ss.at(s1)) -Sh;
                }
            }
            double x = getTSvalue(Ss.at(s), i);
            Res[s][i-P.last() -Sh] = x-x_mine;
        }
    }

    // возвращаем результат
    Residuals = Res;

    return 0;
}

double CommonMathTools::calcStdDeviation(const QVector<double> & ts) {
    double summ1 = 0;
    double summ2 = 0;
    int N = ts.count();
    for (int i = 0; i < N; i++) {
        summ1 += ts.at(i);
        summ2 += ts.at(i)*ts.at(i);
    }
    summ1 /= N;
    summ2 /= N;
    return sqrt(summ2 - summ1*summ1);
}
