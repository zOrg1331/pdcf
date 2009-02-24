#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui>

#include "common_math_tools.h"

#include "ls.h"
#include "pdcf.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
    Q_DISABLE_COPY(MainWindow)
public:
    explicit MainWindow(QWidget *parent = 0);
    virtual ~MainWindow();

private:
    void printResult();
    void printReport();

    Ui::MainWindow *m_ui;

    QButtonGroup *buttonGroup;
    QList<QLabel *> openFileLabels;
    QStringList timeseries;
    QString dirName;

    CommonMathTools *cmtObj;
    LS *lsObj;
    PDCF *pdcfObj;
    QVector<QList<boost::numeric::ublas::matrix<double> > > Ar;
    QVector<QVector<QVector<QVector<double> > > > pdcfResult;

    QTime t;

private slots:
    void on_estPDCButton_clicked();
    void on_estLSButton_clicked();
    void manageOpenButtClick(int);
    void setStatusMsg(QString);
    void estLSFinished();
    void estPDCFinished();
};

#endif // MAINWINDOW_H
