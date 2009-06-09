#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui>

#include "pdcf_shell.h"

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

    Ui::MainWindow *m_ui;

    QButtonGroup *buttonGroup;
    QList<QLabel *> openFileLabels;
    QStringList timeseries;

    PdcfShell *pdcfShell;

private slots:
    void on_estPDCButton_clicked();
    void manageOpenButtClick(int);
    void setStatusMsg(QString);

    void estPDCFinished();

};

#endif // MAINWINDOW_H
