#ifndef PDCF_SHELL_H
#define PDCF_SHELL_H

#include <QtCore>

class PdcfShell : public QObject {
    Q_OBJECT

public:
    PdcfShell(QStringList filesWithData = QStringList(),
              int dimFrom = 1, int dimTo = 0, int dimStep = 0,
              int shiftFrom = 0, int shiftTo = 0, int shiftStep = 0,
              int window = 0,
              int dataStartFrom = 0, int dataStartTo = 0, int dataStartStep = 0,
              int dataEndFrom = 0, int dataEndTo = 0, int dataEndStep = 0,
              int dataNorm = 0,
              int cpuCount = 1,
              bool fromGUI = false,
              int bonf = 0);
    ~PdcfShell();

    void startCalc();

    void setParams(QStringList filesWithData,
                   int dimFrom, int dimTo, int dimStep,
                   int shiftFrom, int shiftTo, int shiftStep,
                   int window,
                   int dataStartFrom, int dataStartTo, int dataStartStep,
                   int dataEndFrom, int dataEndTo, int dataEndStep,
                   int cpuCount,
                   bool fromGUI);

signals:
    void allFinishedSignal();
    void infoMsg(QString);

private:
    void printReport();
    int incDataInterval();
    void prepareCalc();
    void finishWork();
    void printFinalReport();

    QString dirName;

    QTime t;

    QStringList filesWithData;
    int dimFrom;
    int dimTo;
    int dimStep;
    int shiftFrom;
    int shiftTo;
    int shiftStep;
    int window;
    int dataStartFrom;
    int dataStartTo;
    int dataStartStep;
    int dataEndFrom;
    int dataEndTo;
    int dataEndStep;
    int dataNorm;
    int cpuCount;
    bool fromGUI;
    int bonf;

    int currDataFrom;
    int currDataTo;

    QFile log_out;
    QTextStream out_log;

    QString currWindowDir;

    int TSLenAbs;

private slots:
    void setStatusMsg(QString);

};

#endif
