#include "pdcf_shell.h"

#include "pdcfcalc.h"

PdcfShell::PdcfShell(QStringList filesWithData,
                     int dimFrom, int dimTo, int dimStep,
                     int shiftFrom, int shiftTo, int shiftStep,
                     int window,
                     int dataStartFrom, int dataStartTo, int dataStartStep,
                     int dataEndFrom, int dataEndTo, int dataEndStep,
                     int dataNorm,
                     int cpuCount,
                     bool fromGUI,
                     int bonf) :
filesWithData(filesWithData),
dimFrom(dimFrom), dimTo(dimTo), dimStep(dimStep),
shiftFrom(shiftFrom), shiftTo(shiftTo), shiftStep(shiftStep),
window(window),
dataStartFrom(dataStartFrom), dataStartTo(dataStartTo), dataStartStep(dataStartStep),
dataEndFrom(dataEndFrom), dataEndTo(dataEndTo), dataEndStep(dataEndStep),
dataNorm(dataNorm),
cpuCount(cpuCount), fromGUI(fromGUI), bonf(bonf)
{
    dirName = QDateTime::currentDateTime().toString("yyyy-MM-dd_hh-mm-ss");
    QDir dir;
    dir.mkdir(dirName);

    log_out.setFileName(QString("./%1/log.txt").arg(dirName));
    if (!log_out.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    out_log.setDevice(&log_out);

    prepareCalc();
}

PdcfShell::~PdcfShell()
{
    log_out.close();
}

void PdcfShell::prepareCalc()
{
    t.start();

    CommonMathTools *cmtObj = new CommonMathTools;

    if (cmtObj->loadDataFromFiles(filesWithData) == -1) {
        setStatusMsg(trUtf8("Ошибки в файлах данных"));
        return;
    }
    TSLenAbs = cmtObj->getTSlenAbs();
    delete cmtObj;

    if (dataStartTo < dataStartFrom) dataStartTo = dataStartFrom;
    if (dataEndFrom <= dataStartFrom) dataEndFrom = TSLenAbs;
    if (dataEndTo < dataEndFrom) dataEndTo = dataEndFrom;

    currDataFrom = dataStartFrom;
    if (window != 0) currDataTo = dataStartFrom + window;
    else currDataTo = dataEndFrom;

    if (dimTo < dimFrom) dimTo = dimFrom;
    if (dimStep == 0) dimStep = 1;

    if (shiftStep == 0) {
        shiftStep = 1;
        shiftTo = shiftFrom;
    }
}

void PdcfShell::startCalc()
{
    QVector<PDCFcalc *> pdcfCalcObjs;

    int k = 0;
    do {
        for (int dim = dimFrom; dim <= dimTo; dim += dimStep) {
            for (int shift = shiftFrom; shift <= shiftTo; shift += shiftStep) {
                k++;

                CommonMathTools *cmtObjTmp = new CommonMathTools;
                connect(cmtObjTmp, SIGNAL(infoMsg(QString)), this, SLOT(setStatusMsg(QString)), Qt::DirectConnection);

                if (cmtObjTmp->loadDataFromFiles(filesWithData, currDataFrom, currDataTo, dataNorm) == -1) {
                    setStatusMsg(trUtf8("Ошибки в файлах данных"));
                    return;
                }
//                cmtObjTmp->setDataWindow(currDataFrom, currDataTo);

                QString windowDir = QString("./%1/%2-%3").arg(dirName).arg(currDataFrom).arg(currDataTo);
                QDir dir;
                dir.mkdir(windowDir);

                PDCFcalc *pdcfCalcObj = new PDCFcalc(filesWithData, dirName, currDataFrom, currDataTo,
                                                     cmtObjTmp, dim, shift, ((bonf == 1) ? window : 0));
                connect(pdcfCalcObj, SIGNAL(infoMsg(QString)), this, SLOT(setStatusMsg(QString)), Qt::DirectConnection);
                pdcfCalcObjs.append(pdcfCalcObj);

                printReport();

                setStatusMsg(QString("New thread started. dim: %1, shift: %2, dataFrom: %3, dataTo: %4")
                             .arg(dim).arg(shift).arg(currDataFrom).arg(currDataTo));
                pdcfCalcObj->start();

                if (k == cpuCount) {
                    for (int i = 0; i < pdcfCalcObjs.size(); i++) {
                        pdcfCalcObjs.at(i)->wait();
                        setStatusMsg(QString("Thread stopped."));
                    }
                    for (int i = 0; i < pdcfCalcObjs.size(); i++) {
                        delete pdcfCalcObjs[i];
                    }
                    pdcfCalcObjs.resize(0);
                    k = 0;
                }
            }
        }
    } while (incDataInterval() != -1);

    for (int i = 0; i < pdcfCalcObjs.size(); i++) {
        pdcfCalcObjs.at(i)->wait();
    }
    for (int i = 0; i < pdcfCalcObjs.size(); i++) {
        delete pdcfCalcObjs[i];
    }
    pdcfCalcObjs.resize(0);

    finishWork();
}

void PdcfShell::setStatusMsg(QString str)
{
    out_log << t.elapsed() << " : " << currDataFrom << "-" << currDataTo << " " << str << "\n";
    out_log.flush();
    log_out.flush();
    if (fromGUI) {
        emit infoMsg(str);
    }
}

int PdcfShell::incDataInterval()
{
    if (window == 0) {
        if ((dataStartStep == 0) && (dataEndStep == 0)) {
            return -1;
        }

        currDataFrom += dataStartStep;
        if (currDataFrom > dataStartTo) {
            if (dataEndStep != 0) {
                currDataFrom = dataStartFrom;
                currDataTo += dataEndStep;
            } else {
                return -1;
            }
        } else {
            currDataTo = dataEndFrom;
        }
    } else {
        currDataFrom += dataStartStep;
        currDataTo = currDataFrom + window;

        if (currDataTo > TSLenAbs) {
            return -1;
        }
    }

    return 0;
}

void PdcfShell::setParams(QStringList filesWithData_,
                          int dimFrom_, int dimTo_, int dimStep_,
                          int shiftFrom_, int shiftTo_, int shiftStep_,
                          int window_,
                          int dataStartFrom_, int dataStartTo_, int dataStartStep_,
                          int dataEndFrom_, int dataEndTo_, int dataEndStep_,
                          int cpuCount_,
                          bool fromGUI_)
{
    filesWithData = filesWithData_;
    dimFrom = dimFrom_;
    dimTo = dimTo_;
    dimStep = dimStep_;
    shiftFrom = shiftFrom_;
    shiftTo = shiftTo_;
    shiftStep = shiftStep_;
    window = window_;
    dataStartFrom = dataStartFrom_;
    dataStartTo = dataStartTo_;
    dataStartStep = dataStartStep_;
    dataEndFrom = dataEndFrom_;
    dataEndTo = dataEndTo_;
    dataEndStep = dataEndStep_;
    cpuCount = cpuCount_;
    fromGUI = fromGUI_;

    prepareCalc();
}

void PdcfShell::finishWork()
{
    printFinalReport();
    if (fromGUI) emit allFinishedSignal();
    else exit(0);
}

void PdcfShell::printFinalReport()
{
    QFile rep_out(QString("./%1/report.txt").arg(dirName));
    if (!rep_out.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out_rep(&rep_out);

    out_rep << "Analyzed files:\n";
    for (int i = 0; i < filesWithData.count(); i++) {
        out_rep << filesWithData.at(i) << "\n";
    }
    out_rep << "\nParameters:\n";
    out_rep << "data_start_from: " << dataStartFrom << "\n";
    out_rep << "data_start_to:   " << dataStartTo << "\n";
    out_rep << "data_start_step:  " << dataStartStep << "\n";
    out_rep << "data_end_from: " << dataEndFrom << "\n";
    out_rep << "data_end_to:   " << dataEndTo << "\n";
    out_rep << "data_end_step:  " << dataEndStep << "\n";
    out_rep << "dim_from: " << dimFrom << "\n";
    out_rep << "dim_to:   " << dimTo << "\n";
    out_rep << "dim_inc:  " << dimStep << "\n";
    out_rep << "shift_from: " << shiftFrom << "\n";
    out_rep << "shift_to:   " << shiftTo << "\n";
    out_rep << "shift_inc:  " << shiftStep << "\n";
    out_rep << "bonferonni:  " << bonf << "\n";
    out_rep << "time consumed: " << t.elapsed() << "\n";

    rep_out.close();

    int pi = 0;
    for (int p = dimFrom; p <= dimTo; p += dimStep) {
        for (int s = shiftFrom; s != (shiftTo+shiftStep); s += shiftStep) {
            int TSCount = filesWithData.count();

            // автоматически формируем файл для гнуплота
            int current_systems = 0;
            for (int i = 0; i < TSCount; i++) {
                for (int j = 0; j < TSCount; j++) {
                    if (i != j) {
                        QFile cpl_in(QString("./%1/pdc_cumul_p=%2_s=%3_%4-%5_cpl-cnt.txt")
                                       .arg(dirName)
                                       .arg(p).arg(s)
                                       .arg(i+1).arg(j+1));
                        if (!cpl_in.open(QIODevice::ReadOnly | QIODevice::Text))
                            return;

                        QString line = cpl_in.readLine();
                        QString cpl_perc = line.split(' ', QString::SkipEmptyParts).at(2);

                        QFile plot_out(QString("./%1/pdc_cumul_p=%2_s=%3_%4-%5.plt")
                                       .arg(dirName)
                                       .arg(p).arg(s)
                                       .arg(i+1).arg(j+1));
                        if (!plot_out.open(QIODevice::WriteOnly | QIODevice::Text))
                            return;

                        QTextStream out_plot(&plot_out);

                        out_plot << "#!/usr/bin/gnuplot -persist\n";
                        out_plot << "\nreset\n";
                        out_plot << "\n#set terminal png size 800,600\n";
                        out_plot << "\nset terminal postscript eps font \"Arial,16\"\n";
                        out_plot << "\n#set output \"pdc_cumul_p=" << p << "_s=" << s << "_" << i+1 << "-" << j+1 << ".png\"\n";
                        out_plot << "\nset output \"pdc_cumul_p=" << p << "_s=" << s << "_" << i+1 << "-" << j+1 << ".eps\"\n";
                        out_plot << "\nset title \"data: ";
                        for (int k = 0; k < filesWithData.size(); k++) {
                            out_plot << filesWithData.at(k) << " ";
                        }
                        out_plot << "\\ncpl: " << i+1 << "->" << j+1 << " " << cpl_perc << "%\"";
                        out_plot << "\nset parametric";
                        out_plot << "\nset hidden3d\n";
                        out_plot << "\nset xrange [0:*]";
                        out_plot << "\nset yrange [0:0.5]";
                        out_plot << "\nset zrange [0:1]\n";
                        out_plot << "\nunset key\n";
                        out_plot << "\nset ylabel \"frequency, a.u.\"";
                        out_plot << "\nset xlabel \"data window, a.u.\"";
                        out_plot << "\nset zlabel \"P, a.u.\"\n";
                        out_plot << "\nset xtics out";
                        out_plot << "\nset ytics out\n";
                        out_plot << "\nset pm3d map\n";
                        out_plot << "\nset cbrange [0:1]\n";
                        out_plot << "\nset palette model RGB";
                        out_plot << "\nset palette model RGB defined (0 \"white\", 0.01 \"white\", 0.01 \"cyan\", 0.45 \"dark-blue\", 0.45 \"dark-green\", 0.8 \"yellow\", 1 \"red\")\n";
                        out_plot << "\nsplot \"pdc_cumul_p=" << p << "_s=" << s << "_" << i+1 << "-" << j+1 << ".txt\" using ($2):($1):5\n";

                        plot_out.close();
                    }
                    current_systems++;
                }
            }
            pi++;
        }
    }
}

void PdcfShell::printReport()
{
    QFile rep_out(QString("./%1/%2-%3/report.txt").arg(dirName).arg(currDataFrom).arg(currDataTo));
    if (!rep_out.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out_rep(&rep_out);

    out_rep << "Analyzed files:\n";
    for (int i = 0; i < filesWithData.count(); i++) {
        out_rep << filesWithData.at(i) << "\n";
    }
    out_rep << "\nParameters:\n";
    out_rep << "dim_from: " << dimFrom << "\n";
    out_rep << "dim_to:   " << dimTo << "\n";
    out_rep << "dim_inc:  " << dimStep << "\n";
    out_rep << "shift_from: " << shiftFrom << "\n";
    out_rep << "shift_to:   " << shiftTo << "\n";
    out_rep << "shift_inc:  " << shiftStep << "\n";
    out_rep << "bonferonni:  " << bonf << "\n";

    rep_out.close();
}
