#include "pdcf_shell.h"

#include "ls.h"
#include "pdcf.h"

PdcfShell::PdcfShell(QStringList filesWithData,
                     bool calcOnlyAr,
                     int dimFrom,
                     int dimTo,
                     int dimStep,
                     int shiftFrom,
                     int shiftTo,
                     int shiftStep,
                     int window,
                     int dataFrom,
                     int dataTo,
                     int dataStep,
                     int cpuCount) :
filesWithData(filesWithData), calcOnlyAr(calcOnlyAr),
dimFrom(dimFrom), dimTo(dimTo), dimStep(dimStep),
shiftFrom(shiftFrom), shiftTo(shiftTo), shiftStep(shiftStep),
window(window), dataFrom(dataFrom), dataTo(dataTo), dataStep(dataStep),
cpuCount(cpuCount)

{
    cmtObj = new CommonMathTools;
    connect(cmtObj, SIGNAL(infoMsg(QString)), this, SLOT(setStatusMsg(QString)));

    lsObj = new LS;
    connect(lsObj, SIGNAL(infoMsg(QString)), this, SLOT(setStatusMsg(QString)));
    connect(lsObj, SIGNAL(finished()), this, SLOT(estLSFinished()));

    pdcfObj = new PDCF;
    connect(pdcfObj, SIGNAL(infoMsg(QString)), this, SLOT(setStatusMsg(QString)));
    connect(pdcfObj, SIGNAL(finished()), this, SLOT(estPDCFinished()));

    dirName = QDateTime::currentDateTime().toString("yyyy-MM-dd_hh-mm");
    QDir dir;
    dir.mkdir(dirName);

    log_out.setFileName(QString("./%1/log.txt").arg(dirName));
    if (!log_out.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    out_log.setDevice(&log_out);

    currDataFrom = dataFrom;
    currDataTo = dataTo;
}

PdcfShell::~PdcfShell()
{
}

void PdcfShell::startCalc() {
    t.start();

    if (cmtObj->loadTS(filesWithData) == -1) {
        setStatusMsg(trUtf8("Ошибки в файлах данных"));
        return;
    }

    if (currDataTo == 0) {
        currDataTo = cmtObj->getTSlenAbs();
        if (window != 0) currDataTo = currDataFrom + window;
    }

    cmtObj->setDataWindow(currDataFrom, currDataTo);

    QString windowDir = QString("./%1/%2-%3").arg(dirName).arg(currDataFrom).arg(currDataTo);
    QDir dir;
    dir.mkdir(windowDir);

    estLS();
}

void PdcfShell::estLS() {
    int TSCount = filesWithData.count();

    using namespace boost::numeric::ublas;

    matrix<double> Lags(TSCount, TSCount);
    for (int i = 0; i < TSCount; i++) {
        for (int j = 0; j < TSCount; j++) {
            Lags(i, j) = 1;
        }
    }

    matrix<double> Shifts(TSCount, TSCount);
    for (int i = 0; i < TSCount; i++) {
        for (int j = 0; j < TSCount; j++) {
            Shifts(i, j) = 0;
        }
    }

    if (dimTo < dimFrom) dimTo = dimFrom;
    if (dimStep == 0) dimStep = 1;

    if (shiftStep == 0) {
        shiftStep = 1;
        shiftTo = shiftFrom;
    }

    Ar.resize(0);
    Ar.resize(((dimTo-dimFrom)/dimStep + 1)*((shiftTo-shiftFrom)/shiftStep + 1));

    QString dir = QString("%1/%2-%3").arg(dirName).arg(currDataFrom).arg(currDataTo);
    lsObj->setParams(dir,
                     cmtObj,
                     dimFrom, dimTo, dimStep,
                     shiftFrom, shiftTo, shiftStep,
                     Lags, Shifts, &Ar, cpuCount);

    lsObj->start();
}

void PdcfShell::estLSFinished() {

    int pi = 0;
    for (int p = dimFrom; p <= dimTo; p += dimStep) {
        for (int s = shiftFrom; s != (shiftTo+shiftStep); s += shiftStep) {
            QFile file_out(QString("./%1/%2-%3/ar_p=%4_s=%5.txt")
                           .arg(dirName)
                           .arg(currDataFrom).arg(currDataTo)
                           .arg(p).arg(s));
            if (!file_out.open(QIODevice::WriteOnly | QIODevice::Text))
                return;

            QTextStream out(&file_out);
            for (int i = 0; i < Ar.at(pi).count(); i++) {
                for (unsigned int j = 0; j < Ar.at(pi).at(i).size1(); j++) {
                    QString str;
                    for (unsigned int k = 0; k < Ar.at(pi).at(i).size2(); k++) {
                        str += QString("%1 ").arg(Ar.at(pi).at(i)(j, k), 13, 'E', 6, ' ');
                    }
                    out << str << "\n";
                }
                out << QString("\n");
            }

            file_out.close();
            pi++;
        }
    }
    if (!calcOnlyAr) {
        estPDCF();
    } else {
        incDataInterval();
    }
}

void PdcfShell::estPDCF() {

    int TSCount = filesWithData.count();

    using namespace boost::numeric::ublas;

    matrix<double> Lags(TSCount, TSCount);
    for (int i = 0; i < TSCount; i++) {
        for (int j = 0; j < TSCount; j++) {
            Lags(i, j) = 1;
        }
    }

    matrix<double> Shifts(TSCount, TSCount);
    for (int i = 0; i < TSCount; i++) {
        for (int j = 0; j < TSCount; j++) {
            Shifts(i, j) = 0;
        }
    }

    pdcfResult.resize(0);
    pdcfResult.resize(Ar.count());

    double freqFrom = 0;
    double freqTo = 0;
    if ((freqFrom == 0) && (freqTo == 0)) freqTo = 0.5;

    pdcfObj->setParams(cmtObj, Lags, Shifts,
                       dimFrom, dimTo, dimStep,
                       shiftFrom, shiftTo, shiftStep,
                       freqFrom, freqTo, Ar, &pdcfResult, cpuCount);
    pdcfObj->start();
}

void PdcfShell::estPDCFinished() {
    printResult();

    time_elapsed = t.elapsed();
    setStatusMsg(QString("Time elapsed: %1 ms").arg(time_elapsed));

    printReport();

    incDataInterval();
}

void PdcfShell::printResult() {

    int pi = 0;
    for (int p = dimFrom; p <= dimTo; p += dimStep) {
        for (int s = shiftFrom; s != (shiftTo+shiftStep); s += shiftStep) {
            QFile file_out(QString("./%1/%2-%3/pdc_p=%4_s=%5.txt")
                           .arg(dirName)
                           .arg(currDataFrom).arg(currDataTo)
                           .arg(p)
                           .arg(s));
            if (!file_out.open(QIODevice::WriteOnly | QIODevice::Text))
                return;

            QTextStream out(&file_out);

            int TSCount = filesWithData.count();

            // сначала распечатываем в выходной файл заголовок
            out << "freq" << "\t\t";
            for (int ni = 0; ni < TSCount; ni++) {
                for (int nj = 0; nj < TSCount; nj++) {
                    out << ni+1 << "->" << nj+1 << "\t\t";
                    out << "p_level\t\t";
                }
            }
            out << "\n";

            // проходимся по всем частотам и всем компонентам, считаем PDCF и выводим результат в файл
            double freqFrom = 0;
            double freqTo = 0;
            if ((freqFrom == 0) && (freqTo == 0)) freqTo = 0.5;
            int freqI = 0;
            for (double freq = freqFrom; freq < freqTo; freq += 0.5/FREQ_RES, freqI++) {
                QVector<double> line;
                line << freq;
                for (int i = 0; i < TSCount*TSCount; i++) {
                    line << pdcfResult.at(pi).at(freqI).at(i).at(0) << pdcfResult.at(pi).at(freqI).at(i).at(1);
                }
                for (int i = 0; i < line.size(); i++) {
                    QString num;
                    num = QString("%1").arg(line.at(i), 0, 'E', 6);
                    out << num << "\t";
                    if (num == "NAN") out << "\t";
                }
                out << "\n";
            }

            file_out.close();

            // автоматически формируем файл для гнуплота
            QFile plot_out(QString("./%1/%2-%3/pdc_p=%4_s=%5.plt")
                           .arg(dirName)
                           .arg(currDataFrom).arg(currDataTo)
                           .arg(p).arg(s));
            if (!plot_out.open(QIODevice::WriteOnly | QIODevice::Text))
                return;

            QTextStream out_plot(&plot_out);

            out_plot << "#!/usr/bin/gnuplot -persist\n";
            out_plot << "\nreset\n";
            out_plot << "\nset terminal png size 800,600\n";
            out_plot << "\nset output \"pdc" << "_p=" << p << "_s=" << s << ".png\"\n";
            out_plot << "\nset style data lines\n";
            out_plot << "\nset multiplot\n";
            out_plot << "\nset xrange [" << freqFrom << ":" << freqTo << "]\n";
            out_plot << "\nset yrange [0:]\n\n";

            // номер столбца значений PDCF
            int k = 2;
            for (int i = 0; i < TSCount; i++) {
                for (int j = 0; j < TSCount; j++) {
                    out_plot << "set size " << 1.0/TSCount << "," << 1.0/TSCount << "\n";
                    out_plot << "set origin " << (i*1.0/TSCount) << "," << (1-(j+1)*1.0/TSCount) << "\n";
                    out_plot << "set rmargin 0.0\n";
                    if (i == TSCount-1) {
                        out_plot << "set rmargin 1.0\n";
                    }
                    out_plot << "set ylabel \"" << i+1 << "->" << j+1 << "\" offset 0.5,0\n";
                    out_plot << "set format x \"%g\"\n";
                    out_plot << "set format y \"%g\"\n";
                    out_plot << "#unset xtics\n";
                    out_plot << "#unset ytics\n";
                    if (i != 0) {
                        out_plot << "#set format y \"\"\n";
                    } else {
                        out_plot << "#set ytics (\"0\" 0, \"1\" 1.0)\n";
                    }
                    if (j != TSCount-1) {
                        out_plot << "#set format x \"\"\n";
                    } else {
                        QString str = "#set xtics (";
                        for (int tics = 0; tics < 5; tics++) {
                            str += QString("\"%1\" %2").arg(freqFrom+tics*(freqTo-freqFrom)/5.0).arg(freqFrom+tics*(freqTo-freqFrom)/5.0);
                            str += ", ";
                        }
                        str.chop(2);
                        str += ")\n";
                        out_plot << str;
                    }
                    out_plot << "plot 'pdc" << "_p=" << p << "_s=" << s << ".txt' using 1:" << (k  ) << " lc 3 notitle,\\\n";
                    out_plot << "     'pdc" << "_p=" << p << "_s=" << s << ".txt' using 1:" << (k+1) << " lc 0 notitle\n\n";
                    k += 2;
                }
            }

            out_plot << "set nomultiplot\n";

            plot_out.close();

            k = 2;
            for (int system_src = 0; system_src < Ar.at(pi).at(0).size1(); system_src++) {
                for (int system_dst = 0; system_dst < Ar.at(pi).at(0).size1(); system_dst++) {
                    // автоматически формируем файлы для каждого сочетания систем для гнуплота
                    QFile ploti_out(QString("./%1/%6-%7/pdc_p=%2_s=%5_%3to%4.plt")
                                    .arg(dirName)
                                    .arg(p)
                                    .arg(system_src+1)
                                    .arg(system_dst+1)
                                    .arg(s)
                                    .arg(currDataFrom).arg(currDataTo));
                    if (!ploti_out.open(QIODevice::WriteOnly | QIODevice::Text))
                        return;

                    QTextStream outi_plot(&ploti_out);

                    outi_plot << "#!/usr/bin/gnuplot -persist\n";
                    outi_plot << "\nreset\n";
                    outi_plot << "\nset terminal png size 800,600\n";
                    outi_plot << "\nset output \"pdc" << "_p=" << p << "_s=" << s << "_" << system_src+1 << "to" << system_dst+1 << ".png\"\n";
                    outi_plot << "\nset style data lines\n";
                    outi_plot << "\nset xrange [" << 0.9*freqFrom << ":" << 1.1*freqTo << "]\n";
                    QFileInfo fi_src(filesWithData.at(system_src));
                    QFileInfo fi_dst(filesWithData.at(system_dst));
                    outi_plot << "set title \"" << fi_src.fileName() << "->" << fi_dst.fileName() << "\"\n";
                    outi_plot << "set format x \"%g\"\n";
                    outi_plot << "set format y \"%g\"\n";
                    outi_plot << "plot 'pdc" << "_p=" << p << "_s=" << s << ".txt' using 1:" << (k  ) << " lc 3 notitle,\\\n";
                    outi_plot << "     'pdc" << "_p=" << p << "_s=" << s << ".txt' using 1:" << (k+1) << " lc 0 notitle\n\n";

                    ploti_out.close();

                    k += 2;
                }
            }
            pi++;
        }
    }
}

void PdcfShell::setStatusMsg(QString str) {
    out_log << t.elapsed() << " : " << currDataFrom << "-" << currDataTo << " " << str << "\n";
    out_log.flush();
    log_out.flush();
}

void PdcfShell::printReport() {
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
    out_rep << "freq_from: " << 0 << "\n";
    out_rep << "freq_to:   " << 0 << "\n";
    out_rep << "freq_vib:  " << 0 << "\n";
    out_rep << "Lags:\n";
    for (int i = 0; i < filesWithData.count(); i++) {
        for (int j = 0; j < filesWithData.count(); j++) {
            out_rep << 1 << " ";
        }
        out_rep << "\n";
    }
    out_rep << "Shifts:\n";
    for (int i = 0; i < filesWithData.count(); i++) {
        for (int j = 0; j < filesWithData.count(); j++) {
            out_rep << 0 << " ";
        }
        out_rep << "\n";
    }
    out_rep << "time consumed: " << time_elapsed << "\n";

    rep_out.close();

    log_out.close();
}

void PdcfShell::incDataInterval() {
    if (window == 0) {
        if (dataStep == 0) currDataFrom = currDataTo;
        else currDataFrom += dataStep;
    } else {
        if (dataStep == 0) {
            currDataFrom += window;
            currDataTo = currDataFrom+window;
            if (currDataTo > cmtObj->getTSlenAbs()) currDataTo = cmtObj->getTSlenAbs();
        } else {
            currDataFrom += dataStep;
            currDataTo = currDataFrom+window;
            if (currDataTo > cmtObj->getTSlenAbs()) currDataTo = cmtObj->getTSlenAbs();
        }
    }

    if (currDataFrom < currDataTo) {
        cmtObj->setDataWindow(currDataFrom, currDataTo);

        QString windowDir = QString("./%1/%2-%3").arg(dirName).arg(currDataFrom).arg(currDataTo);
        QDir dir;
        dir.mkdir(windowDir);

        estLS();
    } else {
        exit(0);
    }
}
