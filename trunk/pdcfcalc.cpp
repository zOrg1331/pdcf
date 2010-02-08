#include "pdcfcalc.h"

#include "ls.h"
#include "pdcf.h"

PDCFcalc::PDCFcalc(const QStringList & filesWithData,
                   const QString dirName,
                   const int currDataFrom, const int currDataTo,
                   CommonMathTools *cmtObj,
                   const int dim,
                   const int shift,
                   const int window) :
filesWithData(filesWithData), dirName(dirName),
currDataFrom(currDataFrom), currDataTo(currDataTo),
cmtObj(cmtObj),
dim(dim), shift(shift), window(window)
{

}

PDCFcalc::~PDCFcalc()
{
    delete cmtObj;
}

void PDCFcalc::run()
{
    estLS();
}

void PDCFcalc::estLS()
{
    LS lsObj;
    connect(&lsObj, SIGNAL(infoMsg(QString)), this, SIGNAL(infoMsg(QString)), Qt::DirectConnection);
    QString windowDir = QString("./%1/%2-%3").arg(dirName).arg(currDataFrom).arg(currDataTo);
    lsObj.setParams(windowDir,
                    cmtObj,
                    dim,
                    shift,
                    &Ar);
    lsObj.startCalc();

    estLSFinished();
}

void PDCFcalc::estLSFinished()
{
    QFile file_out(QString("./%1/%2-%3/ar_p=%4_s=%5.txt")
                   .arg(dirName)
                   .arg(currDataFrom).arg(currDataTo)
                   .arg(dim).arg(shift));
    if (!file_out.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&file_out);
    for (unsigned int i = 0; i < Ar.size(); i++) {
        for (unsigned int j = 0; j < Ar[i].size1(); j++) {
            QString str;
            for (unsigned int k = 0; k < Ar[i].size2(); k++) {
                str += QString("%1 ").arg(Ar[i](j, k), 13, 'E', 6, ' ');
            }
            out << str << "\n";
        }
        out << QString("\n");
    }

    file_out.close();

    estPDCF();
}

void PDCFcalc::estPDCF()
{
    PDCF pdcfObj;
    connect(&pdcfObj, SIGNAL(infoMsg(QString)), this, SIGNAL(infoMsg(QString)), Qt::DirectConnection);
    pdcfObj.setParams(cmtObj,
                      dim,
                      shift,
                      Ar,
                      window,
                      &pdcfResult);
    pdcfObj.startCalc();

    estPDCFinished();
}

void PDCFcalc::estPDCFinished()
{
    printResult();
}

void PDCFcalc::printResult()
{
    QFile lock(QString("./%1/pdc_cumul_p=%2_s=%3.lock")
               .arg(dirName)
               .arg(dim)
               .arg(shift));
    while (lock.exists()) {}
    lock.open(QIODevice::WriteOnly);

    QFile file_out(QString("./%1/%2-%3/pdc_p=%4_s=%5.txt")
                   .arg(dirName)
                   .arg(currDataFrom).arg(currDataTo)
                   .arg(dim)
                   .arg(shift));
    if (!file_out.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&file_out);

    int TSCount = Ar[0].size1();

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
            line << pdcfResult[freqI][i][0] << pdcfResult[freqI][i][1];
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

    // добавляем информацию об этом окне данных в общий файл отчета по окнам
    int curr_systems = 0;
    for (int i = 0; i < TSCount; i++) {
        for (int j = 0; j < TSCount; j++) {
            if (i != j) {
                QFile file_out_i(QString("./%1/pdc_cumul_p=%2_s=%3_%4-%5.txt")
                                 .arg(dirName)
                                 .arg(dim).arg(shift)
                                 .arg(i+1).arg(j+1));
                if (!file_out_i.open(QIODevice::Append | QIODevice::Text))
                    return;

                QTextStream out_i(&file_out_i);

                int cpl_new = 0;
                int cpl_total_new = 0;
                int freqI = 0;
                for (double freq = freqFrom; freq < freqTo; freq += 0.5/FREQ_RES, freqI++) {
                    QString num0 = QString("%1\t").arg(freq, 0, 'E', 6);
                    QString num1 = QString("%1\t").arg(currDataFrom);
                    QString num2 = QString("%1\t").arg(pdcfResult[freqI][curr_systems][0], 0, 'E', 6);
                    QString num3 = QString("%1\t").arg(pdcfResult[freqI][curr_systems][1], 0, 'E', 6);
                    double cpl = (pdcfResult[freqI][curr_systems][0] > pdcfResult[freqI][curr_systems][1])
                                 ? pdcfResult[freqI][curr_systems][0] : -1;
                    QString num4 = QString("%1\n").arg(cpl, 0, 'E', 6);
                    out_i << num0 << num1 << num2 << num3 << num4;

                    if (cpl != -1) cpl_new++;
                    cpl_total_new++;
                }
                out_i << "\n";
                freqI = 0;
                for (double freq = freqFrom; freq < freqTo; freq += 0.5/FREQ_RES, freqI++) {
                    QString num0 = QString("%1\t").arg(freq, 0, 'E', 6);
                    QString num1 = QString("%1\t").arg(currDataTo);
                    QString num2 = QString("%1\t").arg(pdcfResult[freqI][curr_systems][0], 0, 'E', 6);
                    QString num3 = QString("%1\t").arg(pdcfResult[freqI][curr_systems][1], 0, 'E', 6);
                    double cpl = (pdcfResult[freqI][curr_systems][0] > pdcfResult[freqI][curr_systems][1])
                                 ? pdcfResult[freqI][curr_systems][0] : -1;
                    QString num4 = QString("%1\n").arg(cpl, 0, 'E', 6);
                    out_i << num0 << num1 << num2 << num3 << num4;
                }
                out_i << "\n";

                file_out_i.close();

                int cpl_old = 0;
                int cpl_total_old = 0;
                QFile file_cpl_i(QString("./%1/pdc_cumul_p=%2_s=%3_%4-%5_cpl-cnt.txt")
                                 .arg(dirName)
                                 .arg(dim).arg(shift)
                                 .arg(i+1).arg(j+1));
                if (file_cpl_i.open(QIODevice::ReadOnly | QIODevice::Text)) {
                    QString line = file_cpl_i.readLine();
                    cpl_old = line.split(' ', QString::SkipEmptyParts).at(0).toInt();
                    cpl_total_old = line.split(' ', QString::SkipEmptyParts).at(1).toInt();
                    file_cpl_i.close();
                }

                if (!file_cpl_i.open(QIODevice::WriteOnly | QIODevice::Text))
                    return;
                QTextStream out_cpl_i(&file_cpl_i);
                QString num0 = QString("%1 ").arg(cpl_old+cpl_new);
                QString num1 = QString("%1 ").arg(cpl_total_old+cpl_total_new);
                QString num2 = QString("%1").arg((cpl_old+cpl_new)/(double)(cpl_total_old+cpl_total_new) * 100);
                out_cpl_i << num0 << num1 << num2;
                file_cpl_i.close();
            }
            curr_systems++;
        }
    }


    // автоматически формируем файл для гнуплота
    QFile plot_out(QString("./%1/%2-%3/pdc_p=%4_s=%5.plt")
                   .arg(dirName)
                   .arg(currDataFrom).arg(currDataTo)
                   .arg(dim).arg(shift));
    if (!plot_out.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out_plot(&plot_out);

    out_plot << "#!/usr/bin/gnuplot -persist\n";
    out_plot << "\nreset\n";
    out_plot << "\nset terminal png size 800,600\n";
    out_plot << "\nset output \"pdc" << "_p=" << dim << "_s=" << shift << ".png\"\n";
    out_plot << "\nset style data lines\n";
    out_plot << "\nset multiplot\n";
    out_plot << "\nset xrange [" << freqFrom << ":" << freqTo << "]\n";
    out_plot << "\nset yrange [0:1]\n\n";

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
            out_plot << "plot 'pdc" << "_p=" << dim << "_s=" << shift << ".txt' using 1:" << (k  ) << " lc 3 notitle,\\\n";
            out_plot << "     'pdc" << "_p=" << dim << "_s=" << shift << ".txt' using 1:" << (k+1) << " lc 0 notitle\n\n";
            k += 2;
        }
    }

    out_plot << "set nomultiplot\n";

    plot_out.close();

    k = 2;
    for (unsigned int system_src = 0; system_src < Ar[0].size1(); system_src++) {
        for (unsigned int system_dst = 0; system_dst < Ar[0].size1(); system_dst++) {
            // автоматически формируем файлы для каждого сочетания систем для гнуплота
            QFile ploti_out(QString("./%1/%6-%7/pdc_p=%2_s=%5_%3to%4.plt")
                            .arg(dirName)
                            .arg(dim)
                            .arg(system_src+1)
                            .arg(system_dst+1)
                            .arg(shift)
                            .arg(currDataFrom).arg(currDataTo));
            if (!ploti_out.open(QIODevice::WriteOnly | QIODevice::Text))
                return;

            QTextStream outi_plot(&ploti_out);

            outi_plot << "#!/usr/bin/gnuplot -persist\n";
            outi_plot << "\nreset\n";
            outi_plot << "\nset terminal png size 800,600\n";
            outi_plot << "\nset output \"pdc" << "_p=" << dim << "_s=" << shift << "_" << system_src+1 << "to" << system_dst+1 << ".png\"\n";
            outi_plot << "\nset style data lines\n";
            outi_plot << "\nset xrange [" << 0.9*freqFrom << ":" << 1.1*freqTo << "]\n";
            QFileInfo fi_src(filesWithData.at(system_src));
            QFileInfo fi_dst(filesWithData.at(system_dst));
            outi_plot << "set title \"" << fi_src.fileName() << "->" << fi_dst.fileName() << "\"\n";
            outi_plot << "set format x \"%g\"\n";
            outi_plot << "set format y \"%g\"\n";
            outi_plot << "plot 'pdc" << "_p=" << dim << "_s=" << shift << ".txt' using 1:" << (k  ) << " lc 3 notitle,\\\n";
            outi_plot << "     'pdc" << "_p=" << dim << "_s=" << shift << ".txt' using 1:" << (k+1) << " lc 0 notitle\n\n";

            ploti_out.close();

            k += 2;
        }
    }

    lock.close();
    lock.remove();
}
