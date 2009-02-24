#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "ls.h"
#include "pdcf.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    m_ui(new Ui::MainWindow)
{
    m_ui->setupUi(this);

    timeseries << "./";

    buttonGroup = new QButtonGroup;
    connect(buttonGroup, SIGNAL(buttonClicked(int)),
            this, SLOT(manageOpenButtClick(int)));
    buttonGroup->addButton(m_ui->openFile1Button, 0);
    openFileLabels.append(m_ui->file1NameLabel);

    m_ui->estPDCButton->setDisabled(true);
    m_ui->lagsTableWidget->setRowCount(0);
    m_ui->lagsTableWidget->setColumnCount(0);
    m_ui->lagsTableWidget->setColumnWidth(0, 20);
    m_ui->lagsTableWidget->setItem(0, 0, new QTableWidgetItem(QString("1")));
    m_ui->shiftsTableWidget->setRowCount(0);
    m_ui->shiftsTableWidget->setColumnCount(0);
    m_ui->shiftsTableWidget->setColumnWidth(0, 20);
    m_ui->shiftsTableWidget->setItem(0, 0, new QTableWidgetItem(QString("0")));
    m_ui->arCoeffsTabs->removeTab(1);
    m_ui->arCoeffsTabs->removeTab(0);

    m_ui->estLSButton->setDisabled(true);
    m_ui->estPDCButton->setDisabled(true);

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
}

MainWindow::~MainWindow()
{
    delete m_ui;
}

void MainWindow::manageOpenButtClick(int num) {
    m_ui->estLSButton->setDisabled(false);

    QString fileName = QFileDialog::getOpenFileName(this, trUtf8("Выберите файл..."),
                                                    "./",
                                                    trUtf8("Временные ряды (*.*)"));
    QFile file(fileName);
    QFileInfo fi(fileName);
    if (file.exists()) {
        openFileLabels.at(num)->setText(fi.fileName());
        timeseries[num] = fileName;
        if ((num+1) == buttonGroup->buttons().size()) {
            QPushButton *newButt = new QPushButton(QString(trUtf8("Открыть файл...")));
            QLabel *newLabel = new QLabel(QString(trUtf8("")));
            buttonGroup->addButton(newButt, num+1);
            openFileLabels.append(newLabel);
            m_ui->filesGridLayout->addWidget(newLabel, 2*(num+1), 0, 1, 1);
            m_ui->filesGridLayout->addWidget(newButt, 2*(num+1)+1, 0, 1, 1);
            timeseries << "./";
            m_ui->lagsTableWidget->setRowCount(m_ui->lagsTableWidget->rowCount()+1);
            m_ui->lagsTableWidget->setColumnCount(m_ui->lagsTableWidget->columnCount()+1);
            QStringList hlabels;
            for (int i = 0; i < (num+1); i++) {
                hlabels << QString("%1").arg(i+1);
            }
            m_ui->lagsTableWidget->setVerticalHeaderLabels(hlabels);
            m_ui->lagsTableWidget->setHorizontalHeaderLabels(hlabels);
            m_ui->lagsTableWidget->setColumnWidth(num+1, 20);
            for (int i = 0; i <= (num+1); i++) {
                m_ui->lagsTableWidget->setItem(i, num, new QTableWidgetItem(QString("1")));
                m_ui->lagsTableWidget->setItem(num, i, new QTableWidgetItem(QString("1")));
            }
            m_ui->shiftsTableWidget->setRowCount(m_ui->shiftsTableWidget->rowCount()+1);
            m_ui->shiftsTableWidget->setColumnCount(m_ui->shiftsTableWidget->columnCount()+1);
            m_ui->shiftsTableWidget->setVerticalHeaderLabels(hlabels);
            m_ui->shiftsTableWidget->setHorizontalHeaderLabels(hlabels);
            m_ui->shiftsTableWidget->setColumnWidth(num+1, 20);
            for (int i = 0; i <= (num+1); i++) {
                m_ui->shiftsTableWidget->setItem(i, num, new QTableWidgetItem(QString("0")));
                m_ui->shiftsTableWidget->setItem(num, i, new QTableWidgetItem(QString("0")));
            }
        }
    }
}

void MainWindow::on_estLSButton_clicked() {
    t.start();

    m_ui->estLSButton->setDisabled(true);
    for (int i = 0; i < buttonGroup->buttons().size(); i++) {
        buttonGroup->button(i)->setDisabled(true);
    }
    m_ui->estPDCButton->setDisabled(true);
    m_ui->lagsTableWidget->setDisabled(true);
    m_ui->shiftsTableWidget->setDisabled(true);
    m_ui->dimensionFromEdit->setDisabled(true);
    m_ui->dimensionToEdit->setDisabled(true);
    m_ui->dimensionIncEdit->setDisabled(true);
    m_ui->freqFromEdit->setDisabled(true);
    m_ui->freqToEdit->setDisabled(true);
    m_ui->vibEdit->setDisabled(true);
    m_ui->shiftFromEdit->setDisabled(true);
    m_ui->shiftToEdit->setDisabled(true);
    m_ui->shiftIncEdit->setDisabled(true);
    m_ui->estChBox->setDisabled(true);


    QStringList tst = timeseries;
    if (tst.last() == QString("./")) tst.removeLast();

    if (cmtObj->loadTS(tst) == -1) {
        m_ui->statusLabel2->setText(trUtf8("Ошибки в файлах данных"));
        return;
    }

    int TSCount = tst.count();

    using namespace boost::numeric::ublas;

    matrix<double> Lags(TSCount, TSCount);
    for (int i = 0; i < TSCount; i++) {
        for (int j = 0; j < TSCount; j++) {
            Lags(i, j) = m_ui->lagsTableWidget->itemAt(i, j)->text().toInt();
        }
    }

    matrix<double> Shifts(TSCount, TSCount);
    for (int i = 0; i < TSCount; i++) {
        for (int j = 0; j < TSCount; j++) {
            Shifts(i, j) = m_ui->shiftsTableWidget->itemAt(i, j)->text().toInt();
        }
    }

    int p_from = m_ui->dimensionFromEdit->text().toInt();
    int p_to = m_ui->dimensionToEdit->text().toInt() <= p_from ? p_from : m_ui->dimensionToEdit->text().toInt();
    int p_inc = m_ui->dimensionIncEdit->text().toInt() == 0 ? 1 : m_ui->dimensionIncEdit->text().toInt();

    int s_from = m_ui->shiftFromEdit->text().toInt();
    int s_to = m_ui->shiftToEdit->text().toInt() <= s_from ? s_from : m_ui->shiftToEdit->text().toInt();
    int s_inc = m_ui->shiftIncEdit->text().toInt() == 0 ? 1 : m_ui->shiftIncEdit->text().toInt();

    Ar.resize(0);
    Ar.resize(((p_to-p_from)/p_inc + 1)*((s_to-s_from)/s_inc + 1));

    lsObj->setParams(dirName, cmtObj, p_from, p_to, p_inc, s_from, s_to, s_inc, Lags, Shifts, &Ar);

    lsObj->start();
}

void MainWindow::estLSFinished() {

    int p_from = m_ui->dimensionFromEdit->text().toInt();
    int p_to = m_ui->dimensionToEdit->text().toInt() <= p_from ? p_from : m_ui->dimensionToEdit->text().toInt();
    int p_inc = m_ui->dimensionIncEdit->text().toInt() == 0 ? 1 : m_ui->dimensionIncEdit->text().toInt();

    int s_from = m_ui->shiftFromEdit->text().toInt();
    int s_to = m_ui->shiftToEdit->text().toInt() <= s_from ? s_from : m_ui->shiftToEdit->text().toInt();
    int s_inc = m_ui->shiftIncEdit->text().toInt() == 0 ? 1 : m_ui->shiftIncEdit->text().toInt();

    int pi = 0;
    for (int p = p_from; p <= p_to; p += p_inc) {
        for (int s = s_from; s <= s_to; s += s_inc) {
            QTableWidget *arCoeffsTable = new QTableWidget;
            m_ui->arCoeffsTabs->addTab(arCoeffsTable, QString("p=%1, s=%2").arg(p).arg(s));
            arCoeffsTable->setColumnCount(Ar.at(pi).count()*Ar.at(pi).at(0).size1());
            arCoeffsTable->setRowCount(Ar.at(pi).at(0).size1());
            QStringList header;
            for (int i = 0; i < Ar.at(pi).count(); i++) {
                for (unsigned int j = 0; j < Ar.at(pi).at(i).size1(); j++) {
                    header << QString("p=%1, %2").arg(i+1).arg(j+1);
                    for (unsigned int k = 0; k < Ar.at(pi).at(i).size2(); k++) {
                        arCoeffsTable->setItem(j,
                                               i*(Ar.at(pi).at(0).size1())+k,
                                               new QTableWidgetItem(QString("%1")
                                                                    .arg(Ar.at(pi).at(i)(j, k), 13, 'E', 6, ' ')));
                    }
                }
            }
            arCoeffsTable->setHorizontalHeaderLabels(header);
            pi++;
        }
    }

    pi = 0;
    for (int p = p_from; p <= p_to; p += p_inc) {
        for (int s = s_from; s <= s_to; s += s_inc) {
            QFile file_out(QString("./%1/ar_p=%2_s=%3.txt").arg(dirName).arg(p).arg(s));
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

    m_ui->estPDCButton->setDisabled(false);
    m_ui->estLSButton->setDisabled(false);

    if (m_ui->estChBox->isChecked()) {
        m_ui->estPDCButton->click();
    }
}

void MainWindow::on_estPDCButton_clicked() {

    m_ui->estLSButton->setDisabled(true);
    m_ui->estPDCButton->setDisabled(true);

    QStringList tst = timeseries;
    if (tst.last() == QString("./")) tst.removeLast();
    int TSCount = tst.count();

    using namespace boost::numeric::ublas;

    matrix<double> Lags(TSCount, TSCount);
    for (int i = 0; i < TSCount; i++) {
        for (int j = 0; j < TSCount; j++) {
            Lags(i, j) = m_ui->lagsTableWidget->itemAt(i, j)->text().toInt();
        }
    }

    matrix<double> Shifts(TSCount, TSCount);
    for (int i = 0; i < TSCount; i++) {
        for (int j = 0; j < TSCount; j++) {
            Shifts(i, j) = m_ui->shiftsTableWidget->itemAt(i, j)->text().toInt();
        }
    }

    pdcfResult.resize(0);
    pdcfResult.resize(Ar.count());

    int p_from = m_ui->dimensionFromEdit->text().toInt();
    int p_to = m_ui->dimensionToEdit->text().toInt() <= p_from ? p_from : m_ui->dimensionToEdit->text().toInt();
    int p_inc = m_ui->dimensionIncEdit->text().toInt() == 0 ? 1 : m_ui->dimensionIncEdit->text().toInt();

    int s_from = m_ui->shiftFromEdit->text().toInt();
    int s_to = m_ui->shiftToEdit->text().toInt() <= s_from ? s_from : m_ui->shiftToEdit->text().toInt();
    int s_inc = m_ui->shiftIncEdit->text().toInt() == 0 ? 1 : m_ui->shiftIncEdit->text().toInt();

    double freqFrom = m_ui->freqFromEdit->text().toInt()/(double)m_ui->vibEdit->text().toInt();
    double freqTo = m_ui->freqToEdit->text().toInt()/(double)m_ui->vibEdit->text().toInt();
    if ((freqFrom == 0) && (freqTo == 0)) freqTo = 0.5;

    pdcfObj->setParams(cmtObj, Lags, Shifts, p_from, p_to, p_inc, s_from, s_to, s_inc, freqFrom, freqTo, Ar, &pdcfResult);
    pdcfObj->start();
}

void MainWindow::estPDCFinished() {
    printResult();

    for (int i = 0; i < buttonGroup->buttons().size(); i++) {
        buttonGroup->button(i)->setDisabled(false);
    }
    m_ui->estPDCButton->setDisabled(false);
    m_ui->estLSButton->setDisabled(false);
    m_ui->lagsTableWidget->setDisabled(false);
    m_ui->shiftsTableWidget->setDisabled(false);
    m_ui->dimensionFromEdit->setDisabled(false);
    m_ui->dimensionToEdit->setDisabled(false);
    m_ui->dimensionIncEdit->setDisabled(false);
    m_ui->freqFromEdit->setDisabled(false);
    m_ui->freqToEdit->setDisabled(false);
    m_ui->vibEdit->setDisabled(false);
    m_ui->shiftFromEdit->setDisabled(false);
    m_ui->shiftToEdit->setDisabled(false);
    m_ui->shiftIncEdit->setDisabled(false);
    m_ui->estChBox->setDisabled(false);

    setStatusMsg(QString("Time elapsed: %1 ms").arg(t.elapsed()));

    printReport();
}

void MainWindow::printResult() {

    int p_from = m_ui->dimensionFromEdit->text().toInt();
    int p_to = m_ui->dimensionToEdit->text().toInt() <= p_from ? p_from : m_ui->dimensionToEdit->text().toInt();
    int p_inc = m_ui->dimensionIncEdit->text().toInt() == 0 ? 1 : m_ui->dimensionIncEdit->text().toInt();

    int s_from = m_ui->shiftFromEdit->text().toInt();
    int s_to = m_ui->shiftToEdit->text().toInt() <= s_from ? s_from : m_ui->shiftToEdit->text().toInt();
    int s_inc = m_ui->shiftIncEdit->text().toInt() == 0 ? 1 : m_ui->shiftIncEdit->text().toInt();

    int pi = 0;
    for (int p = p_from; p <= p_to; p += p_inc) {
        for (int s = s_from; s <= s_to; s += s_inc) {

            QFile file_out(QString("./%1/pdc_p=%2_s=%3.txt").arg(dirName).arg(p).arg(s));
            if (!file_out.open(QIODevice::WriteOnly | QIODevice::Text))
                return;

            QTextStream out(&file_out);

            QStringList tst = timeseries;
            if (tst.last() == QString("./")) tst.removeLast();
            int TSCount = tst.count();

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
            double freqFrom = m_ui->freqFromEdit->text().toInt()/(double)m_ui->vibEdit->text().toInt();
            double freqTo = m_ui->freqToEdit->text().toInt()/(double)m_ui->vibEdit->text().toInt();
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
            QFile plot_out(QString("./%1/pdc_p=%2_s=%3.plt").arg(dirName).arg(p).arg(s));
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
            out_plot << "\nset yrange [0:1.1]\n\n";

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
                    out_plot << "unset xtics\n";
                    out_plot << "unset ytics\n";
                    if (i != 0) {
                        out_plot << "set format y \"\"\n";
                    } else {
                        out_plot << "set ytics (\"0\" 0, \"1\" 1.0)\n";
                    }
                    if (j != TSCount-1) {
                        out_plot << "set format x \"\"\n";
                    } else {
                        QString str = "set xtics (";
                        for (int tics = 0; tics < 5; tics++) {
                            str += QString("\"%1\" %2").arg(freqFrom+tics*(freqTo-freqFrom)/5.0).arg(freqFrom+tics*(freqTo-freqFrom)/5.0);
                            str += ", ";
                        }
                        str.chop(2);
                        str += ")\n";
                        out_plot << str;
                    }
                    out_plot << "plot 'pdc" << "_p=" << p << "_s=" << s << ".txt' using 1:" << (k  ) << " lc 3 notitle\n";
                    out_plot << "plot 'pdc" << "_p=" << p << "_s=" << s << ".txt' using 1:" << (k+1) << " lc 0 notitle\n\n";
                    k += 2;
                }
            }

            out_plot << "set nomultiplot\n";

            plot_out.close();

            k = 2;
            for (int system_src = 0; system_src < Ar.at(pi).at(0).size1(); system_src++) {
                for (int system_dst = 0; system_dst < Ar.at(pi).at(0).size1(); system_dst++) {
                    // автоматически формируем файлы для каждого сочетания систем для гнуплота
                    QFile ploti_out(QString("./%1/pdc_p=%2_s=%5_%3to%4.plt")
                                    .arg(dirName)
                                    .arg(p)
                                    .arg(system_src+1)
                                    .arg(system_dst+1)
                                    .arg(s));
                    if (!ploti_out.open(QIODevice::WriteOnly | QIODevice::Text))
                        return;

                    QTextStream outi_plot(&ploti_out);

                    outi_plot << "#!/usr/bin/gnuplot -persist\n";
                    outi_plot << "\nreset\n";
                    outi_plot << "\nset terminal png size 800,600\n";
                    outi_plot << "\nset output \"pdc" << "_p=" << p << "_s=" << s << "_" << system_src+1 << "to" << system_dst+1 << ".png\"\n";
                    outi_plot << "\nset style data lines\n";
                    outi_plot << "\nset xrange [" << 0.9*freqFrom << ":" << 1.1*freqTo << "]\n";
                    QFileInfo fi_src(timeseries.at(system_src));
                    QFileInfo fi_dst(timeseries.at(system_dst));
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

void MainWindow::setStatusMsg(QString str) {
    m_ui->statusLabel2->setText(str);
}

void MainWindow::printReport() {
    QFile rep_out(QString("./%1/report.txt").arg(dirName));
    if (!rep_out.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out_rep(&rep_out);

    out_rep << "Analyzed files:\n";
    for (int i = 0; i < timeseries.count()-1; i++) {
        out_rep << timeseries.at(i) << "\n";
    }
    out_rep << "\nParameters:\n";
    out_rep << "dim_from: " << m_ui->dimensionFromEdit->text() << "\n";
    out_rep << "dim_to:   " << m_ui->dimensionToEdit->text() << "\n";
    out_rep << "dim_inc:  " << m_ui->dimensionIncEdit->text() << "\n";
    out_rep << "shift_from: " << m_ui->shiftFromEdit->text() << "\n";
    out_rep << "shift_to:   " << m_ui->shiftToEdit->text() << "\n";
    out_rep << "shift_inc:  " << m_ui->shiftIncEdit->text() << "\n";
    out_rep << "freq_from: " << m_ui->freqFromEdit->text() << "\n";
    out_rep << "freq_to:   " << m_ui->freqToEdit->text() << "\n";
    out_rep << "freq_vib:  " << m_ui->vibEdit->text() << "\n";
    out_rep << "Lags:\n";
    for (int i = 0; i < timeseries.count()-1; i++) {
        for (int j = 0; j < timeseries.count()-1; j++) {
            out_rep << m_ui->lagsTableWidget->item(i, j)->text() << " ";
        }
        out_rep << "\n";
    }
    out_rep << "Shifts:\n";
    for (int i = 0; i < timeseries.count()-1; i++) {
        for (int j = 0; j < timeseries.count()-1; j++) {
            out_rep << m_ui->shiftsTableWidget->item(i, j)->text() << " ";
        }
        out_rep << "\n";
    }
    out_rep << "time consumed: " << m_ui->statusLabel2->text() << "\n";

    rep_out.close();
}
