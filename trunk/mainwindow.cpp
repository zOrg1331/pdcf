#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <boost/numeric/ublas/matrix.hpp>

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

    m_ui->estPDCButton->setDisabled(true);

    pdcfShell = new PdcfShell();

    connect(pdcfShell, SIGNAL(allFinishedSignal()), this, SLOT(estPDCFinished()));
    connect(pdcfShell, SIGNAL(infoMsg(QString)), this, SLOT(setStatusMsg(QString)));
}

MainWindow::~MainWindow()
{
    delete m_ui;
    delete pdcfShell;
}

void MainWindow::manageOpenButtClick(int num) {
    m_ui->estPDCButton->setDisabled(false);

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

void MainWindow::on_estPDCButton_clicked() {

    m_ui->estPDCButton->setDisabled(true);
    m_ui->lagsTableWidget->setDisabled(true);
    m_ui->shiftsTableWidget->setDisabled(true);
    m_ui->dimensionFromEdit->setDisabled(true);
    m_ui->dimensionToEdit->setDisabled(true);
    m_ui->dimensionIncEdit->setDisabled(true);
    m_ui->shiftFromEdit->setDisabled(true);
    m_ui->shiftToEdit->setDisabled(true);
    m_ui->shiftIncEdit->setDisabled(true);
    m_ui->dataFromEdit->setDisabled(true);
    m_ui->dataToEdit->setDisabled(true);
    m_ui->dataStepEdit->setDisabled(true);
    m_ui->windowEdit->setDisabled(true);
    m_ui->cpuCountEdit->setDisabled(true);
    QList<QAbstractButton *> buttons = buttonGroup->buttons();
    for (int i = 0; i < buttons.size(); i++) {
        buttons.at(i)->setDisabled(true);
    }

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

    int p_from = m_ui->dimensionFromEdit->text().toInt();
    int p_to = m_ui->dimensionToEdit->text().toInt();
    int p_inc = m_ui->dimensionIncEdit->text().toInt();

    int s_from = m_ui->shiftFromEdit->text().toInt();
    int s_to = m_ui->shiftToEdit->text().toInt();
    int s_inc = m_ui->shiftIncEdit->text().toInt();

    int window = m_ui->windowEdit->text().toInt();

    int dataFrom = m_ui->dataFromEdit->text().toInt();
    int dataTo = m_ui->dataToEdit->text().toInt();
    int dataStep = m_ui->dataStepEdit->text().toInt();

    int cpuCount = m_ui->cpuCountEdit->text().toInt();

    pdcfShell->setParams(tst,
                         p_from,
                         p_to,
                         p_inc,
                         s_from,
                         s_to,
                         s_inc,
                         window,
                         dataFrom,
                         dataTo,
                         dataStep,
                         cpuCount,
                         true);

    pdcfShell->startCalc();

}

void MainWindow::estPDCFinished() {
    for (int i = 0; i < buttonGroup->buttons().size(); i++) {
        buttonGroup->button(i)->setDisabled(false);
    }
    m_ui->estPDCButton->setDisabled(false);
    m_ui->lagsTableWidget->setDisabled(false);
    m_ui->shiftsTableWidget->setDisabled(false);
    m_ui->dimensionFromEdit->setDisabled(false);
    m_ui->dimensionToEdit->setDisabled(false);
    m_ui->dimensionIncEdit->setDisabled(false);
    m_ui->shiftFromEdit->setDisabled(false);
    m_ui->shiftToEdit->setDisabled(false);
    m_ui->shiftIncEdit->setDisabled(false);
    m_ui->dataFromEdit->setDisabled(false);
    m_ui->dataToEdit->setDisabled(false);
    m_ui->dataStepEdit->setDisabled(false);
    m_ui->windowEdit->setDisabled(false);
    m_ui->cpuCountEdit->setDisabled(false);
    QList<QAbstractButton *> buttons = buttonGroup->buttons();
    for (int i = 0; i < buttons.size(); i++) {
        buttons.at(i)->setDisabled(false);
    }

}

void MainWindow::setStatusMsg(QString str) {
    m_ui->statusLabel2->setText(str);
}
