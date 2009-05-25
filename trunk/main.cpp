#include <QtCore>

#ifndef QT_GUI_LIB
#include <QCoreApplication>
#endif

#ifdef QT_GUI_LIB
#include <QApplication>

#include "mainwindow.h"
#endif
#include "pdcf_shell.h"

#include <boost/program_options.hpp>
using namespace boost;
namespace po = boost::program_options;

#include <iostream>
using namespace std;

// for daemon mode
#include <sys/time.h>
#include <sys/resource.h>


#include <boost/random.hpp>
void generateData() {
    boost::rand48 rng(0);
    boost::rand48 rng1(1);
    boost::rand48 rng2(2);
    boost::rand48 rng3(3);
    boost::rand48 rng4(4);
    boost::normal_distribution<> one_norm(0.0, 1.0);
    boost::variate_generator<boost::rand48&, boost::normal_distribution<> > norm_rand(rng, one_norm);
    boost::variate_generator<boost::rand48&, boost::normal_distribution<> > norm_rand1(rng1, one_norm);
    boost::variate_generator<boost::rand48&, boost::normal_distribution<> > norm_rand2(rng2, one_norm);
    boost::variate_generator<boost::rand48&, boost::normal_distribution<> > norm_rand3(rng3, one_norm);
    boost::variate_generator<boost::rand48&, boost::normal_distribution<> > norm_rand4(rng4, one_norm);

    QStringList tst;
    tst << "./ts_test1.txt";
    tst << "./ts_test2.txt";
    tst << "./ts_test3.txt";
    tst << "./ts_test4.txt";
    tst << "./ts_test5.txt";

    QFile tst_file1("./ts_test1.txt");
    tst_file1.open(QIODevice::WriteOnly);
    QTextStream tst_stream1(&tst_file1);
    QFile tst_file2("./ts_test2.txt");
    tst_file2.open(QIODevice::WriteOnly);
    QTextStream tst_stream2(&tst_file2);
    QFile tst_file3("./ts_test3.txt");
    tst_file3.open(QIODevice::WriteOnly);
    QTextStream tst_stream3(&tst_file3);
    QFile tst_file4("./ts_test4.txt");
    tst_file4.open(QIODevice::WriteOnly);
    QTextStream tst_stream4(&tst_file4);
    QFile tst_file5("./ts_test5.txt");
    tst_file5.open(QIODevice::WriteOnly);
    QTextStream tst_stream5(&tst_file5);
    QVector<double> x1;
    int P = 20;
    for (int i = 0; i < P; i++) x1 << 1.0*norm_rand();
    x1.resize(50000);
    QVector<double> x2;
    for (int i = 0; i < P; i++) x2 << 1.0*norm_rand();
    x2.resize(50000);
    QVector<double> x3;
    for (int i = 0; i < P; i++) x3 << 1.0*norm_rand();
    x3.resize(50000);
    QVector<double> x4;
    for (int i = 0; i < P; i++) x4 << 1.0*norm_rand();
    x4.resize(50000);
    QVector<double> x5;
    for (int i = 0; i < P; i++) x5 << 1.0*norm_rand();
    x5.resize(50000);
    for (int i = 0; i < P; i++) tst_stream1 << x1.at(i) << "\n";
    for (int i = 0; i < P; i++) tst_stream2 << x2.at(i) << "\n";
    for (int i = 0; i < P; i++) tst_stream3 << x3.at(i) << "\n";
    for (int i = 0; i < P; i++) tst_stream4 << x4.at(i) << "\n";
    for (int i = 0; i < P; i++) tst_stream5 << x5.at(i) << "\n";
    for (int i = 20; i < 50000; i++) {
//        x1[i] = 0.5*x1[i-1] + 0.3*x2[i-1] + 0.4*x3[i-1] + 1.0*norm_rand();
//        x2[i] = -0.5*x1[i-1] + 0.3*x2[i-1] + 1.0*x3[i-1] + 0.0*norm_rand1();
//        x3[i] = -0.3*x2[i-1] - 0.2*x3[i-1] + 1.0*norm_rand2();

        x1[i] = (1.4959)*x1[i-1]+(-0.67032)*x1[i-2]
                +1.0*norm_rand();
        x2[i] = (1.4959)*x2[i-1]+(-0.67032)*x2[i-2]
                +1.0*norm_rand1();

//        x1[i] =  ( 0.6)*x1[i-1]
//                +(0.65)*x2[i-2]
//                + 1.0  *norm_rand();
//        x2[i] =  ( 0.5)*x2[i-1]
//                +(-0.3)*x2[i-2]
//                +(-0.3)*x3[i-4]
//                +( 0.6)*x4[i-1]
//                + 1.0  *norm_rand1();
//        x3[i] =  ( 0.8)*x3[i-1]
//                +(-0.7)*x3[i-2]
//                +(-0.1)*x5[i-3]
//                + 1.0  *norm_rand2();
//        x4[i] =  ( 0.5)*x4[i-1]
//                +( 0.9)*x3[i-2]
//                +( 0.4)*x5[i-2]
//                + 1.0  *norm_rand3();
//        x5[i] =  ( 0.7)*x5[i-1]
//                +(-0.5)*x5[i-2]
//                +(-0.2)*x3[i-1]
//                + 1.0  *norm_rand4();
        tst_stream1 << QString("%1").arg(x1[i], 15, 'E', 8, ' ') << "\n";
        tst_stream2 << QString("%1").arg(x2[i], 15, 'E', 8, ' ') << "\n";
//        tst_stream3 << QString("%1").arg(x3[i], 15, 'E', 8, ' ') << "\n";
//        tst_stream4 << QString("%1").arg(x4[i], 15, 'E', 8, ' ') << "\n";
//        tst_stream5 << QString("%1").arg(x5[i], 15, 'E', 8, ' ') << "\n";
    }
    tst_file1.close();
    tst_file2.close();
    tst_file3.close();
    tst_file4.close();
    tst_file5.close();
}

/********************************************************
  * Full parameter list:
  * - files to analyse
  * - calc only AR models
  * - dimension from, dimenstion to, dimension step
  * - shift from, shift to, shift step
  * - frequency from, frequency to, frequency
  *******************************************************/

int main(int argc, char *argv[]) {

#ifdef QT_GUI_LIB
    QApplication a(argc, argv);
#else
    QCoreApplication a(argc, argv);
#endif

//    generateData();

    if (a.arguments().count() == 1) {
#ifdef QT_GUI_LIB
        MainWindow mw;
        mw.show();
        return a.exec();
#endif
    }
    try {
        std::vector<string> filesWithData;
        bool calcOnlyAr = false;
        int dimFrom = 1;
        int dimTo = 0;
        int dimStep = 0;
        int shiftFrom = 0;
        int shiftTo = 0;
        int shiftStep = 0;
        int window = 0;
        int dataFrom = 0;
        int dataTo = 0;
        int dataStep = 0;

        po::options_description desc("Allowed options");
        desc.add_options()
                ("help", "produce help message")
                ("input-files,I", po::value< std::vector<string> >(&filesWithData),
                 "files with data to analyse")
                ("ar-only", po::value<bool>(&calcOnlyAr)->default_value(0),
                 "calculate only AR models")
                ("d-from", po::value<int>(&dimFrom)->default_value(1),
                 "dimension of the AR models (start value)")
                ("d-to", po::value<int>(&dimTo)->default_value(0),
                 "dimension of the AR models (end value)")
                ("d-step", po::value<int>(&dimStep)->default_value(0),
                 "dimension of the AR models (step value)")
                ("shift-from", po::value<int>(&shiftFrom)->default_value(0),
                 "time shift of the timeseries (start value)")
                ("shift-to", po::value<int>(&shiftTo)->default_value(0),
                 "time shift of the timeseries (end value)")
                ("shift-step", po::value<int>(&shiftStep)->default_value(0),
                 "time shift of the timeseries (step value)")
                ("window", po::value<int>(&window)->default_value(0),
                 "window in witch analyse input data files")
                ("data-from", po::value<int>(&dataFrom)->default_value(0),
                 "data point in timeseries (start value)")
                ("data-to", po::value<int>(&dataTo)->default_value(0),
                 "data point in timeseries (end value)")
                ("data-step", po::value<int>(&dataStep)->default_value(0),
                 "data point in timeseries (step value)")
                ("daemon-mode", "go into daemon mode")
                ;

        po::positional_options_description p;
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << "Usage: pdcf [options]\n";
            cout << desc;
            return 0;
        }

        if (vm.count("daemon-mode")) {
            pid_t pid, sid;
            unsigned int fd;
            struct rlimit flim;
            pid = fork();
            if (pid < 0) {
                exit(-1);
            }
            if (pid > 0) {
                exit(0);
            }
            sid = setsid();
            if (sid < 0) {
                exit(-1);
            }
            getrlimit( RLIMIT_NOFILE, &flim );
            for ( fd = 0; fd < flim.rlim_max; fd++ )
                close( fd );
        }

        QStringList filesWithDataSL;
        for (int i = 0; i < filesWithData.size(); i++) {
            filesWithDataSL << QString::fromStdString(filesWithData[i]);
        }
        PdcfShell pdcfShell(filesWithDataSL,
                            calcOnlyAr,
                            dimFrom,
                            dimTo,
                            dimStep,
                            shiftFrom,
                            shiftTo,
                            shiftStep,
                            window,
                            dataFrom,
                            dataTo,
                            dataStep);

        pdcfShell.startCalc();
        return a.exec();
    }
    catch(std::exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }
}
