#include <QtGui>

#include "mainwindow.h"

#include <boost/random.hpp>
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

//int loadCoeffs(QString fileName) {
//    // на входе ожидается текстовый файл с набором матриц из
//    // p коэффициентов моделей, полученных при анализе связей
//    // между N компонентами.
//
//    QFile coeffsFile(fileName);
//    if (!coeffsFile.exists()) {
//        qDebug() << "unable to load coeffs file";
//        return -1;
//    }
//
//    // по первой строке файла определяем число компонент
//    coeffsFile.open(QIODevice::ReadOnly);
//    QString str = coeffsFile.readLine();
//    QStringList coeffs = str.split(" ", QString::SkipEmptyParts);
//    N = coeffs.count();
//
//    // заполняем строчку первой матрицы уже считанными коэффициентами
//    matrix<double> Ar_tmp(N, N);
//    for (int i = 0; i < coeffs.size(); i++) {
//        bool ok;
//        Ar_tmp(0, i) = coeffs.at(i).toDouble(&ok);
//        if (!ok) {
//            qDebug() << "an ugly symbol " << coeffs.at(i) << " at line " << 0 << " at " << i << " column";
//            return -1;
//        }
//    }
//
//    // продолжаем считывать файл, заполняя матрицы и список матриц
//    int currentRow = 1;
//    while (!coeffsFile.atEnd()) {
//        str = coeffsFile.readLine();
//        coeffs = str.split(" ", QString::SkipEmptyParts);
//        if (coeffs.count() < N) {
//            // дошли до пустой строки, матрица кончилась, увеличиваем счетчик
//            p++;
//            // и добавляем уже полную матрицу в конец списка
//            Ar.append(Ar_tmp);
//            currentRow = 0;
//            continue;
//        }
//        // заполняем в текущей матрице новую строку с коэффициентами
//        for (int i = 0; i < N; i++) {
//            bool ok;
//            Ar_tmp(currentRow, i) = coeffs.at(i).toDouble(&ok);
//            if (!ok) {
//                qDebug() << "an ugly symbol " << coeffs.at(i) << " at line " << (currentRow*(1+p)+p) << " at " << i << " column";
//                return -1;
//            }
//        }
//        currentRow++;
//    }
//    return 0;
//}

void generateData() {
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
    x1 << 1.0*norm_rand() << 1.0*norm_rand() << 1.0*norm_rand() << 1.0*norm_rand();
    x1.resize(50000);
    QVector<double> x2;
    x2 << 1.0*norm_rand1() << 1.0*norm_rand1() << 1.0*norm_rand1() << 1.0*norm_rand1();
    x2.resize(50000);
    QVector<double> x3;
    x3 << 1.0*norm_rand2() << 1.0*norm_rand2() << 1.0*norm_rand2() << 1.0*norm_rand2();
    x3.resize(50000);
    QVector<double> x4;
    x4 << 1.0*norm_rand3() << 1.0*norm_rand3() << 1.0*norm_rand3() << 1.0*norm_rand3();
    x4.resize(50000);
    QVector<double> x5;
    x5 << 0.0 << 0.0 << 0.0 << 0.0;
    x5.resize(50000);
    tst_stream1 << x1.at(0) << "\n";
    tst_stream1 << x1.at(1) << "\n";
    tst_stream1 << x1.at(2) << "\n";
    tst_stream1 << x1.at(3) << "\n";
    tst_stream2 << x2.at(0) << "\n";
    tst_stream2 << x2.at(1) << "\n";
    tst_stream2 << x2.at(2) << "\n";
    tst_stream2 << x2.at(3) << "\n";
    tst_stream3 << x3.at(0) << "\n";
    tst_stream3 << x3.at(1) << "\n";
    tst_stream3 << x3.at(2) << "\n";
    tst_stream3 << x3.at(3) << "\n";
    tst_stream4 << x4.at(0) << "\n";
    tst_stream4 << x4.at(1) << "\n";
    tst_stream4 << x4.at(2) << "\n";
    tst_stream4 << x4.at(3) << "\n";
    tst_stream5 << x5.at(0) << "\n";
    tst_stream5 << x5.at(1) << "\n";
    tst_stream5 << x5.at(2) << "\n";
    tst_stream5 << x5.at(3) << "\n";
    for (int i = 4; i < 50000; i++) {
//        x1[i] = 0.5*x1[i-1] + 0.3*x2[i-1] + 0.4*x3[i-1] + 1.0*norm_rand();
//        x2[i] = -0.5*x1[i-1] + 0.3*x2[i-1] + 1.0*x3[i-1] + 0.0*norm_rand1();
//        x3[i] = -0.3*x2[i-1] - 0.2*x3[i-1] + 1.0*norm_rand2();

//        x1[i] = (1.4959)*x1[i-1]+(-0.67032)*x1[i-2]
//                +1.0*norm_rand();
//        x2[i] = (1.4959)*x2[i-1]+(-0.67032)*x2[i-2]
//                +1.0*norm_rand1();

        x1[i] = (0.6)*x1[i-1]
                +(0.65)*x2[i-2]
                +1.0*norm_rand();
        x2[i] = (0.5)*x2[i-1]
                +(-0.3)*x2[i-2]
                +(-0.3)*x3[i-4]
                +(0.6)*x4[i-1]
                +1.0*norm_rand1();
        x3[i] = (0.8)*x3[i-1]
                +(-0.7)*x3[i-2]
                +(-0.1)*x5[i-3]
                +1.0*norm_rand2();
        x4[i] = (0.5)*x4[i-1]
                +(0.9)*x3[i-2]
                +(0.4)*x5[i-2]
                +1.0*norm_rand3();
        x5[i] = (0.7)*x5[i-1]
                +(-0.5)*x5[i-2]
                +(-0.2)*x3[i-1]
                +1.0*norm_rand4();
        tst_stream1 << QString("%1").arg(x1[i], 15, 'E', 8, ' ') << "\n";
        tst_stream2 << QString("%1").arg(x2[i], 15, 'E', 8, ' ') << "\n";
        tst_stream3 << QString("%1").arg(x3[i], 15, 'E', 8, ' ') << "\n";
        tst_stream4 << QString("%1").arg(x4[i], 15, 'E', 8, ' ') << "\n";
        tst_stream5 << QString("%1").arg(x5[i], 15, 'E', 8, ' ') << "\n";
    }
    tst_file1.close();
    tst_file2.close();
    tst_file3.close();
    tst_file4.close();
    tst_file5.close();
}

int main(int argc, char *argv[]) {

    QApplication a(argc, argv);

    generateData();

    MainWindow mw;
    mw.show();

    return a.exec();
}
