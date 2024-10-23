#include "opt_alg.h"

void lab0();

void lab1();

void lab2();

void lab3();

void lab4();

void lab5();

void lab6();

int main() {
    try {
        lab1();
    }
    catch (string EX_INFO) {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    system("pause");
    return 0;
}

void lab0() {
    //Funkcja testowa
    double epsilon = 1e-2;
    int Nmax = 10000;
    matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
    solution opt;
    a(0) = -1;
    a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Wahadlo
    Nmax = 1000;
    epsilon = 1e-2;
    lb = 0;
    ub = 5;
    double teta_opt = 1;
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Zapis symulacji do pliku csv
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(opt.x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
    ofstream Sout("symulacja_lab0.csv");
    Sout << hcat(Y[0], Y[1]);
    Sout.close();
    Y[0].~matrix();
    Y[1].~matrix();
}

void lab1() {
    try {
        const string FILE_PATH = R"(C:\Users\Dell\Desktop\stuff\studia\semestr 5\Optymalizacja\project\data\data1\data1.csv)";
        fstream csvFile;

        const int numOfElements = 100;
        double minVal = -1, maxVal = 1;
        set<double> uniqueNums;

        default_random_engine randomEngine(random_device{}());
        uniform_real_distribution<double> realDistribution(minVal, maxVal);

        while (uniqueNums.size() < numOfElements) {
            uniqueNums.insert(realDistribution(randomEngine));
        }

        for (double x0 : uniqueNums) {

            // Expansion
            double d = 0.75, alpha = 1.92;
            int Nmax = 1000;
            unique_ptr<double[]> p(expansion(ff1T, x0, d, alpha, Nmax));

            if (!p)
                continue;

            try {
                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << x0 << "," << p[0] << "," << p[1] << "," << solution::f_calls << ",";
                csvFile.close();
            } catch (const exception& e) {
                cout << "Error writing expansion results: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();

            // Fibonacci
            try {
                solution Fibonacci = fib(ff1T, p[0], p[1], 0.0001);
                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << m2d(Fibonacci.x) << "," << m2d(Fibonacci.y) << "," << solution::f_calls << ",";
                csvFile.close();
            } catch (const exception& e) {
                cout << "Error in Fibonacci: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();

            // Lagrange
            try {
                solution Lagrange = lag(ff1T, p[0], p[1], 0.0001, 1e-09, Nmax);

                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << m2d(Lagrange.x) << "," << m2d(Lagrange.y) << "," << solution::f_calls << "\n";
                csvFile.close();
            } catch (const exception& e) {
                cout << "Error in Lagrange: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();
        }

        cout << "Processing completed successfully" << endl;

    } catch (const exception& e) {
        cout << "Fatal error in lab1: " << e.what() << endl;
    } catch (...) {
        cout << "Unknown error in lab1" << endl;
    }
    exit(0);
}

void lab2() {

}

void lab3() {

}

void lab4() {

}

void lab5() {

}

void lab6() {

}
