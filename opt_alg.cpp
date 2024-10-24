#include "opt_alg.h"

solution
MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        while (true) {
            Xopt = rand_mat(N);
            for (int i = 0; i < N; ++i)
                Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
            Xopt.fit_fun(ff, ud1, ud2);
            if (Xopt.y < epsilon) {
                Xopt.flag = 1;
                break;
            }
            if (solution::f_calls > Nmax) {
                Xopt.flag = 0;
                break;
            }
        }
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution MC(...):\n" + ex_info);
    }
}

double *
expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2) {
    try {
        auto *p = new double[2]{0, 0};
        int i = 0;

        solution X0(x0), X1(x0 + d);
        X0.fit_fun(ff, ud1, ud2);
        X1.fit_fun(ff, ud1, ud2);

        if (X0.y == X1.y) {
            p[0] = m2d(X0.x);
            p[1] = m2d(X1.x);
            return p;
        }

        if (X1.y > X0.y) {
            d *= -1;
            X1.x = X0.x + d;
            X1.fit_fun(ff, ud1, ud2);
            if (X1.y >= X0.y) {
                p[0] = m2d(X1.x);
                p[1] = m2d(X0.x - d);
                return p;
            }
        }

        solution X2;
        do {
            X2.x = X0.x + pow(alpha, i) * d;
            X2.fit_fun(ff, ud1, ud2);

            if (solution::f_calls > Nmax) {
                p[0] = m2d(X0.x);
                p[1] = m2d(X2.x);
                return p;
            }

            if (X2.y >= X1.y)
                break;

            X0 = X1;
            X1 = X2;
            i++;
        } while (true);

        if (d > 0) {
            p[0] = m2d(X0.x);
            p[1] = m2d(X2.x);
        } else {
            p[0] = m2d(X2.x);
            p[1] = m2d(X0.x);
        }

        return p;
    }
    catch (string ex_info) {
        throw ("double* expansion(...):\n" + ex_info);
    }
}

//TODO
// dodac obliczanie koncowego parametru b-a (dlugosc przedzialu) dla metody fibonacciego i lagrange'a
// po kazdej iteracji

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        unsigned int k = 1;

        vector<double> fibonacci;
        fibonacci.push_back(1);
        fibonacci.push_back(1);

        while (fibonacci[k] < (b - a) / epsilon) {
            fibonacci.push_back(fibonacci[k - 1] + fibonacci[k]);
            k++;
        }

        double c = b - fibonacci[k - 1] / fibonacci[k] * (b - a);
        solution C(c);
        C.fit_fun(ff, ud1, ud2);

        double d = a + b - c;
        solution D(d);
        D.fit_fun(ff, ud1, ud2);

        for (int i = 0; i < k - 3; i++) {
            if (C.y < D.y)
                b = d;
            else
                a = c;

            c = b - fibonacci[k - i - 2] / fibonacci[k - i - 1] * (b - a);
            C = solution(c);
            C.fit_fun(ff, ud1, ud2);

            d = a + b - c;
            D = solution(d);
            D.fit_fun(ff, ud1, ud2);
        }

        Xopt = solution(c);
        Xopt.fit_fun(ff, ud1, ud2);
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution fib(...):\n" + ex_info);
    }

}

solution
lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1,
    matrix ud2) {
    try {
        solution Xopt;

        int i = 0;

        solution A(a), B(b);
        double c = (a + b) / 2;
        solution C(c);

        A.fit_fun(ff, ud1, ud2);
        B.fit_fun(ff, ud1, ud2);
        C.fit_fun(ff, ud1, ud2);

        double d_prev = INFINITY;
        int iter = 0;

        while (true) {
            double l = m2d(A.y * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y * (pow(A.x(0), 2) - pow(B.x(0), 2)));
            double m = m2d(A.y * (m2d(B.x) - m2d(C.x)) + B.y * (m2d(C.x) - m2d(A.x)) + C.y * (m2d(A.x) - m2d(B.x)));

            if (abs(m) <= epsilon) {
                Xopt = C;
                Xopt.flag = -1;
                return Xopt;
            }

            double d = 0.5 * (l / m);

            if (d < min(m2d(A.x), m2d(B.x)) || d > max(m2d(A.x), m2d(B.x))) {
                Xopt.flag = -2;
                return Xopt;
            }

            solution D(d);
            D.fit_fun(ff, ud1, ud2);

            if (solution::f_calls > Nmax) {
                Xopt = D;
                Xopt.flag = 0;
                return Xopt;
            }

            if (abs(d - d_prev) < gamma) {
                Xopt = D;
                Xopt.flag = 1;
                return Xopt;
            }

            if (d < m2d(C.x)) {
                if (m2d(D.y) < m2d(C.y)) {
                    B = C;
                    C = D;
                } else {
                    A = D;
                }
            } else {
                if (m2d(D.y) < m2d(C.y)) {
                    A = C;
                    C = D;
                } else {
                    B = D;
                }
            }

            double new_interval = fabs(m2d(B.x) - m2d(A.x));
            Xopt.ud.add_row(new_interval);
            cout << "i: " << i << ", interval:" << new_interval << ", ";

            if (new_interval < epsilon) {
                Xopt = C;
                Xopt.flag = 1;
                return Xopt;
            }

            d_prev = d;
            iter++;

            solution::f_calls++;
            i++;
        }
    } catch (...) {
        throw "Error in lag";
    }
}

solution
HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1,
   matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution HJ(...):\n" + ex_info);
    }
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
    try {
        //Tu wpisz kod funkcji

        return XB;
    }
    catch (string ex_info) {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}

solution
Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax,
      matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Rosen(...):\n" + ex_info);
    }
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1,
             matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution pen(...):\n" + ex_info);
    }
}

solution
sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta,
       double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution sym_NM(...):\n" + ex_info);
    }
}

solution
SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution
CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
                matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Newton(...):\n" + ex_info);
    }
}

solution
golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution golden(...):\n" + ex_info);
    }
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Powell(...):\n" + ex_info);
    }
}

solution
EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution EA(...):\n" + ex_info);
    }
}
