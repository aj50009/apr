#include <SearchAlgorithms.hxx>
#include <Matrix.hxx>
#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

#define SEPARATOR "----------------------------------------"
#define START_BLOCK(x) { std::cout << SEPARATOR << std::endl << " " << (x) << std::endl << SEPARATOR << std::endl; }
#define END_BLOCK() { std::cout << std::endl; }

using namespace apr;

/* Rosenbrock banana function, initial point (-1.9, 2), minimum f(1, 1) = 0 */
static double f1(const std::vector<double>& x) {
    return 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
}

/* Nameless, initial point (0.1, 0.3), minimum f(4, 2) = 0 */
static double f2(const std::vector<double>& x) {
    return (x[0] - 4.0) * (x[0] - 4.0) + 4.0 * (x[1] - 2.0) * (x[1] - 2.0);
}

/* Nameless, initial point (0, 0, 0, ..., 0), minimum f(1, 2, 3, ..., n) = 0 */
static double f3(const std::vector<double>& x) {
    double sum = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
        sum += (x[i] - i - 1.0) * (x[i] - i - 1.0);
    return sum;
}

/* Jakobovic function, initial point (5.1, 1.1), minimum f(0, 0) = 0 */
static double f4(const std::vector<double>& x) {
    double p = x[0] * x[0] - x[1] * x[1];
    return ((p > 0.0) ? p : -p) + std::sqrt(x[0] * x[0] + x[1] * x[1]);
}

/* Schaffer function, minimum f(0, 0, 0, ..., 0) = 0 */
static double f5(const std::vector<double>& x) {
    double s = 0.0;
    for (double xi : x)
        s += xi * xi;
    double p = std::sin(std::sqrt(s));
    double p2 = 1.0 + 0.001 * s;
    return 0.5 + (p * p - 0.5) / (p2 * p2);
}

template<typename T, typename K, typename V>
static void print_table(const T* rows, int nrows, const K* columns, int ncolumns, const V* table, const std::string& rowprefix = "", const std::string& columnprefix = "", const std::string& title = "", int width = 15) {
    std::stringstream wss;
    for (std::size_t k = 0; k < width; ++k)
        wss << "-";
    std::string w = wss.str();
    std::stringstream sepss;
    sepss << "---" << w;
    for (std::size_t j = 0; j < ncolumns; ++j)
        sepss << "--" << w;
    sepss << std::endl;
    std::string sep = sepss.str();
    std::cout << sep;
    std::cout << "|" << std::setw(width) << title << " |";
    for (std::size_t j = 0; j < ncolumns; ++j) {
        std::stringstream ss;
        ss << columnprefix << columns[j];
        std::cout << std::setw(width) << ss.str() << " |";
    }
    std::cout << std::endl;
    std::cout << sep;
    for (std::size_t i = 0; i < nrows; ++i) {
        std::stringstream ss;
        ss << rowprefix << rows[i];
        std::cout << "|" << std::setw(width) << ss.str() << " |";
        for (std::size_t j = 0; j < ncolumns; ++j) {
            std::stringstream sss;
            sss << table[i * ncolumns + j];
            std::cout << std::setw(width) << sss.str() << " |";
        }
        std::cout << std::endl;
        std::cout << sep;
    }
}

static Matrix tomat(const std::vector<double> x) {
    Matrix m(1, x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        m[0][i] = x[i];
    return m;
}

int main(int argc, char** argv) {
    
    START_BLOCK("TASK 1") {
        int counter;
        double y;
        std::function<double(double)> f = [&counter](double x) { counter++; return (x - 3) * (x - 3); };
        std::function<double(const std::vector<double>& x)> fv = [&counter](const std::vector<double>& x) { counter++; return (x[0] - 3) * (x[0] - 3); };
        double rows[] = { 10.0, 25.0, 100.0, 500.0 };
        std::string columns[] = { "GSS", "DSM", "PS" };
        std::vector<double> minimums;
        std::vector<double> fminimums;
        std::vector<int> ncalls;
        minimums.reserve((sizeof(rows) / sizeof(*rows)) * (sizeof(columns) / sizeof(*columns)));
        fminimums.reserve((sizeof(rows) / sizeof(*rows)) * (sizeof(columns) / sizeof(*columns)));
        ncalls.reserve((sizeof(rows) / sizeof(*rows)) * (sizeof(columns) / sizeof(*columns)));
        for (double x : rows) {
            counter = 0;
            y = SearchAlgorithms::GoldenSectionSearchPoint(f, x, 1.0);
            minimums.push_back(y);
            fminimums.push_back(f(y));
            ncalls.push_back(counter);
            counter = 0;
            y = SearchAlgorithms::DownhillSimplexMethod(fv, std::vector<double>(1, x))[0];
            minimums.push_back(y);
            fminimums.push_back(f(y));
            ncalls.push_back(counter);
            counter = 0;
            y = SearchAlgorithms::PatternSearch(fv, std::vector<double>(1, x))[0];
            minimums.push_back(y);
            fminimums.push_back(f(y));
            ncalls.push_back(counter);
        }
        std::cout << std::endl;
        print_table(rows, sizeof(rows) / sizeof(*rows), columns, sizeof(columns) / sizeof(*columns), minimums.data(), "x = ", "", "minimums", 9);
        std::cout << std::endl;
        print_table(rows, sizeof(rows) / sizeof(*rows), columns, sizeof(columns) / sizeof(*columns), fminimums.data(), "x = ", "", "fminimums", 12);
        std::cout << std::endl;
        print_table(rows, sizeof(rows) / sizeof(*rows), columns, sizeof(columns) / sizeof(*columns), ncalls.data(), "x = ", "", "ncalls", 8);
    } END_BLOCK()

    START_BLOCK("TASK 2") {
        int counter;
        std::vector<double> y;
        std::function<double(const std::vector<double>&)> fs[] = {
            [&counter](const std::vector<double>& x) { counter++; return f1(x); },
            [&counter](const std::vector<double>& x) { counter++; return f2(x); },
            [&counter](const std::vector<double>& x) { counter++; return f3(x); },
            [&counter](const std::vector<double>& x) { counter++; return f4(x); }
        };
        std::vector<double> ps[] = {
            std::vector<double>({ -1.9, 2.0 }),
            std::vector<double>({ 0.1, 0.3 }),
            std::vector<double>({ 0.0, 0.0, 0.0, 0.0, 0.0 }),
            std::vector<double>({ 5.1, 1.1 })
        };
        std::string rows[] = { "f1", "f2", "f3", "f4" };
        std::string columns[] = { "CAS", "DSM", "PS" };
        std::vector<Matrix> minimums;
        std::vector<double> fminimums;
        std::vector<int> ncalls;
        minimums.reserve((sizeof(rows) / sizeof(*rows)) * (sizeof(columns) / sizeof(*columns)));
        fminimums.reserve((sizeof(rows) / sizeof(*rows)) * (sizeof(columns) / sizeof(*columns)));
        ncalls.reserve((sizeof(rows) / sizeof(*rows)) * (sizeof(columns) / sizeof(*columns)));
        for (std::size_t i = 0; i < sizeof(rows) / sizeof(*rows); ++i) {
            counter = 0;
            y = SearchAlgorithms::CoordinateAxesSearch(fs[i], ps[i], 2.5);
            minimums.push_back(tomat(y));
            fminimums.push_back(fs[i](y));
            ncalls.push_back(counter);
            counter = 0;
            y = SearchAlgorithms::DownhillSimplexMethod(fs[i], ps[i]);
            minimums.push_back(tomat(y));
            fminimums.push_back(fs[i](y));
            ncalls.push_back(counter);
            counter = 0;
            y = SearchAlgorithms::PatternSearch(fs[i], ps[i]);
            minimums.push_back(tomat(y));
            fminimums.push_back(fs[i](y));
            ncalls.push_back(counter);
        }
        std::cout << std::endl;
        print_table(rows, sizeof(rows) / sizeof(*rows), columns, sizeof(columns) / sizeof(*columns), minimums.data(), "", "", "minimums", 40);
        std::cout << std::endl;
        print_table(rows, sizeof(rows) / sizeof(*rows), columns, sizeof(columns) / sizeof(*columns), fminimums.data(), "", "", "fminimums", 12);
        std::cout << std::endl;
        print_table(rows, sizeof(rows) / sizeof(*rows), columns, sizeof(columns) / sizeof(*columns), ncalls.data(), "", "", "ncalls", 7);
    } END_BLOCK()
    
    START_BLOCK("TASK 3") {
        std::vector<double> x = SearchAlgorithms::DownhillSimplexMethod(f4, std::vector<double>({ 5.0, 5.0 }));
        std::cout << "DSM: f(" << tomat(x) << ") = " << f4(x) << std::endl;
        x = SearchAlgorithms::PatternSearch(f4, std::vector<double>({ 5.0, 5.0 }));
        std::cout << "PS:  f(" << tomat(x) << ") = " << f4(x) << std::endl;
    } END_BLOCK()

    START_BLOCK("TASK 4") {
        int counter;
        std::function<double(const std::vector<double>&)> f = [&counter](const std::vector<double>& x) { counter++; return f1(x); };
        std::cout << std::endl << "~~~ POINT (0.5, 0.5) ~~~" << std::endl;
        for (std::size_t i = 1; i <= 20; ++i) {
            counter = 0;
            SearchAlgorithms::DownhillSimplexMethod(f, std::vector<double>({ 0.5, 0.5 }), 1.0, 0.5, 2.0, 0.5, i);
            std::cout << "scale = " << std::setw(2) << i << ", " << "ncalls = " << counter << std::endl;
        }
        std::cout << std::endl << "~~~ POINT (20, 20) ~~~" << std::endl;
        for (std::size_t i = 1; i <= 20; ++i) {
            counter = 0;
            SearchAlgorithms::DownhillSimplexMethod(f, std::vector<double>({ 20.0, 20.0 }), 1.0, 0.5, 2.0, 0.5, i);
            std::cout << "scale = " << std::setw(2) << i << ", " << "ncalls = " << counter << std::endl;
        }
    } END_BLOCK()
    
    START_BLOCK("TASK 5") {
        int num = 100000;
        std::srand(std::time(nullptr));
        std::size_t count = 0;
        for (std::size_t i = 0; i < num; ++i) {
            std::vector<double> point = std::vector<double>({ (std::rand() % 10000) / 100.0 - 50.0, (std::rand() % 10000) / 100.0 - 50.0 });
            try {
                std::vector<double> x = SearchAlgorithms::DownhillSimplexMethod(f5, point, 1.0, 0.5, 2.0, 0.5, 1.0, 1e-12);
                if (f5(x) < 1e-4)
                    count++;
            } catch (const std::exception& ex) {
                num--;
            }
        }
        std::cout << "Probability of global minimum: " << static_cast<double>(100 * count) / num << "%" << std::endl;
    } END_BLOCK()

    return 0;

}
