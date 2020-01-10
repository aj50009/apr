#include <NumInt.hxx>
#include <iostream>
#include <sstream>

#define SEPARATOR "----------------------------------------"
#define START_BLOCK(x) { std::cout << SEPARATOR << std::endl << " " << (x) << std::endl << SEPARATOR << std::endl; }
#define END_BLOCK() { std::cout << std::endl; }

using namespace apr;

static void solve_all(std::size_t printNum, const std::string& destFile, const std::string& appendFile, const Matrix& A, const Matrix& B, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP, const TimeFunction& r = NumInt::DefaultTimeFunction) {
    std::cout << "Runge-Kutta (RK4):" << std::endl;
    NumInt::EnablePrint(printNum);
    NumInt::SetDestFile(destFile, "_rk4", false);
    Matrix x = NumInt::RungeKuttaIter(A, B, x0, tMax, T, r);
    std::cout << "   SOLUTION: x(10) = [ " << ~x << " ]" << std::endl << std::endl;

    std::cout << "Trapezoidal rule:" << std::endl;
    NumInt::EnablePrint(printNum);
    NumInt::SetDestFile(destFile, "_trap");
    x = NumInt::TrapezoidalRuleIter(A, B, x0, tMax, T, r);
    std::cout << "   SOLUTION: x(10) = [ " << ~x << " ]" << std::endl << std::endl;

    std::cout << "Euler method:" << std::endl;
    NumInt::EnablePrint(printNum);
    NumInt::SetDestFile(destFile, "_euler");
    x = NumInt::EulerIter(A, B, x0, tMax, T, r);
    std::cout << "   SOLUTION: x(10) = [ " << ~x << " ]" << std::endl << std::endl;

    std::cout << "Backward Euler method:" << std::endl;
    NumInt::EnablePrint(printNum);
    NumInt::SetDestFile(destFile, "_inveuler");
    x = NumInt::BackwardEulerIter(A, B, x0, tMax, T, r);
    std::cout << "   SOLUTION: x(10) = [ " << ~x << " ]" << std::endl << std::endl;

    std::cout << "PE(CE)2 method:" << std::endl;
    NumInt::EnablePrint(printNum);
    NumInt::SetDestFile(destFile, "_pece2");
    NumInt::PECE2Iter(A, B, x0, tMax, T, r);
    std::cout << "   SOLUTION: x(10) = [ " << ~x << " ]" << std::endl << std::endl;

    std::cout << "PECE method:" << std::endl;
    NumInt::EnablePrint(printNum);
    NumInt::SetDestFile(destFile, "_pece");
    NumInt::PECEIter(A, B, x0, tMax, T, r);
    std::cout << "   SOLUTION: x(10) = [ " << ~x << " ]" << std::endl;

    NumInt::RedirFileContent(destFile, appendFile);
}

static void solve_all(std::size_t printNum, const std::string& destFile, const std::string& appendFile, const Matrix& A, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP) {
    solve_all(printNum, destFile, appendFile, A, Matrix(2, 2), x0, tMax, T, NumInt::ZeroTimeFunction);
}

int main(int argc, char** argv) {

    START_BLOCK("TASK 1") {
        Matrix A(2, 2);
        A[0][0] = 0.0; A[0][1] = 1.0;
        A[1][0] = -1.0; A[1][1] = 0.0;
        Matrix x0(2, 1);
        x0[0][0] = 1.0;
        x0[1][0] = 1.0;
        solve_all(200, "Lab5_1.m", "Lab5_1.m.base", A, x0, 10.0);
    } END_BLOCK()

    START_BLOCK("TASK 2") {
        Matrix A(2, 2);
        A[0][0] = 0.0; A[0][1] = 1.0;
        A[1][0] = -200.0; A[1][1] = -102.0;
        Matrix x0(2, 1);
        x0[0][0] = 1.0;
        x0[1][0] = -2.0;
        solve_all(2, "Lab5_2.m", "Graphs.m.base", A, x0, 1.0, 0.1);
        { std::ifstream ifs("Lab5_2_RK4.m", std::ios::out | std::ios::trunc); }
        for (std::size_t i = 1; i <= 10; ++i) {
            double T = 0.001 * i;
            std::stringstream ss;
            ss << "_" << i;
            NumInt::SetDestFile("Lab5_2_RK4.m", ss.str());
            NumInt::RungeKuttaIter(A, x0, 1.0, T);
        }
        NumInt::RedirFileContent("Lab5_2_RK4.m", "Lab5_2_RK4.m.base");
    } END_BLOCK()

    START_BLOCK("TASK 3") {
        Matrix A(2, 2);
        A[0][0] = 0.0; A[0][1] = -2.0;
        A[1][0] = 1.0; A[1][1] = -3.0;
        Matrix B(2, 2);
        B[0][0] = 2.0; B[0][1] = 0.0;
        B[1][0] = 0.0; B[1][1] = 3.0;
        Matrix x0(2, 1);
        x0[0][0] = 1.0;
        x0[1][0] = 3.0;
        solve_all(200, "Lab5_3.m", "Graphs.m.base", A, B, x0, 10.0, 0.01, NumInt::OneTimeFunction);
    } END_BLOCK()

    START_BLOCK("TASK 4") {
        Matrix A(2, 2);
        A[0][0] = 1.0; A[0][1] = -5.0;
        A[1][0] = 1.0; A[1][1] = -7.0;
        Matrix B(2, 2);
        B[0][0] = 5.0; B[0][1] = 0.0;
        B[1][0] = 0.0; B[1][1] = 3.0;
        Matrix x0(2, 1);
        x0[0][0] = -1.0;
        x0[1][0] = 3.0;
        solve_all(20, "Lab5_4.m", "Graphs.m.base", A, B, x0, 1.0);
    } END_BLOCK()

    return 0;

}
