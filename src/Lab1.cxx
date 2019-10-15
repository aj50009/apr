#include <Matrix.hxx>
#include <iostream>
#include <fstream>
#include <cassert>

#define SEPARATOR "----------------------------------------"
#define START_BLOCK(x) { std::cout << SEPARATOR << std::endl << " " << (x) << std::endl << SEPARATOR << std::endl; }
#define END_BLOCK() { std::cout << std::endl; }

using namespace apr;

void try_solve_lu_lup(const Matrix& A, const Matrix& b) {
    std::cout << "LU:" << std::endl;
    try {
        Matrix decomposition = A.DecompositionLU();
        std::cout << ~(decomposition.BackSupstitution(decomposition.ForwardSupstitution(b))) << std::endl;
    } catch (const std::exception& ex) {
        std::cout << "EXCEPTION [" << ex.what() << "]" << std::endl;
    }
    std::cout << "LUP:" << std::endl;
    try {
        std::size_t permutationsCount;
        Matrix permutationsMatrix;
        Matrix decomposition = A.DecompositionLUP(permutationsCount, permutationsMatrix);
        decomposition = permutationsMatrix * decomposition;
        std::cout << ~(decomposition.BackSupstitution(decomposition.ForwardSupstitution(permutationsMatrix * b))) << std::endl;
    } catch (const std::exception& ex) {
        std::cout << "EXCEPTION [" << ex.what() << "]" << std::endl;
    }
}

int main(int argc, char** argv) {

    START_BLOCK("TASK 1") {
        const double scalar = 1.21e-12;
        Matrix A(2, 2);
        A[0][0] = 3.51252123e-5;
        A[0][1] = 1.33213521e-5;
        A[1][0] = -4.21123514e-5;
        A[1][1] = -1.74123411e-5;
        Matrix B = A * scalar;
        Matrix C = B / scalar;
        std::cout << "[C - A]:" << std::endl << (C - A) << std::endl;
        std::cout << "[C == A]: " << ((C == A) ? "TRUE" : "FALSE") << std::endl;
    } END_BLOCK()

    START_BLOCK("TASK 2") {
        Matrix A, b;
        std::ifstream ifs("Lab1_Task2.txt");
        assert(ifs.good());
        ifs >> A >> b;
        try_solve_lu_lup(A, b);
    } END_BLOCK()

    START_BLOCK("TASK 3") {
        Matrix A, b;
        std::ifstream ifs("Lab1_Task3.txt");
        assert(ifs.good());
        ifs >> A >> b;
        try_solve_lu_lup(A, b);
    } END_BLOCK()

    START_BLOCK("TASK 4") {
        Matrix A, b;
        std::ifstream ifs("Lab1_Task4.txt");
        assert(ifs.good());
        ifs >> A >> b;
        try_solve_lu_lup(A, b);
    } END_BLOCK()

    START_BLOCK("TASK 5") {
        Matrix A, b;
        std::ifstream ifs("Lab1_Task5.txt");
        assert(ifs.good());
        ifs >> A >> b;
        try_solve_lu_lup(A, b);
    } END_BLOCK()

    START_BLOCK("TASK 6") {
        const double scalar = 1e10;
        Matrix A, b;
        std::ifstream ifs("Lab1_Task6.txt");
        assert(ifs.good());
        ifs >> A >> b;
        try_solve_lu_lup(A, b);
        std::cout << "{ A *= " << scalar << "; b *= " << scalar << "; }" << std::endl;
        A *= scalar;
        b *= scalar;
        try_solve_lu_lup(A, b);
    } END_BLOCK()

    START_BLOCK("TASK 7") {
        Matrix A("Lab1_Task7.txt");
        try {
            std::cout << "A^-1 = " << std::endl << A.Inverse() << std::endl;
        } catch (const std::exception& ex) {
            std::cout << "EXCEPTION [" << ex.what() << "]" << std::endl;
        }
    } END_BLOCK()

    START_BLOCK("TASK 8") {
        Matrix A("Lab1_Task8.txt");
        try {
            std::cout << "A^-1 = " << std::endl << A.Inverse() << std::endl;
        } catch (const std::exception& ex) {
            std::cout << "EXCEPTION [" << ex.what() << "]" << std::endl;
        }
    } END_BLOCK()

    START_BLOCK("TASK 9") {
        Matrix A("Lab1_Task9.txt");
        try {
            std::cout << "det(A) = " << A.Determinant() << std::endl;
        } catch (const std::exception& ex) {
            std::cout << "EXCEPTION [" << ex.what() << "]" << std::endl;
        }
    } END_BLOCK()

    START_BLOCK("TASK 10") {
        Matrix A("Lab1_Task10.txt");
        try {
            std::cout << "det(A) = " << A.Determinant() << std::endl;
        } catch (const std::exception& ex) {
            std::cout << "EXCEPTION [" << ex.what() << "]" << std::endl;
        }
    } END_BLOCK()

    return 0;
}
