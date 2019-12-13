#include <GeneticAlgorithm.hxx>
#include <limits>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>
#include <cstddef>
#include <cassert>
#include <cmath>
#include <ctime>

#define SEPARATOR "----------------------------------------"
#define START_BLOCK(x) { std::cout << SEPARATOR << std::endl << " " << (x) << std::endl << SEPARATOR << std::endl; }
#define END_BLOCK() { std::cout << std::endl; }
#define SECTION(x) { std::cout << (x) << ":" << std::endl; }
#define SUBSECTION(x) { std::cout << "- " << (x) << ":" << std::endl; }

/* Rosenbrock's 'banana' function; x0 = (-1.9, 2); xmin = (1, 1); fmin = 0 */
static double f1(const std::vector<double>& x) {
    assert(x.size() == 2);
    double x1 = x[0], x2 = x[1];
    double x1_2 = x1 * x1;
    return 100.0 * (x2 - x1_2) * (x2 - x1_2) + (1.0 - x1) * (1.0 - x1);
}

/* x0 = (0, 0, 0, ..., 0); xmin = (1, 2, 3, ..., n); fmin = 0 */
static double f3(const std::vector<double>& x) {
    assert(x.size() > 0);
    double sum = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        double expr = x[i] - static_cast<double>(i) - 1.0;
        sum += expr * expr;
    }
    return sum;
}

/* Schaffer's f6 function; xmin = (0, 0, 0, ..., 0); fmin = 0  */
static double f6(const std::vector<double>& x) {
    assert(x.size() > 0);
    double sum = 0.0;
    for (double xi : x)
        sum += xi * xi;
    double sin = std::sin(std::sqrt(sum));
    double expr = 1.0 + 0.001 * sum;
    return 0.5 + (sin * sin - 0.5) / (expr * expr);
}

/* similar to Schaffer's f7 function; xmin = (0, 0, 0, ..., 0); fmin = 0  */
static double f7(const std::vector<double>& x) {
    assert(x.size() > 0);
    double sum = 0.0;
    for (double xi : x)
        sum += xi * xi;
    double sin = std::sin(50.0 * std::pow(sum, 0.1));
    return std::pow(sum, 0.25) * (1.0 + sin * sin);
}

static std::ostream& operator<<(std::ostream& outputStream, const std::vector<double>& x) {
    std::stringstream ss;
    ss << "(";
    for (double xi : x)
        ss << " " << xi;
    ss << " )";
    outputStream << ss.str();
    return outputStream;
}

static bool operator==(const std::vector<double>& x1, const std::vector<double>& x2) {
    std::size_t len = x1.size();
    assert(x2.size() == len);
    for (std::size_t index = 0; index < len; ++index)
        if (std::abs(x1[index] - x2[index]) > EPSILON)
            return false;
    return true;
}

static bool operator!=(const std::vector<double>& x1, const std::vector<double>& x2) {
    return !(x1 == x2);
}

using namespace apr;
using GenAlg = GeneticAlgorithm;
using PresP = GeneticAlgorithm::AbstractPresentation::Ptr;
using UnitP = GeneticAlgorithm::AbstractUnit::Ptr;
using BinPres = GeneticAlgorithm::BinaryPresentation;
using BinUnit = GeneticAlgorithm::BinaryUnit;
using FloatPres = GeneticAlgorithm::FloatingPointPresentation;
using FloatUnit = GeneticAlgorithm::FloatingPointUnit;

//#define NO_STEPS

int main(int argc, char** argv) {
    START_BLOCK("TASK 1") {
        std::uint8_t numberOfBitsPerGene = BinPres::GetMinimumNumberOfBitsPerGene(1, { -50.0 }, { 150.0 }, 7);
        SECTION("f1") {
            GenAlg::LogFunction logFun = [](const std::vector<double>& x, const std::vector<std::pair<UnitP, double>>&, std::size_t bestIndex, std::size_t numCalls) {
                static std::vector<double> xPrev(2, std::numeric_limits<double>::min());
                if (xPrev != x) {
#ifndef NO_STEPS
                    std::cout << "       x = " << std::setw(30) << x << "      f(x) = " << std::setw(12) << f1(x) << "   numCalls = " << std::setw(10) << numCalls << std::endl;
#endif
                    xPrev = x;
                }
            };
            SUBSECTION("bin") {
                BinPres* binPres = new BinPres(2, { -50.0, -50.0 }, { 150.0, 150.0 }, numberOfBitsPerGene);
                PresP pres(binPres);
                GenAlg genAlg(5000, pres);
                std::vector<double> xmin = genAlg.Solve(f1, logFun);
                std::cout << "    xmin = " << std::setw(30) << xmin << "   f(xmin) = " << std::setw(12) << f1(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
            }
            SUBSECTION("float") {
                FloatPres* floatPres = new FloatPres(2, { -50.0, -50.0 }, { 150.0, 150.0 });
                PresP pres(floatPres);
                GenAlg genAlg(100, pres, 0.8);
                std::vector<double> xmin = genAlg.Solve(f1, logFun);
                std::cout << "    xmin = " << std::setw(30) << xmin << "   f(xmin) = " << std::setw(12) << f1(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
            }
        }
        SECTION("f3") {
            GenAlg::LogFunction logFun = [](const std::vector<double>& x, const std::vector<std::pair<UnitP, double>>&, std::size_t bestIndex, std::size_t numCalls) {
                static std::vector<double> xPrev(5, std::numeric_limits<double>::min());
                if (xPrev != x) {
#ifndef NO_STEPS
                    std::cout << "       x = " << std::setw(50) << x << "      f(x) = " << std::setw(12) << f3(x) << "   numCalls = " << std::setw(10) << numCalls << std::endl;
#endif
                    xPrev = x;
                }
            };
            SUBSECTION("bin") {
                BinPres* binPres = new BinPres(5, { -50.0, -50.0, -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0, 150.0, 150.0 }, numberOfBitsPerGene);
                PresP pres(binPres);
                GenAlg genAlg(5000, pres);
                std::vector<double> xmin = genAlg.Solve(f3, logFun);
                std::cout << "    xmin = " << std::setw(50) << xmin << "   f(xmin) = " << std::setw(12) << f3(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
            }
            SUBSECTION("float") {
                FloatPres* floatPres = new FloatPres(5, { -50.0, -50.0, -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0, 150.0, 150.0 });
                PresP pres(floatPres);
                GenAlg genAlg(100, pres, 0.8);
                std::vector<double> xmin = genAlg.Solve(f3, logFun);
                std::cout << "    xmin = " << std::setw(50) << xmin << "   f(xmin) = " << std::setw(12) << f3(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
            }
        }
        SECTION("f6") {
            GenAlg::LogFunction logFun = [](const std::vector<double>& x, const std::vector<std::pair<UnitP, double>>&, std::size_t bestIndex, std::size_t numCalls) {
                static std::vector<double> xPrev(2, std::numeric_limits<double>::min());
                if (xPrev != x) {
#ifndef NO_STEPS
                    std::cout << "       x = " << std::setw(30) << x << "      f(x) = " << std::setw(12) << f6(x) << "   numCalls = " << std::setw(10) << numCalls << std::endl;
#endif
                    xPrev = x;
                }
            };
            SUBSECTION("bin") {
                BinPres* binPres = new BinPres(2, { -50.0, -50.0 }, { 150.0, 150.0 }, numberOfBitsPerGene);
                PresP pres(binPres);
                GenAlg genAlg(100, pres);
                std::vector<double> xmin = genAlg.Solve(f6, logFun);
                std::cout << "    xmin = " << std::setw(30) << xmin << "   f(xmin) = " << std::setw(12) << f6(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
            }
            SUBSECTION("float") {
                FloatPres* floatPres = new FloatPres(2, { -50.0, -50.0 }, { 150.0, 150.0 });
                PresP pres(floatPres);
                GenAlg genAlg(100, pres, 0.8);
                std::vector<double> xmin = genAlg.Solve(f6, logFun);
                std::cout << "    xmin = " << std::setw(30) << xmin << "   f(xmin) = " << std::setw(12) << f6(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
            }
        }
        SECTION("f7") {
            GenAlg::LogFunction logFun = [](const std::vector<double>& x, const std::vector<std::pair<UnitP, double>>&, std::size_t bestIndex, std::size_t numCalls) {
                static std::vector<double> xPrev(2, std::numeric_limits<double>::min());
                if (xPrev != x) {
#ifndef NO_STEPS
                    std::cout << "       x = " << std::setw(30) << x << "      f(x) = " << std::setw(12) << f7(x) << "   numCalls = " << std::setw(10) << numCalls << std::endl;
#endif
                    xPrev = x;
                }
            };
            SUBSECTION("bin") {
                BinPres* binPres = new BinPres(2, { -50.0, -50.0 }, { 150.0, 150.0 }, numberOfBitsPerGene);
                binPres->SetSegmentedCrossoverFunction();
                PresP pres(binPres);
                GenAlg genAlg(1000, pres, 0.9);
                std::vector<double> xmin = genAlg.Solve(f7, logFun);
                std::cout << "    xmin = " << std::setw(30) << xmin << "   f(xmin) = " << std::setw(12) << f7(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
            }
            SUBSECTION("float") {
                FloatPres* floatPres = new FloatPres(2, { -50.0, -50.0 }, { 150.0, 150.0 });
                PresP pres(floatPres);
                GenAlg genAlg(250, pres, 0.8);
                std::vector<double> xmin = genAlg.Solve(f7, logFun);
                std::cout << "    xmin = " << std::setw(30) << xmin << "   f(xmin) = " << std::setw(12) << f7(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
            }
        }
    } END_BLOCK()
    START_BLOCK("TASK 2") {
        int numVars[] = { 1, 3, 6, 10 };
        std::initializer_list<double> lowerBounds[] = {
            { -50.0 },
            { -50.0, -50.0, -50.0 },
            { -50.0, -50.0, -50.0, -50.0, -50.0, -50.0 },
            { -50.0, -50.0, -50.0, -50.0, -50.0, -50.0, -50.0, -50.0, -50.0, -50.0 }
        };
        std::initializer_list<double> upperBounds[] = {
            { 150.0 },
            { 150.0, 150.0, 150.0 },
            { 150.0, 150.0, 150.0, 150.0, 150.0, 150.0 },
            { 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0 }
        };
        int spacingabs[] = { 12, 12, 18, 22 };
        int spacingrel[] = { 10, 10, 10, 8 };
        for (int i = 0; i < sizeof(numVars) / sizeof(*numVars); ++i) {
            std::stringstream ss;
            ss << "Number of variables: " << numVars[i];
            SECTION(ss.str()) {
                SUBSECTION("f6") {
                    std::vector<double> xPrev(numVars[i], std::numeric_limits<double>::min());
                    GenAlg::LogFunction logFun = [&xPrev, spacingabs, spacingrel, numVars, i](const std::vector<double>& x, const std::vector<std::pair<UnitP, double>>&, std::size_t bestIndex, std::size_t numCalls) {
                        if (xPrev != x) {
#ifndef NO_STEPS
                            std::cout << "       x = " << std::setw(spacingabs[i] + spacingrel[i] * numVars[i]) << x << "      f(x) = " << std::setw(12) << f6(x) << "   numCalls = " << std::setw(10) << numCalls << std::endl;
#endif
                            xPrev = x;
                        }
                    };
                    FloatPres* floatPres = new FloatPres(numVars[i], lowerBounds[i], upperBounds[i]);
                    PresP pres(floatPres);
                    GenAlg genAlg(250, pres, 0.8);
                    std::vector<double> xmin = genAlg.Solve(f6, logFun);
                    std::cout << "    xmin = " << std::setw(spacingabs[i] + spacingrel[i] * numVars[i]) << xmin << "   f(xmin) = " << std::setw(12) << f6(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
                }
                SUBSECTION("f7") {
                    std::vector<double> xPrev(numVars[i], std::numeric_limits<double>::min());
                    GenAlg::LogFunction logFun = [&xPrev, spacingabs, spacingrel, numVars, i](const std::vector<double>& x, const std::vector<std::pair<UnitP, double>>&, std::size_t bestIndex, std::size_t numCalls) {
                        if (xPrev != x) {
#ifndef NO_STEPS
                            std::cout << "       x = " << std::setw(spacingabs[i] + spacingrel[i] * numVars[i]) << x << "      f(x) = " << std::setw(12) << f7(x) << "   numCalls = " << std::setw(10) << numCalls << std::endl;
#endif
                            xPrev = x;
                        }
                    };
                    FloatPres* floatPres = new FloatPres(numVars[i], lowerBounds[i], upperBounds[i]);
                    PresP pres(floatPres);
                    GenAlg genAlg(250, pres, 0.8);
                    std::vector<double> xmin = genAlg.Solve(f7, logFun);
                    std::cout << "    xmin = " << std::setw(spacingabs[i] + spacingrel[i] * numVars[i]) << xmin << "   f(xmin) = " << std::setw(12) << f7(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
                }
            }
        }
    } END_BLOCK()
    return 0;
}
