#include <GeneticAlgorithm.hxx>
#include <algorithm>
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

static std::vector<double>& operator+=(std::vector<double>& x1, const std::vector<double>& x2) {
    std::size_t len = x1.size();
    assert(x2.size() == len);
    for (std::size_t i = 0; i < len; ++i) {
        x1[i] += x2[i];
    }
    return x1;
}

static std::vector<double>& operator/=(std::vector<double>& x, double k) {
    std::size_t len = x.size();
    for (std::size_t i = 0; i < len; ++i) {
        x[i] /= k;
    }
    return x;
}

using namespace apr;
using GenAlg = GeneticAlgorithm;
using PresP = GeneticAlgorithm::AbstractPresentation::Ptr;
using UnitP = GeneticAlgorithm::AbstractUnit::Ptr;
using BinPres = GeneticAlgorithm::BinaryPresentation;
using BinUnit = GeneticAlgorithm::BinaryUnit;
using FloatPres = GeneticAlgorithm::FloatingPointPresentation;
using FloatUnit = GeneticAlgorithm::FloatingPointUnit;

#define NUM_TRIES 10
//#define NO_STEPS

int main(int argc, char** argv) {
    START_BLOCK("TASK 4") {
        FloatPres* floatPres = new FloatPres(2, { -50.0, -50.0 }, { 150.0, 150.0 });
        PresP pres(floatPres);
        std::size_t bestPopulation;
        double bestMutationChance;
        SECTION("Optimizing population size") {
            double bestFitness = std::numeric_limits<double>::max();
            std::size_t populations[] = { 50, 100, 250, 500, 1000 };
            for (int i = 0; i < sizeof(populations) / sizeof(*populations); ++i) {
                GenAlg genAlg(populations[i], pres, 0.5, 1e5);
                std::vector<double> fitnesses;
                for (int j = 0; j < NUM_TRIES; ++j) {
                    std::vector<double> x = genAlg.Solve(f6);
                    fitnesses.push_back(f6(x));
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                if (medianFitness < bestFitness) {
                    bestPopulation = populations[i];
                    bestFitness = medianFitness;
                }
            }
            std::cout << "    best population: " << bestPopulation << std::endl;
        }
        SECTION("Optimizing mutation chance") {
            double bestFitness = std::numeric_limits<double>::max();
            double chances[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
            for (int i = 0; i < sizeof(chances) / sizeof(*chances); ++i) {
                GenAlg genAlg(bestPopulation, pres, chances[i], 1e5);
                std::vector<double> fitnesses;
                for (int j = 0; j < NUM_TRIES; ++j) {
                    std::vector<double> x = genAlg.Solve(f6);
                    fitnesses.push_back(f6(x));
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                if (medianFitness < bestFitness) {
                    bestMutationChance = chances[i];
                    bestFitness = medianFitness;
                }
            }
            std::cout << "    best mutation chance: " << (100 * bestMutationChance) << "%" << std::endl;
        }
        SECTION("Solving for found parameters") {
            GenAlg::LogFunction logFun = [](const std::vector<double>& x, const std::vector<std::pair<UnitP, double>>&, std::size_t bestIndex, std::size_t numCalls) {
                static std::vector<double> xPrev(2, std::numeric_limits<double>::min());
                if (xPrev != x) {
#ifndef NO_STEPS
                    std::cout << "       x = " << std::setw(30) << x << "      f(x) = " << std::setw(12) << f6(x) << "   numCalls = " << std::setw(10) << numCalls << std::endl;
#endif
                    xPrev = x;
                }
            };
            GenAlg genAlg(bestPopulation, pres, bestMutationChance, 1e5);
            std::vector<double> xmin = genAlg.Solve(f6, logFun);
            std::cout << "    xmin = " << std::setw(30) << xmin << "   f(xmin) = " << std::setw(12) << f6(xmin) << "   numCalls = " << std::setw(10) << genAlg.GetLastSolveNumCalls() << std::endl;
        }
    }
    return 0;
}
