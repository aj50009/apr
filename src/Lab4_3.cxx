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

int main(int argc, char** argv) {
    START_BLOCK("TASK 3") {
        std::uint8_t numberOfBitsPerGene = BinPres::GetMinimumNumberOfBitsPerGene(1, { -50.0 }, { 150.0 }, 4);
        SECTION("Number of variables: 3") {
            SUBSECTION("f6, bin") {
                BinPres* binPres = new BinPres(3, { -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0 }, numberOfBitsPerGene);
                PresP pres(binPres);
                GenAlg genAlg(500, pres, 0.5, 1e5);
                int total = 0;
                std::vector<double> fitnesses;
                for (int i = 0; i < NUM_TRIES; ++i) {
                    std::vector<double> x = genAlg.Solve(f6);
                    double fitness = f6(x);
                    fitnesses.push_back(fitness);
                    if (fitness < EPSILON)
                        total++;
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                std::cout << "    hitrate: " << (100 * total) / NUM_TRIES << "%" << std::endl;
                std::cout << "    median fitness: " << medianFitness << std::endl;
            }
            SUBSECTION("f6, float") {
                FloatPres* floatPres = new FloatPres(3, { -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0 });
                PresP pres(floatPres);
                GenAlg genAlg(500, pres, 0.5, 1e5);
                int total = 0;
                std::vector<double> fitnesses;
                for (int i = 0; i < NUM_TRIES; ++i) {
                    std::vector<double> x = genAlg.Solve(f6);
                    double fitness = f6(x);
                    fitnesses.push_back(fitness);
                    if (fitness < EPSILON)
                        total++;
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                std::cout << "    hitrate: " << (100 * total) / NUM_TRIES << "%" << std::endl;
                std::cout << "    median fitness: " << medianFitness << std::endl;
            }
            SUBSECTION("f7, bin") {
                BinPres* binPres = new BinPres(3, { -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0 }, numberOfBitsPerGene);
                PresP pres(binPres);
                GenAlg genAlg(500, pres, 0.5, 1e5);
                int total = 0;
                std::vector<double> fitnesses;
                for (int i = 0; i < NUM_TRIES; ++i) {
                    std::vector<double> x = genAlg.Solve(f7);
                    double fitness = f7(x);
                    fitnesses.push_back(fitness);
                    if (fitness < EPSILON)
                        total++;
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                std::cout << "    hitrate: " << (100 * total) / NUM_TRIES << "%" << std::endl;
                std::cout << "    median fitness: " << medianFitness << std::endl;
            }
            SUBSECTION("f7, float") {
                FloatPres* floatPres = new FloatPres(3, { -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0 });
                floatPres->SetArithmeticCrossoverFunction();
                PresP pres(floatPres);
                GenAlg genAlg(500, pres, 0.5, 1e5);
                int total = 0;
                std::vector<double> fitnesses;
                for (int i = 0; i < NUM_TRIES; ++i) {
                    std::vector<double> x = genAlg.Solve(f7);
                    double fitness = f7(x);
                    fitnesses.push_back(fitness);
                    if (fitness < EPSILON)
                        total++;
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                std::cout << "    hitrate: " << (100 * total) / NUM_TRIES << "%" << std::endl;
                std::cout << "    median fitness: " << medianFitness << std::endl;
            }
        }
        SECTION("Number of variables: 6") {
            SUBSECTION("f6, bin") {
                BinPres* binPres = new BinPres(6, { -50.0, -50.0, -50.0, -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0, 150.0, 150.0, 150.0 }, numberOfBitsPerGene);
                PresP pres(binPres);
                GenAlg genAlg(500, pres, 0.5, 1e5);
                int total = 0;
                std::vector<double> fitnesses;
                for (int i = 0; i < NUM_TRIES; ++i) {
                    std::vector<double> x = genAlg.Solve(f6);
                    double fitness = f6(x);
                    fitnesses.push_back(fitness);
                    if (fitness < EPSILON)
                        total++;
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                std::cout << "    hitrate: " << (100 * total) / NUM_TRIES << "%" << std::endl;
                std::cout << "    median fitness: " << medianFitness << std::endl;
            }
            SUBSECTION("f6, float") {
                FloatPres* floatPres = new FloatPres(6, { -50.0, -50.0, -50.0, -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0, 150.0, 150.0, 150.0 });
                PresP pres(floatPres);
                GenAlg genAlg(500, pres, 0.5, 1e5);
                int total = 0;
                std::vector<double> fitnesses;
                for (int i = 0; i < NUM_TRIES; ++i) {
                    std::vector<double> x = genAlg.Solve(f6);
                    double fitness = f6(x);
                    fitnesses.push_back(fitness);
                    if (fitness < EPSILON)
                        total++;
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                std::cout << "    hitrate: " << (100 * total) / NUM_TRIES << "%" << std::endl;
                std::cout << "    median fitness: " << medianFitness << std::endl;
            }
            SUBSECTION("f7, bin") {
                BinPres* binPres = new BinPres(6, { -50.0, -50.0, -50.0, -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0, 150.0, 150.0, 150.0 }, numberOfBitsPerGene);
                PresP pres(binPres);
                GenAlg genAlg(500, pres, 0.5, 1e5);
                int total = 0;
                std::vector<double> fitnesses;
                for (int i = 0; i < NUM_TRIES; ++i) {
                    std::vector<double> x = genAlg.Solve(f7);
                    double fitness = f7(x);
                    fitnesses.push_back(fitness);
                    if (fitness < EPSILON)
                        total++;
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                std::cout << "    hitrate: " << (100 * total) / NUM_TRIES << "%" << std::endl;
                std::cout << "    median fitness: " << medianFitness << std::endl;
            }
            SUBSECTION("f7, float") {
                FloatPres* floatPres = new FloatPres(6, { -50.0, -50.0, -50.0, -50.0, -50.0, -50.0 }, { 150.0, 150.0, 150.0, 150.0, 150.0, 150.0 });
                PresP pres(floatPres);
                GenAlg genAlg(500, pres, 0.5, 1e5);
                int total = 0;
                std::vector<double> fitnesses;
                for (int i = 0; i < NUM_TRIES; ++i) {
                    std::vector<double> x = genAlg.Solve(f7);
                    double fitness = f7(x);
                    fitnesses.push_back(fitness);
                    if (fitness < EPSILON)
                        total++;
                }
                std::sort(fitnesses.begin(), fitnesses.end());
                double medianFitness = fitnesses[fitnesses.size()];
                if (fitnesses.size() % 2 == 0)
                    medianFitness = (fitnesses[fitnesses.size() / 2 - 1] + medianFitness) / 2;
                std::cout << "    hitrate: " << (100 * total) / NUM_TRIES << "%" << std::endl;
                std::cout << "    median fitness: " << medianFitness << std::endl;
            }
        }
    }
    return 0;
}
