#include <GradientAlgorithms.hxx>
#include <iostream>
#include <iomanip>

#define SEPARATOR "----------------------------------------"
#define START_BLOCK(x) { std::cout << SEPARATOR << std::endl << " " << (x) << std::endl << SEPARATOR << std::endl; }
#define END_BLOCK() { std::cout << std::endl; }

#define TRY(fn, text, expr) {                                                                                                                                                                                                \
    ncalls = 0;                                                                                                                                                                                                              \
    GradientAlgorithms::ResetComputationCounters();                                                                                                                                                                          \
    std::vector<double> xmin;                                                                                                                                                                                                \
    std::cout << "(~> " #fn " <~) " << (text) << " -->  ";                                                                                                                                                                   \
    try {                                                                                                                                                                                                                    \
        xmin = (expr);                                                                                                                                                                                                       \
        std::cout << "xmin = (" << xmin[0];                                                                                                                                                                                  \
        for (std::size_t i = 1; i < xmin.size(); ++i)                                                                                                                                                                        \
            std::cout << ", " << xmin[i];                                                                                                                                                                                    \
        std::cout << "); fmin = " << (fn)(xmin) << "; ncalls = " << ncalls << "; ngrad = " << GradientAlgorithms::NumGradientComputations() << "; nhesse = " << GradientAlgorithms::NumHessianComputations() << std::endl;   \
    } catch (const std::exception& ex) {                                                                                                                                                                                     \
        std::cout << "FAILED (REASON: " << ex.what() << ")" << std::endl;                                                                                                                                                    \
    }                                                                                                                                                                                                                        \
}

static std::size_t ncalls = 0;

/* Rosenbrock 'banana' function; x0 = (-1.9, 2); xmin = (1, 1); fmin = 0 */
static double f1(const std::vector<double>& x) {
    ncalls++;
    return 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1.0 - x[0]) * (1.0 - x[0]);
}

/* x0 = (0.1, 0.3); xmin = (4, 2); fmin = 0 */
static double f2(const std::vector<double>& x) {
    ncalls++;
    return (x[0] - 4.0) * (x[0] - 4.0) + 4.0 * (x[1] - 2.0) * (x[1] - 2.0);
}

/* x0 = (0, 0); xmin = (2, -3); fmin = 0 */
static double f3(const std::vector<double>& x) {
    ncalls++;
    return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] + 3.0) * (x[1] + 3.0);
}

/* x0 = (0, 0); xmin = (3, 0); fmin = 0 */
static double f4(const std::vector<double>& x) {
    ncalls++;
    return (x[0] - 3.0) * (x[0] - 3.0) + x[1] * x[1];
}

using namespace apr;

int main(int argc, char** argv) {

    START_BLOCK("TASK 1") {
        std::vector<double> x0({ 0.0, 0.0 });
        TRY(f3, "Gradient descent [without GSS]", GradientAlgorithms::GradientDescent(f3, x0))
        TRY(f3, "Gradient descent [with GSS]   ", GradientAlgorithms::GradientDescentGSS(f3, x0))
    } END_BLOCK()

    START_BLOCK("TASK 2") {
        std::vector<double> x0({ -1.9, 2.0 });
        TRY(f1, "Gradient descent      [with GSS] ", GradientAlgorithms::GradientDescentGSS(f1, x0, 1.0, DOUBLE_DELTA, 3e-3, 1e4))
        TRY(f1, "Newton-Raphson method [with GSS] ", GradientAlgorithms::NewtonRaphsonMethodGSS(f1, x0, 1.0, DOUBLE_DELTA, 3e-3, 1e4))

        x0 = std::vector<double>({ 0.1, 0.3 });
        TRY(f2, "Gradient descent      [with GSS] ", GradientAlgorithms::GradientDescentGSS(f2, x0))
        TRY(f2, "Newton-Raphson method [with GSS] ", GradientAlgorithms::NewtonRaphsonMethodGSS(f2, x0))
    } END_BLOCK()

    START_BLOCK("TASK 3") {
        std::vector<std::function<bool(const std::vector<double>&)>> constraints({
            [](const std::vector<double>& x) { return x[1] - x[0] >= 0; },
            [](const std::vector<double>& x) { return 2.0 - x[0] >= 0; }
        });
        std::vector<double> min({ -100.0, -100.0 }), max({ 100.0, 100.0 });

        std::vector<double> x0({ -1.9, 2.0 });
        TRY(f1, "Constrained Box method", GradientAlgorithms::ConstrainedBoxMethod(f1, x0, constraints, min, max))

        x0 = std::vector<double>({ 0.1, 0.3 });
        TRY(f2, "Constrained Box method", GradientAlgorithms::ConstrainedBoxMethod(f2, x0, constraints, min, max))
    } END_BLOCK()

    START_BLOCK("TASK 4") {
        std::vector<std::function<double(const std::vector<double>&)>> cequalities;
        std::vector<std::function<double(const std::vector<double>&)>> cinequalities({
            [](const std::vector<double>& x) { return x[1] - x[0]; },
            [](const std::vector<double>& x) { return 2.0 - x[0]; }
        });

        std::vector<double> x0({ -1.9, 2.0 });
        TRY(f1, "Constrained transformation method", GradientAlgorithms::ConstrainedTransformationMethod(f1, x0, cequalities, cinequalities))

        x0 = std::vector<double>({ 0.1, 0.3 });
        TRY(f2, "Constrained transformation method", GradientAlgorithms::ConstrainedTransformationMethod(f2, x0, cequalities, cinequalities))

        std::cout << "~~~ New starting point for f1 = (0, 10) ~~~" << std::endl;
        x0 = std::vector<double>({ 0.0, 10.0 });
        TRY(f1, "Constrained transformation method", GradientAlgorithms::ConstrainedTransformationMethod(f1, x0, cequalities, cinequalities))
        /* POSSIBLE! Global minimum (1, 1) is within the constrained area. */

        std::cout << "~~~ New starting point for f2 = (3, 4) ~~~" << std::endl;
        x0 = std::vector<double>({ 3.0, 4.0 });
        TRY(f2, "Constrained transformation method", GradientAlgorithms::ConstrainedTransformationMethod(f2, x0, cequalities, cinequalities))
        /* POSSIBLE! Global minimum (4, 2) is outside of the constrained area, therefore (2, 2) is the actual global minimum within the constrained area. */
    } END_BLOCK()

    START_BLOCK("TASK 5") {
        std::vector<std::function<double(const std::vector<double>&)>> cequalities({
            [](const std::vector<double>& x) { return x[1] - 1.0; }
        });
        std::vector<std::function<double(const std::vector<double>&)>> cinequalities({
            [](const std::vector<double>& x) { return 3.0 - x[0] - x[1]; },
            [](const std::vector<double>& x) { return 3.0 + 1.5 * x[0] - x[1]; }
        });

        std::cout << "~~~ New starting point for f4 = (5, 5) ~~~" << std::endl;
        std::vector<double> x0({ 5.0, 5.0 });
        TRY(f4, "Constrained transformation method", GradientAlgorithms::ConstrainedTransformationMethod(f4, x0, cequalities, cinequalities))
    } END_BLOCK()

    return 0;

}
