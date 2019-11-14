#ifndef __APR_GRADIENT_ALGORITHMS_HXX__
#define __APR_GRADIENT_ALGORITHMS_HXX__

#include <Matrix.hxx>
#include <vector>
#include <functional>

#define DOUBLE_DELTA 1e-6
#define DOUBLE_EPSILON 1e-6
#define MAXITER_GA 1000

namespace apr {

    /**
     * @brief Encapsulates several Gradient-based algorithms.
     */
    class GradientAlgorithms {

    public:

        /**
         * @brief Compute approximated gradient at a point of a given function.
         * @param fn Function.
         * @param x0 Point.
         * @param delta Delta offset.
         * @return std::vector<double> Approximated gradient.
         */
        static std::vector<double> ComputeGradient(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double delta = DOUBLE_DELTA);

        /**
         * @brief Compute approximated Hessian matrix at a point of a given function.
         * @param fn Function.
         * @param x0 Point.
         * @param delta Delta offset.
         * @return Matrix Approximated Hessian matrix.
         */
        static Matrix ComputeHessian(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double delta = DOUBLE_DELTA);

        /**
         * @brief Find local function minimum from a given starting point using gradient descent. Moving alongside the gradient vector is in steps of movement ratio size, which is multiplied by the movement decay factor each step.
         * @param fn Function.
         * @param x0 Starting point.
         * @param moveratio Gradient movement ratio.
         * @param movedecay Gradient movement decay.
         * @param delta Delta offset.
         * @param eps Error tolerance.
         * @param maxiter Maximum iterations.
         * @return std::vector<double> Local function minimum.
         */
        static std::vector<double> GradientDescent(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double moveratio = 1.0, double movedecay = 1.0, double delta = DOUBLE_DELTA, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAXITER_GA);

        /**
         * @brief Find local function minimum from a given starting point using gradient descent. Golden section search is used to find the minimum alongside the gradient vector.
         * @param fn Function.
         * @param x0 Starting point.
         * @param offset Golden section search offset.
         * @param delta Delta offset.
         * @param eps Error tolerance.
         * @param maxiter Maximum iterations.
         * @return std::vector<double> Local function minimum.
         */
        static std::vector<double> GradientDescentGSS(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double offset = 1.0, double delta = DOUBLE_DELTA, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAXITER_GA);

        /**
         * @brief Find local function minimum from a given starting point using Newton-Raphson method. Moving alongside the minimization direction vector is in steps of movement ratio size, which is multiplied by the movement decay factor each step.
         * @param fn Function.
         * @param x0 Starting point.
         * @param moveratio Minimization direction movement ratio.
         * @param movedecay Minimization direction movement decay.
         * @param delta Delta offset.
         * @param eps Error tolerance.
         * @param maxiter Maximum iterations.
         * @return std::vector<double> Local function minimum.
         */
        static std::vector<double> NewtonRaphsonMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double moveratio = 1.0, double movedecay = 1.0, double delta = DOUBLE_DELTA, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAXITER_GA);

        /**
         * @brief Find local function minimum from a given starting point using Newton-Raphson method. Golden section search is used to find the minimum alongside the minimization direction vector.
         * @param fn Function.
         * @param x0 Starting point.
         * @param offset Golden section search offset.
         * @param delta Delta offset.
         * @param eps Error tolerance.
         * @param maxiter Maximum iterations.
         * @return std::vector<double> Local function minimum.
         */
        static std::vector<double> NewtonRaphsonMethodGSS(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double offset = 1.0, double delta = DOUBLE_DELTA, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAXITER_GA);

        /**
         * @brief Find local function minimum from a given starting point using Box method with constraints.
         * @param fn Function.
         * @param x0 Starting point.
         * @param nconstraints Number of implicit constraints (inequalities).
         * @param constraints Array of implicit constraints (inequalities).
         * @param min Explicit minimum constraint.
         * @param max Explicit maximum constraint.
         * @param alpha Reflection factor.
         * @param eps Errror tolerance.
         * @param maxiter Maximum iterations.
         * @return std::vector<double> Local function minimum.
         */
        static std::vector<double> ConstrainedBoxMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, std::size_t nconstraints, const std::function<bool(const std::vector<double>&)>* constraints, const std::vector<double>& min, const std::vector<double>& max, double alpha = 1.3, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAXITER_GA);
        
        /**
         * @brief Find local function minimum from a given starting point using Box method with constraints.
         * @param fn Function.
         * @param x0 Starting point.
         * @param constraints Vector of implicit constraints (inequalities).
         * @param min Explicit minimum constraint.
         * @param max Explicit maximum constraint.
         * @param alpha Reflection factor.
         * @param eps Errror tolerance.
         * @param maxiter Maximum iterations.
         * @return std::vector<double> Local function minimum.
         */
        inline static std::vector<double> ConstrainedBoxMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, const std::vector<std::function<bool(const std::vector<double>&)>>& constraints, const std::vector<double>& min, const std::vector<double>& max, double alpha = 1.3, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAXITER_GA) {
            return ConstrainedBoxMethod(fn, x0, constraints.size(), constraints.data(), min, max, alpha, eps, maxiter);
        }

        /**
         * @brief Find local function minimum from a given starting point using Hooke-Jeeves pattern search method with transformed constraints.
         * @param fn Function.
         * @param x0 Starting point.
         * @param ncequalities Number of constraints (equalities).
         * @param cequalities Array of constraints (equalities).
         * @param ncinequalities Number of constraints (inequalities, greater than or equal).
         * @param cinequalities Array of constraints (inequalities, greater than or equal).
         * @param t Transformation factor.
         * @param ts Transformation factor scaler.
         * @param dxscale Step scale.
         * @param decayfactor Decay factor.
         * @param eps Error tolerance.
         * @param maxiter Maximum iterations.
         * @return std::vector<double> Local function minimum.
         */
        static std::vector<double> ConstrainedTransformationMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, std::size_t ncequalities, const std::function<double(const std::vector<double>&)>* cequalities, std::size_t ncinequalities, const std::function<double(const std::vector<double>&)>* cinequalities, double t = 1.0, double ts = 10.0, double dxscale = 0.5, double decayfactor = 0.5, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAXITER_GA);
        
        /**
         * @brief Find local function minimum from a given starting point using Hooke-Jeeves pattern search method with transformed constraints.
         * @param fn Function.
         * @param x0 Starting point.
         * @param cequalities Vector of constraints (equalities).
         * @param cinequalities Vector of constraints (inequalities, greater than or equal).
         * @param t Transformation factor.
         * @param ts Transformation factor scaler.
         * @param dxscale Step scale.
         * @param decayfactor Decay factor.
         * @param eps Error tolerance.
         * @param maxiter Maximum iterations.
         * @return std::vector<double> Local function minimum.
         */
        inline static std::vector<double> ConstrainedTransformationMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, const std::vector<std::function<double(const std::vector<double>&)>>& cequalities, const std::vector<std::function<double(const std::vector<double>&)>>& cinequalities, double t = 1.0, double ts = 10.0, double dxscale = 0.5, double decayfactor = 0.5, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAXITER_GA) {
            return ConstrainedTransformationMethod(fn, x0, cequalities.size(), cequalities.data(), cinequalities.size(), cinequalities.data(), t, ts, dxscale, decayfactor, eps, maxiter);
        }

        /**
         * @brief Get total number of gradient computations.
         * @return std::size_t Total number of gradient computations.
         */
        inline static std::size_t NumGradientComputations() { return s_NumGradientComputations; }

        /**
         * @brief Get total number of Hessian matrix computations.
         * @return std::size_t Total number of Hessian matrix computations.
         */
        inline static std::size_t NumHessianComputations() { return s_NumHessianComputations; }

        /**
         * @brief Reset all computation counters.
         */
        inline static void ResetComputationCounters() { s_NumGradientComputations = s_NumHessianComputations = 0; }

    private:

        /* Emulate static class by deleting the constructor. */
        GradientAlgorithms() = delete;

        /* Computation counters. */
        static std::size_t s_NumGradientComputations, s_NumHessianComputations;

    };

}

#endif
