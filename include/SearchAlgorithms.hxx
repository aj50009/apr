#ifndef __SEARCH_ALGORITHMS_HXX__
#define __SEARCH_ALGORITHMS_HXX__

#include <vector>
#include <functional>

// #define PRINT_STEPS
#define DOUBLE_EPSILON 1e-6
#define MAX_ITERATIONS 1000000

namespace apr {

    /**
     * @brief Encapsulates several search algorithms.
     */
    class SearchAlgorithms {

    public:

        /**
         * @brief Find minimum of an unimodal function interval.
         * @param fn Mathematical function.
         * @param start Start of an unimodal interval.
         * @param end End of an unimodal interval.
         * @param eps Comparison tolerance.
         * @return double Function minimum.
         */
        static double GoldenSectionSearch(const std::function<double(double)>& fn, double start, double end, double eps = DOUBLE_EPSILON);

        /**
         * @brief Find minimum of a function around the given point.
         * @param fn Mathematical function. 
         * @param point Function point.
         * @param offs Search offset.
         * @param eps Comparison tolerance.
         * @return double Function minimum.
         */
        static double GoldenSectionSearchPoint(const std::function<double(double)>& fn, double point, double offs, double eps = DOUBLE_EPSILON);

        /**
         * @brief Find unimodal function interval around the given point.
         * @param fn Mathematical function.
         * @param point Function point.
         * @param offs Search offset.
         * @param start Returned start of an unimodal interval.
         * @param end Returned end of an unimodal interval.
         */
        static void UnimodalIntervalSearch(const std::function<double(double)>& fn, double point, double offs, double& start, double& end);

        /**
         * @brief Find minimum of a multivariable unimodal function around the given point.
         * @param fn Mathematical function.
         * @param point Function point.
         * @param offs Search offset.
         * @param eps Comparison tolerance.
         * @return std::vector<double> Function minimum.
         */
        inline static std::vector<double> CoordinateAxesSearch(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& point, double offs, double eps = DOUBLE_EPSILON) {
            return CoordinateAxesSearch(fn, point, std::vector<double>(point.size(), offs), eps);
        }

        /**
         * @brief Find minimum of a multivariable unimodal function around the given point.
         * @param fn Mathematical function.
         * @param point Function point.
         * @param offs Search offset vector.
         * @param eps Comparison tolerance.
         * @return std::vector<double> Function minimum.
         */
        static std::vector<double> CoordinateAxesSearch(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& point, const std::vector<double>& offs, double eps = DOUBLE_EPSILON);

        /**
         * @brief Find minimum of a multivariable unimodal function around the given point using Nelder-Mead method.
         * @param fn Mathematical function.
         * @param point Function point.
         * @param reflection Reflection scale.
         * @param contraction Contraction scale.
         * @param expansion Expansion scale.
         * @param shrink Shrink scale.
         * @param initscale Scale of the initial right polytope.
         * @param eps Comparison tolerance.
         * @param maxiter Maximum number of iterations.
         * @return std::vector<double> Function minimum.
         */
        static std::vector<double> DownhillSimplexMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& point, double reflection = 1.0, double contraction = 0.5, double expansion = 2.0, double shrink = 0.5, double initscale = 1.0, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAX_ITERATIONS);

        /**
         * @brief Find minimum of a multivariable unimodal function around the given point using Hooke-Jeeves method.
         * @param fn Mathematical function.
         * @param point Function point.
         * @param dxscale Step scale.
         * @param decayfactor Decay factor.
         * @param eps Comparison tolerance.
         * @param maxiter Maximum number of iterations.
         * @return std::vector<double> Function minimum.
         */
        inline static std::vector<double> PatternSearch(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& point, double dxscale = 0.5, double decayfactor = 0.5, double eps = DOUBLE_EPSILON, std::size_t maxiter = MAX_ITERATIONS) {
            return PatternSearch(fn, point, std::vector<double>(point.size(), dxscale), std::vector<double>(point.size(), decayfactor), std::vector<double>(point.size(), eps), maxiter);
        }

        /**
         * @brief Find minimum of a multivariable unimodal function around the given point using Hooke-Jeeves method.
         * @param fn Mathematical function.
         * @param point Function point.
         * @param dx Step vector.
         * @param decay Decay vector.
         * @param epsvec Comparison tolerance vector.
         * @param maxiter Maximum number of iterations.
         * @return std::vector<double> Function minimum.
         */
        static std::vector<double> PatternSearch(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& point, const std::vector<double>& dx, const std::vector<double>& decay, const std::vector<double>& epsvec, std::size_t maxiter = MAX_ITERATIONS);

    private:

        /* Emulate static class by deleting the constructor. */
        SearchAlgorithms() = delete;

    };

}

#endif
