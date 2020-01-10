#ifndef __APR_NUM_INT_HXX__
#define __APR_NUM_INT_HXX__

#include <Matrix.hxx>
#include <functional>
#include <fstream>
#include <cstddef>

#define DEFAULT_TIME_STEP 1e-2

namespace apr {

    /**
     * @brief Time function.
     */
    typedef std::function<Matrix(double)> TimeFunction;

    /**
     * @brief Integration method.
     */
    typedef std::function<Matrix(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)> IntMethod;

    /**
     * @brief Integration method with predicted x matrix.
     */
    typedef std::function<Matrix(const Matrix&, const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)> IntMethod2;

    /**
     * @brief Encapsulates several numerical integration methods.
     */
    class NumInt {

    public:

        /**
         * @brief Find x(t + T) solution of equation x' = Ax + Br(t) using Runge-Kutta method.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x Initial x matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(t + T).
         */
        static Matrix RungeKuttaMethod(const Matrix& A, const Matrix& B, const Matrix& x, double t, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction);

        /**
         * @brief Find x(t + T) solution of equation x' = Ax using Runge-Kutta method.
         * @param A Matrix A.
         * @param x Initial x matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @return Matrix Solution for x(t + T).
         */
        inline static Matrix RungeKuttaMethod(const Matrix& A, const Matrix& x, double t, double T = DEFAULT_TIME_STEP) {
            return RungeKuttaMethod(A, Matrix(2, 2), x, t, T, ZeroTimeFunction);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax + Br(t) iteratively with time step T using Runge-Kutta method.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix RungeKuttaIter(const Matrix& A, const Matrix& B, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction) {
            return ArbitIntIter(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(RungeKuttaMethod), A, B, x0, tMax, T, r);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax iteratively with time step T using Runge-Kutta method.
         * @param A Matrix A.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix RungeKuttaIter(const Matrix& A, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP) {
            return ArbitIntIter(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(RungeKuttaMethod), A, x0, tMax, T);
        }

        /**
         * @brief Find x(t + T) solution of equation x' = Ax + Br(t) using trapezoidal rule.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x Initial x matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(t + T).
         */
        static Matrix TrapezoidalRuleMethod(const Matrix& A, const Matrix& B, const Matrix& x, double t, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction);

        /**
         * @brief Find x(t + T) solution of equation x' = Ax + Br(t) using trapezoidal rule.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x Initial x(t) matrix.
         * @param xp Predicted x(t + T) matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(t + T).
         */
        static Matrix TrapezoidalRuleMethod(const Matrix& A, const Matrix& B, const Matrix& x, const Matrix& xp, double t, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction);

        /**
         * @brief Find x(t + T) solution of equation x' = Ax using trapezoidal rule.
         * @param A Matrix A.
         * @param x Initial x matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @return Matrix Solution for x(t + T).
         */
        inline static Matrix TrapezoidalRuleMethod(const Matrix& A, const Matrix& x, double t, double T = DEFAULT_TIME_STEP) {
            return TrapezoidalRuleMethod(A, Matrix(2, 2), x, t, T, ZeroTimeFunction);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax + Br(t) iteratively with time step T using trapezoidal rule.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix TrapezoidalRuleIter(const Matrix& A, const Matrix& B, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction) {
            return ArbitIntIter(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(TrapezoidalRuleMethod), A, B, x0, tMax, T, r);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax iteratively with time step T using trapezoidal rule.
         * @param A Matrix A.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix TrapezoidalRuleIter(const Matrix& A, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP) {
            return ArbitIntIter(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(TrapezoidalRuleMethod), A, x0, tMax, T);
        }

        /**
         * @brief Find x(t + T) solution of equation x' = Ax + Br(t) using Euler method.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x Initial x matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(t + T).
         */
        static Matrix EulerMethod(const Matrix& A, const Matrix& B, const Matrix& x, double t, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction);

        /**
         * @brief Find x(t + T) solution of equation x' = Ax using Euler method.
         * @param A Matrix A.
         * @param x Initial x matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @return Matrix Solution for x(t + T).
         */
        inline static Matrix EulerMethod(const Matrix& A, const Matrix& x, double t, double T = DEFAULT_TIME_STEP) {
            return EulerMethod(A, Matrix(2, 2), x, t, T, ZeroTimeFunction);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax + Br(t) iteratively with time step T using Euler method.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix EulerIter(const Matrix& A, const Matrix& B, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction) {
            return ArbitIntIter(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(EulerMethod), A, B, x0, tMax, T, r);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax iteratively with time step T using Euler method.
         * @param A Matrix A.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix EulerIter(const Matrix& A, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP) {
            return ArbitIntIter(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(EulerMethod), A, x0, tMax, T);
        }

        /**
         * @brief Find x(t + T) solution of equation x' = Ax + Br(t) using backward Euler method.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x Initial x matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(t + T).
         */
        static Matrix BackwardEulerMethod(const Matrix& A, const Matrix& B, const Matrix& x, double t, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction);
        
        /**
         * @brief Find x(t + T) solution of equation x' = Ax + Br(t) using backward Euler method.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x Initial x(t) matrix.
         * @param xp Predicted x(t + T) matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(t + T).
         */
        static Matrix BackwardEulerMethod(const Matrix& A, const Matrix& B, const Matrix& x, const Matrix& xp, double t, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction);

        /**
         * @brief Find x(t + T) solution of equation x' = Ax using backward Euler method.
         * @param A Matrix A.
         * @param x Initial x matrix.
         * @param t Initial time value.
         * @param T Time step.
         * @return Matrix Solution for x(t + T).
         */
        inline static Matrix BackwardEulerMethod(const Matrix& A, const Matrix& x, double t, double T = DEFAULT_TIME_STEP) {
            return BackwardEulerMethod(A, Matrix(2, 2), x, t, T, ZeroTimeFunction);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax + Br(t) iteratively with time step T using backward Euler method.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix BackwardEulerIter(const Matrix& A, const Matrix& B, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction) {
            return ArbitIntIter(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(BackwardEulerMethod), A, B, x0, tMax, T, r);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax iteratively with time step T using backward Euler method.
         * @param A Matrix A.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix BackwardEulerIter(const Matrix& A, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP) {
            return ArbitIntIter(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(BackwardEulerMethod), A, x0, tMax, T);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax + Br(t) iteratively with time step T using PE(CE)2 method (predictor - Euler, corrector - backward Euler).
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix PECE2Iter(const Matrix& A, const Matrix& B, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction) {
            return ArbitIntIter(s_MethodPECE2, A, B, x0, tMax, T, r);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax iteratively with time step T using PE(CE)2 method (predictor - Euler, corrector - backward Euler).
         * @param A Matrix A.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix PECE2Iter(const Matrix& A, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP) {
            return ArbitIntIter(s_MethodPECE2, A, x0, tMax, T);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax + Br(t) iteratively with time step T using PECE method (predictor - Euler, corrector - trapezoidal rule).
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix PECEIter(const Matrix& A, const Matrix& B, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction) {
            return ArbitIntIter(s_MethodPECE, A, B, x0, tMax, T, r);
        }

        /**
         * @brief Find x(tMax) solution of equation x' = Ax iteratively with time step T using PECE method (predictor - Euler, corrector - trapezoidal rule).
         * @param A Matrix A.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix PECEIter(const Matrix& A, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP) {
            return ArbitIntIter(s_MethodPECE, A, x0, tMax, T);
        }

        /**
         * @brief Arbitrary iterative integration procedure.
         * @param method Arbitrary method.
         * @param A Matrix A.
         * @param B Matrix B.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @param r Time function.
         * @return Matrix Solution for x(tMax).
         */
        static Matrix ArbitIntIter(const IntMethod& method, const Matrix& A, const Matrix& B, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP, const TimeFunction& r = DefaultTimeFunction);

        /**
         * @brief Arbitrary iterative integration procedure.
         * @param method Arbitrary method.
         * @param A Matrix A.
         * @param x0 Initial x(0) condition.
         * @param tMax Final time.
         * @param T Time step.
         * @return Matrix Solution for x(tMax).
         */
        inline static Matrix ArbitIntIter(const IntMethod& method, const Matrix& A, const Matrix& x0, double tMax, double T = DEFAULT_TIME_STEP) {
            return ArbitIntIter(method, A, Matrix(2, 2), x0, tMax, T, ZeroTimeFunction);
        }

        /**
         * @brief Make predictor-corrector integration method.
         * @param pred Predictor method.
         * @param corr Corrector method.
         * @param numCorr Number of times to apply corrector.
         * @return IntMethod Predictor-corrector integration method
         */
        static IntMethod BuildPredCorrMethod(const IntMethod& pred, const IntMethod2& corr, std::size_t numCorr = 1);

        /**
         * @brief Constructs [t t]T time matrix from given t time value.
         * @param t Time value.
         * @return Matrix Time matrix.
         */
        inline static Matrix DefaultTimeFunction(double t) {
            return Matrix(2, 1, t);
        }

        /**
         * @brief Always returns [0 0]T matrix.
         * @param t Time value (ignored).
         * @return Matrix [0 0]T matrix.
         */
        inline static Matrix ZeroTimeFunction(double t = 0.0) {
            return Matrix(2, 1);
        }

        /**
         * @brief Always returns [1 1]T matrix.
         * @param t Time value (ignored).
         * @return Matrix [1 1]T matrix.
         */
        inline static Matrix OneTimeFunction(double t = 0.0) {
            return Matrix(2, 1, 1.0);
        }

        /**
         * @brief Set destination file.
         * @param file File name.
         * @param suffix Array name suffix.
         * @param append Appends on true, else truncates.
         */
        static void SetDestFile(const std::string& file, const std::string& suffix = "", bool append = true);

        /**
         * @brief Unset destination file.
         */
        inline static void UnsetDestFile() {
            s_DestFileStream.close();
            s_Suffix = "";
        }

        /**
         * @brief Redirect file contents from source to destination file.
         * @param destfile Destination file name.
         * @param srcfile Source file name.
         * @param append Appends on true, else truncates.
         */
        static void RedirFileContent(const std::string& destfile, const  std::string& srcfile, bool append = true);

        /**
         * @brief Enable printout every N iterations.
         * @param N Number of iterations. 
         */
        inline static void EnablePrint(std::size_t N) {
            s_PrintCounter = N;
            s_EnablePrint = true;
        }

        /**
         * @brief Disable printout.
         */
        inline static void DisablePrint() {
            s_PrintCounter = 0;
            s_EnablePrint = false;
        }

    private:

        NumInt() = delete;

        static std::ofstream s_DestFileStream;
        static std::string s_Suffix;
        static std::size_t s_PrintCounter;
        static bool s_EnablePrint;
        static IntMethod s_MethodPECE2;
        static IntMethod s_MethodPECE;

    };

}

#endif
