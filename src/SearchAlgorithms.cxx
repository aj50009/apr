#include <SearchAlgorithms.hxx>
#include <Matrix.hxx>
#include <cstddef>
#include <cmath>
#include <stdexcept>

#ifdef PRINT_STEPS
#include <iostream>
#define PRINT(s) { std::cout << s << std::endl; }
#else
#define PRINT(s) { }
#endif

#define GOLDEN_RATIO 1.61803398875
#define MAX_ITER_MSG "Maximum number of iterations exceeded."

namespace apr {

#define PRINT_GSS_INTRO() { PRINT("~~~~~ Golden Section Search for [" << start << ", " << end << "] ~~~~~") }
#define PRINT_GSS_STEP() { PRINT("f(" << start << ") = " << fn(start) << std::endl << "f(" << c << ") = " << fn(c) << std::endl << "f(" << d << ") = " << fn(d) << std::endl << "f(" << end << ") = " << fn(end) << std::endl << "------------------------------------------------------------") }
#define PRINT_GSS_RESULT() { PRINT("Found minimum f(" << (start + end) / 2.0 << ") = " << fn((start + end) / 2.0) << std::endl << "------------------------------------------------------------") }

    double SearchAlgorithms::GoldenSectionSearch(const std::function<double(double)>& fn, double start, double end, double eps) {
        PRINT_GSS_INTRO()
        double c = end - (end - start) / GOLDEN_RATIO;
        double d = start + (end - start) / GOLDEN_RATIO;
        double fc = fn(c);
        double fd = fn(d);
        PRINT_GSS_STEP()
        while ((end - start) > eps) {
            if (fc < fd) {
                end = d;
                d = c;
                c = end - (end - start) / GOLDEN_RATIO;
                fd = fc;
                fc = fn(c);
            } else {
                start = c;
                c = d;
                d = start + (end - start) / GOLDEN_RATIO;
                fc = fd;
                fd = fn(d);
            }
            PRINT_GSS_STEP()
        }
        PRINT_GSS_RESULT()
        return (start + end) / 2.0;
    }

    double SearchAlgorithms::GoldenSectionSearchPoint(const std::function<double(double)>& fn, double point, double offs, double eps) {
        double start, end;
        UnimodalIntervalSearch(fn, point, offs, start, end);
        return GoldenSectionSearch(fn, start, end, eps);
    }

#define PRINT_UIS_INTRO() { PRINT("~~~~~ Unimodal Interval Search around point " << point << " ~~~~~") }
#define PRINT_UIS_STEP() { PRINT("f(" << start << ") = " << fn(start) << std::endl << "f(" << m << ") = " << fn(m) << std::endl << "f(" << end << ") = " << fn(end) << std::endl << "------------------------------------------------------------") }
#define PRINT_UIS_RESULT() { PRINT("Found unimodal interval [" << start << ", " << end << "]" << std::endl << "------------------------------------------------------------") }

    void SearchAlgorithms::UnimodalIntervalSearch(const std::function<double(double)>& fn, double point, double offs, double& start, double& end) {
        PRINT_UIS_INTRO()
        start = point - offs;
        end = point + offs;
        double m = point, fl = fn(start), fr = fn(end), fm = fn(point);
        unsigned int step = 1;
        PRINT_UIS_STEP()
        if ((fm < fr) && (fm < fl)) {
            PRINT_UIS_RESULT()
            return;
        } else if (fm > fr)
            do {
                start = m;
                m = end;
                fm = fr;
                end = point + offs * (step *= 2);
                fr = fn(end);
                PRINT_UIS_STEP()
            } while(fm > fr);
        else 
            do {
                end = m;
                m = start;
                fm = fl;
                start = point - offs * (step *= 2);
                fl = fn(start);
                PRINT_UIS_STEP()
            } while(fm > fl);
        PRINT_UIS_RESULT()
    }

    std::vector<double> SearchAlgorithms::CoordinateAxesSearch(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& point, const std::vector<double>& offs, double eps) {
        std::vector<double> x = point, xs;
        double dist;
        eps *= eps;
        do {
            xs = x;
            dist = 0.0;
            for (std::size_t i = 0; i < point.size(); ++i) {
                x[i] = GoldenSectionSearchPoint([&fn, &x, i](double xi) { x[i] = xi; return fn(x); }, point[i], offs[i], eps);
                dist += (x[i] - xs[i]) * (x[i] - xs[i]);
            }
        } while (dist > eps);
        return x;
    }

#define PRINT_DSM_INTRO() { PRINT("~~~~~ Downhill Simplex Method around point (" << ~(x[0]) << ") ~~~~~") }
#define PRINT_DSM_STEP() { PRINT("f(" << ~xc << ") = " << fn(xc.m_Elements) << std::endl << "------------------------------------------------------------") }
#define PRINT_DSM_RESULT() { PRINT("Found minimum f(" << ~avg << ") = " << fn(avg.m_Elements) << std::endl << "------------------------------------------------------------") }

    std::vector<double> SearchAlgorithms::DownhillSimplexMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& point, double reflection, double contraction, double expansion, double shrink, double initscale, double eps, std::size_t maxiter) {
        Matrix x[point.size() + 1];
        x[0].m_Elements = point;
        for (std::size_t i = 1; i <= point.size(); ++i) {
            x[i].m_Elements = point;
            x[i].m_Elements[i - 1] += initscale;
        }
        PRINT_DSM_INTRO()
        double stdev;
        do {
            if (maxiter-- == 0)
                throw std::runtime_error(MAX_ITER_MSG);
            /* Find min and max value vertices */
            std::size_t h = 0, l = 0;
            double fh, fl;
            fh = fl = fn(x[0].m_Elements);
            for (std::size_t i = 1; i <= point.size(); ++i) {
                double fi = fn(x[i].m_Elements);
                if (fi > fh) {
                    h = i;
                    fh = fi;
                }
                if (fi < fl) {
                    l = i;
                    fl = fi;
                }
            }
            /* Compute centroid */
            Matrix xc(point.size(), 1);
            for (std::size_t i = 0; i <= point.size(); ++i) {
                if (i == h) continue;
                xc += x[i];
            }
            xc /= point.size();
            PRINT_DSM_STEP()
            /* Compute reflection */
            Matrix xr = xc + reflection * (xc - x[h]);
            double fr = fn(xr.m_Elements);
            if (fr < fl) {
                /* Compute expansion */
                Matrix xe = xc + expansion * (xr - xc);
                if (fn(xe.m_Elements) < fl)
                    x[h] = xe;
                else
                    x[h] = xr;
            } else {
                bool is_second_worst = true;
                for (std::size_t i = 0; i <= point.size(); ++i) {
                    if (i == h) continue;
                    if (fr <= fn(x[i].m_Elements)) {
                        is_second_worst = false;
                        break;
                    }
                }
                if (is_second_worst) {
                    if (fr < fh) {
                        x[h] = xr;
                        fh = fr;
                    }
                    /* Compute contraction */
                    Matrix xk = xc + contraction * (x[h] - xc);
                    if (fn(xk.m_Elements) < fh)
                        x[h] = xk;
                    else {
                        /* Perform shrinking */
                        for (std::size_t i = 0; i <= point.size(); ++i) {
                            if (i == l) continue;
                            x[i] = x[l] + shrink * (x[i] - x[l]);
                        }
                    }
                } else
                    x[h] = xr;
            }
            /* Calculate 'standard deviation' */
            double mean = fn(xc.m_Elements);
            stdev = 0.0;
            for (const Matrix& xi : x) {
                double value = fn(xi.m_Elements);
                stdev += (value - mean) * (value - mean);
            }
            stdev /= point.size();
            stdev = std::sqrt(stdev);
        } while (stdev > eps);
        Matrix avg(point.size(), 1);
        for (const Matrix& xi : x)
            avg += xi;
        avg /= point.size() + 1;
        PRINT_DSM_RESULT()
        return avg.m_Elements;
    }

#define PRINT_PS_INTRO() { PRINT("~~~~~ Pattern Search around point (" << ~xb << ") ~~~~~") }
#define PRINT_PS_STEP() { PRINT("xb = (" << ~xb << "), f(xb) = " << fn(xb.m_Elements) << std::endl << "xp = (" << ~xp << "), f(xp) = " << fn(xp.m_Elements) << std::endl << "xn = (" << ~xn << "), f(xn) = " << fn(xn.m_Elements) << std::endl << "------------------------------------------------------------") }
#define PRINT_PS_RESULT() { PRINT("Found minimum f(" << ~xb << ") = " << fn(xb.m_Elements) << std::endl << "------------------------------------------------------------") }
    std::vector<double> SearchAlgorithms::PatternSearch(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& point, const std::vector<double>& dx, const std::vector<double>& decay, const std::vector<double>& epsvec, std::size_t maxiter) {
        Matrix xb, xp, dxm;
        xb.m_Elements = xp.m_Elements = point;
        dxm.m_Elements = dx;
        PRINT_PS_INTRO()
        bool repeat = true;
        do {
            if (maxiter-- == 0)
                throw std::runtime_error(MAX_ITER_MSG);
            Matrix xn = xp;
            for (std::size_t i = 0; i < xn.m_Elements.size(); ++i) {
                double p = fn(xn.m_Elements);
                xn.m_Elements[i] += dxm.m_Elements[i];
                double n = fn(xn.m_Elements);
                if (n > p) {
                    xn.m_Elements[i] -= 2 * dxm.m_Elements[i];
                    n = fn(xn.m_Elements);
                    if (n > p)
                        xn.m_Elements[i] += dxm.m_Elements[i];
                }
            }
            PRINT_PS_STEP()
            if (fn(xn.m_Elements) < fn(xb.m_Elements)) {
                xp = 2 * xn - xb;
                xb = xn;
            } else {
                Matrix ddxm = dxm;
                for (std::size_t i = 0; i < dxm.m_Elements.size(); ++i)
                    dxm.m_Elements[i] *= decay[i];
                ddxm -= dxm;
                /* dxm contains the current dx matrix, while ddxm contains the difference between the previous and the current dx matrix */
                repeat = false;
                for (std::size_t i = 0; i < ddxm.m_Elements.size(); ++i)
                    if (ddxm.m_Elements[i] > epsvec[i]) { /* switch between dxm or ddxm matrices for comparison */
                        repeat = true;
                        break;
                    }
                xp = xb;
            }
        } while (repeat);
        PRINT_PS_RESULT()
        return xb.m_Elements;
    }

}
