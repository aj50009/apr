#include <GradientAlgorithms.hxx>
#include <SearchAlgorithms.hxx>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <iostream>

namespace apr {

    std::vector<double> GradientAlgorithms::ComputeGradient(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double delta) {
        s_NumGradientComputations++;
        double fn0 = fn(x0);
        std::vector<double> xi = x0, grad;
        grad.reserve(xi.size());
        for (std::size_t i = 0; i < xi.size(); ++i) {
            xi[i] += delta;
            grad.push_back((fn(xi) - fn0) / delta);
            xi[i] = x0[i];
        }
        return grad;
    }

    Matrix GradientAlgorithms::ComputeHessian(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double delta) {
        s_NumHessianComputations++;
        double fn0 = fn(x0);
        std::vector<double> fncch;
        fncch.reserve(x0.size());
        Matrix fncch2(x0.size(), x0.size());
        std::vector<double> xij = x0;
        for (std::size_t i = 0; i < x0.size(); ++i) {
            xij[i] += delta;
            fncch.push_back(fn(xij));
            for (std::size_t j = 0; j <= i; ++j) {
                double bkp = xij[j];
                xij[j] += delta;
                fncch2[i][j] = fncch2[j][i] = fn(xij);
                xij[j] = bkp;
            }
            xij[i] = x0[i];
        }
        Matrix hessian(x0.size(), x0.size());
        for (std::size_t i = 0; i < x0.size(); ++i)
            for (std::size_t j = 0; j <= i; ++j)
                hessian[i][j] = hessian[j][i] = (fncch2[i][j] - fncch[i] - fncch[j] + fn0) / (delta * delta);
        return hessian;
    }

    std::vector<double> GradientAlgorithms::GradientDescent(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double moveratio, double movedecay, double delta, double eps, std::size_t maxiter) {
        double eps2 = eps * eps;
        Matrix x, grad;
        x.m_Elements = x0;
        grad.m_Elements = ComputeGradient(fn, x.m_Elements, delta);
        while (grad.LengthSquared() > eps2) {
            if (--maxiter == 0)
                throw std::runtime_error("Max number of iterations exceeded");
            x -= moveratio * grad;
            if ((moveratio *= movedecay) < eps)
                throw std::runtime_error("Move ratio reached zero");
            grad.m_Elements = ComputeGradient(fn, x.m_Elements, delta);
        }
        return x.m_Elements;
    }

    std::vector<double> GradientAlgorithms::GradientDescentGSS(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double offset, double delta, double eps, std::size_t maxiter) {
        double eps2 = eps * eps;
        Matrix x, grad;
        x.m_Elements = x0;
        grad.m_Elements = ComputeGradient(fn, x.m_Elements, delta);
        while (grad.LengthSquared() > eps2) {
            if (--maxiter == 0)
                throw std::runtime_error("Max number of iterations exceeded");
            grad.Normalize();
            double moveratio = SearchAlgorithms::GoldenSectionSearchPoint([&x, &grad, &fn](double k) { Matrix p = x + k * grad; return fn(p.m_Elements); }, 0.0, offset/*, eps*/);
            x += moveratio * grad;
            grad.m_Elements = ComputeGradient(fn, x.m_Elements, delta);
        }
        return x.m_Elements;
    }

    std::vector<double> GradientAlgorithms::NewtonRaphsonMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double moveratio, double movedecay, double delta, double eps, std::size_t maxiter) {
        double eps2 = eps * eps;
        Matrix x, grad;
        x.m_Elements = x0;
        grad.m_Elements = ComputeGradient(fn, x.m_Elements, delta);
        while (grad.LengthSquared() > eps2) {
            if (--maxiter == 0)
                throw std::runtime_error("Max number of iterations exceeded");
            x -= moveratio * ComputeHessian(fn, x.m_Elements, delta).Inverse() * grad;
            if ((moveratio *= movedecay) < eps)
                throw std::runtime_error("Move ratio reached zero");
            grad.m_Elements = ComputeGradient(fn, x.m_Elements, delta);
        }
        return x.m_Elements;
    }

    std::vector<double> GradientAlgorithms::NewtonRaphsonMethodGSS(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, double offset, double delta, double eps, std::size_t maxiter) {
        double eps2 = eps * eps;
        Matrix x, grad;
        x.m_Elements = x0;
        grad.m_Elements = ComputeGradient(fn, x.m_Elements, delta);
        while (grad.LengthSquared() > eps2) {
            if (--maxiter == 0)
                throw std::runtime_error("Max number of iterations exceeded");
            grad = ComputeHessian(fn, x.m_Elements, delta).Inverse() * grad;
            grad.Normalize();
            double moveratio = SearchAlgorithms::GoldenSectionSearchPoint([&x, &grad, &fn](double k) { Matrix p = x + k * grad; return fn(p.m_Elements); }, 0.0, offset/*, eps*/);
            x += moveratio * grad;
            grad.m_Elements = ComputeGradient(fn, x.m_Elements, delta);
        }
        return x.m_Elements;
    }

    std::vector<double> GradientAlgorithms::ConstrainedBoxMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, std::size_t nconstraints, const std::function<bool(const std::vector<double>&)>* constraints, const std::vector<double>& min, const std::vector<double>& max, double alpha, double eps, std::size_t maxiter) {
        if ((x0.size() != min.size()) || (x0.size() != max.size()))
            throw std::logic_error("Mismatched dimensions of starting point and explicit constraints");
        for (std::size_t i = 0; i < x0.size(); ++i)
            if ((x0[i] < min[i]) || (x0[i] > max[i]))
                throw std::logic_error("Starting point does not satisfy explicit constraints");
        for (std::size_t i = 0; i < nconstraints; ++i)
            if (!constraints[i](x0))
                throw std::logic_error("Starting point does not satisfy implicit constraints");
        double eps2 = eps * eps;
        Matrix xc;
        xc.m_Elements = x0;
        std::vector<Matrix> x;
        x.reserve(2 * x0.size());
        for (std::size_t i = 0; i < 2 * x0.size(); ++i) {
            Matrix xi;
            xi.m_Elements.clear();
            xi.m_Elements.reserve(x0.size());
            for (std::size_t j = 0; j < x0.size(); ++j)
                xi.m_Elements.push_back(min[j] + ((std::rand() % 10000) / 10000.0) * (max[j] - min[j]));
            while (true) {
                bool brklp = true;
                for (std::size_t j = 0; j < nconstraints; ++j)
                    if (!constraints[j](xi.m_Elements)) {
                        brklp = false;
                        break;
                    }
                if (brklp)
                    break;
                xi = (xi + xc) / 2.0;
            }
            x.push_back(xi);
            if (i != 2 * x0.size() - 1) {
                xc *= x.size();
                xc += xi;
                xc /= (x.size() + 1);
            }
        }
        double stdev2;
        do {
            if (--maxiter == 0)
                throw std::runtime_error("Max number of iterations exceeded");
            std::size_t h2 = 1, h = 0;
            double fnh2 = fn(x[1].m_Elements), fnh = fn(x[0].m_Elements);
            if (fnh2 > fnh) {
                double bkp = fnh2;
                fnh2 = fnh;
                fnh = bkp;
                h2 = 0;
                h = 1;
            }
            for (std::size_t i = 2; i < x.size(); ++i) {
                double fni = fn(x[i].m_Elements);
                if (fni > fnh) {
                    fnh2 = fnh;
                    fnh = fni;
                    h2 = h;
                    h = i;
                } else if (fni > fnh2) {
                    fnh2 = fni;
                    h2 = i;
                }
            }
            xc = Matrix(x0.size(), 1);
            for (std::size_t i = 0; i < x.size(); ++i)
                if (i != h)
                    xc += x[i];
            xc /= (x.size() - 1);
            Matrix xr = (1.0 + alpha) * xc - alpha * x[h];
            for (std::size_t i = 0; i < x0.size(); ++i)
                if (xr.m_Elements[i] < min[i])
                    xr.m_Elements[i] = min[i];
                else if (xr.m_Elements[i] > max[i])
                    xr.m_Elements[i] = max[i];
            while (true) {
                bool brklp = true;
                for (std::size_t j = 0; j < nconstraints; ++j)
                    if (!constraints[j](xr.m_Elements)) {
                        brklp = false;
                        break;
                    }
                if (brklp)
                    break;
                xr = (xr + xc) / 2.0;
            }
            if (fn(xr.m_Elements) > fnh2)
                xr = (xr + xc) / 2.0;
            x[h] = xr;
            xc *= (x.size() - 1);
            xc += xr;
            xc /= x.size();
            double fnc = fn(xc.m_Elements);
            stdev2 = 0.0;
            for (const Matrix& xi : x) {
                double fni = fn(xi.m_Elements);
                stdev2 += (fni - fnc) * (fni - fnc);
            }
            stdev2 /= x.size();
        } while (stdev2 > eps2);
        return xc.m_Elements;
    }

    std::vector<double> GradientAlgorithms::ConstrainedTransformationMethod(const std::function<double(const std::vector<double>&)>& fn, const std::vector<double>& x0, std::size_t ncequalities, const std::function<double(const std::vector<double>&)>* cequalities, std::size_t ncinequalities, const std::function<double(const std::vector<double>&)>* cinequalities, double t, double ts, double dxscale, double decayfactor, double eps, std::size_t maxiter) {
        double eps2 = eps * eps;
        std::size_t maxiterps = maxiter;
        std::function<double(const std::vector<double>&)> ft = [&fn, ncequalities, cequalities, ncinequalities, cinequalities, &t, eps](const std::vector<double>& x) {
            double gsum = 0.0, hsum = 0.0;
            for (std::size_t i = 0; i < ncinequalities; ++i) {
                double g = cinequalities[i](x);
                if (g < eps)
                    return std::numeric_limits<double>::infinity();
                gsum += std::log(g);
            }
            for (std::size_t i = 0; i < ncequalities; ++i) {
                double h = cequalities[i](x);
                hsum += h * h;
            }
            return fn(x) - (1.0 / t) * gsum + t * hsum;
        };
        std::vector<double> x = x0;
        for (std::size_t i = 0; i < ncinequalities; ++i)
            if (cinequalities[i](x) < -eps) {
                x = SearchAlgorithms::PatternSearch([ncinequalities, cinequalities](const std::vector<double>& x) {
                    double sum = 0.0;
                    for (std::size_t i = 0; i < ncinequalities; ++i) {
                        double value = cinequalities[i](x);
                        if (value < 0.0)
                            sum += value;
                    }
                    return -sum;
                }, x, dxscale, decayfactor, eps, maxiterps);
                break;
            }
        double dist2;
        do {
            if (--maxiter == 0)
                throw std::runtime_error("Max number of iterations exceeded");
            std::vector<double> xn = SearchAlgorithms::PatternSearch(ft, x, dxscale, decayfactor, eps, maxiterps);
            dist2 = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i)
                dist2 += (xn[i] - x[i]) * (xn[i] - x[i]);
            x = xn;
            t *= ts;
        } while (dist2 > eps2);
        return x;
    }

    std::size_t GradientAlgorithms::s_NumGradientComputations, GradientAlgorithms::s_NumHessianComputations = 0;

}
