#include <NumInt.hxx>
#include <iostream>
#include <cmath>
#include <cassert>

namespace apr {

    Matrix NumInt::RungeKuttaMethod(const Matrix& A, const Matrix& B, const Matrix& x, double t, double T, const TimeFunction& r) {
        assert((A.GetHeight() == 2) && (A.GetWidth() == 2) && (B.GetHeight() == 2) && (B.GetWidth() == 2) && (x.GetHeight() == 2) && (x.GetWidth() == 1) && (T != 0.0));
        Matrix k1 = T * (A * x + B * r(t));
        Matrix k2 = T * (A * (x + k1 / 2.0) + B * r(t + T / 2.0));
        Matrix k3 = T * (A * (x + k2 / 2.0) + B * r(t + T / 2.0));
        Matrix k4 = T * (A * (x + k3) + B * r(t + T));
        return x + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }

    Matrix NumInt::TrapezoidalRuleMethod(const Matrix& A, const Matrix& B, const Matrix& x, double t, double T, const TimeFunction& r) {
        assert((A.GetHeight() == 2) && (A.GetWidth() == 2) && (B.GetHeight() == 2) && (B.GetWidth() == 2) && (x.GetHeight() == 2) && (x.GetWidth() == 1) && (T != 0.0));
        return (Matrix(2) - (T / 2.0) * A).Inverse() * (x + (T / 2.0) * (A * x + B * (r(t) + r(t + T))));
    }

    Matrix NumInt::TrapezoidalRuleMethod(const Matrix& A, const Matrix& B, const Matrix& x, const Matrix& xp, double t, double T, const TimeFunction& r) {
        assert((A.GetHeight() == 2) && (A.GetWidth() == 2) && (B.GetHeight() == 2) && (B.GetWidth() == 2) && (x.GetHeight() == 2) && (x.GetWidth() == 1) && (T != 0.0));
        return x + (T / 2.0) * (A * x + B * r(t) + A * xp + B * r(t + T));
    }

    Matrix NumInt::EulerMethod(const Matrix& A, const Matrix& B, const Matrix& x, double t, double T, const TimeFunction& r) {
        assert((A.GetHeight() == 2) && (A.GetWidth() == 2) && (B.GetHeight() == 2) && (B.GetWidth() == 2) && (x.GetHeight() == 2) && (x.GetWidth() == 1) && (T != 0.0));
        return x + T * (A * x + B * r(t));
    }

    Matrix NumInt::BackwardEulerMethod(const Matrix& A, const Matrix& B, const Matrix& x, double t, double T, const TimeFunction& r) {
        assert((A.GetHeight() == 2) && (A.GetWidth() == 2) && (B.GetHeight() == 2) && (B.GetWidth() == 2) && (x.GetHeight() == 2) && (x.GetWidth() == 1) && (T != 0.0));
        return (Matrix(2) - T * A).Inverse() * (x + T * B * r(t + T));
    }

    Matrix NumInt::BackwardEulerMethod(const Matrix& A, const Matrix& B, const Matrix& x, const Matrix& xp, double t, double T, const TimeFunction& r) {
        assert((A.GetHeight() == 2) && (A.GetWidth() == 2) && (B.GetHeight() == 2) && (B.GetWidth() == 2) && (x.GetHeight() == 2) && (x.GetWidth() == 1) && (T != 0.0));
        return x + T * (A * xp + B * r(t + T));
    }

    Matrix NumInt::ArbitIntIter(const IntMethod& method, const Matrix& A, const Matrix& B, const Matrix& x0, double tMax, double T, const TimeFunction& r) {
        assert((A.GetHeight() == 2) && (A.GetWidth() == 2) && (B.GetHeight() == 2) && (B.GetWidth() == 2) && (x0.GetHeight() == 2) && (x0.GetWidth() == 1) && (T != 0.0) && (std::signbit(tMax) == std::signbit(T)));
        Matrix x = x0;
        double t = 0.0;
        std::size_t n = static_cast<std::size_t>(std::ceil(tMax / T)) - 1;
        std::size_t counter = 0;
        if (s_DestFileStream.is_open())
            s_DestFileStream << "xt" << s_Suffix << " = [" << std::endl << "\t" << x[0][0] << ", " << x[1][0] << ", " << t << "; ";
        for (std::size_t i = 0; i < n; ++i) {
            x = method(A, B, x, t, T, r);
            t += T;
            if ((s_EnablePrint) && (++counter >= s_PrintCounter)) {
                counter = 0;
                std::cout << "   x(" << t << ") = [ " << ~x << " ]" << std::endl;
            }
            if (s_DestFileStream.is_open())
                s_DestFileStream << x[0][0] << ", " << x[1][0] << ", " << t << "; "; 
        }
        double Tmod = std::fmod(tMax, T);
        if (Tmod != 0.0)
            T = Tmod;
        x = method(A, B, x, t, T, r);
        t += T;
        if ((s_EnablePrint) && (++counter >= s_PrintCounter)) {
            counter = 0;
            std::cout << "   x(" << t << ") = [ " << ~x << " ]" << std::endl;
        }
        DisablePrint();
        if (s_DestFileStream.is_open()) {
            s_DestFileStream << x[0][0] << ", " << x[1][0] << ", " << t << std::endl << "];" << std::endl;
            UnsetDestFile();
        }
        return x;
    }

    IntMethod NumInt::BuildPredCorrMethod(const IntMethod& pred, const IntMethod2& corr, std::size_t numCorr) {
        assert(numCorr > 0);
        return [pred, corr, numCorr](const Matrix& A, const Matrix& B, const Matrix& x, double t, double T, const TimeFunction& r) {
            Matrix xp = pred(A, B, x, t, T, r);
            for (std::size_t i = 0; i < numCorr; ++i)
                xp = corr(A, B, x, xp, t, T, r);
            return xp;
        };
    }

    void NumInt::SetDestFile(const std::string& file, const std::string& suffix, bool append) {
        UnsetDestFile();
        s_DestFileStream = std::ofstream(file, std::ios::out | (append ? std::ios::app : std::ios::trunc));
        s_Suffix = suffix;
    }

    void NumInt::RedirFileContent(const std::string& destfile, const  std::string& srcfile, bool append) {
        UnsetDestFile();
        std::ofstream ofs(destfile, std::ios::out | (append ? std::ios::app : std::ios::trunc));
        std::ifstream ifs(srcfile, std::ios::in);
        std::string ln;
        while (!ifs.eof()) {
            std::getline(ifs, ln);
            ofs << ln << std::endl;
        }
    }

    std::ofstream NumInt::s_DestFileStream;
    std::string NumInt::s_Suffix = "";
    std::size_t NumInt::s_PrintCounter = 0;
    bool NumInt::s_EnablePrint = false;
    IntMethod NumInt::s_MethodPECE2 = NumInt::BuildPredCorrMethod(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(NumInt::EulerMethod), static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(NumInt::BackwardEulerMethod), 2);
    IntMethod NumInt::s_MethodPECE = NumInt::BuildPredCorrMethod(static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(NumInt::EulerMethod), static_cast<Matrix(*)(const Matrix&, const Matrix&, const Matrix&, const Matrix&, double, double, const TimeFunction&)>(NumInt::TrapezoidalRuleMethod));

}
