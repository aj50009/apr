#include <Matrix.hxx>
#include <stdexcept>
#include <fstream>
#include <cmath>
#include <cctype>
#include <sstream>

#define DOUBLE_EPSILON 1e-9

namespace apr {

    /* Exception messages. */
    static const std::string g_DimensionZeroMessage = "Dimension cannot be zero";
    static const std::string g_HeightZeroMessage = "Height cannot be zero";
    static const std::string g_WidthZeroMessage = "Width cannot be zero";
    static const std::string g_FileStreamBadMessage = "File stream could not be opened";
    static const std::string g_MatricesIncompatibleMessage = "Matrices are incompatible";

    Matrix::Matrix(std::size_t dimension) {
        if (dimension == 0)
            throw std::invalid_argument(g_DimensionZeroMessage);
        m_Width = dimension;
        m_Elements.resize(dimension * dimension, 0.0);
        for (std::size_t i = 0; i < dimension; ++i)
            m_Elements[dimension * i + i] = 1.0;
    }

    Matrix::Matrix(std::size_t height, std::size_t width, double value) {
        if (height == 0)
            throw std::invalid_argument(g_HeightZeroMessage);
        if (width == 0)
            throw std::invalid_argument(g_WidthZeroMessage);
        m_Width = width;
        m_Elements.resize(height * width, value);
    }

    Matrix::Matrix(const std::string& fileName) {
        std::ifstream inputFileStream(fileName);
        if (!inputFileStream.bad())
            throw std::runtime_error(g_FileStreamBadMessage);
        inputFileStream >> *this;
    }

    Matrix Matrix::operator-() const {
        Matrix result = *this;
        for (double& element : result.m_Elements)
            element = -element;
        return result;
    }

    Matrix Matrix::operator~() const {
        std::size_t height = m_Elements.size() / m_Width;
        Matrix result(m_Width, height);
        for (std::size_t i = 0; i < height; ++i)
            for (std::size_t j = 0; j < m_Width; ++j)
                result.m_Elements[height * j + i] = m_Elements[m_Width * i + j];
        return result;
    }

    Matrix& Matrix::operator+=(const Matrix& matrix) {
        if ((m_Width != matrix.m_Width) || (m_Elements.size() != matrix.m_Elements.size()))
            throw std::invalid_argument(g_MatricesIncompatibleMessage);
        for (std::size_t k = 0; k < m_Elements.size(); ++k)
            m_Elements[k] += matrix.m_Elements[k];
        return *this;
    }

    Matrix& Matrix::operator-=(const Matrix& matrix) {
        if ((m_Width != matrix.m_Width) || (m_Elements.size() != matrix.m_Elements.size()))
            throw std::invalid_argument(g_MatricesIncompatibleMessage);
        for (std::size_t k = 0; k < m_Elements.size(); ++k)
            m_Elements[k] -= matrix.m_Elements[k];
        return *this;
    }

    Matrix Matrix::operator*(const Matrix& matrix) const {
        std::size_t midDimension = matrix.m_Elements.size() / matrix.m_Width;
        if (m_Width != midDimension)
            throw std::invalid_argument(g_MatricesIncompatibleMessage);
        std::size_t height = m_Elements.size() / m_Width;
        Matrix result(height, matrix.m_Width);
        for (std::size_t i = 0; i < height; ++i)
            for (std::size_t j = 0; j < result.m_Width; ++j)
                for (std::size_t k = 0; k < midDimension; ++k)
                    result.m_Elements[result.m_Width * i + j] += m_Elements[m_Width * i + k] * matrix.m_Elements[matrix.m_Width * k + j];
        return result;
    }

    Matrix& Matrix::operator*=(double scalar) {
        for (double& element : m_Elements)
            element *= scalar;
        return *this;
    }

    Matrix& Matrix::operator/=(double scalar) {
        for (double& element : m_Elements)
            element /= scalar;
        return *this;
    }

    bool Matrix::operator==(const Matrix& matrix) const {
        if ((m_Width != matrix.m_Width) || (m_Elements.size() != matrix.m_Elements.size()))
            return false;
        for (std::size_t k = 0; k < m_Elements.size(); ++k)
            if (std::fabs(m_Elements[k] - matrix.m_Elements[k]) > DOUBLE_EPSILON)
                return false;
        return true;
    }

    Matrix operator/(double scalar, const Matrix& matrix) {
        Matrix result = matrix;
        for (double& element : result.m_Elements)
            element = scalar / element;
        return result;
    }

    std::istream& operator>>(std::istream& inputStream, Matrix& matrix) {
        if (std::isspace(inputStream.peek()))
            inputStream.get();
        if (inputStream.eof())
            return inputStream;
        matrix.m_Width = 0;
        matrix.m_Elements.clear();
        {
            std::string line;
            std::getline(inputStream, line);
            std::istringstream inputStringStream(line);
            while (!inputStringStream.eof()) {
                double value;
                inputStringStream >> value;
                matrix.m_Width++;
                matrix.m_Elements.push_back(value);
                if (std::isspace(inputStringStream.peek()))
                    inputStringStream.get();
            }
        }
        while (!inputStream.eof()) {
            std::string line;
            std::getline(inputStream, line);
            std::istringstream inputStringStream(line);
            if (std::isspace(inputStringStream.peek()))
                inputStringStream.get();
            if (inputStringStream.eof())
                break;
            for (std::size_t j = 0; j < matrix.m_Width; ++j) {
                double value = 0.0;
                inputStringStream >> value;
                matrix.m_Elements.push_back(value);
            }
        }
        return inputStream;
    }

    std::ostream& operator<<(std::ostream& outputStream, const Matrix& matrix) {
        outputStream << matrix.m_Elements[0];
        for (std::size_t j = 1; j < matrix.m_Width; ++j)
            outputStream << " " << matrix.m_Elements[j];
        for (std::size_t k = matrix.m_Width; k < matrix.m_Elements.size(); k += matrix.m_Width) {
            outputStream << std::endl << matrix.m_Elements[k];
            for (std::size_t j = 1; j < matrix.m_Width; ++j)
                outputStream << " " << matrix.m_Elements[k + j];
        }
        return outputStream;
    }

}
