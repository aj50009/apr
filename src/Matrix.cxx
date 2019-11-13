#include <Matrix.hxx>
#include <stdexcept>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cctype>
#include <sstream>

#define DOUBLE_EPSILON 1e-7

namespace apr {

    /* Exception messages. */
    static const std::string g_DimensionZeroMessage = "Dimension cannot be zero";
    static const std::string g_HeightZeroMessage = "Height cannot be zero";
    static const std::string g_WidthZeroMessage = "Width cannot be zero";
    static const std::string g_FileStreamBadMessage = "File stream could not be opened";
    static const std::string g_MatricesIncompatibleMessage = "Matrices are incompatible";
    static const std::string g_MatrixSquareMessage = "Matrix is not a square matrix";
    static const std::string g_FreeVectorMalformedMessage = "Free vector is malformed";
    static const std::string g_PivotZeroMessage = "Pivot element cannot be zero";
    static const std::string g_IndexRangeMessage = "Index is out of range";
    static const std::string g_MatrixNotVectorMessage = "Matrix is not a vector";
    static const std::string g_NullVectorMessage = "Vector is a null-vector";

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
        if (!inputFileStream.good())
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

    Matrix Matrix::ForwardSupstitution(const Matrix& freeVector) const {
        if (m_Width * m_Width != m_Elements.size())
            throw std::logic_error(g_MatrixSquareMessage);
        if (((freeVector.m_Width != 1) && (freeVector.m_Width != m_Width)) || (freeVector.m_Elements.size() != m_Width))
            throw std::logic_error(g_FreeVectorMalformedMessage);
        Matrix solution = freeVector;
        solution.m_Width = 1;
        for (std::size_t i = 0; i < m_Width - 1; ++i)
            for (std::size_t j = i + 1; j < m_Width; ++j)
                solution.m_Elements[j] -= m_Elements[m_Width * j + i] * solution.m_Elements[i];
        return solution;
    }

    Matrix Matrix::BackSupstitution(const Matrix& freeVector) const {
        if (m_Width * m_Width != m_Elements.size())
            throw std::logic_error(g_MatrixSquareMessage);
        if (((freeVector.m_Width != 1) && (freeVector.m_Width != m_Width)) || (freeVector.m_Elements.size() != m_Width))
            throw std::logic_error(g_FreeVectorMalformedMessage);
        Matrix solution = freeVector;
        solution.m_Width = 1;
        for (std::size_t i = m_Width - 1; i != std::numeric_limits<std::size_t>::max(); --i) {
            if (std::fabs(m_Elements[m_Width * i + i]) <= DOUBLE_EPSILON)
                throw std::runtime_error(g_PivotZeroMessage);
            solution.m_Elements[i] /= m_Elements[m_Width * i + i];
            for (std::size_t j = 0; j < i; ++j)
                solution.m_Elements[j] -= m_Elements[m_Width * j + i] * solution.m_Elements[i];
        }
        return solution;
    }

    Matrix Matrix::DecompositionLU() const {
        if (m_Width * m_Width != m_Elements.size())
            throw std::logic_error(g_MatrixSquareMessage);
        Matrix decomposition = *this;
        for (std::size_t i = 0; i < m_Width - 1; ++i)
            for (std::size_t j = i + 1; j < m_Width; ++j) {
                if (std::fabs(decomposition.m_Elements[m_Width * i + i]) <= DOUBLE_EPSILON)
                    throw std::runtime_error(g_PivotZeroMessage);
                decomposition.m_Elements[m_Width * j + i] /= decomposition.m_Elements[m_Width * i + i];
                for (std::size_t k = i + 1; k < m_Width; ++k)
                    decomposition.m_Elements[m_Width * j + k] -= decomposition.m_Elements[m_Width * j + i] * decomposition.m_Elements[m_Width * i + k];
            }
        return decomposition;
    }

    Matrix Matrix::DecompositionLUP() const {
        std::size_t permutationsCount;
        return DecompositionLUP(permutationsCount);
    }

    Matrix Matrix::DecompositionLUP(std::size_t& permutationsCount) const {
        Matrix permutationsMatrix;
        Matrix decomposition = this->DecompositionLUP(permutationsCount, permutationsMatrix);
        return permutationsMatrix * decomposition;
    }

    Matrix Matrix::DecompositionLUP(std::size_t& permutationsCount, Matrix& permutationsMatrix) const {
        if (m_Width * m_Width != m_Elements.size())
            throw std::logic_error(g_MatrixSquareMessage);
        Matrix decomposition = *this;
        std::vector<std::size_t> indices;
        indices.resize(m_Width);
        for (std::size_t i = 0; i < m_Width; ++i)
            indices[i] = i;
        std::size_t pivot, counter = 0;
        for (std::size_t i = 0; i < m_Width - 1; ++i) {
            pivot = i;
            for (std::size_t j = i + 1; j < m_Width; ++j)
                if (std::fabs(decomposition.m_Elements[m_Width * indices[j] + i]) > std::fabs(decomposition.m_Elements[m_Width * indices[pivot] + i]))
                    pivot = j;
            if (i != pivot) {
                std::swap(indices[i], indices[pivot]);
                counter++;
            }
            for (std::size_t j = i + 1; j < m_Width; ++j) {
                if (std::fabs(decomposition.m_Elements[m_Width * indices[i] + i]) <= DOUBLE_EPSILON)
                    throw std::runtime_error(g_PivotZeroMessage);
                decomposition.m_Elements[m_Width * indices[j] + i] /= decomposition.m_Elements[m_Width * indices[i] + i];
                for (std::size_t k = i + 1; k < m_Width; ++k)
                    decomposition.m_Elements[m_Width * indices[j] + k] -= decomposition.m_Elements[m_Width * indices[j] + i] * decomposition.m_Elements[m_Width * indices[i] + k];
            }
        }
        permutationsCount = counter;
        permutationsMatrix = Matrix(m_Width, m_Width);
        for (std::size_t i = 0; i < m_Width; ++i)
            permutationsMatrix.m_Elements[m_Width * i + indices[i]] = 1.0;
        return decomposition;
    }

    Matrix Matrix::Inverse() const {
        if (m_Width * m_Width != m_Elements.size())
            throw std::logic_error(g_MatrixSquareMessage);
        std::size_t permutationsCount;
        Matrix permutationsMatrix;
        Matrix decomposition = this->DecompositionLUP(permutationsCount, permutationsMatrix);
        decomposition = permutationsMatrix * decomposition;
        Matrix result(m_Width, m_Width);
        for (std::size_t i = 0; i < m_Width; ++i) {
            Matrix freeVector(m_Width, 1);
            freeVector.m_Elements[i] = 1.0;
            Matrix solutionVector = decomposition.BackSupstitution(decomposition.ForwardSupstitution(freeVector));
            std::copy_n(solutionVector.m_Elements.data(), m_Width, result.m_Elements.data() + m_Width * i);
        }
        return ~result * permutationsMatrix;
    }

    double Matrix::Determinant() const {
        if (m_Width * m_Width != m_Elements.size())
            throw std::logic_error(g_MatrixSquareMessage);
        std::size_t permutationsCount;
        Matrix decomposition = this->DecompositionLUP(permutationsCount);
        double determinant = (permutationsCount & 1) ? -1.0 : 1.0;
        for (std::size_t i = 0; i < m_Width; ++i)
            determinant *= decomposition.m_Elements[m_Width * i + i];
        return determinant;
    }

    double Matrix::LengthSquared() const {
        if (m_Width != 1 && m_Width != m_Elements.size())
            throw std::logic_error(g_MatrixNotVectorMessage);
        double norm2 = 0.0;
        for (double x : m_Elements)
            norm2 += x * x;
        return norm2;
    }

    void Matrix::Normalize() {
        double norm = Length();
        if (std::fabs(norm) < DOUBLE_EPSILON)
            throw std::logic_error(g_NullVectorMessage);
        for (double& x : m_Elements)
            x /= norm;
    }

    Matrix Matrix::GetRow(std::size_t index) const {
        if (index >= m_Elements.size() / m_Width)
            throw std::out_of_range(g_IndexRangeMessage);
        Matrix row(1, m_Width);
        std::copy_n(m_Elements.data() + m_Width * index, m_Width, row.m_Elements.data());
        return row;
    }

    Matrix Matrix::GetColumn(std::size_t index) const {
        if (index >= m_Width)
            throw std::out_of_range(g_IndexRangeMessage);
        std::size_t height = m_Elements.size() / m_Width;
        Matrix column(height, 1);
        for (std::size_t i = 0; i < height; ++i)
            column.m_Elements[i] = m_Elements[m_Width * i + index];
        return column;
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
