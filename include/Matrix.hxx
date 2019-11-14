#ifndef __APR_MATRIX_HXX__
#define __APR_MATRIX_HXX__

#include <cstddef>
#include <cmath>
#include <string>
#include <vector>
#include <istream>
#include <ostream>

namespace apr {

    /**
     * @brief Represents a mathematical matrix of arbitrary height and width.
     */
    class Matrix {

    public:

        /**
         * @brief Construct an identity matrix of a given dimension.
         * @param dimension Square matrix dimension.
         */
        Matrix(std::size_t dimension = 1);

        /**
         * @brief Construct a matrix of a given height and width.
         * @param height Matrix height.
         * @param width Matrix width.
         * @param value Value of every element.
         */
        Matrix(std::size_t height, std::size_t width, double value = 0.0);

        /**
         * @brief Construct a matrix from file.
         * @param fileName File name.
         */
        Matrix(const std::string& fileName);

        /**
         * @brief Get matrix height.
         * @return std::size_t Matrix height.
         */
        inline std::size_t GetHeight() const { return m_Elements.size() / m_Width; }

        /**
         * @brief Get matrix width.
         * @return std::size_t Matrix width.
         */
        inline std::size_t GetWidth() const { return m_Width; }

        /**
         * @brief Immutable member access.
         * @param index Row index.
         * @return const double* Immutable row pointer.
         */
        inline const double* operator[](std::size_t index) const { return m_Elements.data() + m_Width * index; }

        /**
         * @brief Mutable member access.
         * @param index Row index.
         * @return double* Mutable row pointer.
         */
        inline double* operator[](std::size_t index) { return m_Elements.data() + m_Width * index; }

        /**
         * @brief Unary plus.
         * @return Matrix Unchanged matrix.
         */
        inline Matrix operator+() const { return *this; }

        /**
         * @brief Unary minus.
         * @return Matrix Matrix with flipped signes.
         */
        Matrix operator-() const;

        /**
         * @brief Transpose matrix.
         * @return Matrix Transposed matrix.
         */
        Matrix operator~() const;

        /**
         * @brief Add two matrices.
         * @param matrix Second matrix.
         * @return Matrix Resulting matrix.
         */
        inline Matrix operator+(const Matrix& matrix) const { Matrix result = *this; result += matrix; return result; }

        /**
         * @brief Add matrix.
         * @param matrix Other matrix.
         * @return Matrix& This matrix.
         */
        Matrix& operator+=(const Matrix& matrix);

        /**
         * @brief Subtract two matrices.
         * @param matrix Second matrix.
         * @return Matrix Resulting matrix.
         */
        inline Matrix operator-(const Matrix& matrix) const { Matrix result = *this; result -= matrix; return result; }

        /**
         * @brief Subtract matrix.
         * @param matrix Other matrix.
         * @return Matrix& This matrix.
         */
        Matrix& operator-=(const Matrix& matrix);

        /**
         * @brief Multiply two matrices.
         * @param matrix Second matrix.
         * @return Matrix Resulting matrix.
         */
        Matrix operator*(const Matrix& matrix) const;

        /**
         * @brief Multiply matrix.
         * @param matrix Other matrix.
         * @return Matrix& This matrix.
         */
        inline Matrix& operator*=(const Matrix& matrix) { return (*this = *this * matrix); }

        /**
         * @brief Multiply matrix with a scalar.
         * @param scalar Scalar value.
         * @return Matrix Resulting matrix.
         */
        inline Matrix operator*(double scalar) const { Matrix result = *this; result *= scalar; return result; };

        /**
         * @brief Multiply matrix with a scalar.
         * @param scalar Scalar value.
         * @return Matrix& This matrix.
         */
        Matrix& operator*=(double scalar);

        /**
         * @brief Divide matrix with a scalar.
         * @param scalar Scalar value.
         * @return Matrix Resulting matrix.
         */
        inline Matrix operator/(double scalar) const { Matrix result = *this; result /= scalar; return result; }

        /**
         * @brief Divide matrix with a scalar.
         * @param scalar Scalar value.
         * @return Matrix& This matrix.
         */
        Matrix& operator/=(double scalar);

        /**
         * @brief Compares two matrices for equality.
         * @param matrix Second matrix.
         * @return true Matrices are equal.
         * @return false Matrices are not equal.
         */
        bool operator==(const Matrix& matrix) const;

        /**
         * @brief Compares two matrics for inequality.
         * @param matrix Second matrix.
         * @return true Matrices are not equal.
         * @return false Matrices are equal.
         */
        inline bool operator!=(const Matrix& matrix) const { return !(*this == matrix); }

        /**
         * @brief Compute solution using forward supstitution.
         * @param freeVector Free vector (matrix).
         * @return Matrix Solution vector (matrix).
         */
        Matrix ForwardSupstitution(const Matrix& freeVector) const;

        /**
         * @brief Compute solution using back supstitution.
         * @param freeVector Free vector (matrix).
         * @return Matrix Solution vector (matrix).
         */
        Matrix BackSupstitution(const Matrix& freeVector) const;

        /**
         * @brief Compute a combined LU matrix using LU decomposition.
         * @return Matrix LU matrix.
         */
        Matrix DecompositionLU() const;

        /**
         * @brief Compute a combined LU matrix using LUP decomposition.
         * @return Matrix LU matrix.
         */
        Matrix DecompositionLUP() const;

        /**
         * @brief Compute a combined LU matrix using LUP decomposition.
         * @param permutationsCount Total number of permutations.
         * @return Matrix LU matrix.
         */
        Matrix DecompositionLUP(std::size_t& permutationsCount) const;

        /**
         * @brief Compute a combined LU matrix using LUP decomposition.
         * @param permutationsCount Total number of permutations.
         * @param permutationsMatrix Permutations matrix.
         * @return Matrix LU matrix (permutations matrix is not premultiplied).
         */
        Matrix DecompositionLUP(std::size_t& permutationsCount, Matrix& permutationsMatrix) const;

        /**
         * @brief Compute inverse using LUP decomposition.
         * @return Matrix Matrix inverse.
         */
        Matrix Inverse() const;

        /**
         * @brief Compute determinant using LUP decomposition.
         * @return double Matrix determinant.
         */
        double Determinant() const;

        /**
         * @brief Compute vector length (Euclidean norm).
         * @note Matrix must be a vector.
         * @return double Vector length (Euclidean norm).
         */
        inline double Length() const {
            return std::sqrt(LengthSquared());
        }

        /**
         * @brief Compute squared vector length (squared Euclidean norm). This is a bit faster than computing an actual vector length (Euclidean norm).
         * @note Matrix must be a vector.
         * @return double Squared vector length (squared Euclidean norm).
         */
        double LengthSquared() const;

        /**
         * @brief Normalize vector by its length (Euclidean norm).
         * @note Matrix must be a vector.
         */
        void Normalize();

        /**
         * @brief Get matrix row of specified index.
         * @param index Specified index.
         * @return Matrix Row vector (matrix).
         */
        Matrix GetRow(std::size_t index) const;

        /**
         * @brief Get matrix column of specified index.
         * @param index Specified index.
         * @return Matrix Column vector (matrix).
         */
        Matrix GetColumn(std::size_t index) const;

    private:

        /* This class and these functions can touch my privates ;) */
        friend class SearchAlgorithms;
        friend class GradientAlgorithms;
        friend std::istream& operator>>(std::istream& inputStream, Matrix& matrix);
        friend std::ostream& operator<<(std::ostream& outputStream, const Matrix& matrix);
        friend Matrix operator/(double scalar, const Matrix& matrix);

        /* Member variables. */
        std::size_t m_Width;
        std::vector<double> m_Elements;

    };

    /**
     * @brief Multiply scalar by a matrix.
     * @param scalar Scalar value.
     * @param matrix Matrix value.
     * @return Matrix Resulting matrix.
     */
    inline Matrix operator*(double scalar, const Matrix& matrix) { return matrix * scalar; }

    /**
     * @brief Divide scalar by a matrix.
     * @param scalar Scalar value.
     * @param matrix Matrix value.
     * @return Matrix Resulting matrix.
     */
    extern Matrix operator/(double scalar, const Matrix& matrix);

    /**
     * @brief Scan matrix elements from an input stream.
     * @param inputStream Input stream.
     * @param matrix Destination matrix.
     * @return std::istream& Input stream.
     */
    extern std::istream& operator>>(std::istream& inputStream, Matrix& matrix);

    /**
     * @brief Print matrix elements to an output stream.
     * @param outputStream Output stream.
     * @param matrix Source matrix.
     * @return std::ostream& Output stream.
     */
    extern std::ostream& operator<<(std::ostream& outputStream, const Matrix& matrix);

}

#endif
