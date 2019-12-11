#include <GeneticAlgorithm.hxx>
#include <sstream>
#include <cassert>
#include <cmath>

namespace apr {

    GeneticAlgorithm::AbstractPresentation::AbstractPresentation(std::uint8_t numberOfGenes, const std::initializer_list<double>& lowerBounds, const std::initializer_list<double>& upperBounds) {
        assert((numberOfGenes > 0) && (lowerBounds.size() == numberOfGenes) && (upperBounds.size() == numberOfGenes));
        m_NumberOfGenes = numberOfGenes;
        m_LowerBounds = lowerBounds;
        const double* upperBoundsIterator = upperBounds.begin();
        for (std::uint8_t index = 0; index < m_NumberOfGenes; ++index)
            SetUpperBound(index, *upperBoundsIterator++);
    }
    GeneticAlgorithm::AbstractPresentation::~AbstractPresentation() { }
    std::uint8_t GeneticAlgorithm::AbstractPresentation::GetNumberOfGenes() const {
        return m_NumberOfGenes;
    }
    double GeneticAlgorithm::AbstractPresentation::GetLowerBound(std::uint8_t index) const {
        assert(index < m_NumberOfGenes);
        return m_LowerBounds[index];
    }
    void GeneticAlgorithm::AbstractPresentation::SetLowerBound(std::uint8_t index, double lowerBound) {
        assert((index < m_NumberOfGenes) && (lowerBound <= m_UpperBounds[index]));
        m_LowerBounds[index] = lowerBound;
    }
    double GeneticAlgorithm::AbstractPresentation::GetUpperBound(std::uint8_t index) const {
        assert(index < m_NumberOfGenes);
        return m_UpperBounds[index];
    }
    void GeneticAlgorithm::AbstractPresentation::SetUpperBound(std::uint8_t index, double upperBound) {
        assert((index < m_NumberOfGenes) && (upperBound >= m_LowerBounds[index]));
        m_UpperBounds[index] = upperBound;
    }

    GeneticAlgorithm::AbstractUnit::AbstractUnit(const AbstractPresentation::Ptr& presentation) {
        m_Presentation = presentation;
        m_Genes = new std::uint64_t[m_Presentation->GetNumberOfGenes()];
    }
    GeneticAlgorithm::AbstractUnit::~AbstractUnit() {
        delete[] m_Genes;
    }
    const GeneticAlgorithm::AbstractPresentation::Ptr& GeneticAlgorithm::AbstractUnit::GetPresentation() const {
        return m_Presentation;
    }
    void GeneticAlgorithm::AbstractUnit::ClampGene(std::uint8_t index) {
        assert(index < m_Presentation->GetNumberOfGenes());
        SetGeneReal(index, std::min(std::max(GetGeneReal(index), m_Presentation->GetLowerBound(index)), m_Presentation->GetUpperBound(index)));
    }
    void GeneticAlgorithm::AbstractUnit::Clamp() {
        std::uint8_t numberOfGenes = m_Presentation->GetNumberOfGenes();
        for (std::uint8_t index = 0; index < numberOfGenes; ++index)
            ClampGene(index);
    }
    std::string GeneticAlgorithm::AbstractUnit::ToString() const {
        std::stringstream stream;
        WriteToOutputStream(stream);
        return stream.str();
    }
    void GeneticAlgorithm::AbstractUnit::FromString(const std::string& str) {
        std::stringstream stream(str);
        ReadFromInputStream(stream);
    }
    std::string GeneticAlgorithm::AbstractUnit::GetGeneString(std::uint8_t index) const {
        assert(index < m_Presentation->GetNumberOfGenes());
        std::stringstream stream;
        stream << GetGeneReal(index);
        return stream.str();
    }
    void GeneticAlgorithm::AbstractUnit::SetGeneString(std::uint8_t index, const std::string& str) {
        assert(index < m_Presentation->GetNumberOfGenes());
        double real = 0.0;
        std::stringstream stream(str);
        stream >> real;
        SetGeneReal(index, real);
    }
    void GeneticAlgorithm::AbstractUnit::WriteToOutputStream(std::ostream& outputStream) const {
        outputStream << "[ ";
        std::uint8_t numberOfGenes = m_Presentation->GetNumberOfGenes();
        for (std::uint8_t index = 0; index < numberOfGenes; ++index)
            outputStream << GetGeneString(index) << " ";
        outputStream << "]";
    }
    void GeneticAlgorithm::AbstractUnit::ReadFromInputStream(std::istream& inputStream) {
        std::string str;
        std::getline(inputStream, str, '[');
        std::getline(inputStream, str, ']');
        std::stringstream stream(str);
        std::vector<std::string> geneStrings;
        while (!stream.eof()) {
            std::string geneString;
            stream >> geneString;
            if (!geneString.empty())
                geneStrings.push_back(geneString);
        }
        std::uint8_t numberOfGenes = m_Presentation->GetNumberOfGenes();
        assert(geneStrings.size() == numberOfGenes);
        for (std::uint8_t index = 0; index < numberOfGenes; ++index)
            SetGeneString(index, geneStrings[index]);
    }
    const std::uint64_t* GeneticAlgorithm::AbstractUnit::GetGenes() const {
        return m_Genes;
    }
    std::uint64_t* GeneticAlgorithm::AbstractUnit::GetGenes() {
        return m_Genes;
    }

}
