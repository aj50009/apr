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

    std::size_t GeneticAlgorithm::BinaryPresentation::GetMinimumNumberOfBitsPerGene(std::uint8_t numberOfGenes, const std::initializer_list<double>& lowerBounds, const std::initializer_list<double>& upperBounds, std::uint8_t numberOfDecimalPlaces) {
        assert((numberOfGenes > 0) && (lowerBounds.size() == numberOfGenes) && (upperBounds.size() == numberOfGenes));
        double largestInterval = 0.0;
        const double* lowerBoundsIterator = lowerBounds.begin();
        const double* upperBoundsIterator = upperBounds.begin();
        for (std::uint8_t index = 0; index < numberOfGenes; ++index)
            largestInterval = std::max(largestInterval, *upperBoundsIterator++ - *lowerBoundsIterator++);
        return static_cast<std::size_t>(std::ceil(std::log2(std::floor(1.0 + largestInterval * std::pow(10.0, numberOfDecimalPlaces)))));
    }
    std::pair<GeneticAlgorithm::AbstractUnit::Ptr, GeneticAlgorithm::AbstractUnit::Ptr> GeneticAlgorithm::BinaryPresentation::SinglePointCrossover(const AbstractUnit::Ptr& firstParentUnit, const AbstractUnit::Ptr& secondParentUnit) {
        BinaryPresentation* binaryPresentation = dynamic_cast<BinaryPresentation*>(firstParentUnit->GetPresentation().get());
        assert(binaryPresentation);
        return SinglePointCrossover(firstParentUnit, secondParentUnit, std::rand() % (binaryPresentation->GetNumberOfBitsPerGene() + 1));
    }
    std::pair<GeneticAlgorithm::AbstractUnit::Ptr, GeneticAlgorithm::AbstractUnit::Ptr> GeneticAlgorithm::BinaryPresentation::SinglePointCrossover(const AbstractUnit::Ptr& firstParentUnit, const AbstractUnit::Ptr& secondParentUnit, std::uint8_t breakingPoint) {
        BinaryPresentation* binaryPresentation = dynamic_cast<BinaryPresentation*>(firstParentUnit->GetPresentation().get());
        BinaryUnit* firstParent = dynamic_cast<BinaryUnit*>(firstParentUnit.get());
        BinaryUnit* secondParent = dynamic_cast<BinaryUnit*>(secondParentUnit.get());
        assert((binaryPresentation) && (firstParent) && (secondParent) && (binaryPresentation == secondParent->GetPresentation().get()) && (breakingPoint <= binaryPresentation->GetNumberOfBitsPerGene()));
        uint64_t bitMask = (1 << (binaryPresentation->GetNumberOfBitsPerGene() - breakingPoint)) - 1;
        BinaryUnit* firstChild = new BinaryUnit(firstParent->GetPresentation());
        BinaryUnit* secondChild = new BinaryUnit(firstParent->GetPresentation());
        uint8_t numberOfGenes = binaryPresentation->GetNumberOfGenes();
        for (std::uint8_t index = 0; index < numberOfGenes; ++index) {
            uint64_t firstParentGene = firstParent->GetGeneBitString(index);
            uint64_t secondParentGene = firstParent->GetGeneBitString(index);
            firstChild->SetGeneBitString(index, (firstParentGene & ~bitMask) | (secondParentGene & bitMask));
            secondChild->SetGeneBitString(index, (firstParentGene & bitMask) | (secondParentGene & ~bitMask));
        }
        return std::pair<AbstractUnit::Ptr, AbstractUnit::Ptr>(AbstractUnit::Ptr(firstChild), AbstractUnit::Ptr(secondChild));
    }
    std::pair<GeneticAlgorithm::AbstractUnit::Ptr, GeneticAlgorithm::AbstractUnit::Ptr> GeneticAlgorithm::BinaryPresentation::SegmentedCrossover(const AbstractUnit::Ptr& firstParentUnit, const AbstractUnit::Ptr& secondParentUnit) {
        return SegmentedCrossover(firstParentUnit, secondParentUnit, (std::rand() % 10001) / 10000.0);
    }
    std::pair<GeneticAlgorithm::AbstractUnit::Ptr, GeneticAlgorithm::AbstractUnit::Ptr> GeneticAlgorithm::BinaryPresentation::SegmentedCrossover(const AbstractUnit::Ptr& firstParentUnit, const AbstractUnit::Ptr& secondParentUnit, double switchChance) {
        BinaryPresentation* binaryPresentation = dynamic_cast<BinaryPresentation*>(firstParentUnit->GetPresentation().get());
        BinaryUnit* firstParent = dynamic_cast<BinaryUnit*>(firstParentUnit.get());
        BinaryUnit* secondParent = dynamic_cast<BinaryUnit*>(secondParentUnit.get());
        assert((binaryPresentation) && (firstParent) && (secondParent) && (binaryPresentation == secondParent->GetPresentation().get()) && (switchChance >= 0.0) && (switchChance <= 1.0));
        uint64_t bitMask = 0;
        bool bit = true;
        uint8_t numberOfBitsPerGene = binaryPresentation->GetNumberOfBitsPerGene();
        for (std::uint8_t index = 0; index < numberOfBitsPerGene; ++index) {
            if (((std::rand() % 10001) / 10000.0) <= switchChance)
                bit = !bit;
            bitMask = (bitMask << 1) | (bit != false);
        }
        BinaryUnit* firstChild = new BinaryUnit(firstParent->GetPresentation());
        BinaryUnit* secondChild = new BinaryUnit(firstParent->GetPresentation());
        uint8_t numberOfGenes = binaryPresentation->GetNumberOfGenes();
        for (std::uint8_t index = 0; index < numberOfGenes; ++index) {
            uint64_t firstParentGene = firstParent->GetGeneBitString(index);
            uint64_t secondParentGene = firstParent->GetGeneBitString(index);
            firstChild->SetGeneBitString(index, (firstParentGene & bitMask) | (secondParentGene & ~bitMask));
            secondChild->SetGeneBitString(index, (firstParentGene & ~bitMask) | (secondParentGene & bitMask));
        }
        return std::pair<AbstractUnit::Ptr, AbstractUnit::Ptr>(AbstractUnit::Ptr(firstChild), AbstractUnit::Ptr(secondChild));
    }
    void GeneticAlgorithm::BinaryPresentation::SimpleMutation(const AbstractUnit::Ptr& unit) {
        BinaryPresentation* binaryPresentation = dynamic_cast<BinaryPresentation*>(unit->GetPresentation().get());
        assert(binaryPresentation);
        SimpleMutation(unit, std::rand() % binaryPresentation->GetNumberOfBitsPerGene());
    }
    void GeneticAlgorithm::BinaryPresentation::SimpleMutation(const AbstractUnit::Ptr& unit, std::uint8_t bitIndex) {
        BinaryPresentation* binaryPresentation = dynamic_cast<BinaryPresentation*>(unit->GetPresentation().get());
        BinaryUnit* binaryUnit = dynamic_cast<BinaryUnit*>(unit.get());
        assert((binaryPresentation) && (binaryUnit) && (bitIndex < binaryPresentation->GetNumberOfBitsPerGene()));
        uint64_t bitMask = 1 << bitIndex;
        uint8_t numberOfGenes = binaryPresentation->GetNumberOfGenes();
        for (std::uint8_t index = 0; index < numberOfGenes; ++index) {
            binaryUnit->SetGeneBitString(index, binaryUnit->GetGeneBitString(index) ^ bitMask);
        }
    }
    GeneticAlgorithm::BinaryPresentation::BinaryPresentation(std::uint8_t numberOfGenes, const std::initializer_list<double>& lowerBounds, const std::initializer_list<double>& upperBounds, std::uint8_t numberOfBitsPerGene) : AbstractPresentation(numberOfGenes, lowerBounds, upperBounds) {
        assert((numberOfBitsPerGene > 0) && (numberOfBitsPerGene <= 64));
        m_NumberOfBitsPerGene = numberOfBitsPerGene;
        SetSinglePointCrossoverFunction();
        SetSimpleMutationFunction();
    }
    GeneticAlgorithm::BinaryPresentation::~BinaryPresentation() { }
    std::uint8_t GeneticAlgorithm::BinaryPresentation::GetNumberOfBitsPerGene() const {
        return m_NumberOfBitsPerGene;
    }
    void GeneticAlgorithm::BinaryPresentation::SetSinglePointCrossoverFunction() {
        m_CrossoverFunction = static_cast<std::pair<AbstractUnit::Ptr, AbstractUnit::Ptr>(*)(const AbstractUnit::Ptr&, const AbstractUnit::Ptr&)>(SinglePointCrossover);
    }
    void GeneticAlgorithm::BinaryPresentation::SetSegmentedCrossoverFunction() {
        m_CrossoverFunction = static_cast<std::pair<AbstractUnit::Ptr, AbstractUnit::Ptr>(*)(const AbstractUnit::Ptr&, const AbstractUnit::Ptr&)>(SegmentedCrossover);
    }
    void GeneticAlgorithm::BinaryPresentation::SetCustomCrossoverFunction(const CrossoverFunction& crossoverFunction) {
        m_CrossoverFunction = crossoverFunction;
    }
    void GeneticAlgorithm::BinaryPresentation::SetSimpleMutationFunction() {
        m_MutationFunction = static_cast<void(*)(const AbstractUnit::Ptr&)>(SimpleMutation);
    }
    void GeneticAlgorithm::BinaryPresentation::SetCustomMutationFunction(const MutationFunction& mutationFunction) {
        m_MutationFunction = mutationFunction;
    }
    std::pair<GeneticAlgorithm::AbstractUnit::Ptr, GeneticAlgorithm::AbstractUnit::Ptr> GeneticAlgorithm::BinaryPresentation::Crossover(const AbstractUnit::Ptr& firstParentUnit, const AbstractUnit::Ptr& secondParentUnit) const {
        return m_CrossoverFunction(firstParentUnit, secondParentUnit);
    }
    void GeneticAlgorithm::BinaryPresentation::Mutate(const AbstractUnit::Ptr& unit) const {
        m_MutationFunction(unit);
    }

    GeneticAlgorithm::BinaryUnit::BinaryUnit(const AbstractPresentation::Ptr& presentation) : AbstractUnit(presentation) {
        assert(dynamic_cast<BinaryPresentation*>(presentation.get()));
    }
    GeneticAlgorithm::BinaryUnit::~BinaryUnit() { }
    std::uint64_t GeneticAlgorithm::BinaryUnit::GetGeneBitString(std::uint8_t index) const {
        assert(index < GetPresentation()->GetNumberOfGenes());
        return GetGenes()[index];
    }
    void GeneticAlgorithm::BinaryUnit::SetGeneBitString(std::uint8_t index, std::uint64_t bitString) {
        assert(index < GetPresentation()->GetNumberOfGenes());
        GetGenes()[index] = bitString;
    }
    double GeneticAlgorithm::BinaryUnit::GetGeneReal(std::uint8_t index) const {
        BinaryPresentation* presentation = dynamic_cast<BinaryPresentation*>(GetPresentation().get());
        assert((presentation) && (index < presentation->GetNumberOfGenes()));
        double lowerBound = presentation->GetLowerBound(index), upperBound = presentation->GetUpperBound(index);
        return lowerBound + GetGenes()[index] * (upperBound - lowerBound) / ((1 << presentation->GetNumberOfBitsPerGene()) - 1);
    }
    void GeneticAlgorithm::BinaryUnit::SetGeneReal(std::uint8_t index, double real) {
        BinaryPresentation* presentation = dynamic_cast<BinaryPresentation*>(GetPresentation().get());
        assert((presentation) && (index < presentation->GetNumberOfGenes()));
        double lowerBound = presentation->GetLowerBound(index), upperBound = presentation->GetUpperBound(index);
        GetGenes()[index] = ((1 << presentation->GetNumberOfBitsPerGene()) - 1) * (real - lowerBound) / (upperBound - lowerBound);
    }
    void GeneticAlgorithm::BinaryUnit::ClampGene(std::uint8_t index) {
        assert(index < GetPresentation()->GetNumberOfGenes());
    }
    void GeneticAlgorithm::BinaryUnit::Clamp() { }
    std::string GeneticAlgorithm::BinaryUnit::GetGeneString(std::uint8_t index) const {
        BinaryPresentation* presentation = dynamic_cast<BinaryPresentation*>(GetPresentation().get());
        assert((presentation) && (index < presentation->GetNumberOfGenes()));
        std::stringstream stream;
        uint64_t bitString = GetGenes()[index];
        for (uint64_t bitMask = (1 << (presentation->GetNumberOfBitsPerGene() - 1)); bitMask; bitMask >>= 1)
            stream << static_cast<int>((bitString & bitMask) != 0);
        return stream.str();
    }
    void GeneticAlgorithm::BinaryUnit::SetGeneString(std::uint8_t index, const std::string& str) {
        BinaryPresentation* presentation = dynamic_cast<BinaryPresentation*>(GetPresentation().get());
        assert((presentation) && (index < presentation->GetNumberOfGenes()) && (str.size() == presentation->GetNumberOfBitsPerGene()));
        uint64_t bitString = 0;
        for (const char* cstr = str.c_str(); cstr; ++cstr)
            bitString = (bitString << 1) | (*cstr != '0');
        GetGenes()[index] = bitString;
    }

}
