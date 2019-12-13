#include <GeneticAlgorithm.hxx>
#include <sstream>
#include <limits>
#include <cassert>
#include <cmath>

namespace apr {

    GeneticAlgorithm::AbstractPresentation::AbstractPresentation(std::uint8_t numberOfGenes, const std::initializer_list<double>& lowerBounds, const std::initializer_list<double>& upperBounds) {
        assert((numberOfGenes > 0) && (lowerBounds.size() == numberOfGenes) && (upperBounds.size() == numberOfGenes));
        m_NumberOfGenes = numberOfGenes;
        m_LowerBounds = lowerBounds;
        m_UpperBounds.resize(m_NumberOfGenes);
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
            uint64_t secondParentGene = secondParent->GetGeneBitString(index);
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
            uint64_t secondParentGene = secondParent->GetGeneBitString(index);
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
    GeneticAlgorithm::AbstractUnit::Ptr GeneticAlgorithm::BinaryPresentation::NewUnit(const AbstractPresentation::Ptr& presentation) const {
        return AbstractUnit::Ptr(new BinaryUnit(presentation));
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
        for (const char* cstr = str.c_str(); *cstr; ++cstr)
            bitString = (bitString << 1) | (*cstr != '0');
        GetGenes()[index] = bitString;
    }

    std::pair<GeneticAlgorithm::AbstractUnit::Ptr, GeneticAlgorithm::AbstractUnit::Ptr> GeneticAlgorithm::FloatingPointPresentation::AverageCrossover(const AbstractUnit::Ptr& firstParentUnit, const AbstractUnit::Ptr& secondParentUnit) {
        FloatingPointPresentation* presentation = reinterpret_cast<FloatingPointPresentation*>(firstParentUnit->GetPresentation().get());
        assert((presentation) && (presentation == secondParentUnit->GetPresentation().get()));
        FloatingPointUnit* firstChild = new FloatingPointUnit(firstParentUnit->GetPresentation());
        FloatingPointUnit* secondChild = new FloatingPointUnit(firstParentUnit->GetPresentation());
        uint8_t numberOfGenes = presentation->GetNumberOfGenes();
        for (uint8_t index = 0; index < numberOfGenes; ++index) {
            double averageGeneReal = (firstParentUnit->GetGeneReal(index) + secondParentUnit->GetGeneReal(index)) / 2.0;
            firstChild->SetGeneReal(index, averageGeneReal);
            secondChild->SetGeneReal(index, averageGeneReal);
        }
        return std::pair<AbstractUnit::Ptr, AbstractUnit::Ptr>(AbstractUnit::Ptr(firstChild), AbstractUnit::Ptr(secondChild));
    }
    std::pair<GeneticAlgorithm::AbstractUnit::Ptr, GeneticAlgorithm::AbstractUnit::Ptr> GeneticAlgorithm::FloatingPointPresentation::ArithmeticCrossover(const AbstractUnit::Ptr& firstParentUnit, const AbstractUnit::Ptr& secondParentUnit) {
        return ArithmeticCrossover(firstParentUnit, secondParentUnit, (std::rand() % 10001) / 10000.0);
    }
    std::pair<GeneticAlgorithm::AbstractUnit::Ptr, GeneticAlgorithm::AbstractUnit::Ptr> GeneticAlgorithm::FloatingPointPresentation::ArithmeticCrossover(const AbstractUnit::Ptr& firstParentUnit, const AbstractUnit::Ptr& secondParentUnit, double interpolationRatio) {
        FloatingPointPresentation* presentation = reinterpret_cast<FloatingPointPresentation*>(firstParentUnit->GetPresentation().get());
        assert((presentation) && (presentation == secondParentUnit->GetPresentation().get()) && (interpolationRatio >= 0.0) && (interpolationRatio <= 1.0));
        FloatingPointUnit* firstChild = new FloatingPointUnit(firstParentUnit->GetPresentation());
        FloatingPointUnit* secondChild = new FloatingPointUnit(firstParentUnit->GetPresentation());
        uint8_t numberOfGenes = presentation->GetNumberOfGenes();
        for (uint8_t index = 0; index < numberOfGenes; ++index) {
            double firstParentGeneReal = firstParentUnit->GetGeneReal(index);
            double secondParentGeneReal = secondParentUnit->GetGeneReal(index);
            firstChild->SetGeneReal(index, firstParentGeneReal + interpolationRatio * (secondParentGeneReal - firstParentGeneReal));
            secondChild->SetGeneReal(index, secondParentGeneReal + interpolationRatio * (firstParentGeneReal - secondParentGeneReal));
        }
        return std::pair<AbstractUnit::Ptr, AbstractUnit::Ptr>(AbstractUnit::Ptr(firstChild), AbstractUnit::Ptr(secondChild));
    }
    void GeneticAlgorithm::FloatingPointPresentation::UniformMutation(const AbstractUnit::Ptr& unit) {
        FloatingPointPresentation* presentation = reinterpret_cast<FloatingPointPresentation*>(unit->GetPresentation().get());
        assert(presentation);
        uint8_t index = std::rand() % presentation->GetNumberOfGenes();
        unit->SetGeneReal(index, unit->GetGeneReal(index) + (std::rand() % 20001) / 10000.0 - 1.0);
    }
    GeneticAlgorithm::FloatingPointPresentation::FloatingPointPresentation(std::uint8_t numberOfGenes, const std::initializer_list<double>& lowerBounds, const std::initializer_list<double>& upperBounds) : AbstractPresentation(numberOfGenes, lowerBounds, upperBounds) {
        SetAverageCrossoverFunction();
        SetUniformMutationFunction();
    }
    GeneticAlgorithm::FloatingPointPresentation::~FloatingPointPresentation() { }
    void GeneticAlgorithm::FloatingPointPresentation::SetAverageCrossoverFunction() {
        m_CrossoverFunction = AverageCrossover;
    }
    void GeneticAlgorithm::FloatingPointPresentation::SetArithmeticCrossoverFunction() {
        m_CrossoverFunction = static_cast<std::pair<AbstractUnit::Ptr, AbstractUnit::Ptr>(*)(const AbstractUnit::Ptr&, const AbstractUnit::Ptr&)>(ArithmeticCrossover);
    }
    void GeneticAlgorithm::FloatingPointPresentation::SetCustomCrossoverFunction(const CrossoverFunction& crossoverFunction) {
        m_CrossoverFunction = crossoverFunction;
    }
    void GeneticAlgorithm::FloatingPointPresentation::SetUniformMutationFunction() {
        m_MutationFunction = UniformMutation;
    }
    void GeneticAlgorithm::FloatingPointPresentation::SetCustomMutationFunction(const MutationFunction& mutationFunction) {
        m_MutationFunction = mutationFunction;
    }
    std::pair<GeneticAlgorithm::AbstractUnit::Ptr, GeneticAlgorithm::AbstractUnit::Ptr> GeneticAlgorithm::FloatingPointPresentation::Crossover(const AbstractUnit::Ptr& firstParentUnit, const AbstractUnit::Ptr& secondParentUnit) const {
        return m_CrossoverFunction(firstParentUnit, secondParentUnit);
    }
    void GeneticAlgorithm::FloatingPointPresentation::Mutate(const AbstractUnit::Ptr& unit) const {
        m_MutationFunction(unit);
    }
    GeneticAlgorithm::AbstractUnit::Ptr GeneticAlgorithm::FloatingPointPresentation::NewUnit(const AbstractPresentation::Ptr& presentation) const {
        return AbstractUnit::Ptr(new FloatingPointUnit(presentation));
    }

    GeneticAlgorithm::FloatingPointUnit::FloatingPointUnit(const AbstractPresentation::Ptr& presentation) : AbstractUnit(presentation) {
        assert(dynamic_cast<FloatingPointPresentation*>(presentation.get()));
    }
    GeneticAlgorithm::FloatingPointUnit::~FloatingPointUnit() { }
    double GeneticAlgorithm::FloatingPointUnit::GetGeneReal(std::uint8_t index) const {
        assert(index < GetPresentation()->GetNumberOfGenes());
        return reinterpret_cast<const double*>(GetGenes())[index];
    }
    void GeneticAlgorithm::FloatingPointUnit::SetGeneReal(std::uint8_t index, double real) {
        assert(index < GetPresentation()->GetNumberOfGenes());
        reinterpret_cast<double*>(GetGenes())[index] = real;
    }

    GeneticAlgorithm::GeneticAlgorithm(std::size_t populationSize, const AbstractPresentation::Ptr& presentation, double mutationChance, std::size_t maxGoalFunctionEvaluations) {
        SetPopulationSize(populationSize);
        SetPresentation(presentation);
        SetMutationChance(mutationChance);
        SetMaxGoalFunctionEvaluations(maxGoalFunctionEvaluations);
    }
    GeneticAlgorithm::~GeneticAlgorithm() { }
    std::size_t GeneticAlgorithm::GetPopulationSize() const {
        return m_PopulationSize;
    }
    void GeneticAlgorithm::SetPopulationSize(std::size_t populationSize) {
        assert(populationSize >= 3);
        m_PopulationSize = populationSize;
    }
    const GeneticAlgorithm::AbstractPresentation::Ptr& GeneticAlgorithm::GetPresentation() const {
        return m_Presentation;
    }
    void GeneticAlgorithm::SetPresentation(const AbstractPresentation::Ptr& presentation) {
        m_Presentation = presentation;
    }
    double GeneticAlgorithm::GetMutationChance() const {
        return m_MutationChance;
    }
    void GeneticAlgorithm::SetMutationChance(double mutationChance) {
        assert((mutationChance >= 0.0) && (mutationChance <= 1.0));
        m_MutationChance = mutationChance;
    }
    std::size_t GeneticAlgorithm::GetMaxGoalFunctionEvaluations() const {
        return m_MaxGoalFunctionEvaluations;
    }
    void GeneticAlgorithm::SetMaxGoalFunctionEvaluations(std::size_t maxGoalFunctionEvaluations) {
        assert(maxGoalFunctionEvaluations > 0);
        m_MaxGoalFunctionEvaluations = maxGoalFunctionEvaluations;
    }
    std::vector<double> GeneticAlgorithm::Solve(const GoalFunction& goalFunction, const LogFunction& logFunction) {
        double lowestBound = m_Presentation->GetLowerBound(0);
        double highestBound = m_Presentation->GetUpperBound(0);
        std::uint8_t numberOfGenes = m_Presentation->GetNumberOfGenes();
        for (std::uint8_t index = 1; index < numberOfGenes; ++index) {
            lowestBound = std::min(lowestBound, m_Presentation->GetLowerBound(index));
            highestBound = std::max(highestBound, m_Presentation->GetUpperBound(index));
        }
        return Solve(goalFunction, std::vector<double>(numberOfGenes, (lowestBound + highestBound) / 2.0), (highestBound - lowestBound) / 2.0, logFunction);
    }
    std::vector<double> GeneticAlgorithm::Solve(const GoalFunction& goalFunction, const std::vector<double>& startingPoint, double initialOffset, const LogFunction& logFunction) {
        std::uint8_t numberOfGenes = m_Presentation->GetNumberOfGenes();
        assert(startingPoint.size() == numberOfGenes);
        std::vector<std::pair<AbstractUnit::Ptr, double>> population;
        AbstractUnit::Ptr firstUnit = m_Presentation->NewUnit(m_Presentation);
        for (std::uint8_t geneIndex = 0; geneIndex < numberOfGenes; ++geneIndex)
            firstUnit->SetGeneReal(geneIndex, startingPoint[geneIndex]);
        firstUnit->Clamp();
        population.push_back(std::pair<AbstractUnit::Ptr, double>(firstUnit, Fitness(goalFunction, firstUnit)));
        for (std::size_t unitIndex = 1; unitIndex < m_PopulationSize; ++unitIndex) {
            AbstractUnit::Ptr unit = m_Presentation->NewUnit(m_Presentation);
            for (std::uint8_t geneIndex = 0; geneIndex < numberOfGenes; ++geneIndex)
                unit->SetGeneReal(geneIndex, startingPoint[geneIndex] + initialOffset * ((std::rand() % 20001) / 10000.0 - 1.0));
            unit->Clamp();
            population.push_back(std::pair<AbstractUnit::Ptr, double>(unit, Fitness(goalFunction, unit)));
        }
        std::size_t goalFunctionEvaluations = m_PopulationSize;
        std::vector<double> x;
        while (goalFunctionEvaluations < m_MaxGoalFunctionEvaluations) {
            const std::pair<AbstractUnit::Ptr, double>* bestUnitPair = &(population[0]);
            std::size_t bestUnitPairIndex = 0;
            for (std::size_t pairIndex = 1; pairIndex < m_PopulationSize; ++pairIndex)
                if (population[pairIndex].second < bestUnitPair->second) {
                    bestUnitPair = &(population[pairIndex]);
                    bestUnitPairIndex = pairIndex;
                }
            x.clear();
            for (std::uint8_t geneIndex = 0; geneIndex < numberOfGenes; ++geneIndex)
                x.push_back(bestUnitPair->first->GetGeneReal(geneIndex));
            if (bestUnitPair->second < EPSILON)
                return x;
            logFunction(x, population, bestUnitPairIndex, goalFunctionEvaluations);
            std::size_t selected[3];
            selected[0] = std::rand() % m_PopulationSize;
            selected[1] = std::rand() % (m_PopulationSize - 1);
            if (selected[1] == selected[0])
                selected[1]++;
            selected[2] = std::rand() % (m_PopulationSize - 2);
            if ((selected[2] == selected[0]) || (selected[2] == selected[1]))
                selected[2]++;
            if ((selected[2] == selected[0]) || (selected[2] == selected[1]))
                selected[2]++;
            std::pair<AbstractUnit::Ptr, double>* worstSelectedUnitPair;
            const std::pair<AbstractUnit::Ptr, double>* firstParentUnitPair;
            const std::pair<AbstractUnit::Ptr, double>* secondParentUnitPair;
            if ((population[selected[0]].second >= population[selected[1]].second) && (population[selected[0]].second >= population[selected[2]].second)) {
                worstSelectedUnitPair = &(population[selected[0]]);
                firstParentUnitPair = &(population[selected[1]]);
                secondParentUnitPair = &(population[selected[2]]);
            } else if ((population[selected[1]].second >= population[selected[0]].second) && (population[selected[1]].second >= population[selected[2]].second)) {
                worstSelectedUnitPair = &(population[selected[1]]);
                firstParentUnitPair = &(population[selected[0]]);
                secondParentUnitPair = &(population[selected[2]]);
            } else {
                worstSelectedUnitPair = &(population[selected[2]]);
                firstParentUnitPair = &(population[selected[0]]);
                secondParentUnitPair = &(population[selected[1]]);
            }
            std::pair<AbstractUnit::Ptr, AbstractUnit::Ptr> children = m_Presentation->Crossover(firstParentUnitPair->first, secondParentUnitPair->first);
            children.first->Clamp();
            children.second->Clamp();
            double firstChildFitness = Fitness(goalFunction, children.first);
            double secondChildFitness = Fitness(goalFunction, children.second);
            goalFunctionEvaluations += 2;
            worstSelectedUnitPair->first = (firstChildFitness <= secondChildFitness) ? children.first : children.second;
            worstSelectedUnitPair->second = (firstChildFitness <= secondChildFitness) ? firstChildFitness : secondChildFitness;
            if (((std::rand() % 10001) / 10000.0) <= m_MutationChance) {
                m_Presentation->Mutate(worstSelectedUnitPair->first);
                worstSelectedUnitPair->first->Clamp();
                worstSelectedUnitPair->second = Fitness(goalFunction, worstSelectedUnitPair->first);
                goalFunctionEvaluations++;
            }
        }
        return x;
    }
    double GeneticAlgorithm::Fitness(const GoalFunction& goalFunction, const AbstractUnit::Ptr& unit) {
        assert(m_Presentation.get() == unit->GetPresentation().get());
        std::vector<double> x;
        uint8_t numberOfGenes = m_Presentation->GetNumberOfGenes();
        for (uint8_t index = 0; index < numberOfGenes; ++index)
            x.push_back(unit->GetGeneReal(index));
        return goalFunction(x);
    }

}
