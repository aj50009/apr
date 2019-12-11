#ifndef __APR_GENETIC_ALGORITHM_HXX__
#define __APR_GENETIC_ALGORITHM_HXX__

#include <initializer_list>
#include <utility>
#include <memory>
#include <vector>
#include <string>
#include <ostream>
#include <istream>
#include <cstdint>

namespace apr {

    class GeneticAlgorithm {
    public:
        using Ptr = std::shared_ptr<GeneticAlgorithm>;

        class AbstractUnit;
        class AbstractPresentation {
        public:
            using Ptr = std::shared_ptr<AbstractPresentation>;
            AbstractPresentation(std::uint8_t numberOfGenes, const std::initializer_list<double>& lowerBounds, const std::initializer_list<double>& upperBounds);
            virtual ~AbstractPresentation();
            std::uint8_t GetNumberOfGenes() const;
            double GetLowerBound(std::uint8_t index) const;
            void SetLowerBound(std::uint8_t index, double lowerBound);
            double GetUpperBound(std::uint8_t index) const;
            void SetUpperBound(std::uint8_t index, double upperBound);
            virtual std::pair<std::shared_ptr<AbstractUnit>, std::shared_ptr<AbstractUnit>> Crossover(const std::shared_ptr<AbstractUnit>& firstParentUnit, const std::shared_ptr<AbstractUnit>& secondParentUnit) const = 0;
            virtual void Mutate(const std::shared_ptr<AbstractUnit>& unit) const = 0;
        private:
            std::uint8_t m_NumberOfGenes;
            std::vector<double> m_LowerBounds;
            std::vector<double> m_UpperBounds;
        };

        class AbstractUnit {
        public:
            using Ptr = std::shared_ptr<AbstractUnit>;
            AbstractUnit(const AbstractPresentation::Ptr& presentation);
            virtual ~AbstractUnit();
            const AbstractPresentation::Ptr& GetPresentation() const;
            virtual double GetGeneReal(std::uint8_t index) const = 0;
            virtual void SetGeneReal(std::uint8_t index, double real) = 0;
            virtual void ClampGene(std::uint8_t index);
            virtual void Clamp();
            std::string ToString() const;
            void FromString(const std::string& str);
        protected:
            virtual std::string GetGeneString(std::uint8_t index) const;
            virtual void SetGeneString(std::uint8_t index, const std::string& str);
            virtual void WriteToOutputStream(std::ostream& outputStream) const;
            virtual void ReadFromInputStream(std::istream& inputStream);
            const std::uint64_t* GetGenes() const;
            std::uint64_t* GetGenes();
        private:
            friend std::ostream& operator<<(std::ostream&, const AbstractUnit&);
            friend std::ostream& operator<<(std::ostream&, const Ptr&);
            friend std::istream& operator>>(std::istream&, AbstractUnit&);
            friend std::istream& operator>>(std::istream&, const Ptr&);
            AbstractPresentation::Ptr m_Presentation;
            std::uint64_t* m_Genes;
        };
    };

    inline std::ostream& operator<<(std::ostream& outputStream, const GeneticAlgorithm::AbstractUnit& unit) { unit.WriteToOutputStream(outputStream); return outputStream; }
    inline std::ostream& operator<<(std::ostream& outputStream, const GeneticAlgorithm::AbstractUnit::Ptr& unit) { unit->WriteToOutputStream(outputStream); return outputStream; }
    inline std::istream& operator>>(std::istream& inputStream, GeneticAlgorithm::AbstractUnit& unit) { unit.ReadFromInputStream(inputStream); return inputStream; }
    inline std::istream& operator>>(std::istream& inputStream, const GeneticAlgorithm::AbstractUnit::Ptr& unit) { unit->ReadFromInputStream(inputStream); return inputStream; }

}

#endif
