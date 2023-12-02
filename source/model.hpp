#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <cstdint>
#include <string>
#include <vector>

#include "Predefines.hpp"
#include "IRandomGenerator.hpp"

enum DgemBinType
{
    POINT,
    LINE,
    HEX_LINES
};

struct RandomGeneratorDescription
{
    std::int32_t seed_{0};
    std::uint64_t numberOfPoints_{0};
    RandomGeneratorType type_{ RandomGeneratorType::NONE };
    bool vectorPerScattering_{ true };
    std::string inputRandomFile_;
    std::string outputRandomFile_;
};

// model parameters
class Model
{
	public:
		static Model& instance(std::vector<Observer> *observers, std::string const& parametersFile)
		{
			static Model mod(observers, parametersFile) ;
			return mod;
		}
		bool fMonteCarlo() const
		{	return fMonteCarlo_;	}
		bool writeScatterings() const
        {   return writeScatterings_;   }
        bool useHEALPixGrid() const
        {   return useHEALPixGrid_; }
		double taumin() const
		{	return taumin_;	}
        double defaultStarRadius() const
        {	return defaultStarRadius_;	}
        DgemBinType dgemBinType() const
        {   return dgemBinType_;    }
		std::uint32_t nscat() const
		{	return nscat_;	}

		IRandomGenerator* createRandomGenerator() const;
        IRandomGenerator* createDgemRandomGenerator() const;

		std::uint32_t SecondaryDirectionsLevel() const
		{	return SecondaryDirectionsLevel_;	}
		std::uint32_t MonteCarloStart() const
		{	return MonteCarloStart_;	}
		std::uint32_t NumOfPrimaryScatterings() const
		{	return NumOfPrimaryScatterings_;	}
		std::uint32_t NumOfSecondaryScatterings() const
		{	return NumOfSecondaryScatterings_;	}
        IGridCPtr grid() const
        {   return grid_; }
		IDustCPtr dust() const
        {   return dust_; }
        SourcesPtr sources()
        {   return sources_; }
	private:
		bool fMonteCarlo_{ true };
		bool useHEALPixGrid_{ true };
		bool writeScatterings_{ false };
		DgemBinType dgemBinType_{ DgemBinType::LINE };
		double taumin_{ 0.0 };
        double defaultStarRadius_{ 0.0 };

        std::uint32_t nscat_{ 1 };
		std::uint32_t MonteCarloStart_{ 1 };
		std::uint32_t SecondaryDirectionsLevel_{ 1 };
        std::uint32_t NumOfPrimaryScatterings_{ 0 };
		std::uint32_t NumOfSecondaryScatterings_{ 1 };

		IDustCPtr dust_;
        IGridPtr grid_;
        SourcesPtr sources_;

		// Random generator parameters
		RandomGeneratorDescription monteCarloGeneratorDescription_;
		RandomGeneratorDescription dgemStratificationGeneratorDescription_;

		Model(std::vector<Observer> *observers, std::string const& parametersFile);
		Model(Model const &);
		Model& operator=( Model const &);
};

#endif
