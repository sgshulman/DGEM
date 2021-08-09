#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <cstdint>
#include <string>
#include <vector>

#include "Predefines.hpp"

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
        bool useHEALPixGrid() const
        {   return useHEALPixGrid_; }
		double taumin() const
		{	return taumin_;	}
        double defaultStarRadius() const
        {	return defaultStarRadius_;	}
		std::uint32_t nscat() const
		{	return nscat_;	}

		Random createRandomGenerator() const;

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
		bool fMonteCarlo_;
		bool useHEALPixGrid_;
		double taumin_;
        double defaultStarRadius_;

		// Random generator parameters
		int32_t iseed_;
		std::string inputRandomFile_;
        std::string outputRandomFile_;

        std::uint32_t nscat_;
		std::uint32_t MonteCarloStart_;
		std::uint32_t SecondaryDirectionsLevel_;
        std::uint32_t NumOfPrimaryScatterings_;
		std::uint32_t NumOfSecondaryScatterings_;

		IDustCPtr dust_;
        IGridCPtr grid_;
        SourcesPtr sources_;

		Model(std::vector<Observer> *observers, std::string const& parametersFile);
		Model(Model const &);
		Model& operator=( Model const &);
};

#endif
