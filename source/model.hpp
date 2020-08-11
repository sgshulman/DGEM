#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <cstdint>
#include <vector>

#include "Predefines.hpp"

// model parameters
class Model
{
	public:
		static Model& instance (std::vector<Observer> *observers)
		{
			static Model mod(observers) ;
			return mod ;
		}
		bool fMonteCarlo() const
		{	return fMonteCarlo_;	}

		double taumin() const
		{	return taumin_;	}
		std::uint32_t nscat() const
		{	return nscat_;	}
		int32_t iseed() const
		{	return iseed_;	}
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
		DustCPtr dust() const
        {   return dust_; }
        SourcesPtr sources()
        {   return sources_; }
	private:
		bool fMonteCarlo_;
		double taumin_;
		int32_t iseed_;

        std::uint32_t nscat_;
		std::uint32_t MonteCarloStart_;
		std::uint32_t SecondaryDirectionsLevel_;
        std::uint32_t NumOfPrimaryScatterings_;
		std::uint32_t NumOfSecondaryScatterings_;

		DustCPtr dust_;
        IGridCPtr grid_;
        SourcesPtr sources_;

		Model(std::vector<Observer> *observers);
		Model(Model const &);
		Model& operator=( Model const &);
};

#endif
