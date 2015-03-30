#ifndef _GaussPrior3D_
#define _GaussPrior3D_

#include <vector>
#include "Distributions/Distribution.h"

class GaussPrior3D:public Distribution
{
	private:
		// The means
		double mean_logA, mean_logTau, mean_logSkew;

		// The coefficients
		double co_ATau, co_ASkew, co_tauSkew;

		// The standard deviations
		// (conditional! for the latter two anyway)
		double sig_logA, sig_logTau, sig_logSkew;

		double perturb_parameters();

	public:
		GaussPrior3D(double x_min, double x_max);

		void fromPrior();

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;

		static const int weight_parameter = 1;
};

#endif

