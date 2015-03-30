#include "GaussPrior3D.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>
#include <gsl/gsl_cdf.h>

using namespace DNest3;

GaussPrior3D::GaussPrior3D(double t_min, double t_max)
:t_min(t_min)
,t_max(t_max)
,t_range(t_max - t_min)
{

}

void GaussPrior3D::fromPrior()
{
	mean_logA = tan(M_PI*(0.97*randomU() - 0.485));
	mean_logA = exp(mean_logA);

	mean_logDuration = log(1E-3*t_range) + log(1E3)*randomU();
//	mean_logSkew = 
}

double GaussPrior3D::perturb_parameters()
{
	double logH = 0.;

	int which = randInt(4);

	if(which == 0)
	{
		mean_logA = log(mean_logA);
		mean_logA = (atan(mean_logA)/M_PI + 0.485)/0.97;
		mean_logA += pow(10., 1.5 - 6.*randomU())*randn();
		mean_logA = mod(mean_logA, 1.);
		mean_logA = tan(M_PI*(0.97*mean_logA - 0.485));
		mean_logA = exp(mean_logA);
	}

	return logH;
}

/*
* vec[0] = spike time
* vec[1] = log(amplitude)
* vec[2] = log(rise time)
* vec[3] = log(s') where (old skew) = s' - 1
*/

double GaussPrior3D::log_pdf(const std::vector<double>& vec) const
{
	double logP = 0.;

	if(vec[0] < t_min || vec[0] > t_max)
		return -1E250;

	logP += -log(sig_logA) - 0.5*pow((vec[1] - mean_logA)/sig_logA, 2);
	logP += -log(sig_logDuration) - 0.5*pow((vec[2] - mean_logDuration - co_ADuration*vec[1])/sig_logDuration, 2);
	logP += -log(sig_logSkew) - 0.5*pow((vec[3] - mean_logSkew - co_ASkew*vec[1] - co_durationSkew*vec[2])/sig_logSkew, 2);

	return 0.;
}

void GaussPrior3D::from_uniform(std::vector<double>& vec) const
{
	vec[0] = t_min + (t_max - t_min)*vec[0];
	vec[1] = mean_logA + sig_logA*gsl_cdf_ugaussian_Pinv(vec[1]);
	vec[2] = mean_logDuration + co_ADuration*vec[1]
			+ sig_logDuration*gsl_cdf_ugaussian_Pinv(vec[2]);
	vec[3] = mean_logSkew + co_ASkew*vec[1] + co_durationSkew*vec[2]
			+ sig_logSkew*gsl_cdf_ugaussian_Pinv(vec[3]);
}

void GaussPrior3D::to_uniform(std::vector<double>& vec) const
{
	vec[3] = gsl_cdf_ugaussian_P((vec[3] - mean_logSkew - co_ASkew*vec[1] - co_durationSkew*vec[2])/sig_logSkew);
	vec[2] = gsl_cdf_ugaussian_P((vec[2] - mean_logDuration - co_ADuration*vec[1])/sig_logDuration);
	vec[1] = gsl_cdf_ugaussian_P((vec[1] - mean_logA)/sig_logA);
	vec[0] = (vec[0] - t_min)/(t_max - t_min);
}

void GaussPrior3D::print(std::ostream& out) const
{
	out<<" ";
}

