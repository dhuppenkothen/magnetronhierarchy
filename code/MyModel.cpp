/*
* Copyright (c) 2009, 2010, 2011, 2012 Brendon J. Brewer.
*
* This file is part of DNest3.
*
* DNest3 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DNest3 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DNest3. If not, see <http://www.gnu.org/licenses/>.
*/

#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>
#include <gsl/gsl_sf_gamma.h>

using namespace std;
using namespace DNest3;

const Data& MyModel::data = Data::get_instance();

MyModel::MyModel()
//:bursts(4, 100, false, ClassicMassInf1D(data.get_t_min(), data.get_t_max(),
//				1E-3*data.get_y_mean(), 1E3*data.get_y_mean()))
:bursts(4, 100, false, GaussPrior3D(data.get_t_min(), data.get_t_max()))
,mu(data.get_t().size())
{

}

void MyModel::calculate_mu()
{
	const vector<double>& t = data.get_t();

	// Update or from scratch?
	bool update = (bursts.get_added().size() < bursts.get_components().size());

	// Get the components
	const vector< vector<double> >& components = (update)?(bursts.get_added()):
				(bursts.get_components());

	// Set the background level
	if(!update)
		mu.assign(mu.size(), background);

	double amplitude, duration, skew;
	double rise, fall, scale;

	for(size_t j=0; j<components.size(); j++)
	{
		amplitude = exp(components[j][1]);
		duration = exp(components[j][2]);
		skew = exp(components[j][3]);

		rise = duration/(1. + skew);
		fall = rise*skew;

		for(size_t i=0; i<mu.size(); i++)
		{
			scale = (t[i] > components[j][0])?(fall):(rise);

			mu[i] += amplitude
					*exp(-fabs(t[i] - components[j][0])/
						scale);
		}
	}
}

void MyModel::fromPrior()
{
//	background = exp(log(1E-3) + log(1E3)*randomU())*data.get_y_mean();
    background = tan(M_PI*(0.97*randomU() - 0.485));
    background = exp(background);
    bursts.fromPrior();
	calculate_mu();
}

double MyModel::perturb()
{
	double logH = 0.;

	if(randomU() <= 0.2)
	{
		for(size_t i=0; i<mu.size(); i++)
			mu[i] -= background;
 
		background = log(background);
		background = (atan(background)/M_PI + 0.485)/0.97;
		background += pow(10., 1.5 - 6.*randomU())*randn();
		background = mod(background, 1.);
		background = tan(M_PI*(0.97*background - 0.485));
		background = exp(background);

		for(size_t i=0; i<mu.size(); i++)
			mu[i] += background;
	}
	else
	{
		logH += bursts.perturb();
		bursts.consolidate_diff();
		calculate_mu();
	}

	return logH;
}

double MyModel::logLikelihood() const
{
        const vector<double>& t = data.get_t();

	double t_start = data.get_t_min();
	double t_end = data.get_t_max();

	double time;
	double amp;
	double duration;
	double rise;
	double skew;
	double y0, y1;
	double mu_int = 0.;

	for(size_t j=0; j<bursts.get_components().size(); j++)
	{
		time = bursts.get_components()[j][0];
		amp = exp(bursts.get_components()[j][1]);
		duration = exp(bursts.get_components()[j][2]);
		skew = exp(bursts.get_components()[j][3]);

		rise = duration/(1. + skew);

		y0 = 1.0 - exp((t_start - time)/rise);
		y1 = skew - skew*exp(-(t_end - time)/(rise*skew));
		mu_int += amp*rise*(y0 + y1);
	}
	mu_int += background*(t_end - t_start);
 
	//const vector<int>& t = data.get_t();

	//double logl = 0.;
	//for(size_t i=0; i<t.size(); i++)
	//	logl += -mu[i] + y[i]*log(mu[i]) - gsl_sf_lngamma(y[i] + 1.);
        double logl = -mu_int;
        for(size_t i=0; i<t.size(); i++)
                logl += log(mu[i]);
 
	return logl;
}

void MyModel::print(std::ostream& out) const
{
	out<<background<<' ';
	bursts.print(out);

	for(size_t i=0; i<mu.size(); i++)
		out<<mu[i]<<' ';
}

string MyModel::description() const
{
	return string("");
}

