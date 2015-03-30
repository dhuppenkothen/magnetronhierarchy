#include "Data.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

Data Data::instance;

Data::Data()
{

}

void Data::load(const char* filename)
{
	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# Failed to open file "<<filename<<"."<<endl;
		return;
	}

	t.clear();

	double temp1, temp2;
	while(fin>>temp1 && fin>>temp2)
		t.push_back(temp1);

	fin.close();
	cout<<"# Found "<<t.size()<<" points in file "<<filename<<"."<<endl;

	compute_summaries();
}

void Data::compute_summaries()
{
	t_min = *min_element(t.begin(), t.end());
	t_max = *max_element(t.begin(), t.end());
	t_range = t_max - t_min;
}

