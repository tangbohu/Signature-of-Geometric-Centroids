
#ifndef _SGC_BIN_H
#define _SGC_BIN_H

#include "math3D.h"
#include "HelpDefines.h"
#include<iostream>

using namespace std;

class SGCBin
{
public:
	double x;
	double y;
	double z;
	int numberOfPoints;

public:
	SGCBin()
	{
		x = y = z = 0;
		numberOfPoints = 0;
	}
	SGCBin(double _x, double _y, double _z, int n)
	{
		x = _x;
		y = _y;
		z = _z;
		numberOfPoints = n;
	}

	void setValue(double _x, double _y, double _z, int n)
	{
		x = _x;
		y = _y;
		z = _z;
		numberOfPoints = n;
	}

	void reset()
	{
		x = y = z = 0;
		numberOfPoints = 0;

	}
	void copy(const SGCBin& that)
	{
		this->x = that.x; this->y = that.y; this->z = that.z;
		this->numberOfPoints = that.numberOfPoints;
	}
	double EclidDistWith(const SGCBin& that)
	{
		return dist3DusePoints_2(this->x, this->y, this->z, that.x, that.y, that.z);
	}
	double distWith(const SGCBin& that)
	{
		return -1 * similarWith(that);
	}
	double similarWith(const SGCBin& that)
	{
		double similarity = 0;
		if (numberOfPoints>0 && that.numberOfPoints>0)
		{
			similarity = log10(numberOfPoints) + log10(that.numberOfPoints) - log10(_MAX(dist3DusePoints_2(this->x, this->y, this->z, that.x, that.y, that.z), 1e-6));
		//	similarity =  - log10(_MAX(dist3DusePoints_2(this->x, this->y, this->z, that.x, that.y, that.z), 1e-6));

			if (similarity < 0)
				cout << "------------------BiG Error------------" << similarity << endl;
		}
		return similarity;
	}

};

#endif