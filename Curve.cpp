//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
// Copyright (c) 2015 Mahyar Khayatkhoei
//
//this cpp file is rewirten by Yuyang Liu
//all the functions are filled by Yuyang Liu

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI

	if (!checkRobust()) {
		return ;
	}


	int size = controlPoints.size();
	float lastPointTime = controlPoints[size - 1].time;
	for (float i = 0; i < lastPointTime; i = i+ window){
		Point p1;
		Point p2;
		calculatePoint(p1, i);
		if (i + window < lastPointTime){
			calculatePoint(p2, i + window);
		}
		else{
			calculatePoint(p2, lastPointTime);
		}
		DrawLib::drawLine(p1, p2, curveColor, curveThickness);
	}

	// Robustness: make sure there is at least two control point: start and end points
	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	// Note that you must draw the whole curve at each frame, that means connecting line segments between each two points on the curve
	
	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
/**
this function is used to sort all the control points in a accending order
the data structure of 'controlpoints' is vector<points>
so we can use C++ sort function
**/
bool ComparePoint(CurvePoint p1, CurvePoint p2)
{
	if (p1.time <= p2.time)
	{
		return true;
	}
	else
		return false;
}
void Curve::sortControlPoints()
{
	std::sort(controlPoints.begin(), controlPoints.end(), ComparePoint);
	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
// Note that this function should return false if the end of the curve is reached, or no next point can be found
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	// Note that nextPoint is an integer containing the index of the next control point
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve given the next control point (nextPoint)
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
/**
I guess this is used to check if the controlpoints are suffcient
then for the HermiteCurve there shold be at least 2 points
for the CatmullCurve there should be at least 3
**/
bool Curve::checkRobust()
{
	if (type == hermiteCurve)
	{
		if (controlPoints.size() >= 2)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else if (type == catmullCurve)
	{
		if (controlPoints.size() >= 3)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else

	return false;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
/**
this is used to find out the index of next poin in vector
**/
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	//used sorted array to find the next
	sortControlPoints();
	for (int i = 0; i < controlPoints.size(); i++)
	{
		if (controlPoints[i].time > time) {
			nextPoint = i;
			return true;
		}
	}
	// run all the t, but still not find next
	//end 
	return false;

}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;

	float T; //the total time from last to next
	float t; // the t

	int last = nextPoint - 1;		//the last point

	T = controlPoints[nextPoint].time - controlPoints[last].time;
	t = (time - controlPoints[last].time) / T;

	// Calculate position at t = time on Hermite curve
	newPosition.x = (controlPoints[last].position.x * (2 * t*t*t - 3 * t*t + 1))
		+ (controlPoints[last].tangent.x * (t*t*t - 2 * t*t + t) * T)
		+ (controlPoints[nextPoint].position.x * (-2 * t*t*t + 3 * t*t))
		+ (controlPoints[nextPoint].tangent.x * (t*t*t - t*t) * T);

	newPosition.y = (controlPoints[last].position.y * (2 * t*t*t - 3 * t*t + 1))
		+ (controlPoints[last].tangent.y * (t*t*t - 2 * t*t + t)* T)
		+ (controlPoints[nextPoint].position.y * (-2 * t*t*t + 3 * t*t))
		+ (controlPoints[nextPoint].tangent.y * (t*t*t - t*t)* T);

	newPosition.z = (controlPoints[last].position.z * (2 * t*t*t - 3 * t*t + 1))
		+ (controlPoints[last].tangent.z * (t*t*t - 2 * t*t + t)* T)
		+ (controlPoints[nextPoint].position.z * (-2 * t*t*t + 3 * t*t))
		+ (controlPoints[nextPoint].tangent.z * (t*t*t - t*t)* T);
	// Return result
	return newPosition;

}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{

	Point newPosition;

	Point R0;
	Point R1;
	if (nextPoint  == 1) //it is the second CP
	{
		Point p0 = controlPoints[0].position;
		Point p1 = controlPoints[1].position;

		R0.x = p1.x - p0.x;
		R0.y = p1.y - p0.y;
		R0.z = p1.z - p0.z;
	}
	else {
		Point p0 = controlPoints[nextPoint - 2].position;
		Point p2 = controlPoints[nextPoint ].position;
		R0.x = (p2.x - p0.x)/2;
		R0.y =(p2.y - p0.y)/2;
		R0.z = (p2.z - p0.z)/2;
	}
	if (nextPoint == 0) //it is the first
	{
		Point p0 = controlPoints[0].position;
		Point p1 = controlPoints[1].position;

		R1.x = p1.x - p0.x;
		R1.y = p1.y - p0.y;
		R1.z = p1.z - p0.z;
	}
	else if (nextPoint == controlPoints.size() - 1) {//it is the last
		Point p1 = controlPoints[nextPoint].position;
		Point p0 = controlPoints[nextPoint - 1].position;

		R1.x = p1.x - p0.x;
		R1.y = p1.y - p0.y;
		R1.z = p1.z - p0.z;
	}
	else {//if p1 is point in the middle
		Point p0 = controlPoints[nextPoint - 1].position;
		Point p2 = controlPoints[nextPoint + 1].position;

		R1.x = (p2.x - p0.x)/2;
		R1.y = (p2.y - p0.y)/2;
		R1.z = (p2.z - p0.z)/2;
	}


	//2.calculate equation H(t)
	float TimeInterval = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	float t = (time - controlPoints[nextPoint - 1].time) / TimeInterval;

	Point p0 = controlPoints[nextPoint - 1].position;
	Point p1 = controlPoints[nextPoint].position;
	newPosition.x = (2.0f * t * t * t - 3.0f * t * t + 1.0f) * p0.x
		+ (t * t * t - 2.0f * t * t + t) * R0.x
		+ (-2.0f * t * t * t + 3.0f * t * t) * p1.x
		+ (t * t * t - t * t) * R1.x;
	newPosition.y = (2.0f * t * t * t - 3.0f * t * t + 1.0f) * p0.y
		+ (t * t * t - 2.0f * t * t + t) * R0.y
		+ (-2.0f * t * t * t + 3.0f * t * t) * p1.y
		+ (t * t * t - t * t) * R1.y;
	newPosition.z = (2.0f * t * t * t - 3.0f * t * t + 1.0f) * p0.z
		+ (t * t * t - 2.0f * t * t + t) * R0.z
		+ (-2.0f * t * t * t + 3.0f * t * t) * p1.z
		+ (t * t * t - t * t) * R1.z;
	// Return result
	return newPosition;
	

}