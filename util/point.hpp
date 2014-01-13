/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef POINT_H_
#define POINT_H_
#include "util.hpp"
template<class T>
class Point
{
	friend class PartFilter;
	friend class PartFilterEmbedded;
	friend class HOGParts;
	friend class Filter;
public:
	Point():x(0),y(0){}
	Point(const T&,const T&);
	Point(const Point<T> &);
	void Initialize(const T&,const T&);
	void operator=(const Point<T>&);
	Point<T> operator+(const Point<T>&);
	bool operator<(const Point<T>&);
	bool operator==(const Point<T>&);
	bool operator>(const Point<T>&);
	Point<T> operator*(const T&); // Scalar * Product...
	Point<T> operator/(const T&); // Point Wise Division...
	T operator*(const Point<T>&); // Dot Product...
	Point operator-(const Point<T>&);
	void GetCoord(T&,T&);
	void SetCoord(T,T);
	void Round(const T& ,const T&);
	void Display();
	void WriteText(ofstream &);
	bool Write(ofstream &);
	bool Read(ifstream &);
	T GetX(){return x;}
	T GetY(){return y;}
	static bool CompareDescendingX(const Point<T>&,const Point<T>&);
	static bool CompareDescendingY(const Point<T>&,const Point<T>&);
	static bool CompareAscendingY(const Point<T>&,const Point<T>&);
private:
	T x,y;
};
//#include "point.h"
template<class T>
Point<T>::Point(const T & x_,const T & y_)
{
	x = x_;
	y =y_;
}
template<class T>
Point<T> :: Point(const Point<T> & pobj)
{
	*this = pobj;
}

template<class T>
void Point<T>::operator=(const Point<T> & pobj)
{
	x = pobj.x;
	y = pobj.y;
}
template<class T>
bool Point<T>::operator<(const Point<T> & pobj)
{
	if(pobj.x < x && pobj.y < y)
		return true;
	return false;
}
template<class T>
bool Point<T>::operator>(const Point<T> & pobj)
{
	if(pobj.x > x && pobj.y > y)
		return true;
	return false;
}
template<class T>
Point<T> Point<T>::operator+(const Point<T> & pobj)
{
	Point tobj;
	tobj.x = x + pobj.x;
	tobj.y = y + pobj.y;

	return tobj;
}
template<class T>
Point<T> Point<T>::operator-(const Point<T> & pobj)
{
	Point tobj;
	tobj.x = x - pobj.x;
	tobj.y = y - pobj.y;

	return tobj;
}
template<class T>
Point<T> Point<T>::operator*(const T & scale)
{
	Point tobj;
	tobj.x = scale * x;
	tobj.y = scale * y;

	return tobj;
}
template<class T>
Point<T> Point<T>::operator/(const T & scale)
{
	Point tobj;
	tobj.x = x/scale;
	tobj.y = y/scale;

	return tobj;
}

template<class T>
T Point<T>::operator*(const Point<T> & pobj)
{
	T tobj = x * pobj.x + y * pobj.y ;
	return tobj;
}

template<class T>
void Point<T>::GetCoord(T &x_,T &y_)
{
	x_ = x;
	y_ = y;
}
template<class T>
void Point<T>::SetCoord(T x_,T y_)
{
	x = x_;
	y = y_;
}

template<class T>
void Point<T>::Round(const T &xcell,const T& ycell)
{
	x = x % xcell == 0 ? x :( x % xcell < xcell/2 ?  x - x % xcell: x + (xcell - x % xcell) ) ;
	y = y % ycell == 0 ? y :( y % ycell < ycell/2 ?  y - y % ycell: y + (ycell - y % ycell) ) ;
}
template<class T>
void Point<T>::Initialize(const T& x_,  const T&y_)
{
	x = x_;
	y = y_;
}
template<class T>
void Point<T>::Display()
{
	cout<< " X = " << x <<" Y= " << y <<endl;
}
template<class T>
bool Point<T>::operator==(const Point<T> & pobj)
{
	return pobj.x == x && pobj.y == y ;
}
template<class T>
void Point<T>::WriteText(ofstream &ofile)
{
	ofile<< " X = "<<x << " Y = " << y<<endl;
}
template<class T>
bool Point<T>::CompareDescendingX(const Point<T>& p1,const Point<T>&p2)
{
	return p1.x > p2.x;
}
template<class T>
bool Point<T>::CompareDescendingY(const Point<T>& p1,const Point<T>&p2)
{
	return p1.y > p2.y;
}
template<class T>
bool Point<T>::CompareAscendingY(const Point<T>& p1,const Point<T>&p2)
{
	return p1.y < p2.y;
}
template<class T>
bool Point<T>::Write(ofstream &ofile)
{
	ofile.write((char*)&x,sizeof(T));
	ofile.write((char*)&y,sizeof(T));
	if(ofile)
		return true;
	return false;
}
template<class T>
bool Point<T>::Read(ifstream &ifile)
{
	if(ifile)
	{
		ifile.read((char*)&x,sizeof(T));
		ifile.read((char*)&y,sizeof(T));
		return true;
	}
	return false;
}

#endif /*POINT_H_*/

