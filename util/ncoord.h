/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef COORD_H_
#define COORD_H_
#include "util.hpp" // to cater for the signed coordinates...
class Coord {
public:
	Coord(int = 0, int = 0, int = 0, int = 0);
	Coord(const Coord&);
	Coord & operator=(const Coord &);
	int GetHeight() const {
		return (ymax - ymin);// TODO add +1 a bug
	}
	int GetWidth() const {
		return (xmax - xmin); // TODO add +1 a bug
	}
	void SetCoord(int, int, int, int);
	void GetCoord(int&, int&, int&, int&) const;
	void GetCoord(Coord &obj) {
		obj = *this;
	}
	const Coord &GetCoord() {
		return *this;
	}
	void SetCoord(const Coord &obj) {
		*this = obj;
	}
	REAL FindOverlap(const Coord &, REAL) const;
	REAL FindOverlap(int, int, int, int, REAL) const;
	REAL FindOverlap(int, int, int, int) const;
	REAL FindNMSOverlap(int, int, int, int) const;
	REAL FindOverlap(const Coord &) const;
	static REAL FindOverlap(const Coord &, const Coord&);
	UINT GetArea() {
		return GetWidth() * GetHeight();
	}
	void Write(ofstream &);
	void WriteText(ofstream &);
	void Read(ifstream &);
	REAL GetAspectRatio() {
		return ((REAL) GetWidth()) / GetHeight();
	}
	void Display();

protected:
	int xmin;
	int ymin, xmax, ymax;
};
inline void Coord::GetCoord(int& xmin_, int& ymin_, int& xmax_, int &ymax_) const {
	xmin_ = xmin;
	ymin_ = ymin;
	xmax_ = xmax;
	ymax_ = ymax;
}
inline Coord::Coord(int xmin_, int ymin_, int xmax_, int ymax_) {
	xmin = xmin_;
	ymin = ymin_;
	xmax = xmax_;
	ymax = ymax_;
}
inline Coord::Coord(const Coord& obj) {
	*this = obj;
}
inline Coord & Coord::operator=(const Coord & obj) {
	xmin = obj.xmin;
	ymin = obj.ymin;
	xmax = obj.xmax;
	ymax = obj.ymax;
	return *this;
}
inline void Coord::SetCoord(int xmin_, int ymin_, int xmax_, int ymax_) {
	xmin = xmin_;
	ymin = ymin_;
	xmax = xmax_;
	ymax = ymax_;
}

#endif /*COORD_H_*/
