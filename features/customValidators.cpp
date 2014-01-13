/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#include "customValidators.h"
#include<iostream>
#include<vector>
#include<string>
void validate(boost::any& v, const std::vector<std::string>& values,
		PyramidType* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(PyramidType(boost::lexical_cast<int>(s)));
}
void validate(boost::any& v, const std::vector<std::string>& values,
		NORM* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(NORM(boost::lexical_cast<int>(s)));
}

void validate(boost::any& v, const std::vector<std::string>& values,
		GridType* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(GridType(boost::lexical_cast<int>(s)));
}

void validate(boost::any& v, const std::vector<std::string>& values,
		HistogramMethod* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(HistogramMethod(boost::lexical_cast<int>(s)));
}

void validate(boost::any& v, const std::vector<std::string>& values,
		PatchType* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(PatchType(boost::lexical_cast<int>(s)));
}

void validate(boost::any& v, const std::vector<std::string>& values,
		PatchEncodingType* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(PatchEncodingType(boost::lexical_cast<int>(s)));
}
void validate(boost::any& v, const std::vector<std::string>& values,
		CBType* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(CBType(boost::lexical_cast<int>(s)));
}
void validate(boost::any& v, const std::vector<std::string>& values,
		WindowFeatureTypes* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(WindowFeatureTypes(boost::lexical_cast<int>(s)));
}
void validate(boost::any& v, const std::vector<std::string>& values,
		HistFeatures* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(HistFeatures(boost::lexical_cast<int>(s)));
}
void validate(boost::any& v, const std::vector<std::string>& values,
		FeatureTypes* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(FeatureTypes(boost::lexical_cast<int>(s)));
}
void validate(boost::any& v, const std::vector<std::string>& values,
		PoolingMethod* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(PoolingMethod(boost::lexical_cast<int>(s)));
}
void validate(boost::any& v, const std::vector<std::string>& values,
		EncodedHistogramNormType* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(EncodedHistogramNormType(boost::lexical_cast<int>(s)));
}
void validate(boost::any& v, const std::vector<std::string>& values,
		QuantizationType* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(QuantizationType(boost::lexical_cast<int>(s)));
}
void validate(boost::any& v, const std::vector<std::string>& values,
		FeatVectorNormalType* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(FeatVectorNormalType(boost::lexical_cast<int>(s)));
}
