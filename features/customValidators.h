/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef CUSTOMVALIDATORS_H_
#define CUSTOMVALIDATORS_H_
#include <boost/program_options.hpp>
#include "../util/util.hpp"
#include "definitions.h"
// define a custom validators to be used with program_options
//These custom validators can be avoided by parsing the values explicitly afterwards like in parse Image-Processor function..
//Specialized templates can be used but are avoided
void validate(boost::any& v, const std::vector<std::string>& values,
		PyramidType* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		NORM* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		GridType* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		HistogramMethod* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		PatchType* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		PatchEncodingType* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		CBType* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		WindowFeatureTypes* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		HistFeatures* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		FeatureTypes* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		PoolingMethod* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		EncodedHistogramNormType* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		QuantizationType* target_type, int);
void validate(boost::any& v, const std::vector<std::string>& values,
		FeatVectorNormalType* target_type, int);
#endif /* DEFINITIONS_H_ */
