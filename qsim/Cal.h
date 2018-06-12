#pragma once

#include <string>
#include <map>
#include <complex>

typedef std::complex<double> CCom;

struct CalExpr;
struct Cal {

	Cal(char const *x);

	void SetVarVal(std::string const &name, CCom const &v);
	CCom GetVarVal(std::string const &name);
	CCom Val();

private:
	std::map<std::string, CCom> fVarVals;
	CalExpr *fExpr;
	std::string fStr;
};


void testCal();