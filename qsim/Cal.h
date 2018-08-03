#pragma once

#include <string>
#include <map>
#include <complex>
#include <vector>

typedef std::complex<double> CCom;

struct CalExpr;
struct Cal {

	Cal(char const *x);
	~Cal();

	void SetVarVal(std::string const &name, CCom const &v);
	CCom &GetVarVal(std::string const &name);
	CCom Val();
	void GenPseudoCode();
	CCom RunPseudoCode();
private:
	friend struct CalExpr;
	std::map<std::string, CCom> fVarVals;
	std::shared_ptr<CalExpr> fExpr;
	std::string fStr;

	std::vector<int> fPseudoCode;
	std::vector<CCom> fStack;
	std::vector<CCom> fConstants;
	std::vector<CCom const*> fVariables;

};


void testCal();