#pragma once

#include <string>
#include <map>
#include <complex>
#include <vector>
#include <memory>

typedef std::complex<double> CCom;

struct CalExpr;

// Usage:
// 
// 
//   Cal cal("1+x");               // step 1: set expresion
//   cal.SetVarVal("x", 1);        // step 2: set variables
//   CCom result = cal.Val();      // step 3: evalute value
//
//   cal.SetVarVal("x", 2);        // you can reset variable value
//   result = cal.Val();           // then re-evaluate
//
//   CCom &x = cal.GetVarVal("x"); // you can take address of variable
//   for(int i = 0; i < 10; ++i) { 
//       x = 1.0*i;                // then reset variable value
//       result = cal.Val()        // and re-evaluate
//   }
//
//   // `SetVarVal` and `GetVarVal` itself is expensive, However reseting value by the reference from GetVarVal is cheap
//   // The first call of `Val` after `SetVarVal` or `constructor` may be expensive
//   // you must call SetVarVal before you call GetVarVal

//  support functions or operators:
//  a ? b : c  =  if a != 0 then a else b)
//     a || b  =  if a != 0 || b != 0 then 1 else 0
//     a && b  =  if a != 0 && b != 0 then 1 else 0
//      a > b  =  if real(a) > real(b) then 1 else 0
//     a >= b  =  if real(a) >= real(b) then 1 else 0
//      a < b  =  if real(a) < real(b) then 1 else 0
//     a <= b  =  if real(a) <= real(b) then 1 else 0
//     a == b  =  if a == b then 1 else 0
//     + - * /
//      a ^ b  =  power(a, b)
//      exp tanh sinh log atanh
//      sin cos tan asin acos atan
//       sign  =  if Re x > 0  then 1 else 0
//       abs   =  sqrt((Re x)^2 + (Im x)^2)   
//         -a (unitary)
//       sqrt
//     guass(x, mu, sigma)

struct Cal {

	// Codes are stil not avaliable
	Cal(char const *x);

	// Alloc and set variable
	// This operation will make codes not avaliable
	void SetVarVal(std::string const &name, CCom const &v);
	// Get variable (address), you must alloc and set variable before you get variable
	CCom &GetVarVal(std::string const &name);
	// Generate (if codes not avaliable) and run codes
	CCom Val();

	~Cal();


	// Generate codes
	// Codes will be avaliable
	void GenPseudoCode();
	// Run codes, Codes must be avaliable
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