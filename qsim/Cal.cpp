#define _CRT_SECURE_NO_WARNINGS
#include "Cal.h"
#include <vector>
#include <complex>
#include <cmath>
#include <math.h>

enum CalExprType {
	CALE_UNK,
	CALE_ADD,
	CALE_SUB,
	CALE_MUL,
	CALE_DIV,
	CALE_POW,
	CALE_NEG,
	CALE_EXP,
	CALE_SIN,
	CALE_COS,
	CALE_ABS,
	CALE_SQRT,
	CALE_GAUS,
	CALE_CONST,
	CALE_VARIABLE,
};

struct CalExpr {

	CalExpr(CalExprType type, std::vector<CalExpr*> const &x) : fSubExprs(x), fType(type)
	{		
	}

	CalExpr(CalExprType type, double a) : fSubExprs({}), fType(type), fConstVal(a, 0)
	{
	}
	CalExpr(CalExprType type, char const *a) : fSubExprs({}), fType(type), fVarName(a)
	{
	}

	CalExprType fType;
	std::vector<CalExpr*> fSubExprs;

	CCom fConstVal;

	std::string fVarName;
	CCom Val(Cal *cal)
	{

		switch (fType) {
		case CALE_ADD:
			return fSubExprs.at(0)->Val(cal) + fSubExprs.at(1)->Val(cal);
			break;
		case CALE_SUB:
			return fSubExprs.at(0)->Val(cal) - fSubExprs.at(1)->Val(cal);
			break;
		case CALE_MUL:
			return fSubExprs.at(0)->Val(cal) * fSubExprs.at(1)->Val(cal);
			break;
		case CALE_DIV:
			return fSubExprs.at(0)->Val(cal) / fSubExprs.at(1)->Val(cal);
			break;
		case CALE_POW:
			return pow(fSubExprs.at(0)->Val(cal), fSubExprs.at(1)->Val(cal));
			break;
		case CALE_EXP:
			return exp(fSubExprs.at(0)->Val(cal));
			break;
		case CALE_SIN:
			return sin(fSubExprs.at(0)->Val(cal));
			break;
		case CALE_COS:
			return cos(fSubExprs.at(0)->Val(cal));
			break;
		case CALE_ABS:
			return abs(fSubExprs.at(0)->Val(cal));
			break;
		case CALE_NEG:
			return (-fSubExprs.at(0)->Val(cal));
			break;
		case CALE_SQRT:
			return sqrt(fSubExprs.at(0)->Val(cal));
			break;
		case CALE_GAUS:
		{
			CCom x = fSubExprs.at(0)->Val(cal);
			CCom mu = fSubExprs.at(1)->Val(cal);
			CCom sigma = fSubExprs.at(2)->Val(cal);
			CCom t = (x - mu) / sigma;
			return 1.0 / (sqrt(2 * 3.14159265358979323846) * sigma) *exp(-0.5*t*t);
			break;
		}
		case CALE_CONST:
		{
			return fConstVal;
			break;
		}
		case CALE_VARIABLE:
		{
			return cal->GetVarVal(fVarName);
			break;
		}
		default:
			throw std::runtime_error("not implemented!");
		}

	}
};

bool isEmpy(char c)
{
	return (c == ' ' || c == '\t' || c == '\n' || c == '\r');
}
void skipEmpty(char const *&s)
{
	for (; *s != 0 && isEmpy(*s); ++s) {
	}
}

bool isNum(char c)
{
	return (c >= '0' && c <= '9');
}

bool isAlpha(char c)
{
	return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c == '_');
}

bool startWith(char const *s, char const *s1)
{
	for (; *s && *s1 && *s == *s1; ++s, ++s1) {
	}
	return *s1 == '\0';
}
CalExpr *parseAddExpr(char const *&s);
CalExpr *parseFunc(char const *&s)
{

	if (startWith(s, "-") || startWith(s, "+")) {
		char c = *s;
		++s;
		skipEmpty(s);
		auto e =  parseFunc(s);
		if (c == '-') {
			return new CalExpr(CALE_NEG, {e});
		} else if (c == '+') {
			return e;
		}
	}

	if (isNum(*s)) {
		double a;
		int n;
		int c = sscanf(s, "%lf%n", &a, &n);
		if (c <= 0) throw std::runtime_error("falied to parse double");
		s += n;
		skipEmpty(s);
		return new CalExpr(CALE_CONST, a);
	}

	if (startWith(s, "(")) {
		++s;
		skipEmpty(s);
		auto e = parseAddExpr(s);
		if (*s == ')') ++s;
		else throw std::runtime_error("expect ')'");
		skipEmpty(s);
		return e;
	}

	std::string name;
	if (isAlpha(*s)) {
		for (; isAlpha(*s) || isNum(*s); ++s) {
			name += *s;
		}
		skipEmpty(s);


		if (*s == '(') {

			CalExprType type = CALE_UNK;
			if (name == "sin") {
				type = CALE_SIN;
			} else if (name == "cos") {
				type = CALE_COS;
			} else if (name == "exp") {
				type = CALE_EXP;
			} else if (name == "abs") {
				type = CALE_ABS;
			} else if (name == "sqrt") {
				type = CALE_SQRT;
			} else if (name == "gauss") {
				type = CALE_GAUS;
			} else {
				throw std::runtime_error("unkown function");
			}

			if (type == CALE_SIN || type == CALE_COS || type == CALE_EXP || type == CALE_ABS || type == CALE_SQRT) {
				++s;
				skipEmpty(s);

				skipEmpty(s);
				auto e1 = parseAddExpr(s);
				if (*s == ')') ++s;
				else throw std::runtime_error("expect ')'");
				skipEmpty(s);

				auto e = new CalExpr(type, { e1 });

				return e;
			} else if (type == CALE_GAUS) {
				++s;
				skipEmpty(s);

				auto e1 = parseAddExpr(s);
				skipEmpty(s);
				if (*s == ',') ++s;
				else throw std::runtime_error("expect ','");
				skipEmpty(s);

				auto e2 = parseAddExpr(s);
				skipEmpty(s);
				if (*s == ',') ++s;
				else throw std::runtime_error("expect ','");
				skipEmpty(s);

				auto e3 = parseAddExpr(s);
				if (*s == ')') ++s;
				else throw std::runtime_error("expect ')'");
				skipEmpty(s);

				auto e = new CalExpr(type, { e1, e2, e3 });

				return e;

			}

		} else {
			return new CalExpr(CALE_VARIABLE, name.c_str());
		}

	}

	throw std::runtime_error("not valid expression!");
}

CalExpr *parsePowExpr(char const *&s)
{
	auto e1 = parseFunc(s);

	if (*s == '^') {
		++s;
		skipEmpty(s);
		auto e2 = parsePowExpr(s);
		e1 = new CalExpr(CALE_POW, { e1, e2 });
	}

	return e1;
}


CalExpr *parseMulExpr(char const *&s)
{
	auto e1 = parsePowExpr(s);
	for (;;) {

		char c = *s;
		if (*s == '*' || *s == '/') {
			++s;
			skipEmpty(s);
			auto e2 = parsePowExpr(s);
			if (c == '*') {
				e1 = new CalExpr(CALE_MUL, { e1, e2 });
			} else if (c == '/') {
				e1 = new CalExpr(CALE_DIV, { e1, e2 });
			}
		} else {
			break;
		}

	}

	return e1;
}

CalExpr *parseAddExpr(char const *&s)
{
	auto e1 = parseMulExpr(s);
	for (;;) {

		char c = *s;
		if (*s == '+' || *s == '-') {
			++s;
			skipEmpty(s);
			auto e2 = parseMulExpr(s);
			if (c == '+') {
				e1 = new CalExpr(CALE_ADD, { e1, e2 });
			} else if (c == '-') {
				e1 = new CalExpr(CALE_SUB, { e1, e2 });
			}
		} else {
			break;
		}
	}

	return e1;
}

CalExpr *parseExpr(char const *s)
{
	skipEmpty(s);
	auto e =  parseAddExpr(s);
	if (*s != '\0') {
		throw std::runtime_error("expect end of file");
	}
	return e;
}

Cal::Cal(char const * x) : fStr(x)
{
	char const * s = fStr.c_str();
	fExpr = parseExpr(s);
	SetVarVal("Pi", 3.14159265358979323846);
	SetVarVal("I", CCom(0, 1));
}

void Cal::SetVarVal(std::string const & name, CCom const & v)
{
	fVarVals[name] = v;
}

CCom Cal::GetVarVal(std::string const & name)
{
	auto it = fVarVals.find(name);
	if (it == fVarVals.end()) {
		throw std::runtime_error("var not found!");
	}
	return fVarVals[name];
}

CCom Cal::Val()
{
	return fExpr->Val(this);
}
