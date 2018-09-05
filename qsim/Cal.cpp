#define _CRT_SECURE_NO_WARNINGS
#include "Cal.h"
#include <vector>
#include <complex>
#include <cmath>
#include <math.h>

enum CalExprType {
	CALE_UNK,

	CALE_CONST,
	CALE_VARIABLE,

	CALE_QUEST,
	CALE_OROR,
	CALE_ANDAND,

	CALE_GT,
	CALE_LT,
	CALE_GE,
	CALE_LE,
	CALE_EQ,

	CALE_ADD,
	CALE_SUB,
	CALE_MUL,
	CALE_DIV,
	CALE_NEG,

	CALE_SIN,
	CALE_COS,
	CALE_TAN,
	CALE_EXP,
	CALE_TANH,
	CALE_SINH,

	CALE_ASIN,
	CALE_ACOS,
	CALE_ATAN,
	CALE_LOG,
	CALE_ATANH,

	CALE_SIGN,
	CALE_ABS,
	CALE_SQRT,
	CALE_POW,
	CALE_GAUS,

};

struct StackUsage {
	int usage;
	int maxu;
	StackUsage()
	{
		usage = 0;
		maxu = 0;
	}
	StackUsage &operator+=(int n)
	{
		usage += n;
		if (usage >= maxu) maxu = usage;
		return *this;
	}
};

struct CalExpr {

	CalExpr(CalExprType type, std::vector<CalExpr*> const &x) : fSubExprs(x), fType(type), fVarValue(nullptr)
	{		
	}

	CalExpr(CalExprType type, double a) : fSubExprs({}), fType(type), fConstVal(a, 0), fVarValue(nullptr)
	{
	}
	CalExpr(CalExprType type, char const *a) : fSubExprs({}), fType(type), fVarName(a), fVarValue(nullptr)
	{
	}

	~CalExpr()
	{
		for (auto x : fSubExprs) {
			delete x;
		}
	}
	CalExprType fType;
	std::vector<CalExpr*> fSubExprs;

	CCom fConstVal;

	std::string fVarName;
	CCom const *fVarValue;

	void Gen(Cal *cal, StackUsage &nstack);

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
CalExpr *parseQuestionExpr(char const *&s);
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
		if (c <= 0) throw std::runtime_error(std::string("falied to parse double for: ") + s);
		s += n;
		skipEmpty(s);
		return new CalExpr(CALE_CONST, a);
	}

	if (startWith(s, "(")) {
		++s;
		skipEmpty(s);
		auto e = parseQuestionExpr(s);
		if (*s == ')') ++s;
		else throw std::runtime_error(std::string("expect ')' before: ") + s);
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
			} else if (name == "tan") {
				type = CALE_TAN;
			} else if (name == "asin") {
				type = CALE_ASIN;
			} else if (name == "acos") {
				type = CALE_ACOS;
			} else if (name == "atan") {
				type = CALE_ATAN;
			} else if (name == "sign") {
				type = CALE_SIGN;
			} else if (name == "exp") {
				type = CALE_EXP;
			} else if (name == "tanh") {
				type = CALE_TANH;
			} else if (name == "sinh") {
				type = CALE_SINH;
			} else if (name == "log") {
				type = CALE_LOG;
			} else if (name == "atanh") {
				type = CALE_ATANH;
			} else if (name == "abs") {
				type = CALE_ABS;
			} else if (name == "sqrt") {
				type = CALE_SQRT;
			} else if (name == "gauss") {
				type = CALE_GAUS;
			} else {
				throw std::runtime_error("unkown function");
			}

			if (type == CALE_SIN || type == CALE_COS || type == CALE_TAN ||
				type == CALE_ASIN || type == CALE_ACOS || type == CALE_ATAN ||
				type == CALE_EXP || type == CALE_TANH || type == CALE_SINH ||
				type == CALE_LOG || type == CALE_ATANH ||
				type == CALE_ABS || type == CALE_SQRT) {
				++s;
				skipEmpty(s);

				skipEmpty(s);
				auto e1 = parseQuestionExpr(s);
				if (*s == ')') ++s;
				else throw std::runtime_error(std::string("expect ')' before: ") + s);
				skipEmpty(s);

				auto e = new CalExpr(type, { e1 });

				return e;
			} else if (type == CALE_GAUS) {
				++s;
				skipEmpty(s);

				auto e1 = parseQuestionExpr(s);
				skipEmpty(s);
				if (*s == ',') ++s;
				else throw std::runtime_error(std::string("expect ',' before: ") + s);
				skipEmpty(s);

				auto e2 = parseQuestionExpr(s);
				skipEmpty(s);
				if (*s == ',') ++s;
				else throw std::runtime_error(std::string("expect ',' before: ") + s);
				skipEmpty(s);

				auto e3 = parseQuestionExpr(s);
				if (*s == ')') ++s;
				else throw std::runtime_error(std::string("expect ')' before: ") + s);
				skipEmpty(s);

				auto e = new CalExpr(type, { e1, e2, e3 });

				return e;

			}

		} else {
			return new CalExpr(CALE_VARIABLE, name.c_str());
		}

	}

	if(*s == '\0') 	throw std::runtime_error(std::string("unexpected: ") + "end of file");
	else throw std::runtime_error(std::string("unexpected: ") + s);
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


CalExpr *parseCmpExpr(char const *&s)
{
	auto e1 = parseAddExpr(s);
	for (;;) {

		char const *ss = s;
		if (startWith(s, "<") || startWith(s, "<=") ||
			startWith(s, ">") || startWith(s, ">=") ||
			startWith(s, "==")) {
			
			if (startWith(ss, "<=") || startWith(ss, ">=") || startWith(ss, ">=")) s += 2;
			else s += 1;

			skipEmpty(s);
			auto e2 = parseAddExpr(s);
			if (startWith(ss, "<=")) {
				e1 = new CalExpr(CALE_LE, { e1, e2 });
			} else if (startWith(ss, ">=")) {
				e1 = new CalExpr(CALE_GE, { e1, e2 });
			} else if (startWith(ss, "==")) {
				e1 = new CalExpr(CALE_EQ, { e1, e2 });
			} else if (startWith(ss, "<")) {
				e1 = new CalExpr(CALE_LT, { e1, e2 });
			} else if (startWith(ss, ">")) {
				e1 = new CalExpr(CALE_GT, { e1, e2 });
			}
		} else {
			break;
		}
	}

	return e1;

}

CalExpr *parseAndAndExpr(char const *&s)
{
	auto e1 = parseCmpExpr(s);
	for (;;) {

		char c = *s;
		if (startWith(s, "&&")) {
			s += 2;
			skipEmpty(s);
			auto e2 = parseCmpExpr(s);
			e1 = new CalExpr(CALE_ANDAND, { e1, e2 });
		} else {
			break;
		}
	}

	return e1;

}

CalExpr *parseOrOrExpr(char const *&s)
{
	auto e1 = parseAndAndExpr(s);
	for (;;) {

		if (startWith(s, "||")) {
			s += 2;
			skipEmpty(s);
			auto e2 = parseAndAndExpr(s);
			e1 = new CalExpr(CALE_OROR, { e1, e2 });
		} else {
			break;
		}
	}

	return e1;

}

CalExpr *parseQuestionExpr(char const *&s)
{
	auto e1 = parseOrOrExpr(s);
	for (;;) {

		if (startWith(s, "?")) {
			s += 1;
			skipEmpty(s);

			auto e2 = parseQuestionExpr(s);

			if (*s != ':') throw std::runtime_error(std::string("expect ':' before: ") + s);
			s += 1;
			skipEmpty(s);

			auto e3 = parseOrOrExpr(s);
			skipEmpty(s);

			e1 = new CalExpr(CALE_QUEST, { e1, e2, e3 });
		} else {
			break;
		}
	}

	return e1;
}

CalExpr *parseExpr(char const *s)
{
	skipEmpty(s);
	auto e = parseQuestionExpr(s);
	if (*s != '\0') {
		throw std::runtime_error(std::string("expect end of input, but find: ") + s);
	}
	return e;
}

Cal::Cal(char const * x) : fStr(x)
{
	char const * s = fStr.c_str();
	fExpr.reset(parseExpr(s));
	SetVarVal("Pi", 3.14159265358979323846);
	SetVarVal("I", CCom(0, 1));
}

Cal::~Cal()
{
}

void Cal::SetVarVal(std::string const & name, CCom const & v)
{
	if (fVarVals.find(name) == fVarVals.end()) {
		// not alloc ever ?
		fPseudoCode.clear();
	}
	fVarVals[name] = v;
}

CCom &Cal::GetVarVal(std::string const & name)
{
	auto it = fVarVals.find(name);
	if (it == fVarVals.end()) {
		throw std::runtime_error(std::string("var not defined: ") + name);
	}
	return it->second;
}


















void CalExpr::Gen(Cal *cal, StackUsage &nstack)
{
	for (size_t i = 0; i < fSubExprs.size(); ++i) {
		fSubExprs.at(i)->Gen(cal, nstack);
	}

	cal->fPseudoCode.push_back(fType);
	switch (fType) {
	case CALE_QUEST:
		nstack += -2;
		break;
	case CALE_ANDAND:
	case CALE_OROR:
	case CALE_GT:
	case CALE_GE:
	case CALE_LT:
	case CALE_LE:
	case CALE_EQ:
	case CALE_ADD:
	case CALE_SUB:
	case CALE_MUL:
	case CALE_DIV:
	case CALE_POW:
		nstack += -1;
		break;
	case CALE_EXP:
	case CALE_TANH:
	case CALE_SINH:
	case CALE_LOG:
	case CALE_ATANH:
	case CALE_SIN:
	case CALE_COS:
	case CALE_TAN:
	case CALE_ASIN:
	case CALE_ACOS:
	case CALE_ATAN:
	case CALE_SIGN:
	case CALE_ABS:
	case CALE_NEG:
	case CALE_SQRT:
		break;
	case CALE_GAUS:
	{
		nstack += -2;
		break;
	}
	case CALE_CONST:
	{
		cal->fConstants.push_back(fConstVal);
		nstack += 1;
		break;
	}
	case CALE_VARIABLE:
	{
		if (!fVarValue) fVarValue = &cal->GetVarVal(fVarName);
		cal->fVariables.push_back(fVarValue);
		nstack += 1;
		break;
	}
	default:
		throw std::runtime_error("not implemented!");
	}

}

CCom Cal::Val()
{
	if (fPseudoCode.size() == 0) GenPseudoCode();
	return RunPseudoCode();
}

void Cal::GenPseudoCode()
{
	fPseudoCode.clear();
	fStack.clear();
	fConstants.clear();
	fVariables.clear();

	StackUsage su;
	fExpr->Gen(this, su);
	fPseudoCode.push_back(0);
	fStack.resize(su.maxu);
}

CCom Cal::RunPseudoCode()
{
	int *ip = fPseudoCode.data();
	CCom *stack = fStack.data() - 1;
	CCom *constant = fConstants.data();
	CCom const **variables = fVariables.data();

	for (; *ip; ++ip) {
		CalExprType code = (CalExprType)*ip;

		switch (code) {
		case CALE_CONST:
		{
			*(++stack) = *(constant++);
			break;
		}
		case CALE_VARIABLE:
		{
			*(++stack) = **(variables++);
			break;
		}

		case CALE_QUEST:
			stack -= 2;
			stack[0] = stack[0] != 0.0 ? stack[1] : stack[2];
			break;
		case CALE_ANDAND:
			stack -= 1;
			stack[0] = (stack[0] != 0.0 && stack[1] != 0.0) ? 1 : 0;
			break;
		case CALE_OROR:
			--stack;
			stack[0] = (stack[0] != 0.0 || stack[1] != 0.0) ? 1 : 0;
			break;
		case CALE_GT:
			--stack;
			stack[0] = stack[0].real() > stack[1].real();
			break;
		case CALE_GE:
			--stack;
			stack[0] = stack[0].real() >= stack[1].real();
			break;
		case CALE_LT:
			--stack;
			stack[0] = stack[0].real() < stack[1].real();
			break;
		case CALE_LE:
			--stack;
			stack[0] = stack[0].real() <= stack[1].real();
			break;
		case CALE_EQ:
			--stack;
			stack[0] = stack[0] == stack[1];
			break;
		case CALE_ADD:
			--stack;
			stack[0] = stack[0] + stack[1];
			break;
		case CALE_SUB:
			--stack;
			stack[0] = stack[0] - stack[1];
			break;
		case CALE_MUL:
			--stack;
			stack[0] = stack[0] * stack[1];
			break;
		case CALE_DIV:
			--stack;
			stack[0] = stack[0] / stack[1];
			break;
		case CALE_POW:
			stack -= 1;
			stack[0] = pow(stack[0], stack[1]);
			break;
		case CALE_EXP:
			stack[0] = exp(stack[0]);
			break;
		case CALE_TANH:
			stack[0] = tanh(stack[0]);
			break;
		case CALE_SINH:
			stack[0] = sinh(stack[0]);
			break;
		case CALE_LOG:
			stack[0] = log(stack[0]);
			break;
		case CALE_ATANH:
			stack[0] = atan(stack[0]);
			break;
		case CALE_SIN:
			stack[0] = sin(stack[0]);
			break;
		case CALE_COS:
			stack[0] = cos(stack[0]);
			break;
		case CALE_TAN:
			stack[0] = tan(stack[0]);
			break;
		case CALE_ASIN:
			stack[0] = asin(stack[0]);
			break;
		case CALE_ACOS:
			stack[0] = acos(stack[0]);
			break;
		case CALE_ATAN:
			stack[0] = atan(stack[0]);
			break;
		case CALE_SIGN:
			stack[0] = stack[0].real() > 0 ? 1 : -1;
			break;
		case CALE_ABS:
			stack[0] = abs(stack[0]);
			break;
		case CALE_NEG:
			stack[0] = -stack[0];
			break;
		case CALE_SQRT:
			stack[0] = sqrt(stack[0]);
			break;
		case CALE_GAUS:
		{
			stack -= 2;
			CCom x = stack[0];
			CCom mu = stack[1];
			CCom sigma = stack[2];
			CCom t = (x - mu) / sigma;
			stack[0] = 1.0 / (sqrt(2 * 3.14159265358979323846) * sigma) *exp(-0.5*t*t);
			break;
		}
		default:
			throw std::runtime_error("not implemented!");
		}

	}
	return *stack;
}
