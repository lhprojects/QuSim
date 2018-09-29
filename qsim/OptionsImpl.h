#pragma once

#include <map>
#include <string>
#include "QuSim.h"

enum class OptionType {
	Unkown,
	Bool,
	Int,
	Real,
	Complex,
	String,
};

struct OptionValue {
	OptionType fType;
	std::string fString;
	Complex fComplex;
	Real fReal;
	Int fInt;
	bool fBool;
};

struct OptionsImpl {

	void SetBool(char const *k, bool v);
	void SetInt(char const *k, Int v);
	void SetReal(char const *k, Real v);
	void SetComplex(char const *k, Complex v);
	void SetString(char const *k, char const *v);

	// key must exist and has type of bool, else an exception will be thrown
	bool GetBool(char const *k) const;
	// if key doesn't exist, v will be returned;
	// else if key exist however, however, without a type of bool, an exception will be thrown
	bool GetBool(char const *k, bool v) const;
	Int GetInt(char const *k) const;
	Int GetInt(char const *k, Int v) const;
	Real GetReal(char const *k) const;
	Real GetReal(char const *k, Real v) const;
	Complex GetComplex(char const *k) const;
	Complex GetComplex(char const *k, Complex v) const;
	std::string GetString(char const *k) const;
	std::string GetString(char const *k, std::string const &v) const;

	// if key doesn't exist, v will be returned;
	// else if key exist however, however, without a type of bool, an exception will be thrown
	bool Get(char const *k, bool &) const;
	bool Get(char const *k, std::string &) const;
private:
	std::map<std::string, OptionValue> fOpts;
};
