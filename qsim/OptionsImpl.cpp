#include "OptionsImpl.h"


void OptionsImpl::SetBool(char const *k, bool v)
{
	auto &opt = fOpts[k];
	if (opt.fType != OptionType::Unkown && opt.fType != OptionType::Bool) {
		throw std::runtime_error("the key has set and with different type");
	}
	opt.fType = OptionType::Bool;
	opt.fBool = v;
}

void OptionsImpl::SetInt(char const *k, Int v)
{
	auto &opt = fOpts[k];
	if (opt.fType != OptionType::Unkown && opt.fType != OptionType::Int) {
		throw std::runtime_error("the key has set and with different type");
	}
	opt.fType = OptionType::Int;
	opt.fInt = v;
}

void OptionsImpl::SetReal(char const *k, Real v)
{
	auto &opt = fOpts[k];
	if (opt.fType != OptionType::Unkown && opt.fType != OptionType::Real) {
		throw std::runtime_error("the key has set and with different type");
	}
	opt.fType = OptionType::Real;
	opt.fReal = v;
}

void OptionsImpl::SetComplex(char const *k, Complex v)
{
	auto &opt = fOpts[k];
	if (opt.fType != OptionType::Unkown && opt.fType != OptionType::Complex) {
		throw std::runtime_error("the key has set and with different type");
	}
	opt.fType = OptionType::Complex;
	opt.fComplex = v;
}

void OptionsImpl::SetString(char const *k, char const *v)
{
	auto &opt = fOpts[k];
	if (opt.fType != OptionType::Unkown && opt.fType != OptionType::String) {
		throw std::runtime_error("the key has set and with different type");
	}
	opt.fType = OptionType::String;
	opt.fString = v;
}

bool OptionsImpl::GetBool(char const *k) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end()) {
		std::string info = std::string("key not found: ") + k;
		throw std::runtime_error(info.c_str());
	} else if(it->second.fType != OptionType::Bool) {
		throw std::runtime_error("value type error: ");
	}
	return it->second.fBool;
}

bool OptionsImpl::GetBool(char const *k, bool v) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end()) {
		return v;
	} else if (it->second.fType != OptionType::Bool) {
		throw std::runtime_error("value type wrong");
	}
	return it->second.fBool;
}

Int OptionsImpl::GetInt(char const *k) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end() || it->second.fType != OptionType::Int) {
		throw std::runtime_error("value not set");
	}
	return it->second.fInt;
}

Int OptionsImpl::GetInt(char const *k, Int v) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end()) {
		return v;
	} else if (it->second.fType != OptionType::Int) {
		throw std::runtime_error("value type wrong");
	}
	return it->second.fInt;
}

Real OptionsImpl::GetReal(char const *k) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end() || it->second.fType != OptionType::Real) {
		throw std::runtime_error("value not set");
	}
	return it->second.fReal;
}

Real OptionsImpl::GetReal(char const *k, Real v) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end()) {
		return v;
	} else if (it->second.fType != OptionType::Real) {
		throw std::runtime_error("value type wrong");
	}
	return it->second.fReal;
}


Complex OptionsImpl::GetComplex(char const *k) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end() || it->second.fType != OptionType::Complex) {
		throw std::runtime_error("value not set");
	}
	return it->second.fComplex;
}

Complex OptionsImpl::GetComplex(char const *k, Complex v) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end()) {
		return v;
	} else if (it->second.fType != OptionType::Complex) {
		throw std::runtime_error("value type wrong");
	}
	return it->second.fComplex;
}


std::string OptionsImpl::GetString(char const *k) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end() || it->second.fType != OptionType::String) {
		throw std::runtime_error("value not set");
	}
	return it->second.fString;
}

std::string OptionsImpl::GetString(char const *k, std::string const &v) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end()) {
		return v;
	} else if (it->second.fType != OptionType::String) {
		throw std::runtime_error("value type wrong");
	}
	return it->second.fString;
}

bool OptionsImpl::Get(char const * k, bool &v) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end()) {
		return false;
	} else if (it->second.fType != OptionType::Bool) {
		throw std::runtime_error("value type wrong");
	}
	v = it->second.fBool;
	return true;
}

bool OptionsImpl::Get(char const * k, std::string &v) const
{
	auto it = fOpts.find(k);
	if (it == fOpts.end()) {
		return false;
	} else if (it->second.fType != OptionType::String) {
		throw std::runtime_error("value type wrong");
	}
	v = it->second.fString;
	return true;
}


