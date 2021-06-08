#include "OptionsImpl.h"
#include <stdio.h>

static const std::string Unkown_ = "Unkown";
static const std::string Bool_ = "Bool";
static const std::string Int_ = "Int";
static const std::string Real_ = "Real";
static const std::string Complex_ = "Complex";
static const std::string String_ = "String";

std::string to_string(Complex c)
{
    return "(" + std::to_string(c.real()) + "," + std::to_string(c.imag()) + ")";
}
std::string to_string(std::string s)
{
    return std::move(s);
}

std::string GetOptionTypeName(OptionType type)
{
    switch (type) {
    case OptionType::Bool:
        return Bool_;
    case OptionType::Int:
        return Int_;
    case OptionType::Real:
        return Real_;
    case OptionType::Complex:
        return Complex_;
    case OptionType::String:
        return String_;
    default:
        return "";
    }
}

template<class V>
void OptionsImpl::Set_(std::string k, V v)
{
    auto type = GetOptionType(v);

    auto r = fOpts.insert(std::make_pair(std::move(k), OptionValue()) );
    auto iter = r.first;
    auto& k_ = iter->first;
    auto& opt = iter->second;
    if (!r.second && opt.fType != type) {
        fprintf(stderr, "key `%s` set with type `%s` not the same the previous set `%s`",
            k_.c_str(), GetOptionTypeName(type).c_str(), GetOptionTypeName(opt.fType).c_str());
    }
    opt.fType = type;
    opt.fValue = std::move(v);
}

void OptionsImpl::SetBool(std::string k, bool v)
{
    this->Set(std::move(k), v);
}

void OptionsImpl::SetInt(std::string k, Int v)
{
    this->Set(std::move(k), v);
}

void OptionsImpl::SetReal(std::string k, Real v)
{
    this->Set(std::move(k), v);
}

void OptionsImpl::SetComplex(std::string k, Complex v)
{
    this->Set(std::move(k), v);
}

void OptionsImpl::SetString(std::string k, std::string v)
{
    this->Set(std::move(k), std::move(v));
}


bool OptionsImpl::GetBool(std::string const &k) const
{
    return Get<bool>(k);
}

bool OptionsImpl::GetBool(std::string const& k, bool v) const
{
    return this->Get(k, v);
}

Int OptionsImpl::GetInt(std::string const& k) const
{
    return this->Get<Int>(k);
}

Int OptionsImpl::GetInt(std::string const& k, Int v) const
{
    return this->GetDefaut(k, v);
}

Real OptionsImpl::GetReal(std::string const& k) const
{
    return Get<Real>(k);
}

Real OptionsImpl::GetReal(std::string const &k, Real v) const
{
    return this->GetDefaut(k, v);
}


Complex OptionsImpl::GetComplex(std::string const &k) const
{
    return Get<Complex>(k);
}

Complex OptionsImpl::GetComplex(std::string const &k, Complex v) const
{
    return this->GetDefaut(k, v);
}


std::string OptionsImpl::GetString(std::string const &k) const
{
    return Get<std::string>(k);
}

std::string OptionsImpl::GetString(std::string const &k, std::string const &v) const
{
    return this->GetDefaut(k, v);
}

std::string OptionsImpl::GetString(std::string const& k, std::string && v) const
{
    return GetDefaut<std::string>(v, std::move(v));
}

bool OptionsImpl::Contains(std::string const& k) const
{
    auto it = fOpts.find(k);
    return it != fOpts.find(k);
}

template<class V>
std::remove_cv_t< std::remove_reference_t<V> > OptionsImpl::Get(std::string const& k) const
{
    using Return = std::remove_const_t< std::remove_cv_t<V> >;
    auto it = fOpts.find(k);
    if (it == fOpts.end()) {
        throw std::runtime_error("key not set" + k);
    } else if (it->second.fType != GetOptionType(V())) {
        throw std::runtime_error("value type wrong");
    }
    return std::get<Return>(it->second.fValue);
}

template<class V>
bool OptionsImpl::Get(std::string const& k, V& v) const
{
    auto type = GetOptionType(v);
    auto it = fOpts.find(k);
    if (it == fOpts.end()) {
        return false;
    } else if (it->second.fType != GetOptionType(V())) {
        return false;
    }
    v = std::get<std::remove_cv_t<V> >(it->second.fValue);
    return true;
}

template bool OptionsImpl::Get<bool>(std::string const& k, bool& v) const;
template bool OptionsImpl::Get<Int>(std::string const& k, Int& v) const;
template bool OptionsImpl::Get<Real>(std::string const& k, Real& v) const;
template bool OptionsImpl::Get<Complex>(std::string const& k, Complex& v) const;
template bool OptionsImpl::Get<std::string>(std::string const& k, std::string& v) const;

template<class V>
V OptionsImpl::GetDefault_(std::string const& k, V v) const
{
    using Return = std::remove_const_t< std::remove_cv_t<V> >;
    auto it = fOpts.find(k);
    if (it == fOpts.end()) {
        return std::move(v);
    } else if (it->second.fType != GetOptionType(V())) {
        throw std::runtime_error("value type wrong");
    }
    return std::get<Return>(it->second.fValue);
}
