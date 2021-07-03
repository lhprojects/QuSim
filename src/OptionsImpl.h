#pragma once

#include <map>
#include <string>
#include "QuSim.h"
#include <variant>

std::string to_string(Complex c);
std::string to_string(std::string s);

enum class OptionType {
    Bool,
    Int,
    Real,
    Complex,
    String,
};

constexpr OptionType  GetOptionType(bool) { return OptionType::Bool; }
constexpr OptionType  GetOptionType(Int) { return OptionType::Int; }
constexpr OptionType  GetOptionType(Real) { return OptionType::Real; }
constexpr OptionType  GetOptionType(Complex) { return OptionType::Complex; }
constexpr OptionType  GetOptionType(std::string const&) { return OptionType::String; }

std::string GetOptionTypeName(OptionType type);

struct OptionValue {
    OptionType fType;
    std::variant<bool, Int, Real, Complex, std::string> fValue;
};

struct OptionsImpl {

    // Set `k` with new value `v`.
    // If the type of `v` is not the same as the previous set (If you have set one),
    // it will emit a warning from stderr, but no exception will be thrown
    void SetBool(std::string k, bool v);
    void SetInt(std::string k, Int v);
    void SetReal(std::string k, Real v);
    void SetComplex(std::string k, Complex v);
    void SetString(std::string k, std::string v);

    template<class String, class V>
    void Set(String &&k, V &&v)
    {
        using PureType = std::remove_cv_t< std::remove_reference_t<V> >;
        this->Set_<PureType>(std::forward<String>(k), std::forward<V>(v));
    }

    // key must exist and has type of bool, else an exception will be thrown
    bool GetBool(std::string const &k) const;
    // if key doesn't exist, default_ will be returned;
    // else if the key exists, however without a type of bool, an exception will be thrown
    bool GetBool(std::string const& k, bool default_) const;
    Int GetInt(std::string const& k) const;
    Int GetInt(std::string const& k, Int default_) const;
    Real GetReal(std::string const& k) const;
    Real GetReal(std::string const& k, Real default_) const;
    Complex GetComplex(std::string const& k) const;
    Complex GetComplex(std::string const& k, Complex default_) const;
    std::string GetString(std::string const& k) const;
    std::string GetString(std::string const& k, std::string const &v) const;
    std::string GetString(std::string const& k, std::string &&v) const;

    template<class V>
    std::remove_cv_t< std::remove_reference_t<V> > Get(std::string const& k) const;


    bool Contains(std::string const& k) const;

    template<class V>
    std::remove_cv_t< std::remove_reference_t<V> > GetDefaut(std::string const& k, V&& v) const
    {
        using PureType = std::remove_cv_t< std::remove_reference_t<V> >;
        return this->GetDefault_<PureType>(k, std::forward<V>(v));
    }

    // return true, if we have the key and the type is right
    // note `int a; Get("", a)` may have linking error, if the `int` is not the same the `Int`
    template<class V>
    bool Get(std::string const& k, V& v) const;
private:

    template<class PureType>
    PureType GetDefault_(std::string const& k, PureType v) const;

    template<class V>
    void Set_(std::string k, V v);

    std::map<std::string, OptionValue> fOpts;
};
