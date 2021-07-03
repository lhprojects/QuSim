#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>

template<class PsiVector>
inline void dump_comp(PsiVector const &v, char const *fn)
{
    FILE *f = fopen(fn, " w");
    for (int i = 0; i < (int)v.size(); ++i) {
        auto x = v[i];
        fprintf(f, "(% .20lf, % .20lf)\n", x.real(), x.imag());
    }
    fclose(f);
};

template<class PsiVector>
inline void dump_real(PsiVector const &v, char const *fn)
{
    FILE *f = fopen(fn, " w");
    for (int i = 0; i < (int)v.size(); ++i) {
        fprintf(f, "% .20lf\n", (double)v(i));
    }
    fclose(f);
};

template<class Matrix>
inline void dump_matrix_real(Matrix const &v, char const *fn)
{
    FILE *f = fopen(fn, " w");
    for (int i = 0; i < (int)v.rows(); ++i) {
        for (int j = 0; j < (int)v.cols(); ++j) {
            fprintf(f, "%+.20lf ", (double)v(i, j));
        }
        fprintf(f, "\n");
    }
    fclose(f);
};

template<class Matrix>
inline void dump_matrix_comp(Matrix const &v, char const *fn)
{
    FILE *f = fopen(fn, " w");
    for (int i = 0; i < (int)v.rows(); ++i) {
        for (int j = 0; j < (int)v.cols(); ++j) {
            fprintf(f, "(% .20lf, % .20lf)", v(i, j).real(), v(i, j).imag());
        }
        fprintf(f, "\n");
    }
    fclose(f);
};
