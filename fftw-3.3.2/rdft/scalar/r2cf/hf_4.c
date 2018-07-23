/*
 * Copyright (c) 2003, 2007-11 Matteo Frigo
 * Copyright (c) 2003, 2007-11 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Sat Apr 28 11:02:42 EDT 2012 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc.native -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 4 -dit -name hf_4 -include hf.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 16 additions, 6 multiplications, 6 fused multiply/add),
 * 31 stack variables, 0 constants, and 16 memory accesses
 */
#include "hf.h"

static void hf_4(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 6); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 6, MAKE_VOLATILE_STRIDE(rs)) {
	       E To, Te, Tm, T8, Ty, Tw, Tq, Tk;
	       {
		    E T1, Tv, Tu, T7, Tg, Tj, Tf, Ti, Tp, Th;
		    T1 = cr[0];
		    Tv = ci[0];
		    {
			 E T3, T6, T2, T5;
			 T3 = cr[WS(rs, 2)];
			 T6 = ci[WS(rs, 2)];
			 T2 = W[2];
			 T5 = W[3];
			 {
			      E Ta, Td, Tc, Tn, Tb, Tt, T4, T9;
			      Ta = cr[WS(rs, 1)];
			      Td = ci[WS(rs, 1)];
			      Tt = T2 * T6;
			      T4 = T2 * T3;
			      T9 = W[0];
			      Tc = W[1];
			      Tu = FNMS(T5, T3, Tt);
			      T7 = FMA(T5, T6, T4);
			      Tn = T9 * Td;
			      Tb = T9 * Ta;
			      Tg = cr[WS(rs, 3)];
			      Tj = ci[WS(rs, 3)];
			      To = FNMS(Tc, Ta, Tn);
			      Te = FMA(Tc, Td, Tb);
			      Tf = W[4];
			      Ti = W[5];
			 }
		    }
		    Tm = T1 - T7;
		    T8 = T1 + T7;
		    Tp = Tf * Tj;
		    Th = Tf * Tg;
		    Ty = Tv - Tu;
		    Tw = Tu + Tv;
		    Tq = FNMS(Ti, Tg, Tp);
		    Tk = FMA(Ti, Tj, Th);
	       }
	       {
		    E Tr, Ts, Tl, Tx;
		    Tr = To - Tq;
		    Ts = To + Tq;
		    Tl = Te + Tk;
		    Tx = Tk - Te;
		    ci[WS(rs, 3)] = Ts + Tw;
		    cr[WS(rs, 2)] = Ts - Tw;
		    cr[WS(rs, 1)] = Tm + Tr;
		    ci[0] = Tm - Tr;
		    ci[WS(rs, 2)] = Tx + Ty;
		    cr[WS(rs, 3)] = Tx - Ty;
		    cr[0] = T8 + Tl;
		    ci[WS(rs, 1)] = T8 - Tl;
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 4},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 4, "hf_4", twinstr, &GENUS, {16, 6, 6, 0} };

void X(codelet_hf_4) (planner *p) {
     X(khc2hc_register) (p, hf_4, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc.native -compact -variables 4 -pipeline-latency 4 -n 4 -dit -name hf_4 -include hf.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 16 additions, 6 multiplications, 6 fused multiply/add),
 * 13 stack variables, 0 constants, and 16 memory accesses
 */
#include "hf.h"

static void hf_4(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 6); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 6, MAKE_VOLATILE_STRIDE(rs)) {
	       E T1, Tp, T6, To, Tc, Tk, Th, Tl;
	       T1 = cr[0];
	       Tp = ci[0];
	       {
		    E T3, T5, T2, T4;
		    T3 = cr[WS(rs, 2)];
		    T5 = ci[WS(rs, 2)];
		    T2 = W[2];
		    T4 = W[3];
		    T6 = FMA(T2, T3, T4 * T5);
		    To = FNMS(T4, T3, T2 * T5);
	       }
	       {
		    E T9, Tb, T8, Ta;
		    T9 = cr[WS(rs, 1)];
		    Tb = ci[WS(rs, 1)];
		    T8 = W[0];
		    Ta = W[1];
		    Tc = FMA(T8, T9, Ta * Tb);
		    Tk = FNMS(Ta, T9, T8 * Tb);
	       }
	       {
		    E Te, Tg, Td, Tf;
		    Te = cr[WS(rs, 3)];
		    Tg = ci[WS(rs, 3)];
		    Td = W[4];
		    Tf = W[5];
		    Th = FMA(Td, Te, Tf * Tg);
		    Tl = FNMS(Tf, Te, Td * Tg);
	       }
	       {
		    E T7, Ti, Tj, Tm;
		    T7 = T1 + T6;
		    Ti = Tc + Th;
		    ci[WS(rs, 1)] = T7 - Ti;
		    cr[0] = T7 + Ti;
		    Tj = T1 - T6;
		    Tm = Tk - Tl;
		    ci[0] = Tj - Tm;
		    cr[WS(rs, 1)] = Tj + Tm;
	       }
	       {
		    E Tn, Tq, Tr, Ts;
		    Tn = Tk + Tl;
		    Tq = To + Tp;
		    cr[WS(rs, 2)] = Tn - Tq;
		    ci[WS(rs, 3)] = Tn + Tq;
		    Tr = Th - Tc;
		    Ts = Tp - To;
		    cr[WS(rs, 3)] = Tr - Ts;
		    ci[WS(rs, 2)] = Tr + Ts;
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 4},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 4, "hf_4", twinstr, &GENUS, {16, 6, 6, 0} };

void X(codelet_hf_4) (planner *p) {
     X(khc2hc_register) (p, hf_4, &desc);
}
#endif				/* HAVE_FMA */
