#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>

typedef double Float;
Float Half(0.5);
Float Zero(0);
Float One(1);
Float Two(2);
Float Three(3);
Float Four(3);

// exp(aV) exp(T) exp((1-a)V)
void solve2()
{
	Float L[3];
	Float M[3];
	Float R[3];

	Float a(0.5);
	L[0] = 1;
	L[1] = a;
	L[2] = Half * a*a;
	M[0] = 1;
	M[1] = 1;
	M[2] = Half;
	R[0] = 1;
	R[1] = One - a;
	R[2] = Half * (One - a)*(One - a);

	Float COne = L[0]*M[0]*R[0];
	Float CV = L[1] + R[1];
	Float CT = M[1];
	Float CVV = L[2] * M[0] * R[0] + L[1] * M[0] * R[1] + L[0] * M[0] * R[2];
	Float CVT = L[1] * M[1] * R[0];
	Float CTV = L[0] * M[1] * R[1];
	Float CTT = L[0] * M[2] * R[0];

	std::cout << "COne " << COne << std::endl;
	std::cout << "CV " << CV << std::endl;
	std::cout << "CT " << CT << std::endl;
	std::cout << "CVV " << CVV << std::endl;
	std::cout << "CVT " << CVT << std::endl;
	std::cout << "CTV " << CTV << std::endl;
	std::cout << "CTT " << CTT << std::endl;
}

// exp(0.5 a V) exp(0.5T) exp((1-a)V) exp(0.5T) exp(0.5 aV)
void solve3()
{
	Float CC[5][5];

	Float a(One / Three);
	CC[0][0] = 1;
	CC[0][1] = 0.5 * a;
	CC[0][2] = Half * (0.5*a)*(0.5*a);
	CC[0][3] = Half / Three * (0.5*a)*(0.5*a)*(0.5*a);
	CC[0][4] = Half / Three / Four * (0.5*a)*(0.5*a)*(0.5*a)*(0.5*a);

	CC[1][0] = 1;
	CC[1][1] = 0.5;
	CC[1][2] = Half * Float(0.5 * 0.5);
	CC[1][3] = (Half / Three) * Float(0.5 * 0.5 * 0.5);
	CC[1][4] = (Half / Three / Four) * Float(0.5*0.5*0.5*0.5);

	CC[2][0] = 1;
	CC[2][1] = 1 - a;
	CC[2][2] = Half * (1 - a)* (1 - a);
	CC[2][3] = Half / Three * (1 - a) * (1 - a)* (1 - a);
	CC[2][4] = Half / Three / Four * (1 - a) * (1 - a)* (1 - a)* (1 - a);

	for (int i = 0; i < 5; ++i) {
		CC[3][i] = CC[1][i];
		CC[4][i] = CC[0][i];
	}

	std::map<std::string, Float> map;
	for (int i1 = 0; i1 < 5; i1++) {
		for (int i2 = 0; i2 < 5; i2++) {
			for (int i3 = 0; i3 < 5; i3++) {
				for (int i4 = 0; i4 < 5; i4++) {
					for (int i5 = 0; i5 < 5; i5++) {

						if (i1 + i2 + i3 + i4 + i5 >= 4) continue;

						std::string name;
						char last = '\0';
						if (i1 != 0) {
							last = 'V';
							name.push_back('V');
							name.push_back('0' + i1);
						}
						if (i2 != 0) {
							if (last == 'T') {
								name.back() += i2;
							} else {
								last = 'T';
								name.push_back('T');
								name.push_back('0' + i2);
							}
						}
						if (i3 != 0) {
							if (last == 'V') {
								name.back() += i3;
							} else {
								last = 'V';
								name.push_back('V');
								name.push_back('0' + i3);
							}
						}
						if (i4 != 0) {
							if (last == 'T') {
								name.back() += i4;
							} else {
								last = 'T';
								name.push_back('T');
								name.push_back('0' + i4);
							}
						}
						if (i5 != 0) {
							if (last == 'V') {
								name.back() += i5;
							} else {
								last = 'V';
								name.push_back('V');
								name.push_back('0' + i5);
							}
						}

						map[name] += CC[0][i1] * CC[1][i2] * CC[2][i3] * CC[3][i4] * CC[4][i5];


					}
				}
			}
		}
	}

	std::map<std::string, Float> map2;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < (1 << i); ++j) {
			std::string name;
			char last = '\0';
			for (int k = 0; k < i; ++k) {
				bool b = (1 << k) & j;
				char v = b ? 'V' : 'T';
				if (last == v) {
					name.back() += 1;
				} else {
					name.push_back(v);
					name.push_back('1');
				}
				last = v;
			}
			Float c = 1;
			for (int j = 1; j <= i; ++j) {
				c /= j;
			}
			map2[name] += c;
		}
	}
	for (auto b = map.begin(), b2 = map2.begin(); b != map.end() && b2 != map2.end(); ++b, ++b2) {
		printf("%10s %f\n", b->first.c_str(), b->second);
		printf("%10s %f\n", b2->first.c_str(), b2->second);
	}
}

// exp(0.5 a V) exp(0.5T) exp((1-a)V) exp(0.5T) exp(0.5 aV)
void solve4()
{
	Float CC[5][5];

	Float a(One/Three);
	CC[0][0] = 1;
	CC[0][1] = 0.5 * a;
	CC[0][2] = Half * (0.5*a)*(0.5*a);
	CC[0][3] = Half / Three * (0.5*a)*(0.5*a)*(0.5*a);
	CC[0][4] = Half / Three / Four * (0.5*a)*(0.5*a)*(0.5*a)*(0.5*a);

	CC[1][0] = 1;
	CC[1][1] = 0.5;
	CC[1][2] = Half * Float(0.5 * 0.5);
	CC[1][3] = (Half / Three) * Float(0.5 * 0.5 * 0.5);
	CC[1][4] = (Half / Three / Four) * Float(0.5*0.5*0.5*0.5);

	CC[2][0] = 1;
	CC[2][1] = 1 - a;
	CC[2][2] = Half * (1 - a)* (1 - a);
	CC[2][3] = Half / Three * (1 - a) * (1 - a)* (1 - a);
	CC[2][4] = Half / Three / Four * (1 - a) * (1 - a)* (1 - a)* (1 - a);

	for (int i = 0; i < 5; ++i) {
		CC[3][i] = CC[1][i];
		CC[4][i] = CC[0][i];
	}

	std::map<std::string, Float> map;
	for (int i1 = 0; i1 < 5; i1++) {
		for (int i2 = 0; i2 < 5; i2++) {
			for (int i3 = 0; i3 < 5; i3++) {
				for (int i4 = 0; i4 < 5; i4++) {
					for (int i5 = 0; i5 < 5; i5++) {

						if (i1 + i2 + i3 + i4 + i5 >= 5) continue;

						std::string name;
						char last = '\0';
						if (i1 != 0) {
							last = 'V';
							name.push_back('V');
							name.push_back('0' + i1);
						}
						if (i2 != 0) {
							if (last == 'T') {
								name.back() += i2;
							} else {
								last = 'T';
								name.push_back('T');
								name.push_back('0' + i2);
							}
						}
						if (i3 != 0) {
							if (last == 'V') {
								name.back() += i3;
							} else {
								last = 'V';
								name.push_back('V');
								name.push_back('0' + i3);
							}
						}
						if (i4 != 0) {
							if (last == 'T') {
								name.back() += i4;
							} else {
								last = 'T';
								name.push_back('T');
								name.push_back('0' + i4);
							}
						}
						if (i5 != 0) {
							if (last == 'V') {
								name.back() += i5;
							} else {
								last = 'V';
								name.push_back('V');
								name.push_back('0' + i5);
							}
						}

						map[name] += CC[0][i1] * CC[1][i2] * CC[2][i3] * CC[3][i4] * CC[4][i5];


					}
				}
			}
		}
	}

	std::map<std::string, Float> map2;
	for(int i = 0; i < 5; ++i)
	{
		for (int j = 0; j < (1 << i); ++j) {
			std::string name;
			char last = '\0';
			for (int k = 0; k < i; ++k) {
				bool b = (1 << k) & j;
				char v = b ? 'V' : 'T';
				if (last == v) {
					name.back() += 1;
				} else {
					name.push_back(v);
					name.push_back('1');
				}
				last = v;
			}
			Float c = 1;
			for (int j = 1; j <= i; ++j) {
				c /= j;
			}
			map2[name] += c;
		}
	}
	for (auto b = map.begin(), b2 = map2.begin(); b != map.end() && b2 != map2.end(); ++b, ++b2) {
		printf("%10s %f\n", b->first.c_str(), b->second);
		printf("%10s %f\n", b2->first.c_str(), b2->second);
	}
}

int main()
{
	//float_precision_ctrl.precision(100);
	solve2();
	solve3();

    return 0;
}

