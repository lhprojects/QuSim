#ifndef GUI_H
#define GUI_H

#include <nana/gui.hpp>
#include <nana/paint/pixel_buffer.hpp>
#include <nana/gui/widgets/button.hpp>
#include <nana/gui/widgets/label.hpp>
#include <nana/gui/widgets/textbox.hpp>
#include <nana/gui/widgets/checkbox.hpp>
#include <nana/gui/widgets/group.hpp>
#include <nana/gui/place.hpp>

#include <type_traits>
#include <locale>
#include <string>

inline std::string & ltrim(std::string & str)
{
	auto it2 = std::find_if(str.begin(), str.end(), [](char ch) { return !std::isspace<char>(ch, std::locale::classic()); });
	str.erase(str.begin(), it2);
	return str;
}

inline std::string & rtrim(std::string & str)
{
	auto it1 = std::find_if(str.rbegin(), str.rend(), [](char ch) { return !std::isspace<char>(ch, std::locale::classic()); });
	str.erase(it1.base(), str.end());
	return str;
}

inline std::string & trim(std::string & str)
{
	ltrim(rtrim(str));
	return str;
}


struct Button : nana::button {
	Button(nana::form &fm, std::string const &caption) : nana::button(fm, caption) {
		edge_effects(false);
	}
};

struct TextBox : nana::textbox {
	TextBox(nana::form &fm, std::string const &caption) : nana::textbox(fm, caption) {
	}
	TextBox(nana::form &fm) : nana::textbox(fm) {
	}

	double GetDouble() {
		std::string text = caption();
		return stod(trim(text));
	}
	int GetInt() {
		std::string text = caption();
		return stoi(trim(text));
	}
};

struct WaveFunction1D : nana::panel<true> {

	WaveFunction1D(nana::form &fm) : nana::panel<true>(fm) {
		bgcolor(nana::colors::black);
		nana::drawing drawing(*this);
		drawing.draw(std::bind(&WaveFunction1D::Draw, this, std::placeholders::_1));
		fShowPotential = true;
		fAddInit = true;
		fShowWaveFunction = true;
	}

	std::shared_ptr<QuScatteringProblemSolver1D> fSolver;
	bool fShowPotential;
	bool fAddInit;
	bool fShowWaveFunction;

	void DrawLine(nana::paint::graphics &g, std::vector<nana::point> const &ls, nana::color c) {
		if (ls.size() >= 2) {
			int lastx = ls[0].x;
			for (int i = 0; i < (int)ls.size(); i += 1) {
				if (i == 0) g.line_begin(ls[i].x, ls[i].y);
				else {
					if (ls[i].x > lastx) {
						g.line_to(ls[i], c);
						lastx = ls[i].x;
					}
				}
			}
		}
	}

	void DrawPotential(nana::paint::graphics &g, VectorView<double> const &pot) {
		if (pot.size() < 2) return;

		double vmin = *std::min_element(pot.begin(), pot.end());
		double vmax = *std::max_element(pot.begin(), pot.end());
		vmax = 1.4*std::max(abs(vmin), abs(vmax));

		int w = g.width();
		int h = g.height();
		using namespace nana;
		std::vector<point> abc;

		for (int i = 0; i < w; i += 1) {
			int xi = (int)((1.0 * i / w * (pot.size() - 1)));
			int v = (int)(0.5*(1 - pot[xi] / vmax)*h);
			abc.push_back(point(i, v));
		}

		DrawLine(g, abc, color(0, 255, 255));

	}

	void DrawPsi(nana::paint::graphics &g, std::vector<Complex> const &psi) {
		if (psi.size() < 2) return;

		int w = g.width();
		int h = g.height();
		using namespace nana;
		std::vector<point> abc;
		std::vector<point> rec;
		std::vector<point> imc;
		std::vector<point> vc;

		for (int i = 0; i < w; i += 1) {
			int xi = (int)((1.0 * i / w * (psi.size() - 1)));
			int ab = (int)(-abs(psi[xi]) * 100 + h / 2);
			int re = (int)(-psi[xi].real() * 100 + h / 2);
			int im = (int)(-psi[xi].imag() * 100 + h / 2);
			abc.push_back(point(i, ab));
			rec.push_back(point(i, re));
			imc.push_back(point(i, im));
		}

		DrawLine(g, abc, color(255, 0, 0));
		DrawLine(g, rec, color(255, 0, 255));
		DrawLine(g, imc, color(255, 255, 0));

	}

	void Draw(nana::paint::graphics &g) {
		using namespace nana;
		if (fSolver) {
			auto txcolor = colors::white;

			nana::paint::font font("Courier New", 14);
			g.typeface(font);

			auto const width = g.width();
			auto const height = g.height();
			int x0 = 20;
			int dy = 20;
			int y = 0;
			char b[1000];

			sprintf(b, "%9s: %12g", "Energy", fSolver->GetEnergy());
			g.string(point(x0, y += dy), b, txcolor);
			sprintf(b, "%9s: %12g", "T", fSolver->GetT());
			g.string(point(x0, y += dy), b, txcolor);
			sprintf(b, "%9s: %12g", "R", fSolver->GetR());
			g.string(point(x0, y += dy), b, txcolor);


			g.line(point(0, height / 2), point(width, height / 2), txcolor);

			if (fShowPotential)
				DrawPotential(g, fSolver->GetV());
			if (fShowWaveFunction) {
				std::vector<Complex> psi(fSolver->GetPsi().begin(),
					fSolver->GetPsi().end());
				if (fAddInit) {
					auto psi0 = fSolver->GetPsi0();
					for (int i = 0; i < (int)psi0.size(); ++i)
						psi[i] += psi0(i);
				}
				DrawPsi(g, psi);
			}
		}
	}

};

#endif
