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
#include <QuSim.h>
#include <exception>

std::string & ltrim(std::string & str)
{
    auto it2 = std::find_if(str.begin(), str.end(), [](char ch) { return !std::isspace<char>(ch, std::locale::classic()); });
    str.erase(str.begin(), it2);
    return str;
}

std::string & rtrim(std::string & str)
{
    auto it1 = std::find_if(str.rbegin(), str.rend(), [](char ch) { return !std::isspace<char>(ch, std::locale::classic()); });
    str.erase(it1.base(), str.end());
    return str;
}

std::string & trim(std::string & str)
{
    ltrim(rtrim(str));
    return str;
}

struct Button : nana::button {
    Button(nana::form &fm, std::string const &caption) : nana::button(fm, caption) {
        //edge_effects(false);
    }
};

struct TextBox : nana::textbox {
    TextBox(nana::form &fm, std::string const &caption) : nana::textbox(fm, caption) {
        multi_lines(false);
    }
    TextBox(nana::form &fm) : nana::textbox(fm) {
        multi_lines(false);
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


    struct Legend {

        int i = 0;
        nana::color text_color;
        nana::paint::graphics& g;

        Legend(nana::paint::graphics& g) : g(g) {}

        void Add(std::string const& str,
            nana::color line_color)
        {
            int w = g.width();
            g.line_begin(w - 100, 20 * i + 20);
            g.line_to(nana::point(w - 80, 20 * i + 20), line_color);
            g.string(nana::point(w - 70, 20 * i - 10 + 20), str, text_color);

            ++i;
        }
    };

    void DrawPotential(nana::paint::graphics &g, 
        VectorView<Complex> const &pot,
        Legend &leg) {
        if (pot.size() < 2) return;

        int w = g.width();
        int h = g.height();
        using namespace nana;
        std::vector<point> pv;
        std::vector<point> al;

        for (int i = 0; i < w; i += 1) {
            int xi = (int)((1.0 * i / w * (pot.size() - 1)));
            int v = (int)(-(+pot[xi].real()) * 100 + h / 2);
            int a = (int)(-(-pot[xi].imag()) * 100 + h / 2);

            pv.push_back(point(i, v));
            al.push_back(point(i, a));
        }

        DrawLine(g, pv, color(0, 255, 0));
        DrawLine(g, al, color(0, 0, 255));
        leg.Add("Pot", color(0, 255, 0));
        leg.Add("Absorb", color(0, 0, 255));

    }

    void DrawPsi(nana::paint::graphics &g,
        std::vector<Complex> const &psi,
        Legend &leg) {
        if (psi.size() < 2) return;

        int w = g.width();
        int h = g.height();
        using namespace nana;
        std::vector<point> abc;
        std::vector<point> rec;
        std::vector<point> imc;

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

        leg.Add("Abs", color(255, 0, 0));
        leg.Add("Real", color(255, 0, 255));
        leg.Add("Imag", color(255, 255, 0));

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

            Legend leg(g);
            leg.text_color = color(255, 255, 255);

            if (fShowPotential)
                DrawPotential(g, fSolver->GetV(), leg);
            if (fShowWaveFunction) {
                std::vector<Complex> psi(fSolver->GetPsi().begin(),
                    fSolver->GetPsi().end());
                if (fAddInit) {
                    auto psi0 = fSolver->GetPsi0();
                    for (int i = 0; i < (int)psi0.size(); ++i)
                        psi[i] += psi0(i);
                }
                DrawPsi(g, psi, leg);
            }
        }
    }

};

#if defined(_WIN32)
#include <Windows.h>
int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
    _In_opt_ HINSTANCE hPrevInstance,
    _In_ LPWSTR    lpCmdLine,
    _In_ int       nCmdShow)
#else
int main()
#endif
{
    nana::form fm;
    fm.caption("Scattering 1D");
    fm.size(nana::size(1000, 600));
    fm.bgcolor(nana::color(240, 240, 240));

    //fm.bgcolor(nana::colors::white);
    nana::label initWaveLabel(fm, "Intial Wave Function:");
    nana::label energyLabel(fm, "Energy");
    TextBox energy(fm);
    nana::label dxLabel(fm, "DirectionX");
    TextBox dx(fm);
    nana::label potentialLabel(fm, "Potential");
    TextBox potential(fm, "f(x)");
    nana::label x0Label(fm, "x0");
    TextBox x0(fm);
    nana::label x1Label(fm, "x1");
    TextBox x1(fm);
    nana::label binsLabel(fm, "Bins");
    TextBox bins(fm);
    nana::label massLabel(fm, "mass");
    TextBox mass(fm);
    nana::label hbarLabel(fm, "hbar");
    TextBox hbar(fm);

    Button compute(fm, "compute");
    Button increase(fm, "increase");
    nana::label solverLabel(fm, "solver");
    nana::radio_group sover;
    nana::checkbox spaceO2(fm, "MatInvSpaceO2");
    nana::checkbox spaceO4(fm, "MatInvSpaceO4");
    nana::checkbox bornser(fm, "BornSeries");
    sover.add(spaceO2);
    sover.add(spaceO4);
    sover.add(bornser);
    nana::label boundaryLabel(fm, "BoundaryCondtion");
    nana::radio_group boundary;
    nana::checkbox period(fm, "Period");
    boundary.add(period);
    nana::checkbox absorb(fm, "AbsorbationWall");
    Button redraw(fm, "Redraw");
    nana::checkbox drawWaveFunction(fm, "Wave Function");
    nana::checkbox addInit(fm, "Add Initial Wave Function");
    nana::checkbox drawPoential(fm, "Potential");
    WaveFunction1D wavefunction(fm);
    //period.enabled(false);
    period.check(true);
    absorb.enabled(false);


    nana::place layout(fm);
    char const *laystr = "vert\
<vert weight=180 gap=10 margin=[10,10,10,10]\
<l1 weight=30 fit gap=5 arrange=[variable,variable,200,variable,100,variable,100]>\
<l2 weight=30 fit gap=5 arrange=[variable,400]>\
<l3 weight=30 fit gap=5 arrange=[variable,100,variable,100,variable,100,variable,100]>\
<l4 weight=30 fit gap=5 arrange=[variable,variable,variable,120,120,110,variable,110,120]>\
<l5 weight=30 fit gap=5 arrange=[variable,120,170,120]>\
>\
<l6>\
";
    try {
        layout.div(laystr);
    }
    catch (std::exception &exp) {
        nana::msgbox(exp.what());
        exit(1);
    }
    layout.field("l1") << initWaveLabel;
    layout.field("l1") << energyLabel << energy;
    layout.field("l1") ;
    layout.field("l1") << dxLabel << dx;
    layout.field("l2") << potentialLabel << potential;
    layout.field("l3") << x0Label << x0;
    layout.field("l3") << x1Label << x1;
    layout.field("l3") << binsLabel << bins;
    layout.field("l4") << compute << increase;
    layout.field("l4") << solverLabel;
    layout.field("l4") << spaceO2;
    layout.field("l4") << spaceO4;
    layout.field("l4") << bornser;
    layout.field("l4") << boundaryLabel;
    layout.field("l4") << period;
    layout.field("l4") << absorb;
    layout.field("l5") << redraw << drawWaveFunction << addInit << drawPoential;
    layout.field("l6") << wavefunction;
    layout.collocate();



    spaceO2.check(true);
    absorb.check(true);
    period.check(true);
    drawWaveFunction.check(true);
    addInit.check(true);
    drawPoential.check(true);
    x0.caption("-150");
    x1.caption("150");
    bins.caption("2000");
    mass.caption("1");
    hbar.caption("1");
    energy.caption("0.5");
    dx.caption("1");

    potential.caption("exp(-x*x)");

    addInit.events().checked([&](nana::arg_checkbox const &check) {
        wavefunction.fAddInit = check.widget->checked();
    });

    drawPoential.events().checked([&](nana::arg_checkbox const &check) {
        wavefunction.fShowPotential = check.widget->checked();
    });

    drawWaveFunction.events().checked([&](nana::arg_checkbox const &check) {
        wavefunction.fShowWaveFunction = check.widget->checked();
    });

    redraw.events().click([&](nana::arg_click const &click) {
        // redraw
        //fm.collocate();
        //layout.collocate();
        //nana::API::update_window(wavefunction);
        nana::size sz = fm.size();
        fm.size(fm.size() + nana::size(1, 0));
        fm.size(sz);
    });

    increase.events().click([&](nana::arg_click const &click) {
        try {
            if (wavefunction.fSolver) {
                if (wavefunction.fSolver->GetMethod() == SolverMethod::BornSerise) {
                    wavefunction.fSolver->Compute();
                    redraw.events().click.emit(nana::arg_click(), fm.handle());
                }
            }
        } catch (std::exception &exp) {
            nana::msgbox error;
            error << exp.what();
            error.show();
        }
    });
    compute.events().click([&](nana::arg_click const &click) {

        try {

            double x0_ = x0.GetDouble();
            double x1_ = x1.GetDouble();
            int bins_ = bins.GetInt();
            double en = energy.GetDouble();
            double dx_ = dx.GetDouble();
            FunctorWrapper pot(potential.caption().c_str());
            SolverMethod met = SolverMethod::ImplicitMidpointMethod;
            double mass_ = mass.GetDouble();
            double hbar_ = hbar.GetDouble();

            Options opts;
            if (bornser.checked()) {
                opts.Order(1);
                opts.Preconditional(true);
                opts.VellekoopPreconditioner();
                met = SolverMethod::BornSerise;
            } else if (spaceO4.checked()) {
                met = SolverMethod::MatrixInverse;
                opts.SpaceOrder(4);
            } else if (spaceO2.checked()) {
                met = SolverMethod::MatrixInverse;
                opts.SpaceOrder(2);
            }

            if (met == SolverMethod::BornSerise) {
                auto solver = std::make_shared<QuPerturbation1D>();
                solver->init(pot, x0_, x1_, bins_, en, 0, dx_, met, mass_, hbar_, opts);
                solver->Compute();
                wavefunction.fSolver = solver;
            } else {
                auto solver = std::make_shared<QuScatteringInverseMatrix1D>();
                solver->init(pot, x0_, x1_, bins_, en, dx_, met, mass_, hbar_, opts);
                solver->Compute();
                wavefunction.fSolver = solver;
            }
            redraw.events().click.emit(nana::arg_click(), fm.handle());
        }
        catch (std::exception &exp) {
            nana::msgbox error;
            error << exp.what();
            error.show();
        }
    });



    fm.show();
    nana::exec();
}
