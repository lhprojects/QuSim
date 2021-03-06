#define _CRT_SECURE_NO_WARNINGS
#include "header.h"
#include "QuSolver.h"
#include <thread>
#include <mutex>
#include <vector>

// core
#include <QuSim.h>
#include "../src/View.h"


// windows theme
#pragma comment(linker,"\"/manifestdependency:type='win32' \
name='Microsoft.Windows.Common-Controls' version='6.0.0.0' \
processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

#pragma comment(lib, "Comctl32.lib")  
// gdi plus
#pragma comment(lib, "gdiplus.lib")  
using namespace Gdiplus;

// user message
#define WM_DATAREADY (WM_USER + 12)


// windows
HWND hPotential;
HWND hInitalPsi;
HWND hPsi;
HWND hPsiPrime;
HWND hE;
HWND hShowPsi;
HWND hShowPotential;
HWND hX0;
HWND hX1;
HWND hBins;
HWND hRun;
HWND hCanvas;
HWND hMainWin;
HWND hHbar;
HWND hMass;
HWND hMidpoint;
HWND hGauss;
HWND hRK4;
HWND hOpt;
HWND hRng;
HWND hFor;

std::vector<HWND> enableControlList;
void enableAllWindows(bool enable)
{
    for (auto hwnd : enableControlList) {
        EnableWindow(hwnd, enable);
    }
}

struct Data {
    bool fReady;
    Real fT;
    Real fR;
    Complex fFinalPsi;
    Complex fFinalPsiPrime;
    std::vector<Complex> psi;
    std::vector<Complex> v;
    Int fN;
};

Data gData;
Solver1D syst;


#define MAX_LOADSTRING 100


HINSTANCE hInst;
WCHAR szTitle[MAX_LOADSTRING];
WCHAR szWindowClass[MAX_LOADSTRING];
WCHAR szCanvasWindowClass[] = { L"Canvas" };


ATOM                MyRegisterClass(HINSTANCE hInstance);
ATOM                RegisterCanvasClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
LRESULT CALLBACK    CanvasWndProc(HWND, UINT, WPARAM, LPARAM);

INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
    _In_opt_ HINSTANCE hPrevInstance,
    _In_ LPWSTR    lpCmdLine,
    _In_ int       nCmdShow)
{
    // high dpi
    SetProcessDPIAware();

    // windows theme
    InitCommonControls();

    ULONG_PTR gdiplusToken;
    GdiplusStartupInput gdiplusStartupInput;
    GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);


    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_QUSOLVER, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);
    RegisterCanvasClass(hInstance);

    if (!InitInstance(hInstance, nCmdShow)) {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_QUSOLVER));

    MSG msg;

    while (GetMessage(&msg, nullptr, 0, 0)) {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg)) {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    GdiplusShutdown(gdiplusToken);

    return (int)msg.wParam;
}



ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc = WndProc;
    wcex.cbClsExtra = 0;
    wcex.cbWndExtra = 0;
    wcex.hInstance = hInstance;
    wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDC_QUSOLVER));
    wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground = (HBRUSH)(CTLCOLOR_STATIC);
    wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_QUSOLVER);
    wcex.lpszClassName = szWindowClass;
    wcex.hIconSm = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}


ATOM RegisterCanvasClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc = CanvasWndProc;
    wcex.cbClsExtra = 0;
    wcex.cbWndExtra = 0;
    wcex.hInstance = hInstance;
    wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDC_QUSOLVER));
    wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_QUSOLVER);;
    wcex.lpszClassName = szCanvasWindowClass;
    wcex.hIconSm = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
    hInst = hInstance;

#if 0
    HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN,
        CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);
#else
    HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);
#endif

    hMainWin = hWnd;

    if (!hWnd) {
        return FALSE;
    }

    ShowWindow(hWnd, nCmdShow);
    UpdateWindow(hWnd);

    return TRUE;
}

#include <algorithm>

void DrawPotential(Gdiplus::Graphics &graphics,
    std::vector<Complex> const &V,
    long left, long top, long w, long h)
{
    std::vector<Point> vc;

    auto it = std::max_element(V.begin(), V.end(), [](Complex u, Complex v){
        return std::abs(u) < std::abs(v);
    });
    double vmax = std::abs(*it);

    // draw potential axis
    Real axis_max = 0;
    Real vmax_log10 = log10(vmax);
    Real int_vmax_log10 = floor(vmax_log10);
    Real remain = vmax_log10 - int_vmax_log10;
    if (remain > log10(5)) {
        axis_max = pow(10, int_vmax_log10 + 1);
    } else {
        axis_max = 5 * pow(10, int_vmax_log10);
    }

    Pen      black_pen(Color(255, 0, 0, 0));
    Font myFont(L"Courier New", 14);
    StringFormat format;
    format.SetAlignment(StringAlignmentNear);
    SolidBrush blackBrush(Color(255, 0, 0, 0));

    for (int i = -10; i <= 10; ++i) {
        Point b, e;

        b.X = w - 10;
        e.X = w;
        b.Y = e.Y = (int)((-1.0*i / 20 + 0.5)*h);
        graphics.DrawLine(&black_pen, b, e);

        wchar_t bf[100];
        swprintf(bf, 100, L"% 5g", (double)(i / 10.0*axis_max));
        PointF layoutRect((float)(b.X - 100), (float)(b.Y - 10));

        graphics.DrawString(bf, lstrlenW(bf), &myFont, layoutRect, &blackBrush);
    }


    // draw potential lines
    for (int i = 0; i < w; i += 1) {
        int xi = (int)(1.0 * i / w * (V.size() - 1));
        int v = 0;
        if (vmax != 0) {
            v = (int)(-std::abs(V[xi]) / axis_max * h / 2 + h / 2);
        }
        vc.push_back(Point(i, v));
    }

    Pen      pen4(Color(255, 0, 255, 0));
    graphics.DrawLines(&pen4, vc.data(), (int)vc.size());

}

void DrawPsi(Gdiplus::Graphics &graphics,
    std::vector<Complex> const &psi,
    long left, long top, long w, long h)
{
    std::vector<Point> abc;
    std::vector<Point> rec;
    std::vector<Point> imc;
    std::vector<Point> vc;


    for (int i = 0; i < w; i += 1) {
        int xi = (int)((1.0 * i / w * (psi.size() - 1)));
        int ab = (int)(-abs(psi[xi]) * 50 + h / 2);
        int re = (int)(-psi[xi].real() * 50 + h / 2);
        int im = (int)(-psi[xi].imag() * 50 + h / 2);
        abc.push_back(Point(i, ab));
        rec.push_back(Point(i, re));
        imc.push_back(Point(i, im));
    }

    Pen      pen1(Color(255, 0, 0, 255));
    Pen      pen2(Color(255, 0, 255, 255));
    Pen      pen3(Color(255, 255, 0, 255));
    graphics.DrawCurve(&pen1, abc.data(), (int)abc.size());
    graphics.DrawCurve(&pen2, rec.data(), (int)rec.size());
    graphics.DrawCurve(&pen3, imc.data(), (int)imc.size());

}

void OnPaint1(Gdiplus::Graphics &graphics, long left, long top, long w, long h)
{
    if (!gData.fReady) return;

    graphics.Clear(Color(255, 255, 255, 255));


    Font myFont(L"Courier New", 14);
    RectF layoutRect((float)(left + 25), (float)(top + 20), 500.0f, 800.0f);
    StringFormat format;
    format.SetAlignment(StringAlignmentNear);
    SolidBrush blackBrush(Color(255, 0, 0, 0));
    static const int len = 300;
    wchar_t text[len];
    wchar_t *b = text;

    b += swprintf(b, text + len - b, L"T           % .16g\n", gData.fT);
    b += swprintf(b, text + len - b, L"R           % .16g\n", gData.fR);
    b += swprintf(b, text + len - b, L"T+R         % .16g\n", gData.fR + gData.fT);
    b += swprintf(b, text + len - b, L"Final Psi  (% .5g, % .5g)\n", gData.fFinalPsi.real(), gData.fFinalPsi.imag());
    b += swprintf(b, text + len - b, L"Final Psi' (% .5g, % .5g)\n", gData.fFinalPsiPrime.real(), gData.fFinalPsiPrime.imag());

    graphics.DrawString(text, lstrlenW(text), &myFont, layoutRect, &format, &blackBrush);

    Pen      blackPen(Color(255, 0, 0, 0));

    graphics.DrawLine(&blackPen, Point(0, h / 2), Point(w, h / 2));

    DrawPotential(graphics, gData.v, left, top, w, h);
    DrawPsi(graphics, gData.psi, left, top, w, h);
}

void GetDouble(HWND hx, double &db, char const *inf)
{
    std::vector<char> psiStr;
    int len = GetWindowTextLength(hx);
    psiStr.resize(len + 1);
    GetWindowTextA(hx, psiStr.data(), len + 1);
    int na = sscanf(psiStr.data(), "%lf", &db);
    if (na < 1) {
        throw std::runtime_error(std::string("can't parse '") + psiStr.data() + "' as floating point for " + inf);
    }
}


struct IntialParameters {
    double x0;
    double x1;
    double mass;
    double hbar;
    double en;
    size_t n;
    bool fn;
    bool fnes;
    Complex initPsi;
    Complex initPsiPrime;
    SolverMethod met;
    std::function<Complex(Real)> v;

    IntialParameters(
        std::function<Complex(Real)> const & v,
        Real x0,
        Real x1,
        size_t n,
        Real en,
        Complex initPsi,
        Complex initPsiPrime,
        SolverMethod met,
        Real mass,
        Real hbar) : v(v), x0(x0), x1(x1), n(n), en(en), initPsi(initPsi),
        initPsiPrime(initPsiPrime), met(met), mass(mass), hbar(hbar)
    {

    }

};

IntialParameters InitParameters();



void SolveEnergy()
{
    double rng;
    GetDouble(hRng, rng, "energy range");
    rng = abs(rng);

    IntialParameters pars = InitParameters();

    std::vector<char> psiStr;
    {
        int len = GetWindowTextLength(hFor);
        psiStr.resize(len + 1);
        GetWindowTextA(hFor, psiStr.data(), len + 1);
    }
    Calculator cal(psiStr.data());

    double e1;
    double v1;
    {
        Solver1D syst1;
        syst1.init(pars.v, pars.x0, pars.x1, pars.n, pars.en - rng, pars.initPsi,
            pars.initPsiPrime, pars.met, pars.mass, pars.hbar, Options());

        syst1.Compute();
        cal.SetVaraible("Psi", syst1.FinalPsi());
        cal.SetVaraible("PsiPrime", syst1.FinalPsiPrime());
        cal.SetVaraible("E", syst1.GetEnergy());
        v1 = cal.Evaluate().real();
        e1 = syst1.GetEnergy();
    }

    double e2 = pars.en + rng;
    double change = 0;
    int i = 0;
    for (; i < 20; ++i) {

        Solver1D syst1;
        syst1.init(pars.v, pars.x0, pars.x1, pars.n, e2, pars.initPsi,
            pars.initPsiPrime, pars.met, pars.mass, pars.hbar, Options());

        syst1.Compute();
        cal.SetVaraible("Psi", syst1.FinalPsi());
        cal.SetVaraible("PsiPrime", syst1.FinalPsiPrime());
        cal.SetVaraible("E", syst1.GetEnergy());
        double v2 = cal.Evaluate().real();
        if (v2 == 0) { break; }
        if (v1 == v2) break;
        double f = v2 / (v2 - v1);
        double e = f*e1 + (1-f) * e2;

        double last_change = change;
        change = abs(e - e2);
        if (i > 5 && change == 0) break;
        if (i > 5 && change > last_change) break;
        v1 = v2;
        e1 = e2;
        e2 = e;
    }

    if (i == 20) {
        throw std::runtime_error("iters: 20, not converge");
    }

    char b[32];
    sprintf(b, "%.17E", e2);
    SetWindowTextA(hE, b);
    SendMessage(hMainWin, WM_COMMAND, MAKELONG(0, BN_CLICKED), (LPARAM)hRun);



}

void InitialASystem1D(Solver1D &syst, IntialParameters&pars)
{
    Options opts;
    opts.ComplexPotential();
    syst.init(pars.v, pars.x0, pars.x1, pars.n, pars.en, pars.initPsi,
        pars.initPsiPrime, pars.met, pars.mass, pars.hbar,
        opts);
}

IntialParameters InitParameters()
{
    std::string pot;
    double x0 = -1;
    double x1 = 1;
    double mass = 1;
    double hbar = 1;
    double e = 1;
    int n = 100;
    bool fn = false;
    bool fnes = false;
    Complex dx;
    Complex psi;
    Complex psiprime;
    SolverMethod sl;


    if (SendMessage(hMidpoint, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        sl = SolverMethod::ImplicitMidpointMethod;
    } else if (SendMessage(hGauss, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        sl = SolverMethod::GaussLegendreO4;
    } else if (SendMessage(hRK4, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        sl = SolverMethod::ExplicitRungeKuttaO4Classical;
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hMass);
        psiStr.resize(len + 1);
        GetWindowTextA(hMass, psiStr.data(), len + 1);
        int na = sscanf(psiStr.data(), "%lf", &mass);
        if (na < 1) {
            throw std::runtime_error("can't parse mass");
        }
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hHbar);
        psiStr.resize(len + 1);
        GetWindowTextA(hHbar, psiStr.data(), len + 1);
        int na = sscanf(psiStr.data(), "%lf", &hbar);
        if (na < 1) {
            throw std::runtime_error("can't parse hbar");
        }
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hX0);
        psiStr.resize(len + 1);
        GetWindowTextA(hX0, psiStr.data(), len + 1);
        int na = sscanf(psiStr.data(), "%lf", &x0);
        if (na < 1) {
            throw std::runtime_error("can't parse x0");
        }
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hX1);
        psiStr.resize(len + 1);
        GetWindowTextA(hX1, psiStr.data(), len + 1);
        int na = sscanf(psiStr.data(), "%lf", &x1);
        if (na < 1) {
            throw std::runtime_error("can't parse x1");
        }
    }


    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hE);
        psiStr.resize(len + 1);
        GetWindowTextA(hE, psiStr.data(), len + 1);
        int na = sscanf(psiStr.data(), "%lf", &e);
        if (na < 1) {
            throw std::runtime_error("can't parse e");
        }
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hPsi);
        psiStr.resize(len + 1);
        GetWindowTextA(hPsi, psiStr.data(), len + 1);
        Calculator cal(psiStr.data());
        cal.SetVaraible("E", e);
        psi = cal.Evaluate();
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hPsiPrime);
        psiStr.resize(len + 1);
        GetWindowTextA(hPsiPrime, psiStr.data(), len + 1);
        Calculator cal(psiStr.data());
        cal.SetVaraible("Psi", psi);
        cal.SetVaraible("E", e);
        psiprime = cal.Evaluate();
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hBins);
        psiStr.resize(len + 1);
        GetWindowTextA(hBins, psiStr.data(), len + 1);
        int na = sscanf(psiStr.data(), "%d", &n);
        if (na < 1) {
            throw std::runtime_error("can't parse bins");
        }
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hPotential);
        psiStr.resize(len + 1);
        GetWindowTextA(hPotential, psiStr.data(), len + 1);
        pot.assign(psiStr.data(), len);
    }


    return IntialParameters(FunctorWrapper(pot.data()), x0, x1, n, e, psi, psiprime, sl, mass, hbar);
}

VOID OnPaint(HDC hdc, long left, long top, long w, long h)
{
    Bitmap bitmap(w, h);
    Gdiplus::Graphics graphics(&bitmap);

    graphics.TranslateTransform(0, 15);
    OnPaint1(graphics, left, top, w, h - 30);

    Graphics gr(hdc);
    CachedBitmap cbitmap(&bitmap, &gr);
    gr.DrawCachedBitmap(&cbitmap, 0, 0);
}

LRESULT CALLBACK CanvasWndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message) {
    case WM_ERASEBKGND:
    {
        // no need to erase background
        return 1;
    }
    break;
    case WM_PAINT:
    {
        RECT rect;
        GetClientRect(hWnd, &rect);

        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hWnd, &ps);
        OnPaint(hdc, rect.left, rect.top, rect.right - rect.left, rect.bottom - rect.top);
        EndPaint(hWnd, &ps);
    }
    break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}


LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message) {
    case WM_CREATE:
    {
        RECT client_rect;
        GetClientRect(hWnd, &client_rect);

        NONCLIENTMETRICS metrics = {};
        metrics.cbSize = sizeof(metrics);
        SystemParametersInfo(SPI_GETNONCLIENTMETRICS, metrics.cbSize, &metrics, 0);

        HFONT guiFont = CreateFontIndirect(&metrics.lfCaptionFont);

        int const x0 = 30;
        int const h = 30;
        int x = x0;
        int y = 30;


        ////////////////////////
        x = x0;
        y += 30;
        HWND hs3 = CreateWindow(
            TEXT("STATIC"), TEXT("Potential"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 100, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hs3, WM_SETFONT, (LPARAM)guiFont, true);
        x += 100;

        hPotential = CreateWindow(TEXT("EDIT"), TEXT("exp(-x*x)"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 500, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hPotential, WM_SETFONT, (LPARAM)guiFont, true);
        x += 500;
        enableControlList.push_back(hPotential);

        hShowPotential = CreateWindow(
            TEXT("BUTTON"), TEXT("Show"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTOCHECKBOX,
            x, y, 100, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hShowPotential, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hShowPotential, BM_SETCHECK, BST_CHECKED, NULL);
        x += 100;
        enableControlList.push_back(hShowPotential);

        ////////////////////////
        x = x0;
        y += 30;

        HWND static_en = CreateWindow(
            TEXT("STATIC"), TEXT("Energy"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 100, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_en, WM_SETFONT, (LPARAM)guiFont, true);
        x += 100;

        hE = CreateWindow(TEXT("EDIT"), TEXT("0.5"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 240, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hE, WM_SETFONT, (LPARAM)guiFont, true);
        x += 240;
        enableControlList.push_back(hE);

        HWND static_psi = CreateWindow(
            TEXT("STATIC"), TEXT("Psi(x0)"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 100, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_psi, WM_SETFONT, (LPARAM)guiFont, true);
        x += 100;

        hPsi = CreateWindow(TEXT("EDIT"), TEXT("1"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 160, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hPsi, WM_SETFONT, (LPARAM)guiFont, true);
        x += 160;
        enableControlList.push_back(hPsi);


        HWND static_psiprime = CreateWindow(
            TEXT("STATIC"), TEXT("Psi'(x0)"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 100, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_psiprime, WM_SETFONT, (LPARAM)guiFont, true);
        x += 100;

        hPsiPrime = CreateWindow(TEXT("EDIT"), TEXT("-I"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 160, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hPsiPrime, WM_SETFONT, (LPARAM)guiFont, true);
        x += 160;
        enableControlList.push_back(hPsiPrime);


        HWND static_x0 = CreateWindow(
            TEXT("STATIC"), TEXT("x0"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 50, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_x0, WM_SETFONT, (LPARAM)guiFont, true);
        x += 50;

        hX0 = CreateWindow(TEXT("EDIT"), TEXT("-10"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hX0, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hX0);

        HWND static_x1 = CreateWindow(
            TEXT("STATIC"), TEXT("x1"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 50, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_x1, WM_SETFONT, (LPARAM)guiFont, true);
        x += 50;

        hX1 = CreateWindow(TEXT("EDIT"), TEXT("10"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hX1, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hX1);

        HWND static_bins = CreateWindow(
            TEXT("STATIC"), TEXT("Bins"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 50, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_bins, WM_SETFONT, (LPARAM)guiFont, true);
        x += 50;

        hBins = CreateWindow(TEXT("EDIT"), TEXT("1000"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hBins, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hBins);
        

        hShowPsi = CreateWindow(
            TEXT("BUTTON"), TEXT("Show"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTOCHECKBOX,
            x, y, 230, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hShowPsi, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hShowPsi, BM_SETCHECK, BST_CHECKED, NULL);
        x += 230;
        enableControlList.push_back(hShowPsi);

        ////////////////////////
        x = x0;
        y += h;

        hRun = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Compute"),
            WS_TABSTOP | WS_VISIBLE | WS_CHILD,
            x, y, 100, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hRun, WM_SETFONT, (LPARAM)guiFont, true);
        x += 100;

        HWND hs1 = CreateWindow(
            TEXT("STATIC"), TEXT("Solver:"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hs1, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;

        hMidpoint = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Midpoint"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
            x, y, 120, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hMidpoint, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hMidpoint, BM_SETCHECK, BST_CHECKED, NULL);

        x += 120;
        enableControlList.push_back(hMidpoint);

        hGauss = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Gauss"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
            x, y, 120, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hGauss, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hGauss);

        hRK4 = CreateWindow(
            TEXT("BUTTON"),
            TEXT("RK4"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
            x, y, 120, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hRK4, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hRK4);

        HWND static_hbar = CreateWindow(
            TEXT("STATIC"), TEXT("hbar"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 60, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_hbar, WM_SETFONT, (LPARAM)guiFont, true);
        x += 60;

        hHbar = CreateWindow(TEXT("EDIT"), TEXT("1"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hHbar, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hHbar);

        HWND static_mass = CreateWindow(
            TEXT("STATIC"), TEXT("mass"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 60, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_mass, WM_SETFONT, (LPARAM)guiFont, true);
        x += 60;

        hMass = CreateWindow(TEXT("EDIT"), TEXT("1"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hMass, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hMass);

        /////////////////////////////////////
        x = 30;
        y += h;
        hOpt = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Solve"),
            WS_TABSTOP | WS_VISIBLE | WS_CHILD,
            x, y, 100, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hOpt, WM_SETFONT, (LPARAM)guiFont, true);
        x += 100;

        HWND static_rng = CreateWindow(
            TEXT("STATIC"), TEXT(" energy in range \u00B1 "),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 150, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_rng, WM_SETFONT, (LPARAM)guiFont, true);
        x += 150;

        hRng = CreateWindow(TEXT("EDIT"), TEXT("0.01"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hRng, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hRng);

        HWND static_opt = CreateWindow(
            TEXT("STATIC"), TEXT(" for "),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 60, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_opt, WM_SETFONT, (LPARAM)guiFont, true);
        x += 60;

        hFor = CreateWindow(TEXT("EDIT"), TEXT("Psi"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hFor, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hFor);

        HWND static_not = CreateWindow(
            TEXT("STATIC"), TEXT("= 0 at end point"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_not, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;



        ////////////////////////
        y += h;
        hCanvas = CreateWindow(
            szCanvasWindowClass, TEXT("Canvas"),
            WS_CHILD | WS_VISIBLE,
            0, y, client_rect.right, client_rect.bottom,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        x += 180;

        // When you're done with the font, don't forget to call
        //DeleteObject(guiFont);
    }
    break;
    case WM_SIZE:
    {
        RECT client_rect;
        GetClientRect(hWnd, &client_rect);

        SetWindowPos(hCanvas, NULL, 0, 180,
            client_rect.right, client_rect.bottom - 180, 0);
    }
    break;
    case WM_COMMAND:
    {

        if (lParam != 0) 
        { // from button
            if ((HWND)lParam == hRun && HIWORD(wParam) == BN_CLICKED) {

                try {
                    Solver1D syst;
                    auto pars = InitParameters();
                    InitialASystem1D(syst, pars);
                    syst.Compute();
                    gData.psi = ToVector(syst.GetPsi());
                    gData.fN = syst.GetNPoints();

                    gData.v.resize(pars.n);
                    for (int i = 0; i < pars.n; ++i) {
                        gData.v[i] = pars.v(pars.x0 + i*(pars.x1 - pars.x0)/pars.n);
                    }

                    gData.fT = syst.GetT();
                    gData.fR = syst.GetR();
                    gData.fFinalPsi = syst.FinalPsi();
                    gData.fFinalPsiPrime = syst.FinalPsiPrime();
                    gData.fReady = true;

                    InvalidateRect(hWnd, NULL, false);

                } catch (std::exception const &exp) {
                    wchar_t w[1000];
                    mbstowcs(w, exp.what(), 100);
                    MessageBox(hWnd, w, L"Error", MB_OK);
                } catch (...) {
                    MessageBox(hWnd, L"err", L"err", MB_OK);
                }
            } else if ((HWND)lParam == hOpt && HIWORD(wParam) == BN_CLICKED) {
                try {
                    SolveEnergy();

                } catch (std::exception const &exp) {
                    wchar_t w[1000];
                    mbstowcs(w, exp.what(), 100);
                    MessageBox(hWnd, w, L"Error", MB_OK);
                } catch (...) {
                    MessageBox(hWnd, L"err", L"err", MB_OK);
                }
            }
        } else {
            int wmId = LOWORD(wParam);

            switch (wmId) {
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }

        }
    }
    break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message) {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL) {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
