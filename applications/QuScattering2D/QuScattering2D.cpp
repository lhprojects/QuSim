#define _CRT_SECURE_NO_WARNINGS
#include "header.h"
#include "QuScattering2D.h"
#include <thread>
#include <mutex>
#include <memory>

// windows theme
#include <Commctrl.h>
// windows theme
#pragma comment(linker,"\"/manifestdependency:type='win32' \
name='Microsoft.Windows.Common-Controls' version='6.0.0.0' \
processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")
#pragma comment(lib, "Comctl32.lib")

// gdi plus
#include <objidl.h>
#include <GdiPlus.h> 
#include <GdiplusGraphics.h>  
// gdi plus
#pragma comment(lib, "gdiplus.lib")  
using namespace Gdiplus;

#include "../QuGui/WinGui.h"

// core
#include <QuSim.h>
#include "../src/View.h"


// windows
HWND hPotential;
HWND hEnergy;
HWND hDirectionX;
HWND hDirectionY;
HWND hShowPsi;
HWND hAddInitPsi;
HWND hShowPotential;
HWND hX0;
HWND hX1;
HWND hXBins;
HWND hY0;
HWND hY1;
HWND hYBins;
HWND hCompute;
HWND hIncrease;
HWND hRedraw;
HWND hBoundaryCondition;
HWND hCanvas;
HWND hMainWin;
HWND hPeriod;
HWND hAbsorbtionLayer;
HWND hHbar;
HWND hMass;
HWND hEigen;
HWND hMidpoint;
HWND hGaussO4;
HWND hPerturbation;

std::vector<HWND> enableControlList;
void enableAllWindows(bool enable)
{
    for (auto hwnd : enableControlList) {
        EnableWindow(hwnd, enable);
    }
}


struct Data {
    Eigen::MatrixXcd psiToDraw;
    Eigen::MatrixXcd v;
    Int fNx = 0;
    Int fNy = 0;
};

struct GUI {
    // gui thread
    bool data_ready = false;
    bool incrase = false;
    int times;
    Data data;
    std::unique_ptr<QuScatteringProblemSolver2D> solver;
};
GUI gui;


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
    LoadStringW(hInstance, IDC_QUSCATTERING2D, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);
    RegisterCanvasClass(hInstance);

    if (!InitInstance(hInstance, nCmdShow)) {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_QUSCATTERING2D));

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
    wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDC_QUSCATTERING2D));
    wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground = (HBRUSH)(CTLCOLOR_STATIC);
    wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_QUSCATTERING2D);
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
    wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDC_QUSCATTERING2D));
    wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_QUSCATTERING2D);;
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

void OnPaint1(Gdiplus::Graphics &graphics, long left, long top, long w, long h)
{

    graphics.Clear(Color(255, 0, 0, 0));

    if (!gui.data_ready) return;
    bool show_psi = SendMessage(hShowPsi, BM_GETCHECK, 0, 0) == BST_CHECKED;
    bool show_pot = SendMessage(hShowPotential, BM_GETCHECK, 0, 0) == BST_CHECKED;

    Pen      blackPen(Color(255, 0, 0, 0));
    graphics.DrawLine(&blackPen, Point(0, h / 2), Point(w, h / 2));

    if (show_pot) DrawPotential2D(graphics, gui.data.v, left, top, w, h);
    if (show_psi) DrawPsi2D(graphics, gui.data.psiToDraw, left, top, w, h);

    Font myFont(L"Courier New", 14);
    RectF layoutRect((float)left + 25, (float)top + 20, 500.0f, 800.0f);
    StringFormat format;
    format.SetAlignment(StringAlignmentNear);
    SolidBrush textBrush(Color(255, 0, 255, 0));
    static const int len = 300;
    wchar_t text[len];
    wchar_t *b = text;
    wsprintf(text, L"");

    graphics.DrawString(text, lstrlenW(text), &myFont, layoutRect, &format, &textBrush);

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

void GetInt(HWND hx, int &i, char const *inf)
{
    std::vector<char> psiStr;
    int len = GetWindowTextLength(hx);
    psiStr.resize(len + 1);
    GetWindowTextA(hx, psiStr.data(), len + 1);
    int na = sscanf(psiStr.data(), "%d", &i);
    if (na < 1) {
        throw std::runtime_error(std::string("can't parse '") + psiStr.data() + "' as integer for " + inf);
    }
}


void InitialASystem2D(GUI &gui)
{
    std::string psi;
    std::string pot;
    double x0 = -1;
    double x1 = 1;
    double y0 = -1;
    double y1 = 1;
    double mass = 1;
    double hbar = 1;
    double en = 1;
    double dx = 1;
    double dy = 0;
    int nx = 100;
    int ny = 100;
    bool fn = false;
    bool fnes = false;
    bool splito4 = false;
    Complex deltaT;
    BoundaryCondition bc;
    Options opts;
    SolverMethod sm;

    if (SendMessage(hMidpoint, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        opts.SpaceOrder(2);
        sm = SolverMethod::MatrixInverse;
    } else if (SendMessage(hGaussO4, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        opts.SpaceOrder(4);
        sm = SolverMethod::MatrixInverse;
    } else if (SendMessage(hPerturbation, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        sm = SolverMethod::BornSerise;
    }

    if (SendMessage(hPeriod, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        bc = BoundaryCondition::Period;
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hMass);
        psiStr.resize(len + 1);
        GetWindowTextA(hMass, psiStr.data(), len + 1);
        int na = sscanf(psiStr.data(), "%lf", &mass);
        if (na < 1) {
            throw std::runtime_error("can't parse x0");
        }
    }

    GetDouble(hX0, x0, "x0");
    GetDouble(hX1, x1, "x1");
    GetDouble(hY0, y0, "y0");
    GetDouble(hY1, y1, "y1");
    GetDouble(hEnergy, en, "Energy");
    GetDouble(hDirectionX, dx, "DirectionX");
    GetDouble(hDirectionY, dy, "DirectionY");
    GetInt(hXBins, nx, "# bins x");
    GetInt(hYBins, ny, "# bins y");

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hPotential);
        psiStr.resize(len + 1);
        GetWindowTextA(hPotential, psiStr.data(), len + 1);
        pot.assign(psiStr.data(), len);
    }



    if (sm == SolverMethod::MatrixInverse) {
        auto ptr = new QuScatteringInverseMatrix2D();
        gui.solver.reset(ptr);
        ptr->init(Functor2DWrapper(pot.c_str()), x0, x1, nx, y0, y1, ny,
            en, dx, dy, SolverMethod::MatrixInverse,
            mass, hbar, opts);
        gui.incrase = false;
    } else if (sm == SolverMethod::BornSerise) {
#ifdef QUSIM_USE_CUDA
        opts.Cuda();
#endif

        opts.Preconditional(true);
        opts.VellekoopPreconditioner();
        opts.Order(100);
        auto ptr = new QuPerturbation2D();
        gui.solver.reset(ptr);
        ptr->init(Functor2DWrapper(pot.c_str()), x0, x1, nx, y0, y1, ny,
            en, 0, dx, dy, SolverMethod::BornSerise,
            mass, hbar, opts);
        gui.incrase = true;
    }

}

void OnPaint2(Gdiplus::Graphics &graphics, long left, long top, long w, long h)
{

    graphics.Clear(Color(255, 255, 255, 255));

    Pen      blackPen(Color(255, 0, 0, 0));
    graphics.DrawLine(&blackPen, Point(0, h / 2), Point(w, h / 2));

    Evolver2D syst;

    bool show_psi = SendMessage(hShowPsi, BM_GETCHECK, 0, 0) == BST_CHECKED;
    bool show_pot = SendMessage(hShowPotential, BM_GETCHECK, 0, 0) == BST_CHECKED;

    std::string psi;
    std::string pot;
    double x0 = -1;
    double x1 = 1;
    double y0 = -1;
    double y1 = 1;
    int ny = 100;
    int nx = 100;
    bool fn = false;

    try {
        GetDouble(hX0, x0, "x0");
        GetDouble(hX1, x1, "x1");
        GetDouble(hY0, y0, "y0");
        GetDouble(hY1, y1, "y1");
        GetInt(hXBins, nx, "# bins x");
        GetInt(hYBins, ny, "# bins y");

        if (show_pot) {
            std::vector<char> psiStr;
            int len = GetWindowTextLength(hPotential);
            psiStr.resize(len + 1);
            GetWindowTextA(hPotential, psiStr.data(), len + 1);
            pot.assign(psiStr.data(), len);
            if (pot.empty()) show_pot = false;
        }
        Eigen::MatrixXd vv;
        Eigen::MatrixXcd psiv;

        if (show_pot) {
            vv.resize(ny, nx);

            Calculator cal(pot.c_str());
            cal.SetVaraible("x", 0);
            cal.SetVaraible("y", 0);

            Complex *px = &cal.GetVaraible("x");
            Complex *py = &cal.GetVaraible("y");

            for_each_2d_idx(x0, (x1 - x0) / nx, nx,
                y0, (y1 - y0) / ny, ny,
                [&](size_t i, size_t j, Real x, Real y) {
                *px = x;
                *py = y;
                vv(i, j) = cal.Evaluate().real();
            });
            DrawPotential2D(graphics, vv, left, top, w, h);
        }

    } catch (std::exception const &exp) {
        wchar_t w[100];
        mbstowcs(w, exp.what(), 100);
        MessageBox(hMainWin, w, L"Error", MB_OK);
    } catch (...) {
        MessageBox(hMainWin, L"err", L"err", MB_OK);
    }


}
VOID OnPaint(HDC hdc, long left, long top, long w, long h)
{
    Bitmap bitmap(w, h);
    Gdiplus::Graphics graphics(&bitmap);

    graphics.TranslateTransform(0, 15);
    if (gui.data_ready) {
        OnPaint1(graphics, left, top, w, h - 30);
    } else {
        OnPaint2(graphics, left, top, w, h - 30);
    }

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

        HWND hs11 = CreateWindow(
            TEXT("STATIC"), TEXT("Initial Wave Function"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 180, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hs11, WM_SETFONT, (LPARAM)guiFont, true);
        x += 180;

        HWND hs2 = CreateWindow(
            TEXT("STATIC"), TEXT("Energy"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hs2, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;

        hEnergy = CreateWindow(TEXT("EDIT"), TEXT("1"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hEnergy, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hEnergy);


        HWND hSDX = CreateWindow(
            TEXT("STATIC"), TEXT("DirectionX"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hSDX, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;

        hDirectionX = CreateWindow(TEXT("EDIT"), TEXT("1"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hDirectionX, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hDirectionX);

        HWND hSDY = CreateWindow(
            TEXT("STATIC"), TEXT("DirectionY"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hSDY, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;

        hDirectionY = CreateWindow(TEXT("EDIT"), TEXT("0"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 120, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hDirectionY, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hDirectionY);

        ////////////////////////
        x = x0;
        y += 30;
        HWND hs3 = CreateWindow(
            TEXT("STATIC"), TEXT("Potential"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 180, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hs3, WM_SETFONT, (LPARAM)guiFont, true);
        x += 180;

        hPotential = CreateWindow(TEXT("EDIT"), TEXT("exp(-x*x-y*y)"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 500, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hPotential, WM_SETFONT, (LPARAM)guiFont, true);
        x += 500;
        enableControlList.push_back(hPotential);


        ////////////////////////
        x = x0;
        y += 30;


        HWND static_x0 = CreateWindow(
            TEXT("STATIC"), TEXT("x0"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 50, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_x0, WM_SETFONT, (LPARAM)guiFont, true);
        x += 50;

        hX0 = CreateWindow(TEXT("EDIT"), TEXT("-100"),
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

        hX1 = CreateWindow(TEXT("EDIT"), TEXT("100"),
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

        hXBins = CreateWindow(TEXT("EDIT"), TEXT("200"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hXBins, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hXBins);

        HWND static_y0 = CreateWindow(
            TEXT("STATIC"), TEXT("y0"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 50, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_y0, WM_SETFONT, (LPARAM)guiFont, true);
        x += 50;

        hY0 = CreateWindow(TEXT("EDIT"), TEXT("-100"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hY0, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hY0);

        HWND static_y1 = CreateWindow(
            TEXT("STATIC"), TEXT("y1"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 50, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_y1, WM_SETFONT, (LPARAM)guiFont, true);
        x += 50;

        hY1 = CreateWindow(TEXT("EDIT"), TEXT("100"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hY1, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hY1);

        HWND static_ybins = CreateWindow(
            TEXT("STATIC"), TEXT("Bins"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 50, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_ybins, WM_SETFONT, (LPARAM)guiFont, true);
        x += 50;

        hYBins = CreateWindow(TEXT("EDIT"), TEXT("200"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hYBins, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hYBins);

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

        ////////////////////////
        x = x0;
        y += h;

        hCompute = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Compute"),
            WS_TABSTOP | WS_VISIBLE | WS_CHILD,
            x, y, 90, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hCompute, WM_SETFONT, (LPARAM)guiFont, true);
        x += 90;

        hIncrease = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Increase"),
            WS_TABSTOP | WS_VISIBLE | WS_CHILD,
            x, y, 90, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hIncrease, WM_SETFONT, (LPARAM)guiFont, true);
        x += 90;

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
            TEXT("SpaceO2"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON | WS_GROUP,
            x, y, 120, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hMidpoint, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hMidpoint, BM_SETCHECK, BST_CHECKED, NULL);
        x += 120;
        enableControlList.push_back(hMidpoint);

        hGaussO4 = CreateWindow(
            TEXT("BUTTON"),
            TEXT("SpaceO4"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
            x, y, 120, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hGaussO4, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hGaussO4);

        hPerturbation = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Perturbation"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
            x, y, 120, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hPerturbation, WM_SETFONT, (LPARAM)guiFont, true);
        x += 120;
        enableControlList.push_back(hPerturbation);


        HWND hs4 = CreateWindow(
            TEXT("STATIC"), TEXT("Boundary Condition:"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 160, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hs4, WM_SETFONT, (LPARAM)guiFont, true);
        x += 160;

        hPeriod = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Periodic"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON | WS_GROUP,
            x, y, 100, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hPeriod, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hPeriod, BM_SETCHECK, BST_CHECKED, NULL);
        x += 100;
        enableControlList.push_back(hPeriod);

        hAbsorbtionLayer = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Absorbtion Layer"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON | WS_GROUP,
            x, y, 160, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hAbsorbtionLayer, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hAbsorbtionLayer, BM_SETCHECK, BST_CHECKED, NULL);
        x += 160;
        enableControlList.push_back(hAbsorbtionLayer);

        ////////////////////////
        x = x0;
        y += h;

        hRedraw = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Redraw"),
            WS_TABSTOP | WS_VISIBLE | WS_CHILD,
            x, y, 180, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hRedraw, WM_SETFONT, (LPARAM)guiFont, true);
        x += 180;

        hShowPsi = CreateWindow(
            TEXT("BUTTON"), TEXT("Show Wave Function"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTOCHECKBOX,
            x, y, 200, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hShowPsi, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hShowPsi, BM_SETCHECK, BST_CHECKED, NULL);
        x += 200;
        enableControlList.push_back(hShowPsi);

        hAddInitPsi = CreateWindow(
            TEXT("BUTTON"), TEXT("Add Initial Wave Function"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTOCHECKBOX,
            x, y, 240, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hAddInitPsi, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hAddInitPsi, BM_SETCHECK, BST_CHECKED, NULL);
        x += 240;
        enableControlList.push_back(hAddInitPsi);

        hShowPotential = CreateWindow(
            TEXT("BUTTON"), TEXT("Show Potential"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTOCHECKBOX,
            x, y, 200, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hShowPotential, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hShowPotential, BM_SETCHECK, BST_CHECKED, NULL);
        x += 200;
        enableControlList.push_back(hShowPotential);

        ////////////////////////
        x = x0;
        y += h;
        hCanvas = CreateWindow(
            szCanvasWindowClass, TEXT("Canvas"),
            WS_CHILD | WS_VISIBLE,
            0, y, client_rect.right, client_rect.bottom,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hs2, WM_SETFONT, (LPARAM)guiFont, true);
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

        if (lParam != 0) { // from button

            if ((HWND)lParam == hCompute && HIWORD(wParam) == BN_CLICKED) {
                try {
                    InitialASystem2D(gui);
                    gui.times = 0;
                    SendMessage(hWnd, WM_COMMAND, MAKELONG(0, BN_CLICKED), (LPARAM)hIncrease);
                }
                catch (std::exception &e) {
                    gui.solver.reset();
                    MessageBoxA(hWnd, e.what(), "Error", MB_OK);
                }
            }  else if ((HWND)lParam == hIncrease && HIWORD(wParam) == BN_CLICKED) {
                try {
                    if (gui.solver.get() && (gui.times == 0 || gui.incrase)) {
                        gui.solver->Compute();
                        gui.data_ready = true;
                        gui.times += 1;
                    }
                }
                catch (std::exception &e) {
                    std::string str;
                    str = e.what();
                    MessageBoxA(hWnd, str.c_str(), "Error", MB_OK);
                }

                SendMessage(hWnd, WM_COMMAND, MAKELONG(0, BN_CLICKED), (LPARAM)hRedraw);
            } else if ((HWND)lParam == hRedraw && HIWORD(wParam) == BN_CLICKED) {

                if (gui.solver.get() && gui.data_ready) {

                    gui.data.v = ToEigen(gui.solver->GetV());
                    auto psi = gui.solver->GetPsi();
                    auto psi0 = gui.solver->GetPsi0();
                    bool add = SendMessage(hAddInitPsi, BM_GETCHECK, 0, 0) == BST_CHECKED;

                    gui.data.psiToDraw.resize(gui.data.v.rows(), gui.data.v.cols());
                    for (size_t i = 0; i < gui.data.v.rows(); ++i) {
                        for (size_t j = 0; j < gui.data.v.cols(); ++j) {
                            if (add)
                                gui.data.psiToDraw(i, j) = psi(i, j) + psi0(i, j);
                            else
                                gui.data.psiToDraw(i, j) = psi(i, j);
                        }
                    }
                }
                InvalidateRect(hWnd, NULL, false);
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
