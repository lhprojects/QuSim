#define _CRT_SECURE_NO_WARNINGS
#include "header.h"
#include "QuSim2D.h"
#include <thread>
#include <mutex>

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


// user message
#define WM_DATAREADY (WM_USER + 12)


// windows
HWND hUpdate;
HWND hPotential;
HWND hInitalPsi;
HWND hNormInitPsi;
HWND hNormPsiEachStep;
HWND hShowPsi;
HWND hShowPotential;
HWND hX0;
HWND hX1;
HWND hXBins;
HWND hY0;
HWND hY1;
HWND hYBins;
HWND hStop;
HWND hRun;
HWND hBoundaryCondition;
HWND hCanvas;
HWND hMainWin;
HWND hDeltaT;
HWND hPeriod;
HWND hHbar;
HWND hMass;
HWND hEigen;
HWND hVTV;
HWND hSplitO4;
HWND hMidpoint;
HWND hGaussO4;
HWND hGaussO6;

std::vector<HWND> enableControlList;
void enableAllWindows(bool enable)
{
    for (auto hwnd : enableControlList) {
        EnableWindow(hwnd, enable);
    }
}


enum class WhatToDo {
    Pause,
    InitAndRun,
    Run,
};

struct Messager {
    std::condition_variable cond;
    std::mutex mutex;
    WhatToDo what; // Pause by default
    WhatToDo GetWhat()
    {
        std::lock_guard<std::mutex> lk(mutex);
        return what;
    }
    void SetWhat(WhatToDo w)
    {
        std::lock_guard<std::mutex> lk(mutex);
        what = w;
    }

    WhatToDo GetWhatWait()
    {
        std::unique_lock<std::mutex> lk(mutex);
        cond.wait(lk);
        return what;
    }

    void GetWhatNotify()
    {
        cond.notify_all();
    }
};
Messager messager;

struct DurationEstimate
{
    using duration_t = std::chrono::high_resolution_clock::duration;

    bool init = false;
    double fDuration = 0;
    std::chrono::high_resolution_clock::time_point fLast;

    double duration()
    {
        auto last = fLast;
        fLast = std::chrono::high_resolution_clock::now();
        double dur = 0;
        if (!init) {
            init = true;
        } else {
            dur = std::chrono::duration_cast<std::chrono::milliseconds>(fLast - last).count();
        }
        fDuration = fDuration * 0.9 + dur * 0.1;
        return fDuration;
    }

};

struct Data {
    Eigen::MatrixXcd psi;
    Real norm;
    Real kin;
    Real pot;
    Real time;
    Eigen::MatrixXd v;
    Int fNx;
    Int fNy;
    Real FPS;
};

struct WorkerToGUI {
    // data shard package 1
    Data data;
    std::mutex mutex;
};

struct GUIToWorker {
    // data shard package 2
    std::shared_ptr<Evolver2D> syst;
    std::mutex mutex;
};


int const STATE_RUNNING = 2;
int const STATE_STOPPED = 0;
int const STATE_PAUSED = 1;

struct GUI {
    // gui thread
    bool data_ready;
    Data data;
    int runningState;
    bool showRunning;
};

struct Worker {
    // woker thread
    std::shared_ptr<Evolver2D> syst;
    DurationEstimate fDuration;

};

// the worker may access these data
// we will never free the these memory
WorkerToGUI *workerToGUI = new WorkerToGUI();
GUIToWorker *guiToWorker = new GUIToWorker();
Worker *worker = new Worker();

GUI gui;

void simulate()
{

    WhatToDo what = messager.GetWhat();

    for (;;) {
        if (what == WhatToDo::InitAndRun) {
            messager.SetWhat(WhatToDo::Run);

            worker->syst = guiToWorker->syst;
            {
                std::lock_guard<std::mutex> lk(workerToGUI->mutex);
                workerToGUI->data.psi = ToEigen(worker->syst->GetPsi());
                workerToGUI->data.norm = worker->syst->Norm2();
                workerToGUI->data.time = worker->syst->Time();
                workerToGUI->data.fNx = worker->syst->GetNx();
                workerToGUI->data.fNy = worker->syst->GetNy();
                workerToGUI->data.v = ToEigen(worker->syst->GetV());
                workerToGUI->data.pot = worker->syst->PotEn();
                workerToGUI->data.kin = worker->syst->KinEn();
                workerToGUI->data.FPS = 0;
                worker->fDuration.duration();
            }
            PostMessage(hMainWin, WM_DATAREADY, 0, 0);
            what = messager.GetWhat();
        } else if (what == WhatToDo::Run) {
            //double x = worker->syst.Norm2();
            worker->syst->step();
            {
                std::lock_guard<std::mutex> lk(workerToGUI->mutex);
                workerToGUI->data.psi = ToEigen(worker->syst->GetPsi());
                workerToGUI->data.norm = worker->syst->Norm2();
                workerToGUI->data.time = worker->syst->Time();
                workerToGUI->data.fNx = worker->syst->GetNx();
                workerToGUI->data.fNy = worker->syst->GetNy();
                workerToGUI->data.v = ToEigen(worker->syst->GetV());
                workerToGUI->data.pot = worker->syst->PotEn();
                workerToGUI->data.kin = worker->syst->KinEn();
                workerToGUI->data.FPS = 1000/worker->fDuration.duration();
            }
            PostMessage(hMainWin, WM_DATAREADY, 0, 0);
            what = messager.GetWhat();

        } else if (what == WhatToDo::Pause) {
            what = messager.GetWhatWait();
        }
    }
}
std::thread simulator;


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

    simulator = std::move(std::thread(simulate));

    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);


    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_QUSIM2D, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);
    RegisterCanvasClass(hInstance);

    if (!InitInstance(hInstance, nCmdShow)) {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_QUSIM2D));

    MSG msg;

    while (GetMessage(&msg, nullptr, 0, 0)) {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg)) {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    GdiplusShutdown(gdiplusToken);

    auto native_handler = simulator.native_handle();
    simulator.detach();
    // we will terminate this thread manually
    TerminateThread(native_handler, 0);
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
    wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_QUSIM2D));
    wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground = (HBRUSH)(CTLCOLOR_STATIC);
    wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_QUSIM2D);
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
    wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_QUSIM2D));
    wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_QUSIM2D);;
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

    graphics.Clear(Color(255, 255, 255, 255));

    if (!gui.data_ready) return;

    Pen      blackPen(Color(255, 0, 0, 0));
    graphics.DrawLine(&blackPen, Point(0, h / 2), Point(w, h / 2));

    std::vector<double> v2(gui.data.v.data(), gui.data.v.data() + gui.data.v.size());
    DrawPotential2D(graphics, gui.data.v, left, top, w, h);
    DrawPsi2D(graphics, gui.data.psi, left, top, w, h);

    Font myFont(L"Courier New", 14);
    RectF layoutRect((float)left + 25, (float)top + 20, 500.0f, 800.0f);
    StringFormat format;
    format.SetAlignment(StringAlignmentNear);
    SolidBrush textBrush(Color(255, 0, 255, 0));
    static const int len = 300;
    wchar_t text[len];
    wchar_t *b = text;


    double fps = gui.data.FPS;
    b += swprintf(b, text + len - b, L"FPS        %5.1f\n", fps);
    b += swprintf(b, text + len - b, L"time       % .15g\n", gui.data.time);
    b += swprintf(b, text + len - b, L"norm       % .16g\n", gui.data.norm);
    b += swprintf(b, text + len - b, L"kin        % .15g\n", gui.data.kin);
    b += swprintf(b, text + len - b, L"pot        % .15g\n", gui.data.pot);
    b += swprintf(b, text + len - b, L"energy     % .15g\n", gui.data.pot + gui.data.kin);

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


void InitialASystem2D(Evolver2D &syst)
{
    std::string psi;
    std::string pot;
    double x0 = -1;
    double x1 = 1;
    double y0 = -1;
    double y1 = 1;
    double mass = 1;
    double hbar = 1;
    int nx = 100;
    int ny = 100;
    bool fn = false;
    bool fnes = false;
    bool splito4 = false;
    Complex deltaT;
    BoundaryCondition bc;
    SolverMethod sl;

    fn = SendMessage(hNormInitPsi, BM_GETCHECK, 0, 0) == BST_CHECKED;
    fnes = SendMessage(hNormPsiEachStep, BM_GETCHECK, 0, 0) == BST_CHECKED;

    if (SendMessage(hSplitO4, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        sl = SolverMethod::SplittingMethodO4;
    } else if (SendMessage(hMidpoint, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        sl = SolverMethod::ImplicitMidpointMethod;
    } else if (SendMessage(hGaussO4, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        sl = SolverMethod::GaussLegendreO4;
    } else if (SendMessage(hGaussO6, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        sl = SolverMethod::GaussLegendreO6;
    } else {
        sl = SolverMethod::SplittingMethodO2;
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
    GetInt(hXBins, nx, "# bins x");
    GetInt(hYBins, ny, "# bins y");

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hInitalPsi);
        psiStr.resize(len + 1);
        GetWindowTextA(hInitalPsi, psiStr.data(), len + 1);
        psi.assign(psiStr.data(), len);
        fn = SendMessage(hNormInitPsi, BM_GETCHECK, 0, 0) == BST_CHECKED;
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hPotential);
        psiStr.resize(len + 1);
        GetWindowTextA(hPotential, psiStr.data(), len + 1);
        pot.assign(psiStr.data(), len);
    }

    {
        std::vector<char> psiStr;
        int len = GetWindowTextLength(hDeltaT);
        psiStr.resize(len + 1);
        GetWindowTextA(hDeltaT, psiStr.data(), len + 1);
        Calculator cal(psiStr.data());
        deltaT = cal.Evaluate();
    }
    Options opts;
#ifdef QUSIM_USE_CUDA
    if (sl == SolverMethod::SplittingMethodO2 || sl == SolverMethod::SplittingMethodO4) {
        //opts.Cuda().Batch(10);
        opts.SetString("fft_lib", "FFTW");
    }
#endif
    syst.init(Functor2DWrapper(psi.c_str()), fn, deltaT, fnes, Functor2DWrapper(pot.c_str()),
        x0, x1, nx, y0, y1, ny, bc, sl, mass, hbar, opts);

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

        if (show_psi) {
            std::vector<char> psiStr;
            int len = GetWindowTextLength(hInitalPsi);
            psiStr.resize(len + 1);
            GetWindowTextA(hInitalPsi, psiStr.data(), len + 1);
            psi.assign(psiStr.data(), len);
            fn = SendMessage(hNormInitPsi, BM_GETCHECK, 0, 0) == BST_CHECKED;
            if (psi.empty()) show_psi = false;
        }
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
            vv.resize(nx, ny);

            Calculator cal(pot.c_str());
            cal.SetVaraible("x", 0);
            cal.SetVaraible("y", 0);

            Complex *x = &cal.GetVaraible("x");
            Complex *y = &cal.GetVaraible("y");
            for (int j = 0; j < nx; ++j) {
                *x = x0 + j * (x1 - x0) / nx;
                for (int i = 0; i < ny; ++i) {
                    *y = y0 + i * (y1 - y0) / ny;
                    vv(j, i) = cal.Evaluate().real();
                }
            }
            DrawPotential2D(graphics, vv, left, top, w, h);
        }
        if (show_psi) {
            psiv.resize(nx, ny);

            Calculator cal(psi.c_str());
            cal.SetVaraible("x", 0);
            cal.SetVaraible("y", 0);

            Complex *x = &cal.GetVaraible("x");
            Complex *y = &cal.GetVaraible("y");
            for (int j = 0; j < nx; ++j) {
                *x = x0 + j * (x1 - x0) / nx;
                for (int i = 0; i < ny; ++i) {
                    *y = y0 + i * (y1 - y0) / ny;
                    psiv(j, i) = cal.Evaluate();
                }
            }

            if (fn) {
                double k = 1 / sqrt(psiv.squaredNorm()*(x1 - x0) / nx * (y1 - y0) / ny);
                psiv *= k;
            }
            DrawPsi2D(graphics, psiv, left, top, w, h);
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
    if (true) {
        Bitmap bitmap(w, h);
        Gdiplus::Graphics graphics(&bitmap);

        graphics.TranslateTransform(0, 15);
        if (gui.showRunning) {
            OnPaint1(graphics, left, top, w, h - 30);
        } else {
            OnPaint2(graphics, left, top, w, h - 30);
        }

        Graphics gr(hdc);
        CachedBitmap cbitmap(&bitmap, &gr);
        gr.DrawCachedBitmap(&cbitmap, 0, 0);
    } else {

        Graphics graphics(hdc);
        graphics.TranslateTransform(0, 15);
        if (gui.showRunning) {
            OnPaint1(graphics, left, top, w, h - 30);
        } else {
            OnPaint2(graphics, left, top, w, h - 30);
        }
    }
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

        HWND hs2 = CreateWindow(
            TEXT("STATIC"), TEXT("Intial Wave Function"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 180, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hs2, WM_SETFONT, (LPARAM)guiFont, true);
        x += 180;

        hInitalPsi = CreateWindow(TEXT("EDIT"), TEXT("gauss(x, -50, 10)*gauss(y, 0, 30)*exp(5*I*x)"),
            WS_CHILD | WS_VISIBLE | WS_BORDER /*边框*/ | ES_AUTOHSCROLL /*水平滚动*/,
            x, y, 500, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hInitalPsi, WM_SETFONT, (LPARAM)guiFont, true);
        x += 500;
        enableControlList.push_back(hInitalPsi);

        hNormInitPsi = CreateWindow(
            TEXT("BUTTON"), TEXT("Normalize Initial Wave Function"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTOCHECKBOX,
            x, y, 260, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hNormInitPsi, WM_SETFONT, (LPARAM)guiFont, true);
        x += 260;
        enableControlList.push_back(hNormInitPsi);
        SendMessage(hNormInitPsi, BM_SETCHECK, BST_CHECKED, NULL);

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

        hUpdate = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Update"),
            WS_TABSTOP | WS_VISIBLE | WS_CHILD,
            x, y, 200, 60,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hUpdate, WM_SETFONT, (LPARAM)guiFont, true);
        x += 200;
        enableControlList.push_back(hUpdate);

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

        hPotential = CreateWindow(TEXT("EDIT"), TEXT("50*exp(-x*x/0.25)*(1-exp(-(y-5)*(y-5)))*(1-exp(-(y+5)*(y+5)))"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL ,
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

        hXBins = CreateWindow(TEXT("EDIT"), TEXT("1000"),
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

        hYBins = CreateWindow(TEXT("EDIT"), TEXT("500"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hYBins, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hYBins);

        ////////////////////////
        x = x0;
        y += 30;

        HWND static_dt = CreateWindow(
            TEXT("STATIC"), TEXT("\u0394t(Complex Allowed)"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 180, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(static_dt, WM_SETFONT, (LPARAM)guiFont, true);
        x += 180;

        hDeltaT = CreateWindow(TEXT("EDIT"), TEXT("0.01 + 0*I"),
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
            x, y, 500, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hDeltaT, WM_SETFONT, (LPARAM)guiFont, true);
        x += 500;
        enableControlList.push_back(hDeltaT);

        hNormPsiEachStep = CreateWindow(
            TEXT("BUTTON"), TEXT("Normalize Wave Function After Each Step"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTOCHECKBOX,
            x, y, 440, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hNormPsiEachStep, WM_SETFONT, (LPARAM)guiFont, true);
        x += 440;
        enableControlList.push_back(hNormPsiEachStep);

        ////////////////////////
        x = x0;
        y += h;

        hRun = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Run"),
            WS_TABSTOP | WS_VISIBLE | WS_CHILD,
            x , y , 100, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hRun, WM_SETFONT, (LPARAM)guiFont, true);
        x += 100;

        hStop = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Stop"),
            WS_TABSTOP | WS_VISIBLE | WS_CHILD,
            x , y , 100 , 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hStop, WM_SETFONT, (LPARAM)guiFont, true);
        EnableWindow(hStop, false);

        x += 100;

        HWND hs1 = CreateWindow(
            TEXT("STATIC"), TEXT("Solver:"),
            WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
            x, y, 80, 30,
            hWnd, (HMENU)NULL, hInst, NULL
        );
        SendMessage(hs1, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;


        hVTV = CreateWindow(
            TEXT("BUTTON"),
            TEXT("SplitO2"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON | WS_GROUP,
            x , y , 80, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hVTV, WM_SETFONT, (LPARAM)guiFont, true);
        SendMessage(hVTV, BM_SETCHECK, BST_CHECKED, NULL);
        x += 80;
        enableControlList.push_back(hVTV);

        hSplitO4 = CreateWindow(
            TEXT("BUTTON"),
            TEXT("SplitO4"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
            x, y, 80, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hSplitO4, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hSplitO4);

        hMidpoint = CreateWindow(
            TEXT("BUTTON"),
            TEXT("Midpoint"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
            x, y, 80, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hMidpoint, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hMidpoint);

        hGaussO4 = CreateWindow(
            TEXT("BUTTON"),
            TEXT("GaussO4"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
            x, y, 80, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hGaussO4, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hGaussO4);

        hGaussO6 = CreateWindow(
            TEXT("BUTTON"),
            TEXT("GaussO6"),
            WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
            x, y, 80, 30,
            hWnd, (HMENU)0, hInst, NULL
        );
        SendMessage(hGaussO6, WM_SETFONT, (LPARAM)guiFont, true);
        x += 80;
        enableControlList.push_back(hGaussO6);


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

            if ((HWND)lParam == hRun && HIWORD(wParam) == BN_CLICKED) {
                if (gui.runningState == STATE_STOPPED || gui.runningState == STATE_PAUSED) {
                    try {
                        if (gui.runningState == STATE_STOPPED) {
                            std::shared_ptr<Evolver2D> syst(new Evolver2D);
                            InitialASystem2D(*syst);
                            std::lock_guard<std::mutex> lk(guiToWorker->mutex);
                            guiToWorker->syst = syst;
                            messager.SetWhat(WhatToDo::InitAndRun);
                            messager.GetWhatNotify();

                        } else {
                            messager.SetWhat(WhatToDo::Run);
                            messager.GetWhatNotify();
                        }
                        SendMessage(hRun, WM_SETTEXT, 0, (LPARAM)TEXT("Pause"));
                        if (gui.runningState == STATE_STOPPED) {
                            gui.showRunning = true;
                            enableAllWindows(false);
                            EnableWindow(hStop, true);
                            UpdateWindow(hWnd);
                        }
                        gui.runningState = STATE_RUNNING;

                    } catch (std::exception const &exp) {
                        wchar_t w[100];
                        mbstowcs(w, exp.what(), 100);
                        MessageBox(hWnd, w, L"Intialization Error", MB_OK);
                    } catch (...) {
                        MessageBox(hWnd, L"err", L"err", MB_OK);
                    }

                } else if (gui.runningState == STATE_RUNNING) {
                    messager.SetWhat(WhatToDo::Pause);
                    SendMessage(hRun, WM_SETTEXT, 0, (LPARAM)TEXT("Run"));
                    gui.runningState = STATE_PAUSED;
                }
            } else if ((HWND)lParam == hStop && HIWORD(wParam) == BN_CLICKED) {
                if (gui.runningState == STATE_RUNNING || gui.runningState == STATE_PAUSED) {

                    messager.SetWhat(WhatToDo::Pause);
                    SendMessage(hRun, WM_SETTEXT, 0, (LPARAM)TEXT("Run"));
                    gui.runningState = STATE_STOPPED;

                    // don't set showRunning
                    enableAllWindows(true);
                    EnableWindow(hStop, false);

                }
            } else if ((HWND)lParam == hUpdate && HIWORD(wParam) == BN_CLICKED) {
                if (gui.runningState == STATE_STOPPED) {
                    gui.showRunning = false;
                    InvalidateRect(hWnd, NULL, false);
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
    case WM_DATAREADY:
    {
        {
            std::lock_guard<std::mutex> lc(workerToGUI->mutex);
            gui.data = workerToGUI->data;
            gui.data_ready = true;
        }
        InvalidateRect(hWnd, NULL, false);
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
