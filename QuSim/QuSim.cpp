#define _CRT_SECURE_NO_WARNINGS
#include "header.h"
#include "QuSim.h"
#include <thread>
#include <mutex>

// core
#include "../qsim/System.h"
#include "../qsim/Cal.h"


// windows theme
#pragma comment(linker,"\"/manifestdependency:type='win32' \
name='Microsoft.Windows.Common-Controls' version='6.0.0.0' \
processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

// gdi plus
#pragma comment(lib, "gdiplus.lib")  
using namespace Gdiplus;

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
HWND hBins;
HWND hStop;
HWND hRun;
HWND hBoundaryCondition;
HWND hCanvas;
HWND hMainWin;
HWND hDeltaT;
HWND hInfiniteWall;
HWND hPeriod;
HWND hHbar;
HWND hMass;
HWND hExtend;

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

struct Data {
	PsiVector psi;
	Real norm;
	Real kin;
	Real pot;
	Real normLeft;
	Real normRight;
	Real xavg;
	Real time;
	Real enPartialTime;
	std::vector<Real> v;
	Int fN;
};

struct WorkerToGUI {
	// data shard package 1
	Data data;
	std::mutex mutex;
};

struct GUIToWorker {
	// data shard package 2
	System syst;
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
	System syst;

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
				workerToGUI->data.psi = worker->syst.fPsi;
				workerToGUI->data.norm = worker->syst.Norm2();
				workerToGUI->data.xavg = worker->syst.Xavg();
				workerToGUI->data.time = worker->syst.Time();
				workerToGUI->data.fN = worker->syst.fN;
				workerToGUI->data.v = worker->syst.fV;
				workerToGUI->data.normLeft = worker->syst.NormLeft();
				workerToGUI->data.normRight = worker->syst.NormRight();
				workerToGUI->data.pot = worker->syst.PotEn();
				workerToGUI->data.kin = worker->syst.KinEn();
				workerToGUI->data.enPartialTime = worker->syst.EnPartialT();
			}
			PostMessage(hMainWin, WM_DATAREADY, 0, 0);
			what = messager.GetWhat();
		} else if (what == WhatToDo::Run) {

			worker->syst.step();
			{
				std::lock_guard<std::mutex> lk(workerToGUI->mutex);
				workerToGUI->data.psi = worker->syst.fPsi;
				workerToGUI->data.norm = worker->syst.Norm2();
				workerToGUI->data.xavg = worker->syst.Xavg();
				workerToGUI->data.time = worker->syst.Time();
				workerToGUI->data.fN = worker->syst.fN;
				workerToGUI->data.v = worker->syst.fV;
				workerToGUI->data.normLeft = worker->syst.NormLeft();
				workerToGUI->data.normRight = worker->syst.NormRight();
				workerToGUI->data.pot = worker->syst.PotEn();
				workerToGUI->data.kin = worker->syst.KinEn();
				workerToGUI->data.enPartialTime = worker->syst.EnPartialT();
			}
			PostMessage(hMainWin, WM_DATAREADY, 0, 0);
			what = messager.GetWhat();

		} else if(what == WhatToDo::Pause) {
			what = messager.GetWhatWait();
		}
	}
}
std::thread simulator;


#define MAX_LOADSTRING 100

// 全局变量: 
HINSTANCE hInst;                                // 当前实例
WCHAR szTitle[MAX_LOADSTRING];                  // 标题栏文本
WCHAR szWindowClass[MAX_LOADSTRING];
WCHAR szCanvasWindowClass[] = { L"Canvas" };

// 此代码模块中包含的函数的前向声明: 
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

	// 在第一次使用GDI + 对象前，调用以下代码：
	ULONG_PTR gdiplusToken; // 这个变量需要保存下来
	GdiplusStartupInput gdiplusStartupInput;
	GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

	simulator = std::move(std::thread(simulate));

    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: 在此放置代码。

    // 初始化全局字符串
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_QUSIM, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);
	RegisterCanvasClass(hInstance);

    // 执行应用程序初始化: 
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_QUSIM));

    MSG msg;

    // 主消息循环: 
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

	// 在最后一次使用GDI+对象之后，调用以下代码：
	GdiplusShutdown(gdiplusToken);

	auto native_handler = simulator.native_handle();
	simulator.detach();
	// we will terminate this thread manually
	TerminateThread(native_handler, 0);
	return (int) msg.wParam;
}



//
//  函数: MyRegisterClass()
//
//  目的: 注册窗口类。
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_QUSIM));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(CTLCOLOR_STATIC);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_QUSIM);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}


//
//  函数: MyRegisterClass()
//
//  目的: 注册窗口类。
//
ATOM RegisterCanvasClass(HINSTANCE hInstance)
{
	WNDCLASSEXW wcex;

	wcex.cbSize = sizeof(WNDCLASSEX);

	wcex.style = CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc = CanvasWndProc;
	wcex.cbClsExtra = 0;
	wcex.cbWndExtra = 0;
	wcex.hInstance = hInstance;
	wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_QUSIM));
	wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
	wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
	wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_QUSIM);;
	wcex.lpszClassName = szCanvasWindowClass;
	wcex.hIconSm = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

	return RegisterClassExW(&wcex);
}

//
//   函数: InitInstance(HINSTANCE, int)
//
//   目的: 保存实例句柄并创建主窗口
//
//   注释: 
//
//        在此函数中，我们在全局变量中保存实例句柄并
//        创建和显示主程序窗口。
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // 将实例句柄存储在全局变量中

#if 0
   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN,
      CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);
#else
   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
	   CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);
#endif

   hMainWin = hWnd;

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}

#include <algorithm>

void DrawPotential(Gdiplus::Graphics &graphics,
	std::vector<Real> const &V,
	long left, long top, long w, long h)
{
	std::vector<Point> vc;

	auto it = std::minmax_element(V.begin(), V.end());
	Real vmin = *it.first;
	Real vmax = *it.second;
	vmax = 1.5*max(abs(vmin), abs(vmax));

	for (int i = 0; i < w; i += 1) {
		int xi = (int)(1.0 * i / w * (V.size() - 1));
		int v = 0;
		if (vmax != 0) {
			v = (int)(-V[xi] / vmax * h / 2 + h / 2);
		}
		vc.push_back(Point(i, v));
	}

	Pen      pen4(Color(255, 0, 255, 0));
	graphics.DrawLines(&pen4, vc.data(), vc.size());

}

void DrawPsi(Gdiplus::Graphics &graphics,
	PsiVector const &psi,
	long left, long top, long w, long h)
{
	std::vector<Point> abc;
	std::vector<Point> rec;
	std::vector<Point> imc;
	std::vector<Point> vc;


	for (int i = 0; i < w; i += 1) {
		int xi = (int)(1.0 * i / w * (psi.size() - 1));
		int ab = -abs(psi[xi]) * 100 + h / 2;
		int re = -psi[xi].real() * 100 + h / 2;
		int im = -psi[xi].imag() * 100 + h / 2;
		abc.push_back(Point(i, ab));
		rec.push_back(Point(i, re));
		imc.push_back(Point(i, im));
	}

	Pen      pen1(Color(255, 0, 0, 255));
	Pen      pen2(Color(255, 0, 255, 255));
	Pen      pen3(Color(255, 255, 0, 255));
	graphics.DrawCurve(&pen1, abc.data(), abc.size());
	graphics.DrawCurve(&pen2, rec.data(), rec.size());
	graphics.DrawCurve(&pen3, imc.data(), imc.size());

}

void OnPaint1(Gdiplus::Graphics &graphics, long left, long top, long w, long h)
{

	graphics.Clear(Color(255, 255, 255, 255));

	if (!gui.data_ready) return;

	Font myFont(L"Courier New", 14);
	RectF layoutRect(left+25, top+20, 500.0f, 800.0f);
	StringFormat format;
	format.SetAlignment(StringAlignmentNear);
	SolidBrush blackBrush(Color(255, 0, 0, 0));
	static const int len = 300;
	wchar_t text[len];
	wchar_t *b = text;

	b += swprintf(b, text + len - b, L"time       % .15g\n", gui.data.time);
	b += swprintf(b, text + len - b, L"xavg       % .16g\n", gui.data.xavg);
	b += swprintf(b, text + len - b, L"norm       % .16g\n", gui.data.norm);
	b += swprintf(b, text + len - b, L"norm left  % .16g\n", gui.data.normLeft);
	b += swprintf(b, text + len - b, L"norm right % .15g\n", gui.data.normRight);
	b += swprintf(b, text + len - b, L"kin        % .15g\n", gui.data.kin);
	b += swprintf(b, text + len - b, L"pot        % .15g\n", gui.data.pot);
	b += swprintf(b, text + len - b, L"energy     % .15g\n", gui.data.pot + gui.data.kin);
	b += swprintf(b, text + len - b, L"energy(\u0394T) % .15g\n", gui.data.enPartialTime);

	graphics.DrawString(text, lstrlenW(text), &myFont, layoutRect, &format, &blackBrush);

	Pen      blackPen(Color(255, 0, 0, 0));

	graphics.DrawLine(&blackPen, Point(0, h / 2), Point(w, h / 2));

	DrawPotential(graphics, gui.data.v, left, top, w, h);
	DrawPsi(graphics, gui.data.psi, left, top, w, h);
}

void InitialASystem(System &syst)
{
	std::string psi;
	std::string pot;
	double x0 = -1;
	double x1 = 1;
	double mass = 1;
	double hbar = 1;
	int n = 100;
	bool fn = false;
	bool fnes = false;
	Complex deltaT;
	BoundaryCondition bc;

	fn = SendMessage(hNormInitPsi, BM_GETCHECK, 0, 0) == BST_CHECKED;
	fnes = SendMessage(hNormPsiEachStep, BM_GETCHECK, 0, 0) == BST_CHECKED;

	if (SendMessage(hInfiniteWall, BM_GETCHECK, 0, 0) == BST_CHECKED) {
		bc = BoundaryCondition::InfiniteWall;
	} else if(SendMessage(hPeriod, BM_GETCHECK, 0, 0) == BST_CHECKED){
		bc = BoundaryCondition::Period;
	} else {
		bc = BoundaryCondition::ExtendInfinity;
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

	{
		std::vector<char> psiStr;
		int len = GetWindowTextLength(hHbar);
		psiStr.resize(len + 1);
		GetWindowTextA(hHbar, psiStr.data(), len + 1);
		int na = sscanf(psiStr.data(), "%lf", &hbar);
		if (na < 1) {
			throw std::runtime_error("can't parse x0");
		}
	}

	{
		std::vector<char> psiStr;
		int len = GetWindowTextLength(hInitalPsi);
		psiStr.resize(len + 1);
		GetWindowTextA(hX0, psiStr.data(), len + 1);
		int na = sscanf(psiStr.data(), "%lf", &x0);
		if (na < 1) {
			throw std::runtime_error("can't parse x0");
		}
	}

	{
		std::vector<char> psiStr;
		int len = GetWindowTextLength(hInitalPsi);
		psiStr.resize(len + 1);
		GetWindowTextA(hX1, psiStr.data(), len + 1);
		int na = sscanf(psiStr.data(), "%lf", &x1);
		if (na < 1) {
			throw std::runtime_error("can't parse x1");
		}
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
		Cal cal(psiStr.data());
		deltaT = cal.Val();
	}

	syst.init(psi.c_str(), fn, deltaT, fnes, pot.c_str(),
		x0, x1, n, bc, mass, hbar);

}

void OnPaint2(Gdiplus::Graphics &graphics, long left, long top, long w, long h)
{

	graphics.Clear(Color(255, 255, 255, 255));

	Pen      blackPen(Color(255, 0, 0, 0));
	graphics.DrawLine(&blackPen, Point(0, h / 2), Point(w, h / 2));

	System syst;

	bool show_psi = SendMessage(hShowPsi, BM_GETCHECK, 0, 0) == BST_CHECKED;
	bool show_pot = SendMessage(hShowPotential, BM_GETCHECK, 0, 0) == BST_CHECKED;

	std::string psi;
	std::string pot;
	double x0 = -1;
	double x1 = 1;
	int n = 100;
	bool fn = false;

	{
		std::vector<char> psiStr;
		int len = GetWindowTextLength(hInitalPsi);
		psiStr.resize(len + 1);
		GetWindowTextA(hX0, psiStr.data(), len + 1);
		int na = sscanf(psiStr.data(), "%lf", &x0);
		if (na < 1) {
			show_psi = false;
			show_pot = false;
		}
	}

	{
		std::vector<char> psiStr;
		int len = GetWindowTextLength(hInitalPsi);
		psiStr.resize(len + 1);
		GetWindowTextA(hX1, psiStr.data(), len + 1);
		int na = sscanf(psiStr.data(), "%lf", &x1);
		if (na < 1) {
			show_psi = false;
			show_pot = false;
		}
	}

	{
		std::vector<char> psiStr;
		int len = GetWindowTextLength(hBins);
		psiStr.resize(len + 1);
		GetWindowTextA(hBins, psiStr.data(), len + 1);
		int na = sscanf(psiStr.data(), "%d", &n);
		if (na < 1) {
			show_psi = false;
			show_pot = false;
		}
	}

	if(show_psi) {
		std::vector<char> psiStr;
		int len = GetWindowTextLength(hInitalPsi);
		psiStr.resize(len + 1);
		GetWindowTextA(hInitalPsi, psiStr.data(),len+1);
		psi.assign(psiStr.data(), len);
		fn = SendMessage(hNormInitPsi, BM_GETCHECK, 0, 0) == BST_CHECKED;
		if (psi.empty()) show_psi = false;
	}
	
	if(!show_psi){
		psi = "0";
	}

	if (show_pot) {
		std::vector<char> psiStr;
		int len = GetWindowTextLength(hPotential);
		psiStr.resize(len + 1);
		GetWindowTextA(hPotential, psiStr.data(), len + 1);
		pot.assign(psiStr.data(), len);
		if (pot.empty()) show_pot = false;
	}
	if(!show_pot) {
		pot = "0";
	}

	try {
		syst.init(psi.c_str(), fn, 1, false, pot.c_str(),
			x0, x1, n, BoundaryCondition::Period, 1);
	} catch (...) {
		show_psi = false;
		show_pot = false;
	}
	if (show_psi) {
		DrawPsi(graphics, syst.fPsi, left, top, w, h);
	}
	if (show_pot) {
		DrawPotential(graphics, syst.fV, left, top, w, h);
	}

}
VOID OnPaint(HDC hdc, long left, long top, long w, long h)
{
	Bitmap bitmap(w, h);
	Gdiplus::Graphics graphics(&bitmap);

	if (gui.showRunning) {
		OnPaint1(graphics, left, top, w, h);
	} else {
		OnPaint2(graphics, left, top, w, h);
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

//
//  函数: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  目的:    处理主窗口的消息。
//
//  WM_COMMAND  - 处理应用程序菜单
//  WM_PAINT    - 绘制主窗口
//  WM_DESTROY  - 发送退出消息并返回
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message)
    {
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

		hInitalPsi = CreateWindow(TEXT("EDIT"), TEXT("gauss(x, 0, 2)*exp(I*x)"),
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

		hPotential = CreateWindow(TEXT("EDIT"), TEXT("0.5*x*x"),
			WS_CHILD | WS_VISIBLE | WS_BORDER /*边框*/ | ES_AUTOHSCROLL /*水平滚动*/,
			x, y, 500, 30,
			hWnd, (HMENU)NULL, hInst, NULL
		);
		SendMessage(hPotential, WM_SETFONT, (LPARAM)guiFont, true);
		x += 500;
		enableControlList.push_back(hPotential);

		HWND static_x0 = CreateWindow(
			TEXT("STATIC"), TEXT("x0"),
			WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
			x, y, 50, 30,
			hWnd, (HMENU)NULL, hInst, NULL
		);
		SendMessage(static_x0, WM_SETFONT, (LPARAM)guiFont, true);
		x += 50;

		hX0 = CreateWindow(TEXT("EDIT"), TEXT("-10"),
			WS_CHILD | WS_VISIBLE | WS_BORDER /*边框*/ | ES_AUTOHSCROLL,
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
			WS_CHILD | WS_VISIBLE | WS_BORDER /*边框*/ | ES_AUTOHSCROLL,
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
			x /*X坐标*/, y /*Y坐标*/, 100 /*宽度*/, 30/*高度*/,
			hWnd, (HMENU)0, hInst, NULL
		);
		SendMessage(hRun, WM_SETFONT, (LPARAM)guiFont, true);
		x += 100;

		hStop = CreateWindow(
			TEXT("BUTTON"),
			TEXT("Stop"),
			WS_TABSTOP | WS_VISIBLE | WS_CHILD,
			x /*X坐标*/, y /*Y坐标*/, 100 /*宽度*/, 30/*高度*/,
			hWnd, (HMENU)0, hInst, NULL
		);
		SendMessage(hStop, WM_SETFONT, (LPARAM)guiFont, true);
		EnableWindow(hStop, false);

		x += 100;

		HWND hs1 = CreateWindow(
			TEXT("STATIC"), TEXT("Boundary Condition:"),
			WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
			x, y, 180, 30,
			hWnd, (HMENU)NULL, hInst, NULL
		);
		SendMessage(hs1, WM_SETFONT, (LPARAM)guiFont, true);
		x += 180;

		hInfiniteWall = CreateWindow(
			TEXT("BUTTON"),
			TEXT("Infinite Wall"),
			WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
			x /*X坐标*/, y /*Y坐标*/, 120/*宽度*/, 30/*高度*/,
			hWnd, (HMENU)0, hInst, NULL
		);
		SendMessage(hInfiniteWall, WM_SETFONT, (LPARAM)guiFont, true);
		SendMessage(hInfiniteWall, BM_SETCHECK, BST_CHECKED, NULL);
		x += 120;
		enableControlList.push_back(hInfiniteWall);

		hPeriod = CreateWindow(
			TEXT("BUTTON"),
			TEXT("Periodic"),
			WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
			x /*X坐标*/, y /*Y坐标*/, 100 /*宽度*/, 30/*高度*/,
			hWnd, (HMENU)0, hInst, NULL
		);
		SendMessage(hPeriod, WM_SETFONT, (LPARAM)guiFont, true);
		x += 100;
		enableControlList.push_back(hPeriod);

		hExtend = CreateWindow(
			TEXT("BUTTON"),
			TEXT("Extend Infinity"),
			WS_CHILD | WS_VISIBLE | BS_LEFT | BS_AUTORADIOBUTTON,
			x /*X坐标*/, y /*Y坐标*/, 160 /*宽度*/, 30/*高度*/,
			hWnd, (HMENU)0, hInst, NULL
		);
		SendMessage(hExtend, WM_SETFONT, (LPARAM)guiFont, true);
		x += 160;
		enableControlList.push_back(hExtend);

		HWND static_hbar = CreateWindow(
			TEXT("STATIC"), TEXT("hbar"),
			WS_CHILD | WS_VISIBLE | SS_CENTERIMAGE | SS_CENTER,
			x, y, 60, 30,
			hWnd, (HMENU)NULL, hInst, NULL
		);
		SendMessage(static_hbar, WM_SETFONT, (LPARAM)guiFont, true);
		x += 60;

		hHbar = CreateWindow(TEXT("EDIT"), TEXT("1"),
			WS_CHILD | WS_VISIBLE | WS_BORDER /*边框*/ | ES_AUTOHSCROLL,
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
			WS_CHILD | WS_VISIBLE | WS_BORDER /*边框*/ | ES_AUTOHSCROLL,
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

		SetWindowPos(hCanvas, NULL, 0, 150,
			client_rect.right, client_rect.bottom - 50, 0);
	}
	break;
	case WM_COMMAND:
	{

		if (lParam != 0) { // from button

			if ((HWND)lParam == hRun && HIWORD(wParam) == BN_CLICKED) {
				if (gui.runningState == STATE_STOPPED || gui.runningState == STATE_PAUSED) {
					if (gui.runningState == STATE_STOPPED) {
						try {
							System syst;
							InitialASystem(syst);
							std::lock_guard<std::mutex> lk(guiToWorker->mutex);
							guiToWorker->syst = syst;
							messager.SetWhat(WhatToDo::InitAndRun);
							messager.GetWhatNotify();

						} catch (...) {
							MessageBox(hWnd, L"err", L"err", MB_OK);
						}
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
			// 分析菜单选择: 
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

// “关于”框的消息处理程序。
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
