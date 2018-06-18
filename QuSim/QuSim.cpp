//QuSim.cpp: 定义应用程序的入口点。
//

#include "header.h"
#include "QuSim.h"
#include <thread>
#include <mutex>

// core
#include "../qsim/System.h"


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
HWND hStop;
int stopState;

HWND hMainWin;


enum class WhatToDo {
	Pause,
	Run,
	Init,
	Abort,
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

struct GUI {
	// gui thread
	bool data_ready;
	Data data;
};

struct Worker {
	// woker thread
	System syst;

};

// the worker may access these data
// we will never free the these memory
WorkerToGUI *workerToGUI = new WorkerToGUI();
Worker *worker = new Worker();

GUI gui;

void simulate()
{

	WhatToDo what = messager.GetWhat();

	for (;;) {
		if (what == WhatToDo::Init) {

		} else if (what == WhatToDo::Run) {
			if (worker->syst.fN == 0) {
				int n = 1;
				int m = 1;
				worker->syst.init("gauss(x, -16, 4)*exp(I*x)", true, 1E-3,
					"(x < 0.500001 && x > -0.499999)? 1 : 0", -50*n, 50*n , 2000*n*m,
					BoundaryCondition::InfiniteWall,
					1, 1);

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
			}
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
		} else if (what == WhatToDo::Abort) {
			break;
		}
	}
}
std::thread simulator;


#define MAX_LOADSTRING 100

// 全局变量: 
HINSTANCE hInst;                                // 当前实例
WCHAR szTitle[MAX_LOADSTRING];                  // 标题栏文本
WCHAR szWindowClass[MAX_LOADSTRING];            // 主窗口类名

// 此代码模块中包含的函数的前向声明: 
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
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
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_QUSIM);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

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

   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN,
      CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

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

void OnPaint1(Gdiplus::Graphics &graphics, long left, long top, long w, long h)
{

	graphics.Clear(Color(255, 255, 255, 255));

	if (!gui.data_ready) return;

	Font myFont(L"Courier New", 14);
	RectF layoutRect(left, top, 500.0f, 800.0f);
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
	b += swprintf(b, text + len - b, L"energy(∆T) % .15g\n", gui.data.enPartialTime);

	graphics.DrawString(text, lstrlenW(text), &myFont, layoutRect, &format, &blackBrush);

	Pen      blackPen(Color(255, 0, 0, 0));

	graphics.DrawLine(&blackPen, Point(0, h / 2), Point(w, h / 2));

	std::vector<Point> abc;
	std::vector<Point> rec;
	std::vector<Point> imc;
	std::vector<Point> vc;

	auto it = std::minmax_element(gui.data.v.begin(), gui.data.v.end());
	Real vmin = *it.first;
	Real vmax = *it.second;
	vmax = 1.5*max(abs(vmin), abs(vmax));

	for (int i = 0; i < w; i += 1) {
		int xi = (int)(1.0 * i / w * (gui.data.fN - 1));
		int ab = -abs(gui.data.psi[xi]) * 100 + h / 2;
		int re = -gui.data.psi[xi].real() * 100 + h / 2;
		int im = -gui.data.psi[xi].imag() * 100 + h / 2;
		int v = 0;
		if (vmax != 0) {
			v = -(gui.data.v[xi] / vmax) * h / 2 + h / 2;
		}
		abc.push_back(Point(i, ab));
		rec.push_back(Point(i, re));
		imc.push_back(Point(i, im));
		vc.push_back(Point(i, v));
	}

	Pen      pen1(Color(255, 0, 0, 255));
	Pen      pen2(Color(255, 0, 255, 255));
	Pen      pen3(Color(255, 255, 0, 255));
	Pen      pen4(Color(255, 0, 255, 0));
	graphics.DrawCurve(&pen1, abc.data(), abc.size());
	graphics.DrawCurve(&pen2, rec.data(), rec.size());
	graphics.DrawCurve(&pen3, imc.data(), imc.size());
	graphics.DrawLines(&pen4, vc.data(), vc.size());



}

VOID OnPaint(HDC hdc, long left, long top, long w, long h)
{
	Bitmap bitmap(w, h);
	Gdiplus::Graphics graphics(&bitmap);

	OnPaint1(graphics, left, top + 100, w, h);

	Graphics gr(hdc);
	CachedBitmap cbitmap(&bitmap, &gr);
	gr.DrawCachedBitmap(&cbitmap, 0, 0);
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
		NONCLIENTMETRICS metrics = {};
		metrics.cbSize = sizeof(metrics);
		SystemParametersInfo(SPI_GETNONCLIENTMETRICS, metrics.cbSize, &metrics, 0);

		HFONT guiFont = CreateFontIndirect(&metrics.lfCaptionFont);


		hStop = CreateWindow(
			TEXT("BUTTON"),
			TEXT("Run"),
			WS_TABSTOP | WS_VISIBLE | WS_CHILD ,
			30 /*X坐标*/, 20 /*Y坐标*/, 100 /*宽度*/, 30/*高度*/,
			hWnd, (HMENU)0, hInst, NULL
		);
		SendMessage(hStop, WM_SETFONT, (LPARAM)guiFont, true);

		// When you're done with the font, don't forget to call
		//DeleteObject(guiFont);
	}
	break;
	case WM_COMMAND:
	{

		if (lParam != 0) { // from button

			if ((HWND)lParam == hStop && HIWORD(wParam) == BN_CLICKED) {
				if (stopState == 0) { // paused
					messager.SetWhat(WhatToDo::Run);
					messager.GetWhatNotify();
					SendMessage(hStop, WM_SETTEXT, 0, (LPARAM)TEXT("Pause"));
					stopState = 1;

				} else if (stopState == 1) { // running
					messager.SetWhat(WhatToDo::Pause);
					SendMessage(hStop, WM_SETTEXT, 0, (LPARAM)TEXT("Run"));
					stopState = 0;
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
    case WM_PAINT:
        {
		
			RECT rect;
			GetClientRect( hWnd, &rect);

            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
			OnPaint(hdc, rect.left, rect.top, rect.right - rect.left, rect.bottom - rect.top);
            EndPaint(hWnd, &ps);
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
