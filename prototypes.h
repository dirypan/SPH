// prototypes.h
//

#pragma once

#ifndef _PROTOTYPES_H_

#define _PROTOTYPES_H_

#include <vector>

using namespace std;

extern int Nx = 40; //横向网格个数
extern int Ny = 40; //纵向网格个数
extern int N_x = Nx + 2;
extern int N_y = Nx + 2;
extern int N = N_x*N_y; //总网格个数
extern int N0 = 4; //每个网格单元内初态的粒子数
extern double h = 0.5; //网格宽度
extern double m = 5.0; //粒子质量
extern double Width = 2.0 * h * Nx; //视场宽度
extern double Height = 2.0 * h * Ny; //视场高度
extern double pumpForce = 50.0; //泵力系数
extern double gravity = 30.0; //重力场系数
extern double Gamma = 0.1; //Gruneisen 系数
extern double c0 = 150.0; //声速
extern double s = 0.3; //激波系数
extern double correctionFactor = -0.5; //位置修正系数
extern double alpha = 1.0; //人为粘性可调参数
extern double h_0 = 1.5; //人为粘性可调参数
extern double e = 0.1; //人为粘性可调参数
extern double a = 2.0; //人为粘性可调参数
extern double b = 2.0; //人为粘性可调参数
extern double g1 = 1.0; //人为热流可调参数
extern double g2 = 1.0; //人为热流可调参数
extern const double pi = 3.141592653589793238462643; //圆周率

class vec //矢量类
{
public:
	double x, y; //分量

	vec(const double &X = 0.0, const double &Y = 0.0); //构造函数

	double mold2(); //获取模方

	double mold(); //获取模

	void rotationAntiClock(); //逆时针旋转

	void rotationClock(); //顺时针旋转

	vec operator +(const vec &v); //加法

	vec operator +=(const vec &v); //加法

	vec operator -(const vec &v); //减法

	vec operator -=(const vec &v); //减法

	vec operator -(); //负矢量

	vec operator *(const double &k); //数乘

	vec operator *=(const double &k); //数乘

	vec operator /(const double &k); //数乘

	vec operator /=(const double &k); //数乘

	double operator *(const vec &v); //点乘

	bool operator ==(const vec &v); //判断相等

	bool operator !=(const vec &v); //判断不等
};

vec operator *(const double &k, const vec &v);//数乘

struct particle //粒子结构
{
	vec r, v, a;
	double u, pressure, density, density0;
};

double W(double rx, double ry); //样条核函数

double W(const vec &r); //样条核函数

double W(const vec &r1, const vec &r2); //样条核函数

vec gradW(double rx, double ry = 0.0); //样条核函数梯度

vec gradW(const vec &r); //样条核函数梯度

vec gradW1(const vec &r1, const vec &r2); //样条核函数梯度

vec gradW2(const vec &r1, const vec &r2); //样条核函数梯度

int periodicity(int nx, int ny); //周期性条件

void periodicity(particle &particlein); //周期性条件

void setWall(int nx1, int ny1, int nx2, int ny2); //设置墙

void setWall(int nx[], int ny[], int n, char mod = 'B'); //设置墙

void setPump(int nx1, int ny1, int nx2, int ny2, char mod = 'F', double pump_Force = pumpForce); //设置泵

void setPump(int nx1, int ny1, int nx2, int ny2, vec force); //设置泵

void setGravity(vec &direction, double sgravity = gravity); //设置重力场

void setGravity(double x, double y, double sgravity = gravity); //设置重力场

double viscosity(particle in1, particle in2); //人为粘性

double heatFactor(vector<particle> netin[], int nx, int ny, int np); //人为热流因子

double heat(vector<particle> netin[], int nx, int ny, int np); //人为热流

double pressure(double u, double density, double density0); //压强

void next(vector<particle> netin[], vector<particle> netout[], double &dt);

void divide(vector<particle> netin[], vector<particle> netout[]); //划分网格

void initialization(vector<particle> netin[]); //初始化

extern vec g = { 0.0,0.0 }; //重力场
extern vec v0 = { 0.0,0.0 }; //粒子初速度
extern vec *pump = {}; //泵
extern vec *wall = {}; //墙

#endif // !_PROTOTYPES_H_
