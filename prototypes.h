// prototypes.h
//

#pragma once

#ifndef _PROTOTYPES_H_

#define _PROTOTYPES_H_

#include <vector>

using namespace std;

extern const double h = 0.5;
extern const double m = 5;
extern const int Nx = 30; //横向网格个数
extern const int Ny = 30; //纵向网格个数
extern const int N_x = Nx + 2;
extern const int N_y = Nx + 2;
extern const int N = N_x*N_y; //总网格个数
extern const int N0 = 1; //每个网格单元内初态的粒子数
extern const double Width = 2.0 * h * Nx;
extern const double Height = 2.0 * h * Ny;
extern const double Gamma = 0.1; //Gruneisen 系数
extern const double c0 = 150.0; //声速
extern const double pi = 3.141592653589793238462643; //圆周率
extern const double s = 0.3;
extern const double correction_factor = -0.5; //位置修正系数
extern const double alpha = 1.0; //viscosity中的可调参数
extern const double h_0 = 1.5; //viscosity中的可调参数
extern const double e = 0.1;	 //e, a, b are the factors of viscosity;
extern const double a = 2;
extern const double b = 2.0;
extern const double g1 = 1.0; //g1, g2 are factors of fake heat flow
extern const double g2 = 1.0;
extern const double pumpForce = 50.0;
extern const double gravity = 50.0;

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

int periodicity(int nx, int ny);

void periodicity(particle &particlein);

void setWall(int nx1, int ny1, int nx2, int ny2);

void setWall(int nx[], int ny[], int n, char mod = 'B');

void setPump(int nx1, int ny1, int nx2, int ny2, char mod = 'F', double pump_Force = pumpForce);

void setPump(int nx1, int ny1, int nx2, int ny2, vec force);

void setGravity(vec &direction, double sgravity = gravity);

void setGravity(double x, double y, double sgravity = gravity);

double viscosity(particle in1, particle in2);

double heatFactor(vector<particle> netin[], int nx, int ny, int np);

double heat(vector<particle> netin[], int nx, int ny, int np);

double pressure(double u, double density, double density0); //压强函数

void next(vector<particle> netin[], vector<particle> netout[], double &dt); //函数

void divide(vector<particle> netin[], vector<particle> netout[]); //划分网格

void initialization(vector<particle> netin[]); //初始化

extern vec g = { 0.0,0.0 }, pump[N] = {}, wall[N] = {};

#endif // !_PROTOTYPES_H_
