// functions.h
//

#pragma once

#include "stdafx.h"

#ifndef _FUNCTIONS_H_

#define _FUNCTIONS_H_

#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include "prototypes.h"

using namespace std;

vec::vec(const double &X, const double &Y) //构造函数
{
	x = X, y = Y;
}

double vec::mold2() //获取模方
{
	return x*x + y*y;
}

double vec::mold() //获取模
{
	return sqrt(x*x + y*y);
}

void vec::rotationAntiClock() //逆时针旋转
{
	double temp = x;
	x = -y, y = temp;
}

void vec::rotationClock() //顺时针旋转
{
	double temp = x;
	x = y, y = -temp;
}

vec vec::operator +(const vec &v) //加法
{
	return{ x + v.x, y + v.y };
}

vec vec::operator +=(const vec &v) //加法
{
	*this = *this + v;
	return *this;
}

vec vec::operator -(const vec &v) //减法
{
	return{ x - v.x, y - v.y };
}

vec vec::operator -=(const vec &v) //减法
{
	*this = *this - v;
	return *this;
}

vec vec::operator -() //负矢量
{
	return{ -x, -y };
}

vec vec::operator *(const double &k) //数乘
{
	return{ k*x, k*y };
}

vec operator *(const double &k, const vec &v) //数乘
{
	return{ k*v.x, k*v.y };
}

vec vec::operator *=(const double &k) //数乘
{
	*this = (*this) * k;
	return *this;
}

vec vec::operator /(const double &k) //数乘
{
	return{ x / k, y / k };
}

vec vec::operator /=(const double &k) //数乘
{
	*this = (*this) / k;
	return *this;
}

double vec::operator *(const vec &v) //点乘
{
	return x*v.x + y*v.y;
}

bool vec::operator ==(const vec &v) //判断相等
{
	return (x == v.x && y == v.y);
}

bool vec::operator !=(const vec &v) //判断不等
{
	return (x != v.x || y != v.y);
}

double W(double rx, double ry) //样条核函数
{
	double W;
	double s = sqrt(rx * rx + ry * ry) / h;
	if (s > 2.0)
		W = 0.0;
	else if (s > 1.0)
		W = 0.25 * (2.0 - s) * (2.0 - s) * (2.0 - s);
	else
		W = 1 - 3 * s*s*(2 - s) / 4;
	return W * 10.0 / (7.0 * pi * h * h);
}

double W(const vec &r) //样条核函数
{
	return W(r.x, r.y);
}

double W(const vec &r1, const vec &r2) //样条核函数
{
	return W(r2.x - r1.x, r2.y - r1.y);
}

vec gradW(double rx, double ry) //样条核函数梯度
{
	double GradW;
	double s = sqrt(rx * rx + ry * ry) / h;
	if (s > 2.0)
		GradW = 0.0;
	else if (s > 1.0)
		GradW = -0.75 * (2.0 - s) * (2.0 - s);
	else
		GradW = -3 * s*(4 - 3 * s) / 4;
	return (GradW * 10.0*s / (7.0 * pi * h * h*h*h))*vec(rx, ry);
}

vec gradW(const vec &r) //样条核函数梯度
{
	return gradW(r.x, r.y);
}

vec gradW1(const vec &r1, const vec &r2) //样条核函数梯度
{
	return gradW(r1.x - r2.x, r1.y - r2.y);
}

vec gradW2(const vec &r1, const vec &r2) //样条核函数梯度
{
	return gradW(r2.x - r1.x, r2.y - r1.y);
}

int periodicity(int nx, int ny) //周期性条件
{
	return (nx + Nx - 1) % Nx + 1 + ((ny + Ny - 1) % Ny + 1)*N_x;
}

void periodicity(particle &particlein) //周期性条件
{
	int k = (particlein.r.x < 2.0*h ? 1 : -1);
	while (particlein.r.x<2.0*h || particlein.r.x>Width + 2.0*h)
		particlein.r.x += k*Width;
	k = (particlein.r.y < 2.0*h ? 1 : -1);
	while (particlein.r.y<2.0*h || particlein.r.y>Height + 2.0*h)
		particlein.r.y += k*Height;
}

void setWall(int nx1, int ny1, int nx2, int ny2) //设置墙
{
	if (abs(nx2 - nx1) > abs(ny2 - ny1))
	{
		if (nx2 < nx1)
		{
			int temp_nx = nx1, temp_ny = ny1;
			nx1 = nx2, ny1 = ny2;
			nx2 = temp_nx, ny2 = temp_ny;
		}
		vec walln = { double(ny2 - ny1), double(nx1 - nx2) };
		walln /= walln.mold();
		for (int nx = nx1; nx <= nx2; ++nx)
			wall[nx + ((ny2 - ny1)*(nx - nx1) / (nx2 - nx1) + ny1)*N_x] = walln;
	}
	else if (ny1 != ny2)
	{
		if (ny2 < ny1)
		{
			int temp_nx = nx1, temp_ny = ny1;
			nx1 = nx2, ny1 = ny2;
			nx2 = temp_nx, ny2 = temp_ny;
		}
		vec walln = { double(ny2 - ny1), double(nx1 - nx2) };
		walln /= walln.mold();
		for (int ny = ny1; ny <= ny2; ++ny)
			wall[(nx2 - nx1)*(ny - ny1) / (ny2 - ny1) + nx1 + ny*N_x] = walln;
	}
}

void setWall(int nx[], int ny[], int n, char mod) //设置墙
{
	for (int i = 0; i < n - 1; ++i)
		setWall(nx[i], ny[i], nx[i + 1], ny[i + 1]);
	if (n > 2 && mod == 'C')
		setWall(nx[n - 1], ny[n - 1], nx[0], ny[0]);
}

void setPump(int nx1, int ny1, int nx2, int ny2, char mod, double pump_Force) //设置泵
{
	if (abs(nx2 - nx1) > abs(ny2 - ny1))
	{
		if (nx2 < nx1)
		{
			int tnx = nx1, tny = ny1;
			nx1 = nx2, ny1 = ny2;
			nx2 = tnx, ny2 = tny;
			mod = 'F' + 'R' - mod;
		}
		vec force = { double(ny2 - ny1), double(nx1 - nx2) };
		force *= pump_Force / force.mold();
		double k = (ny2 - ny1) / (nx2 - nx1), d;
		bool *flag = new bool[N];
		if (mod != 'F')
			force = -force;
		for (int nx = nx1, ny; nx <= nx2; ++nx)
		{
			ny = int(ny1 + k*(nx - nx1));
			for (int i = 0; i < 4; ++i)
			{
				for (int tnx = nx - i; tnx <= nx + i; ++tnx)
				{
					for (int tny = ny - i; tny <= ny + i; ++tny)
					{
						if (flag[periodicity(tnx, tny)]) 
						{
							d = abs(ny1 + k*(tnx - nx1) - tny) / sqrt(1 + k*k);
							pump[periodicity(tnx, tny)] += force / (1 + d);
							flag[periodicity(tnx, tny)] = false;
						}
					}
				}
			}
		}
		delete[] flag;
	}
	else if (ny1 != ny2)
	{
		if (ny2 < ny1)
		{
			int tnx = nx1, tny = ny1;
			nx1 = nx2, ny1 = ny2;
			nx2 = tnx, ny2 = tny;
			mod = 'F' + 'R' - mod;
		}
		vec force = { double(ny2 - ny1), double(nx1 - nx2) };
		force *= pump_Force / force.mold();
		double k = (nx2 - nx1) / (ny2 - ny1), d;
		bool *flag = new bool[N];
		if (mod != 'F')
			force = -force;
		for (int ny = ny1, nx; ny <= ny2; ++ny)
		{
			nx = int(nx1 + k*(ny - ny1));
			for (int i = 0; i < 4; ++i)
			{
				for (int tnx = nx - i; tnx <= nx + i; ++tnx)
				{
					for (int tny = ny - i; tny <= ny + i; ++tny)
					{
						if (flag[periodicity(tnx, tny)])
						{
							d = abs(nx1 + k*(tny - ny1) - tnx) / sqrt(1 + k*k);
							pump[periodicity(tnx, tny)] += force / (1 + d);
							flag[periodicity(tnx, tny)] = false;
						}
					}
				}
			}
		}
		delete[] flag;
	}
}

void setPump(int nx1, int ny1, int nx2, int ny2, vec force) //设置泵
{
	if (abs(nx2 - nx1) > abs(ny2 - ny1))
	{
		if (nx2 < nx1)
		{
			int tnx = nx1, tny = ny1;
			nx1 = nx2, ny1 = ny2;
			nx2 = tnx, ny2 = tny;
		}
		double k = (ny2 - ny1) / (nx2 - nx1), d;
		bool *flag = new bool[N];
		for (int nx = nx1, ny; nx <= nx2; ++nx)
		{
			ny = int(ny1 + k*(nx - nx1));
			for (int i = 0; i < 4; ++i)
			{
				for (int tnx = nx - i; tnx <= nx + i; ++tnx)
				{
					for (int tny = ny - i; tny <= ny + i; ++tny)
					{
						if (flag[periodicity(tnx, tny)])
						{
							d = abs(ny1 + k*(tnx - nx1) - tny) / sqrt(1 + k*k);
							pump[periodicity(tnx, tny)] += force / (1 + d);
							flag[periodicity(tnx, tny)] = false;
						}
					}
				}
			}
		}
		delete[] flag;
	}
	else if (ny1 != ny2)
	{
		if (ny2 < ny1)
		{
			int tnx = nx1, tny = ny1;
			nx1 = nx2, ny1 = ny2;
			nx2 = tnx, ny2 = tny;
		}								
		double k = (nx2 - nx1) / (ny2 - ny1), d;
		bool *flag = new bool[N];
		for (int ny = ny1, nx; ny <= ny2; ++ny)
		{
			nx = int(nx1 + k*(ny - ny1));
			for (int i = 0; i < 4; ++i)
			{
				for (int tnx = nx - i; tnx <= nx + i; ++tnx)
				{
					for (int tny = ny - i; tny <= ny + i; ++tny)
					{
						if (flag[periodicity(tnx, tny)])
						{
							d = abs(nx1 + k*(tny - ny1) - tnx) / sqrt(1 + k*k);
							pump[periodicity(tnx, tny)] += force / (1 + d);
							flag[periodicity(tnx, tny)] = false;
						}
					}
				}
			}
		}
		delete[] flag;
	}
}

void setGravity(vec &direction, double sgravity) //设置重力场
{
	if (direction.mold2() > 0.0)
		g = (sgravity / direction.mold())*direction;
	else
		g = vec(0.0, 0.0);
}

void setGravity(double x, double y, double sgravity) //设置重力场
{
	setGravity(vec(x, y), sgravity);
}

double viscosity(particle in1, particle in2) //人为粘性
{
	vec r = in1.r - in2.r;
	vec v = in1.v - in2.v;
	if (r*v >= 0.0)
		return 0.0;
	else
	{
		double density_av = (in1.density + in2.density) / 2.0;
		double factor = h*(v*r) / (r*r + e*h*h);
		return (-a*c0*factor + b*factor*factor) / density_av;
	}
}

double heatFactor(vector<particle> netin[], int nx, int ny, int np) //人为热流因子
{
	particle particle_i = netin[periodicity(nx, ny)][np], particle_j = {};
	double delta_v = 0.0;
	for (int tnx = nx - 1; tnx <= nx + 1; ++tnx)
	{
		for (int tny = ny - 1; tny <= ny + 1; ++tny)
		{
			for (int j = 0; j < netin[periodicity(tnx, tny)].size(); ++j)
			{
				particle_j = netin[periodicity(tnx, tny)][j];
				if (particle_j.density != 0.0)
					delta_v += particle_j.v*gradW(particle_i.r - particle_j.r) / particle_j.density;
			}
		}
	}
	return g1*h*c0 + g2*h*h*(abs(delta_v) - delta_v);
}

double heat(vector<particle> netin[], int nx, int ny, int np) //人为热流
{
	particle particle_i = netin[periodicity(nx, ny)][np], particle_j = {};
	vec rij;
	double heatFactori = heatFactor(netin, nx, ny, np), heatFactorj, heatFactor_av, density_av, heat = 0.0;
	for (int tnx = nx - 1; tnx <= nx + 1; ++tnx)
	{
		for (int tny = ny - 1; tny <= ny + 1; ++tny)
		{
			for (int j = 0; j < netin[periodicity(tnx, tny)].size(); ++j)
			{
				particle_j = netin[periodicity(tnx, tny)][j];
				rij = particle_i.r - particle_j.r;
				heatFactorj = heatFactor(netin, tnx, tny, j);
				heatFactor_av = (heatFactori + heatFactorj) / 2.0;
				density_av = (particle_i.density + particle_j.density) / 2.0;
				if (rij.mold2() != 0 && density_av != 0)
					heat += heatFactor_av*(particle_i.u - particle_j.u)*rij*gradW(rij) / (density_av*rij.mold2());
			}
		}
	}
	return 2.0 * m*heat;
}

double pressure(double u, double density, double density0) //压强
{
	if (density == 0.0 || density0 == 0.0)
		return 0.0;
	else
	{
		double P_H = c0*c0*(1.0 / density0 - 1.0 / density) / ((1.0 / density0 - s / density0 + s / density) * (1.0 / density0 - s / density0 + s / density));
		double u_H = 0.5*P_H*(1.0 / density0 - 1.0 / density);
		return (P_H + Gamma*density0*(u - u_H));
	}
}

void next(vector<particle> netin[], vector<particle> netout[], double &dt)
{
	int net_i, net_i_size, net_j, net_j_size;
	particle particlein_i = {}, particlein_j = {};
	double density = 0.0, du_dt = 0.0, P_d2;
	vec a = { 0.0, 0.0 }, P_d2_grad = { 0.0, 0.0 }, v_correction = { 0.0,0.0 };
	for (int nx = 1; nx <= Nx; ++nx)
	{
		for (int ny = 1; ny <= Ny; ++ny)
		{
			net_i = nx + ny*N_x;
			net_i_size = netin[net_i].size();
			netout[net_i].resize(net_i_size);
			for (int i = 0; i < net_i_size; ++i)
			{
				particlein_i = netin[net_i][i];
				a = { 0.0, 0.0 }, P_d2_grad = { 0.0, 0.0 }, v_correction = { 0.0,0.0 }, density = 0.0, du_dt = 0.0, P_d2 = particlein_i.pressure / (particlein_i.density*particlein_i.density);
				for (int tnx = nx - 1; tnx <= nx + 1; ++tnx)
				{
					for (int tny = ny - 1; tny <= ny + 1; ++tny)
					{
						net_j = periodicity(tnx, tny);
						net_j_size = netin[net_j].size();
						for (int j = 0; j < net_j_size; ++j)
						{
							particlein_j = netin[net_j][j];
							P_d2_grad = (particlein_j.pressure / (particlein_j.density*particlein_j.density) + P_d2 + viscosity(particlein_j, particlein_i))* gradW1(particlein_i.r, particlein_j.r);
							density += W(particlein_j.r - particlein_i.r);
							du_dt += (particlein_i.v - particlein_j.v)*(P_d2 + 0.5*viscosity(particlein_j, particlein_i))* (gradW1(particlein_i.r, particlein_j.r));
							a += P_d2_grad;
							v_correction += ((2.0*W(particlein_j.r - particlein_i.r)) / (particlein_i.density + particlein_j.density))*(particlein_i.v - particlein_j.v);
						}
					}
				}
				density *= m, du_dt *= m / 2.0, a = -m*a + g + pump[net_i];
				if (density > netin[net_i][i].density)
					du_dt += heat(netin, nx, ny, i);
				netout[net_i][i].density = density;
				netout[net_i][i].u = particlein_i.u + du_dt*dt;
				netout[net_i][i].pressure = pressure(netout[net_i][i].u, density, netin[net_i][i].density0);
				netout[net_i][i].a = a;
				netout[net_i][i].v = particlein_i.v + a*dt;
				v_correction *= correctionFactor*m;
				netout[net_i][i].r = particlein_i.r + (2.0*particlein_i.v + a*dt + 2.0*v_correction) *(dt / 2.0);
				periodicity(netout[net_i][i]);
			}
		}
	}
}

void divide(vector<particle> netin[], vector<particle> netout[]) //划分网格
{
	particle particlein_i = {};
	int net_i;
	for (int net_i = 0; net_i < N; ++net_i)
		netout[net_i].clear();
	for (int nx = 1; nx <= Nx; ++nx)
	{
		for (int ny = 1; ny <= Ny; ++ny)
		{
			for (int i = 0; i < netin[nx + ny*N_x].size(); ++i)
			{
				particlein_i = netin[nx + ny*N_x][i];
				net_i = (int((particlein_i.r.x) / (2.0 * h))) + (int((particlein_i.r.y) / (2.0 * h)))*N_x;
				if (wall[net_i].mold2()>0.0)
					particlein_i.v -= (particlein_i.v*wall[net_i])*wall[net_i];
				netout[net_i].push_back(particlein_i);
			}
		}
	}
}

void initialization(vector<particle> netin[]) //初始化
{
	srand(time(NULL));
	for (int nx = 1; nx <= Nx; ++nx)
	{
		for (int ny = 1; ny <= Ny; ++ny)
		{
			netin[nx + ny*N_x].resize(N0);
			for (int i = 0; i < N0; ++i)
				netin[nx + ny*N_x][i] = { { (nx + rand() / (RAND_MAX + 1.0)) * 2.0 * h , (ny + rand() / (RAND_MAX + 1.0)) * 2.0 * h } ,v0,{ 0.0,0.0 }, 0.0, 0.0, 0.0, 0.0 };
		}
	}
	double density = 0.0;
	for (int nx = 1; nx <= Nx; ++nx)
	{
		for (int ny = 1; ny <= Ny; ++ny)
		{
			for (int i = 0; i < N0; ++i)
			{
				density = 0.0;
				for (int x = nx - 1; x <= nx + 1; ++x)
				{
					for (int y = ny - 1; y <= ny + 1; ++y)
					{
						for (int j = 0; j < N0; ++j)
							density += W(netin[periodicity(x, y)][j].r - netin[nx + ny*N_x][i].r);
					}
				}
				density *= m;
				netin[nx + ny*N_x][i].density = density;
				netin[nx + ny*N_x][i].pressure = pressure(0.0, density, density);
				netin[nx + ny*N_x][i].density0 = density;
			}
		}
	}
}

#endif // !_FUNCTIONS_H_
