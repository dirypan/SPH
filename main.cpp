// main.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <vector>
#include "prototypes.h"
#include "functions.h"

using namespace std;

fstream fin, fout;

const time_t t0 = time(NULL);

void set() //配置 
{
	fin >> Nx; //横向网格个数
	fin >> Ny; //纵向网格个数
	N_x = Nx + 2;
	N_y = Nx + 2;
	N = N_x*N_y; //总网格个数
	fin >> N0; //每个网格单元内初态的粒子数
	fin >> h; //网格宽度
	fin >> m; //粒子质量
	Width = 2.0 * h * Nx; //视场宽度
	Height = 2.0 * h * Ny; //视场高度
	fin >> pumpForce; //泵力系数
	fin >> gravity; //重力场系数
	fin >> Gamma; //Gruneisen 系数
	fin >> c0; //声速
	fin >> s; //激波系数
	fin >> correctionFactor; //位置修正系数
	fin >> alpha; //人为粘性可调参数
	fin >> h_0; //人为粘性可调参数
	fin >> e;	 //人为粘性可调参数
	fin >> a; //人为粘性可调参数
	fin >> b; //人为粘性可调参数
	fin >> g1; //人为热流可调参数
	fin >> g2; //人为热流可调参数
	fin >> v0.x >> v0.y; //粒子初速度
	pump = new vec[N], wall = new vec[N];
	cout << "横向网格个数 Nx: " << Nx << endl;
	cout << "纵向网格个数 Ny: " << Ny << endl;
	cout << "每个网格单元内初态的粒子数 N0: " << N0 << endl;
	cout << "网格宽度 h: " << h << endl;
	cout << "视场宽度 Width: " << Width << endl;
	cout << "视场高度 Height: " << Height << endl;
	cout << "粒子质量 m: " << m << endl;
	cout << "泵力系数 pumpForce: " << pumpForce << endl;
	cout << "重力场系数 gravity: " << gravity << endl;
	cout << "Gruneisen 系数 Gamma: " << Gamma << endl;
	cout << "声速 c0: " << c0 << endl;
	cout << "激波系数 s: " << s << endl;
	cout << "位置修正系数 correctionFactor: " << correctionFactor << endl;
	cout << "人为粘性可调参数 alpha: " << alpha << endl;
	cout << "人为粘性可调参数 h_0: " << h_0 << endl;
	cout << "人为粘性可调参数 e: " << e << endl;
	cout << "人为粘性可调参数 a: " << a << endl;
	cout << "人为粘性可调参数 b: " << b << endl;
	cout << "人为热流可调参数 g1: " << g1 << endl;
	cout << "人为热流可调参数 g2: " << g2 << endl << endl;
	cout << "粒子初速度 v0: (" << v0.x << ", " << v0.y << ')' << endl;
	char mod;
	int nx1, ny1, nx2, ny2;
	double gx, gy;
	while (!fin.eof()) 
	{
		cout << endl;
		fin >> mod;
		if (mod == 'G')
		{
			fin >> gx >> gy;
			cout << "Set gravity direction to (" << gx << ", " << gy << ')' << endl;
			setGravity(gx, gy);
		}
		else if (mod == 'P')
		{
			fin >> nx1 >> ny1 >> nx2 >> ny2;
			cout << "Set pump from (" << nx1 << ", " << ny1 << ") to (" << nx2 << ", " << ny2 << ')' << endl;
			setPump(nx1, ny1, nx2, ny2);
		}
		else 
		{
			fin >> nx1 >> ny1 >> nx2 >> ny2;
			cout << "Set wall from (" << nx1 << ", " << ny1 << ") to (" << nx2 << ", " << ny2 << ')' << endl;
			setWall(nx1, ny1, nx2, ny2);
		}
	}
	cout << endl;
}

void output(vector<particle> netin[]) //输出
{
	int tN = 0, net_i;
	for (int net_i = 0; net_i < N; ++net_i)
	{
		tN += netin[net_i].size();
		if (pump[net_i].mold2() > 0.0)
			tN += 5;
		if (wall[net_i].mold2() > 0.0)
			tN += 5;
	}
	fout << tN << endl;
	for (int nx = 0; nx <= Nx + 1; ++nx)
	{
		for (int ny = 0; ny <= Ny + 1; ++ny)
		{
			net_i = nx + ny*N_x;
			for (int i = 0; i < netin[net_i].size(); ++i)
			{
				fout << setw(8) << netin[net_i][i].r.x << " ";
				fout << setw(8) << netin[net_i][i].r.y << " ";
				fout << setw(8) << netin[net_i][i].v.mold() << endl;
			}
			if (pump[net_i].mold2() > 0.0)
			{
				vec pumpt = pump[net_i];
				pumpt.rotationClock();
				if (abs(pumpt.x) > abs(pumpt.y))
					pumpt *= 0.4*h / pumpt.x;
				else
					pumpt *= 0.4*h / pumpt.y;
				for (int i = -2; i <= 2; ++i)
					fout << setw(8) << (2 * nx + 1)*h + i*pumpt.x << " " << setw(8) << (2 * ny + 1)*h + i*pumpt.y << " " << setw(8) << 50.0 << endl;
			}
			if (wall[net_i].mold2() > 0.0)
			{						  
				vec wallt = wall[net_i];
				wallt.rotationClock();
				if (abs(wallt.x) > abs(wallt.y))
					wallt *= 0.4*h / wallt.x;
				else
					wallt *= 0.4*h / wallt.y;
				for (int i = -2; i <= 2; ++i)
					fout << setw(8) << (2 * nx + 1)*h + i*wallt.x << " " << setw(8) << (2 * ny + 1)*h + i*wallt.y << " " << setw(8) << 0.0 << endl; 
			}
		}
	}
}

int main()
{
	fin.open("D:\in.txt", ios_base::in);
	fout.open("D:\out.txt", ios_base::out);
	if (fin.is_open())
		set(); //配置
	else
		pump = new vec[N], wall = new vec[N];
	fout << N_x << ' ' << N_y << endl;
	vector<particle>* net1 = new vector<particle>[N], *net2 = new vector<particle>[N];
	double dt = 0.001;
	initialization(net1); //初始化 
	for (int f = 0; f < 2500; ++f)
	{
		cout << "frame " << f << " : " << time(NULL) - t0 << 's' << endl;
		output(net1); //输出
		next(net1, net2, dt);
		divide(net2, net1);
	}
	return 0;
}
