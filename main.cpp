// main.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "prototypes.h"
#include "functions.h"

using namespace std;

fstream fout;

const time_t t0 = time(NULL);

void output(vector<particle> netin[])
{
	int sum = 0, net_i;
	for (int net_i = 0; net_i < N; ++net_i)
	{
		sum += netin[net_i].size();
		if (pump[net_i].mold2() > 0.0)
			++sum;
		if (wall[net_i].mold2() > 0.0)
			++sum;
	}
	fout << sum << endl;
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
				fout << setw(8) << (2 * nx + 1)*h << " " << setw(8) << (2 * ny + 1)*h << " " << setw(8) << 100.0 << endl;
			if (wall[net_i].mold2() > 0.0)
				fout << setw(8) << (2 * nx + 1)*h << " " << setw(8) << (2 * ny + 1)*h << " " << setw(8) << 0.0 << endl;
		}
	}
}

int main()
{
	fout.open("D:\out.txt", ios_base::out);
	fout << N_x << ' ' << N_y << endl;
	vector<particle> net1[N], net2[N];
	double dt = 0.001;
	initialization(net1); //初始化
	setWall(Nx / 3, Ny / 2, 2 * Nx / 3, Ny / 2); 
	setWall(1, 1, Nx, 1); 
	setPump(Nx / 3, 1, Nx / 3, Ny);
	setPump(2 * Nx / 3, 1, 2 * Nx / 3, Ny); 
	setGravity(0.0, -1.0);
	for (int f = 0; f < 2500; ++f)
	{
		cout << "frame " << f << " : " << time(NULL) - t0 << 's' << endl;
		output(net1);
		next(net1, net2, dt);
		cout << "frame " << f << " after next() : " << time(NULL) - t0 << 's' << endl;
		divide(net2, net1);
	}
	return 0;
}
