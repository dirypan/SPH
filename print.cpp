// print.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <gl\glut.h>
#include <iostream>
#include <fstream>
#include <Windows.h>

using namespace std;

fstream fin;

int N;
int width = 20;
int height = 20;
int pointsize = 5;
const double Vmax = 50.0;

void init()
{
	glClearColor(1, 0, 0, 0);
	//glShadeModel(GL_FLAT);
}

void switchColor(double v)
{
	double r, g, b;
	g = (Vmax - v) / Vmax;
	r = v / Vmax;
	b = 0;
	glColor3f(r, g, b);
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(pointsize);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH, GL_NICEST);
	int i;
	double x, y, v;
	glBegin(GL_POINTS);
	for (i = 1; i <= N; i++)
	{
		fin >> x >> y >> v;
		fin.get();
		switchColor(v);
		glVertex2f(x, y);
	}
	glEnd();
	glFlush();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, width, 0, height);
}

bool Click = false;
int SleepFrequency = 1;

void finN()
{
	fin >> N;
	Sleep(1000 / SleepFrequency);
	glutPostRedisplay();
}

void motion(int x, int y)
{
	if (Click)
	{
		SleepFrequency = x / 5;
		if (SleepFrequency < 1) SleepFrequency = 1;
		if (SleepFrequency > 200) SleepFrequency = 200;
		cout << SleepFrequency << endl;
	}
}

void mouse(int button, int state, int x, int y)
{
	switch (button)
	{
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN)
		{
			Click = true;
			glutIdleFunc(finN);
		};
		break;
	case GLUT_RIGHT_BUTTON:
		if (state == GLUT_DOWN)
		{
			glutIdleFunc(NULL);
			Click = false;
		};
		break;
	default:
		break;
	}
}

int main(int argc, char* argv[])
{
	fin.open("D:\out.txt", ios_base::in);
	fin >> width >> height;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowSize(width * 50, height * 50);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(argv[0]);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutMainLoop();
	return 0;
}
