#include "olcConsoleGameEngine.h"
using namespace std;

struct vec3d
{
	float x, y, z;
};

struct triangle
{
	vec3d p[3];
};

struct mesh
{
	vector<triangle> tris;
};

class Engine3D : public olcConsoleGameEngine
{
public:
	Engine3D()
	{
		m_sAppName = L"3D DEMO";
	}

private:
	mesh meshCube;

public:
	bool OnUserCreate() override
	{
		meshCube.tris = {

		};

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		return true;
	}
};

int main()
{
	Engine3D demo;
	if (demo.ConstructConsole(256, 240, 4, 4))
		demo.Start();

	return 0;
}