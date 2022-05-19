#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Mesh* > meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

    void ModellingTransform(Vec3 **pVec3);

    void ViewingTransform(Camera *pCamera, Vec3 **pVec3);

    void TriangleRasterization();

    void MidPointDrawer(Vec4 &v0, Vec4 &v1, Color & c0, Color & c1);

    bool LineClipper(Vec4 & line_1, Vec4 & line_2, Color & c0, Color &c1);

    void TriangleRasterizer(const Color *c0, const Color *c1, const Color *c2, Vec4 & v0,
                            Vec4 & v1,
                            Vec4 & v2, int nx, int ny);

    double Clamper(double value, double min, double max);

    bool CheckCulling(Vec3 & vec4, Vec3 & vec41, Vec3 & vec42, Vec3 camPos);
};

#endif
