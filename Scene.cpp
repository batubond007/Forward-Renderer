#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function.

    Vec3 **vertex_copy = new Vec3*[vertices.size()];
    for (int i = 0; i < vertices.size(); ++i) {
        vertex_copy[i] = new Vec3(*vertices[i]);
    }
    ModellingTransform(vertex_copy);
    ViewingTransform(camera, vertex_copy);

    delete [] vertex_copy;
}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL) {
		str = pElement->GetText();

		if (strcmp(str, "enabled") == 0) {
			cullingEnabled = true;
		}
		else {
			cullingEnabled = false;
		}
	}

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		// read projection type
		str = pCamera->Attribute("type");

		if (strcmp(str, "orthographic") == 0) {
			cam->projectionType = 0;
		}
		else {
			cam->projectionType = 1;
		}

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read meshes
	pElement = pRoot->FirstChildElement("Meshes");

	XMLElement *pMesh = pElement->FirstChildElement("Mesh");
	XMLElement *meshElement;
	while (pMesh != NULL)
	{
		Mesh *mesh = new Mesh();

		pMesh->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = pMesh->Attribute("type");

		if (strcmp(str, "wireframe") == 0) {
			mesh->type = 0;
		}
		else {
			mesh->type = 1;
		}

		// read mesh transformations
		XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *clone_str;
		int v1, v2, v3;
		XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
		clone_str = strdup(str);

		row = strtok(clone_str, "\n");
		while (row != NULL)
		{
			sscanf(row, "%d %d %d", &v1, &v2, &v3);
			mesh->triangles.push_back(Triangle(v1, v2, v3));
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		meshes.push_back(mesh);

		pMesh = pMesh->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}

bool Contains(vector<int> &vec, int item){
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] == item)
            return true;
    }
    return false;
}

Matrix4 RotationMatrixCreator(Rotation rotation){
    Vec3 u = Vec3(rotation.ux, rotation.uy, rotation.uz, -1);
    Vec3 v, w;
    double abs_x = abs(u.x);
    double abs_y = abs(u.y);
    double abs_z = abs(u.z);
    if(abs_x <= abs_y && abs_x <= abs_z){
        v = Vec3(0, -1 * u.z, u.y, -1);
    }
    else if(abs_y <= abs_x && abs_y <= abs_z){
        v = Vec3(- u.z, 0 , u.x, -1);
    }
    else if(abs_z <= abs_x && abs_z <= abs_y){
        v = Vec3(- u.y, u.x , 0, -1);
    }
    w = crossProductVec3(u, v);
    v = normalizeVec3(v);
    w = normalizeVec3(w);

    double val_M[4][4] = {
            {u.x, u.y, u.z, 0},
            {v.x, v.y, v.z, 0},
            {w.x, w.y, w.z, 0},
            {0, 0, 0, 1}
    };

    Matrix4 matrix_M(val_M);
    double val_M_inverse[4][4] = {
            {u.x, v.x, w.x, 0},
            {u.y, v.y, w.y, 0},
            {u.z, v.z, w.z, 0},
            {0, 0, 0, 1}
    };
    Matrix4 matrix_M_inverse(val_M_inverse);


    double val_x[4][4] = {
            {1, 0, 0, 0},
            {0, cos(rotation.angle * M_PI / 180), -sin(rotation.angle * M_PI / 180), 0},
            {0, sin(rotation.angle * M_PI / 180), cos(rotation.angle * M_PI / 180), 0},
            {0, 0, 0, 1}
    };
    Matrix4 result_x(val_x);

    Matrix4 result(multiplyMatrixWithMatrix(matrix_M_inverse, multiplyMatrixWithMatrix(result_x, matrix_M)));
    return result;

}

void Scene::ModellingTransform(Vec3 **vertex_copy) {
    for (int i = 0; i < meshes.size(); ++i) {
        Mesh *_mesh = meshes[i];
        Matrix4 matrix_transform = getIdentityMatrix();
        vector<int> mesh_vertices;

        for (int j = 0; j < _mesh->numberOfTransformations; ++j) {
            int t_id = _mesh->transformationIds[j];
            char t_type = _mesh->transformationTypes[j];
            if (t_type == 't'){
                Translation* _translation = translations[t_id-1];
                double val[4][4] = {{1, 0, 0, _translation->tx},
                                    {0, 1, 0, _translation->ty},
                                    {0, 0, 1, _translation->tz},
                                    {0, 0, 0, 1}};
                Matrix4 mat(val);

                matrix_transform = multiplyMatrixWithMatrix(mat, matrix_transform);
            }
            else if(t_type == 'r'){
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Rotation *_rotation = rotations[t_id-1];
                Matrix4 rotation_matrix = RotationMatrixCreator(*_rotation);
                matrix_transform = multiplyMatrixWithMatrix(rotation_matrix, matrix_transform);
            }
            else if(t_type == 's'){
                Scaling *_scaling = scalings[t_id-1];
                double val[4][4] = {{_scaling->sx, 0, 0, 0},
                                    {0, _scaling->sy, 0, 0},
                                    {0, 0, _scaling->sz, 0},
                                    {0, 0, 0, 1}};
                Matrix4 mat(val);

                matrix_transform = multiplyMatrixWithMatrix(mat, matrix_transform);
            }
        }


        for (int t = 0; t < _mesh->numberOfTriangles; ++t) {
            Triangle _triangle = _mesh->triangles[t];

            if (!Contains(mesh_vertices, _triangle.getFirstVertexId())){
                mesh_vertices.push_back(_triangle.getFirstVertexId());
            }
            if(!Contains(mesh_vertices, _triangle.getSecondVertexId())){
                mesh_vertices.push_back(_triangle.getSecondVertexId());
            }
            if(!Contains(mesh_vertices, _triangle.getThirdVertexId())){
                mesh_vertices.push_back(_triangle.getThirdVertexId());
            }
        }

        for (int j = 0; j < mesh_vertices.size(); ++j) {
            Vec4 vertex_result;
            Vec3* vertex = vertex_copy[mesh_vertices[j]-1];

            vertex_result.x = vertex->x;
            vertex_result.y = vertex->y;
            vertex_result.z = vertex->z;
            vertex_result.t = 1;

            vertex_result = multiplyMatrixWithVec4(matrix_transform, vertex_result);

            vertex->x = vertex_result.x;
            vertex->y = vertex_result.y;
            vertex->z = vertex_result.z;
        }
    }
}

void Scene::ViewingTransform(Camera *pCamera, Vec3 **vertex_copy) {
    /*** CAMERA TRANSFORMAION ***/

    /** Translation **/
    double camera_translation_val_M1[4][4] = {
            {1, 0, 0, -(pCamera->pos.x)},
            {0, 1, 0, -(pCamera->pos.y)},
            {0, 0, 1, -(pCamera->pos.z)},
            {0, 0, 0, 1}
    };
    Matrix4 camera_translation_matrix_ToCenter(camera_translation_val_M1);

    /** Rotation **/
    Vec3 cam_w = pCamera->w;
    Vec3 cam_u = pCamera->u;
    Vec3 cam_v = pCamera->v;

    double camera_rotation_val[4][4] = {
            {cam_u.x, cam_u.y, cam_u.z, 0},
            {cam_v.x, cam_v.y, cam_v.z, 0},
            {cam_w.x, cam_w.y, cam_w.z, 0},
            {0, 0, 0, 1}
    };
    Matrix4 camera_rotation_matrix(camera_rotation_val);

    /** Camera transformation matrix **/
    Matrix4 camera_transformation_matrix = multiplyMatrixWithMatrix(camera_rotation_matrix, camera_translation_matrix_ToCenter);

    /*** ORTHOGRAPHIC PROJECTION ***/
    double r = pCamera->right;
    double l = pCamera->left;
    double t = pCamera->top;
    double b = pCamera->bottom;
    double f = pCamera->far;
    double n = pCamera->near;
    Matrix4 projection_matrix;
    if (pCamera->projectionType == 0){
        double orthographic_val[4][4] ={
                { 2/ (r-l), 0, 0, - (r+l) / (r-l)},
                {0, 2 / (t-b), 0, -(t+b)/ (t-b)},
                {0, 0, -2 / (f-n), -(f+n)/ (f-n)},
                {0, 0, 0, 1}
        };
        projection_matrix = Matrix4(orthographic_val);
    }
    else if(pCamera->projectionType == 1){
        double perspective_val[4][4] ={
                { 2*n/ (r-l), 0, (r+l) / (r-l), 0},
                { 0, 2 *n/ (t-b), (t+b)/ (t-b), 0},
                {0, 0, -(f +n)/ (f-n), -(2*f*n)/ (f-n)},
                {0, 0, -1, 0}
        };
        projection_matrix = Matrix4(perspective_val);
    }

    /*** VIEWPORT TRANSFORMATION ***/
    double ny = pCamera->verRes;
    double nx = pCamera->horRes;
    double viewport_val[4][4] ={
            {nx/2, 0, 0, (nx-1) / 2.0},
            {0, ny/2.0, 0, (ny-1) / 2.0},
            {0, 0, 0.5, 0.5},
            {0, 0, 0, 1}
    };
    Matrix4 viewport_matrix(viewport_val);

    Matrix4 camera_projection_matrix = multiplyMatrixWithMatrix(projection_matrix, camera_transformation_matrix);;

    for (int i = 0; i < meshes.size(); ++i) {
        Mesh *_mesh = meshes[i];

        for (int j = 0; j < _mesh->numberOfTriangles; ++j) {
            Triangle triangle = _mesh->triangles[j];
            Vec3 *i0 = vertex_copy[triangle.getFirstVertexId()-1];
            Vec3 *i1 = vertex_copy[triangle.getSecondVertexId()-1];
            Vec3 *i2 = vertex_copy[triangle.getThirdVertexId()-1];

            // Check culling in world space before going into projection space
            if (CheckCulling(*i0, *i1, *i2, pCamera->pos)){
                continue;
            }

            Color *c0 = colorsOfVertices[triangle.getFirstVertexId()-1];
            Color *c1 = colorsOfVertices[triangle.getSecondVertexId()-1];
            Color *c2 = colorsOfVertices[triangle.getThirdVertexId()-1];

            Vec4 pro_i0 = multiplyMatrixWithVec4(camera_projection_matrix, Vec4(i0->x, i0->y, i0->z, 1, i0->colorId));
            Vec4 pro_i1 = multiplyMatrixWithVec4(camera_projection_matrix, Vec4(i1->x, i1->y, i1->z, 1, i1->colorId));
            Vec4 pro_i2 = multiplyMatrixWithVec4(camera_projection_matrix, Vec4(i2->x, i2->y, i2->z, 1, i2->colorId));


            /// Wireframe
            if (_mesh->type == 0){

                /** Create Lines **/

                ///Line_01
                Vec4 Line_01_p1 = pro_i0;
                Vec4 Line_01_p2 = pro_i1;
                Color Line_01_c1(*c0);
                Color Line_01_c2(*c1);
                Line_01_p1.x /= Line_01_p1.t;
                Line_01_p1.y /= Line_01_p1.t;
                Line_01_p1.z /= Line_01_p1.t;
                Line_01_p1.t = 1;
                Line_01_p2.x /= Line_01_p2.t;
                Line_01_p2.y /= Line_01_p2.t;
                Line_01_p2.z /= Line_01_p2.t;
                Line_01_p2.t = 1;

                ///Line_12
                Vec4 Line_12_p1 = pro_i1;
                Vec4 Line_12_p2 = pro_i2;
                Color Line_12_c1(*c1);
                Color Line_12_c2(*c2);
                Line_12_p1.x /= Line_12_p1.t;
                Line_12_p1.y /= Line_12_p1.t;
                Line_12_p1.z /= Line_12_p1.t;
                Line_12_p1.t = 1;
                Line_12_p2.x /= Line_12_p2.t;
                Line_12_p2.y /= Line_12_p2.t;
                Line_12_p2.z /= Line_12_p2.t;
                Line_12_p2.t = 1;

                ///Line_20
                Vec4 Line_20_p1 = pro_i2;
                Vec4 Line_20_p2 = pro_i0;
                Color Line_20_c1(*c2);
                Color Line_20_c2(*c0);
                Line_20_p1.x /= Line_20_p1.t;
                Line_20_p1.y /= Line_20_p1.t;
                Line_20_p1.z /= Line_20_p1.t;
                Line_20_p1.t = 1;
                Line_20_p2.x /= Line_20_p2.t;
                Line_20_p2.y /= Line_20_p2.t;
                Line_20_p2.z /= Line_20_p2.t;
                Line_20_p2.t = 1;

                /// Clip Lines and check if they visible
                bool line_01_visible = LineClipper(Line_01_p1, Line_01_p2, Line_01_c1, Line_01_c2);
                bool line_12_visible = LineClipper(Line_12_p1, Line_12_p2, Line_12_c1, Line_12_c2);
                bool line_20_visible = LineClipper(Line_20_p1, Line_20_p2, Line_20_c1, Line_20_c2);


                if (line_01_visible){
                    ///Change to viewport space
                    Line_01_p1 = multiplyMatrixWithVec4(viewport_matrix, Line_01_p1);
                    Line_01_p2 = multiplyMatrixWithVec4(viewport_matrix, Line_01_p2);

                    MidPointDrawer(Line_01_p1, Line_01_p2, Line_01_c1, Line_01_c2);
                }
                if (line_12_visible){
                    ///Change to viewport space
                    Line_12_p1 = multiplyMatrixWithVec4(viewport_matrix, Line_12_p1);
                    Line_12_p2 = multiplyMatrixWithVec4(viewport_matrix, Line_12_p2);

                    MidPointDrawer(Line_12_p1, Line_12_p2, Line_12_c1, Line_12_c2);
                }
                if (line_20_visible){
                    ///Change to viewport space
                    Line_20_p1 = multiplyMatrixWithVec4(viewport_matrix, Line_20_p1);
                    Line_20_p2 = multiplyMatrixWithVec4(viewport_matrix, Line_20_p2);

                    MidPointDrawer(Line_20_p1, Line_20_p2, Line_20_c1, Line_20_c2);
                }

            }
            /// Solid
            else{
                pro_i0.x /= pro_i0.t;
                pro_i0.y /= pro_i0.t;
                pro_i0.z /= pro_i0.t;
                pro_i0.t = 1;

                pro_i1.x /= pro_i1.t;
                pro_i1.y /= pro_i1.t;
                pro_i1.z /= pro_i1.t;
                pro_i1.t = 1;

                pro_i2.x /= pro_i2.t;
                pro_i2.y /= pro_i2.t;
                pro_i2.z /= pro_i2.t;
                pro_i2.t = 1;

                Vec4 viewport_ind0 = multiplyMatrixWithVec4(viewport_matrix, pro_i0);
                Vec4 viewport_ind1 = multiplyMatrixWithVec4(viewport_matrix, pro_i1);
                Vec4 viewport_ind2 = multiplyMatrixWithVec4(viewport_matrix, pro_i2);

                TriangleRasterizer(c0, c1, c2, viewport_ind0, viewport_ind1, viewport_ind2, (int) nx, (int) ny);
            }
        }
    }

}

void Scene::MidPointDrawer(Vec4 &v0, Vec4 &v1, Color & c0, Color & c1) {
    double dx = v1.x - v0.x;
    double dy = v1.y - v0.y;

    double d, sign = 1;
    Color dc, c;

    /// Slope is [0,1]
    if (abs(dy) <= abs(dx)){
        if (v1.x < v0.x){
            swap(v0, v1);
            swap(c0, c1);
        }
        if (v1.y < v0.y)
            sign = -1;

        int x0 = (int) v0.x;
        int x1 = (int) v1.x;
        int y0 = (int) v0.y;
        int y1 = (int) v1.y;

        /// Basic Middle Point
        int y = y0;
        d = (2 * (y0 - y1) + sign * (x1 - x0));
        c = c0;
        dc = (c1 - c0) / (x1 - x0);

        double sum_val_1 = 2*((y0 - y1) + sign * (x1 - x0));
        double sum_val_2 = 2*(y0- y1);

        for (int x = x0; x < x1; ++x) {
            image[x][y] = c.Round();
            if (d*sign < 0){ //choose NE
                y = y + sign;
                d += sum_val_1;
            }
            else{   // Choose E
                d += sum_val_2;
            }
            c = c + dc;
        }
    }
    /// Slope is [1, infinite]
    else if (abs(dy) > abs(dx)){
        if (v1.y < v0.y){
            swap(v0, v1);
            swap(c0, c1);
        }
        if (v1.x < v0.x)
            sign = -1;

        int x0 = (int) v0.x;
        int x1 = (int) v1.x;
        int y0 = (int) v0.y;
        int y1 = (int) v1.y;

        int x = x0;
        c = c0;
        d = ( (x1 - x0) + sign * 0.5 * (y0 - y1) );
        dc = (c1 - c0) / (y1 - y0);

        double sum_val_1 = ((x1 - x0) + sign * (y0 - y1));
        double sum_val_2 = (x1- x0);

        for (int y = y0; y < y1; ++y) {
            image[x][y] = c.Round();
            if (d * sign > 0){
                x = x + sign;
                d += sum_val_1;
            }
            else{
                d+= sum_val_2;
            }
            c = c + dc;
        }
    }
}

double LineEquation(double x, double y, double x_first, double y_first, double x_sec, double y_sec){
    return (x * (y_first - y_sec)) + (y * (x_sec - x_first)) + (x_first * y_sec) - (y_first * x_sec);
}

void Scene::TriangleRasterizer(const Color * c0, const Color * c1, const Color * c2,
                               Vec4 & v0, Vec4 & v1, Vec4 & v2, int nx, int ny) {

    int x_min = (int) Clamper(min(min(v0.x, v1.x), v2.x), 0, nx-1);
    int y_min = (int) Clamper(min(min(v0.y, v1.y), v2.y), 0, ny-1);
    int x_max = (int) Clamper(max(max(v0.x, v1.x), v2.x), 0, nx-1);
    int y_max = (int) Clamper(max(max(v0.y, v1.y), v2.y), 0, ny-1);

    for(int y = y_min; y <= y_max; ++y){
        for(int x = x_min; x <= x_max; ++x){
            double f_12_xy = LineEquation(x,y, v1.x, v1.y, v2.x, v2.y);
            double f_12_x0y0 = LineEquation(v0.x,v0.y, v1.x,v1.y, v2.x,v2.y);

            double f_20_xy = LineEquation(x,y, v2.x, v2.y, v0.x, v0.y);
            double f_20_x1y1 = LineEquation(v1.x,v1.y, v2.x,v2.y, v0.x,v0.y);

            double f_01_xy = LineEquation(x,y, v0.x, v0.y, v1.x, v1.y);
            double f_01_x2y2 = LineEquation(v2.x,v2.y, v0.x,v0.y, v1.x,v1.y);;

            double alpha = f_12_xy / f_12_x0y0;
            double beta = f_20_xy / f_20_x1y1;
            double gamma = f_01_xy / f_01_x2y2;

            if(alpha >= 0 && beta >= 0 && gamma >= 0){
                image[x][y] = (( *c0 * alpha) + ( *c1 * beta) + ( *c2 * gamma)).Round() ;
            }
        }
    }
}

bool Visible(double den, double num, double & pE, double & pL){
    if (den > 0){   // Potentially entering
        double t = num / den;
        if(t > pL)
            return false;
        else if (t > pE)
            pE = t;
    }
    else if (den < 0){  // Potentially leaving
        double t = num / den;
        if (t < pE)
            return false;
        else if (t < pL)
            pL = t;
    }
    else if (num > 0)   // Line parallel to edge
        return false;
    return true;
}


bool Scene::LineClipper(Vec4 &line_1, Vec4 &line_2,Color & c1, Color &c2){

    double pE = 0, pL = 1;
    bool visible = false;

    double x0 = line_1.x;
    double x1 = line_2.x;
    double y0 = line_1.y;
    double y1 = line_2.y;
    double z0 = line_1.z;
    double z1 = line_2.z;

    double dx = x1 - x0;
    double dy = y1 - y0;
    double dz = z1 - z0;

    double dcr = c2.r - c1.r;
    double dcg = c2.g - c1.g;
    double dcb = c2.b - c1.b;
    Color dc = Color(dcr, dcg, dcb);

    double x_min = -1, y_min = -1, z_min = -1;
    double x_max = 1, y_max = 1, z_max = 1;

    if (
        Visible(dx, x_min - x0, pE, pL)             // Left
        && Visible(-dx, x0 - x_max, pE, pL)         // Right
        && Visible(dy, y_min - y0, pE, pL)          // Bottom
        && Visible(-dy, y1 - y_max, pE, pL)         // Up
        && Visible(dz, z_min - z0, pE, pL)          // Back
        && Visible(-dz, z0 - z_max, pE, pL))        // Forward
    {
        visible = true;

        if (pL < 1) {
            x1 = x0 + (dx * pL);
            y1 = y0 + (dy * pL);
            z1 = z0 + (dz * pL);
            c2 = c1 + (dc * pL);
        }
        if (pE > 0) {
            x0 = x0 + (dx * pE);
            y0 = y0 + (dy * pE);
            z0 = z0 + (dz * pE);
            c1 = c1 + (dc * pE);
        }
    }

    line_1.x = x0;
    line_1.y = y0;
    line_1.z = z0;

    line_2.x = x1;
    line_2.y = y1;
    line_2.z = z1;

    return visible;
}

double Scene::Clamper(double value, double min, double max) {
    return std::max(std::min(value, max), min);
}

bool Scene::CheckCulling(Vec3 & vert0, Vec3 & vert1, Vec3 & vert2, Vec3 camPos) {
    if (!cullingEnabled)
        return false;

    Vec3 a_b = Vec3(vert1.x - vert0.x,
                    vert1.y - vert0.y,
                    vert1.z - vert0.z, 1);
    Vec3 a_c = Vec3(vert2.x - vert0.x,
                    vert2.y - vert0.y,
                    vert2.z - vert0.z, 1);

    Vec3 triangle_normal = normalizeVec3(crossProductVec3(a_b, a_c));

    Vec3 middle_point;
    middle_point.x = (vert0.x + vert1.x + vert2.x) / 3;
    middle_point.y = (vert0.y + vert1.y + vert2.y) / 3;
    middle_point.z = (vert0.z + vert1.z + vert2.z) / 3;

    Vec3 eye_to_point = subtractVec3(middle_point, camPos);

    double dot_product = dotProductVec3(triangle_normal, eye_to_point);

    return dot_product > 0;
}
