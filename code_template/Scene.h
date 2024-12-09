#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color>> image;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;
    std::vector<std::vector<double>> depthBuffer;

	Scene(const char *xmlPath);

	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName, int osType);
	void forwardRenderingPipeline(Camera *camera);

    Matrix4 translate(Translation *translation, Matrix4 matrix);
    Matrix4 scale(Scaling *scaling, Matrix4 matrix);
    Matrix4 rotate(Rotation *rotation, Matrix4 matrix);
    Matrix4 modeling_transformation(Scene *scene, Mesh *mesh);
    bool is_visible(double diff, double num, double t_e, double t_l);
    bool clipping(Camera *cam, Vec4 point1, Vec4 point2);
    Matrix4 cameraTransformation(Camera *cam);
    Matrix4 projection(Camera *cam);
    Matrix4 viewportTransformation(Camera *cam);
    Vec4 perspectiveDivide(Vec4 vec);
    void lineRasterization(Vec4 vertex0, Vec4 vertex1, Color c0, Color c1, Camera *cam, std::vector<std::vector<double>> depthBuffer);
    void triangleRasterization(Vec4 vertex0 , Vec4 vertex1, Vec4 vertex2, Color c0, Color c1, Color c2,Camera *cam, std::vector<std::vector<double>> depthBuffer);
    bool backfaceCheck(Vec4 vertex0, Vec4 vertex1, Vec4 vertex2);
};

#endif
