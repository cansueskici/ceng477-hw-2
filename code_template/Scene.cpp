#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
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

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
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
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}


void Scene::initializeDepthBuffer(Camera *camera)
{
    if (this->depthBuffer.empty())
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            std::vector<double> rowOfDepths;

            for (int j = 0; j < camera->verRes; j++)
            {
                rowOfDepths.push_back(std::numeric_limits<double>::infinity());
            }

            this->depthBuffer.push_back(rowOfDepths);
        }
    }
    else
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            for (int j = 0; j < camera->verRes; j++)
            {
                this->depthBuffer[i][j] = std::numeric_limits<double>::infinity();
            }
        }
    }

}

Color Scene::clamp(Color &color)
{
    Color result;
    result.r = makeBetweenZeroAnd255(color.r);
    result.g = makeBetweenZeroAnd255(color.g);
    result.b = makeBetweenZeroAnd255(color.b);

    return result;
}


Matrix4 Scene::translate(Translation *translation, Matrix4 &matrix){
    Matrix4 translationMatrix = getIdentityMatrix();

    translationMatrix.values[0][3] = translation->tx;
    translationMatrix.values[1][3] = translation->ty;
    translationMatrix.values[2][3] = translation->tz;

    return multiplyMatrixWithMatrix(translationMatrix, matrix);
}


Matrix4 Scene::scale(Scaling *scaling, Matrix4 &matrix){
    Matrix4 scalingMatrix = getIdentityMatrix();
    scalingMatrix.values[0][0] = scaling->sx;
    scalingMatrix.values[1][1] = scaling->sy;
    scalingMatrix.values[2][2] = scaling->sz;

    return multiplyMatrixWithMatrix(scalingMatrix, matrix);
}


Matrix4 Scene::rotate(Rotation *rotation, Matrix4 &matrix){
    Vec3 u(rotation->ux, rotation->uy, rotation->uz), v, w;
    u = normalizeVec3(u);
    double angle_radians = (rotation->angle * M_PI)/180.0;
    double smallest = min(u.x, min(u.y, u.z));

    if (smallest == abs(u.x)){
        v.x = 0.0;
        v.y = -u.z;
        v.z = u.y;
    }else if(smallest == abs(u.y)){
        v.x = -u.z;
        v.y = 0.0;
        v.z = u.x;
    } else{
        v.x = u.y;
        v.y = -u.x;
        v.z = 0.0;
    }

    w = crossProductVec3(u, v);
    v = normalizeVec3(v);
    w = normalizeVec3(w);

    double m_inverse_values[4][4] = {{u.x, v.x, w.x, 0}, {u.y, v.y, w.y, 0}, {u.z, v.z, w.z, 0}, {0, 0, 0, 1}};
    Matrix4 m_inverse(m_inverse_values);

    double m_values[4][4] = {{u.x, u.y, u.z, 0}, {v.x, v.y, v.z, 0}, {w.x, w.y, w.z, 0}, {0, 0, 0, 1}};
    Matrix4 m(m_values);

    Matrix4 rotationMatrix = getIdentityMatrix();
    rotationMatrix.values[1][1] = cos(angle_radians);
    rotationMatrix.values[1][2] = sin(angle_radians)*(-1.0);
    rotationMatrix.values[2][1] = sin(angle_radians);
    rotationMatrix.values[2][2] = cos(angle_radians);

    return multiplyMatrixWithMatrix(m_inverse, multiplyMatrixWithMatrix(rotationMatrix, multiplyMatrixWithMatrix(m,matrix)));

}


Matrix4 Scene::modeling_transformation(Mesh *mesh){
    Matrix4 result = getIdentityMatrix();

    for (int i = 0; i < mesh->numberOfTransformations; ++i) {
        char type = mesh->transformationTypes[i];
        int id = mesh->transformationIds[i]-1;

        if (type == 's'){
            result = scale(scalings[id], result);
        }else if(type == 't'){
            result = translate(translations[id], result);
        }else{
            result = rotate(rotations[id], result);
        }
    }
    return result;
}


bool Scene::is_visible(double diff, double num, double &t_e, double &t_l){
    double t;

    if(diff > 0){
        t = num/diff;
        if (t > t_l){
            return false;
        }
        if (t > t_e){
            t_e = t;
        }
    } else if(diff < 0){
        t = num/diff;

        if (t < t_e){
            return false;

        } else if(t < t_l){
            t_l = t;
        }
    } else if(num > 0){
        return false;
    }

    return true;
}

bool Scene::clipping(Camera *cam, Vec4 &point1, Vec4 &point2){
    //each line will enter and leave twice,
    // if first L(eave) is before last E(ntrance) -> not visible/False
    // if t_pl < t_pe : return false

    bool visible = false;
    double t_pe = 0, t_pl = 1;
    double dx = point2.x - point1.x;
    double dy = point2.y - point1.y;
    double dz = point2.z - point1.z;
    double xmin = -1, ymin = -1, zmin = -1;
    double xmax = 1, ymax = 1, zmax = 1;

    double what = sqrt(pow(dx ,2) + pow(dy, 2) + pow(dz, 2));

    if(is_visible(dx, xmin-point1.x, t_pe, t_pl)){
        if (is_visible(dx*(-1), (point1.x-xmax), t_pe, t_pl)){
            if (is_visible(dy, (ymin-point1.y), t_pe, t_pl)){
                if (is_visible(dy*(-1), (point1.y-ymax), t_pe, t_pl)){
                    if (is_visible(dz, (zmin -point1.z ),t_pe, t_pl)){
                        if (is_visible(dz*(-1), (point1.z - zmax), t_pe, t_pl)){
                            visible = true;
                            if(t_pl < 1){
                                point2.x = point1.x + dx*t_pl;
                                point2.y = point1.y + dy*t_pl;
                                point2.z = point1.z + dz*t_pl;

                            }

                            if (t_pe > 0){
                                point1.x = point1.x + dx*t_pe;
                                point1.y = point1.y + dy*t_pe;
                                point1.z = point1.z + dz*t_pe;
                            }
                        }
                    }
                }
            }
        }
    }

    return visible;
}

Matrix4 Scene::cameraTransformation(Camera *cam){
    double cam_position_product1 =-1.0*(dotProductVec3(cam->u, cam->position));
    double cam_position_product2 =-1.0*(dotProductVec3(cam->v, cam->position));
    double cam_position_product3 =-1.0*(dotProductVec3(cam->w, cam->position));
    double m_cam_values[4][4] = {{cam->u.x, cam->u.y, cam->u.z,cam_position_product1 },
                                 {cam->v.x, cam->v.y, cam->v.z, cam_position_product2},
                                 {cam->w.x, cam->w.y, cam->w.z, cam_position_product3},
                                 {0, 0, 0, 1}};
    Matrix4 result(m_cam_values);
    return result;
}

Matrix4 Scene::projection(Camera *cam){
    double l = cam->left, r = cam->right, b = cam->bottom, t = cam->top, n = cam->near, f = cam->far;
    Matrix4 p2o;
    Matrix4 m_orth = getIdentityMatrix();
    double r_l = r - l, t_b = t - b, f_n = f - n;

    /*do M_orth
     *
     * ask if perspective, if yes calculate P20 and return the product.*/

    m_orth.values[0][0] = 2/r_l;
    m_orth.values[0][3] = -(r+l)/r_l;
    m_orth.values[1][1] = 2/t_b;
    m_orth.values[1][2] = -(t+b)/t_b;
    m_orth.values[2][2] = -2/f_n;
    m_orth.values[2][3] = -(f+n)/f_n;

    if(cam->projectionType == ORTOGRAPHIC_PROJECTION){
        return m_orth;
    } else{
        p2o.values[0][0] = n;
        p2o.values[1][1] = n;
        p2o.values[2][2] = f+n;
        p2o.values[2][3] = f*n;
        p2o.values[3][2] = -1;

        return m_orth*p2o;
    }

}

Matrix4 Scene::viewportTransformation(Camera *cam){
    double nx = cam->horRes, ny = cam->verRes;
    Matrix4 m_vp = getIdentityMatrix();
    m_vp.values[0][0] = nx/2;
    m_vp.values[0][3] = (nx-1)/2;
    m_vp.values[1][1] = ny/2;
    m_vp.values[1][3] = (ny-1)/2;
    m_vp.values[2][2] = 0.5;
    m_vp.values[2][3] = 0.5;

    return m_vp;
}

Vec4 Scene::perspectiveDivide(Vec4 vec){

    if (vec.t == 1){
        return vec;
    }

    Vec4 result;

    if (vec.t != 0){
        result.x = vec.x/vec.t;
        result.y = vec.y/vec.t;
        result.z = vec.z/vec.t;
        result.t = 1.0;
    }

    return result;
}

void Scene::depthBufferCheck(int x, int y, double depth, Color &c){
    if(x<0 || y<0 || x>=this->image.size() || y>= this->image[0].size()){
        return;
    }

    if(0 <= depth && depth <= 1 && depth < this->depthBuffer[x][y]){
        this->depthBuffer[x][y] = depth;
        this->image[x][y] = c;
    }
}

// deneysel
void Scene::lineRasterization(Vec4 &vertex0, Vec4 &vertex1, Color &c0, Color &c1, Camera *cam, std::vector<std::vector<double> > &buffer) {
    double x0 = vertex0.x, y0 = vertex0.y, z0 = vertex0.z, x1 = vertex1.x, y1 = vertex1.y, z1 = vertex1.z;
    double m, depth, multiplier, inc, dec;
    Color c;

    if (x0 == x1) {
        m = y1 > y0 ? MAXFLOAT : INT32_MIN;
    } else {
        m = (y1 - y0) / (x1 - x0);
    }

    if (x0 > x1) {
        swap(vertex0, vertex1);
        swap(c0, c1);
    }

    if (m > 0 && m <= 1.0) { //0,1
        double y = y0;
        double d = (y0 - y1) + ((x1 - x0) / 2);
        inc = (y0 - y1) + (x1 - x0);
        dec = y0 - y1;

        for (int x = x0; x <= x1; x++) {
            multiplier = (x - x0) / (x1 - x0);
            c = c0 * (1 - multiplier) + c1 * multiplier;
            depth = z0 * (1 - multiplier) + z1 * multiplier;
            depthBufferCheck(x, y, depth, c);
            if (d < 0) {
                y++;
                d += inc;
            } else {
                d += dec;
            }
        }

    } else if (m > 1.0) {
        double x = x0;
        double d = (x0 - x1) + ((y1 - y0) / 2);

        inc = (x0 - x1) + (y1 - y0);
        dec = x0 - x1;

        for (int y = y0; y <= y1; y++) {
            multiplier = (y - y0) / (y1 - y0);
            c = c0 * (1 - multiplier) + c1 * multiplier;
            depth = z0 * (1 - multiplier) + z1 * multiplier;
            depthBufferCheck(x, y, depth, c);

            if (d < 0) {
                x++;
                d += inc;
            } else {
                d += dec;
            }

        }

    } else if (m <= 0 && m >= -1.0) {
        double y = y1;
        double d = -(y0 - y1) + ((x1 - x0) / 2);
        inc = -(y0 - y1) + (x1 - x0);
        dec = y1 - y0;

        for (int x = x1; x >= x0; x--) {
//            if (x >= 0 && y >= 0 && x < cam->horRes && y < cam->verRes) {
//            }
            multiplier = (x - x0) / (x1 - x0);
            c = c0 * (1 - multiplier) + c1 * multiplier;
            depth = z0 * (1 - multiplier) + z1 * multiplier;
            depthBufferCheck(x, y, depth, c);
            if (d < 0) {
                y++;
                d += inc;
            } else {
                d += dec;
            }
        }


    } else if (m < -1.0) {
        double x = x0;
        double d = (x0 - x1) - ((y1 - y0) / 2);

        inc = (x0 - x1) - (y1 - y0);
        dec = x0 - x1;

        for (int y = y0; y >= y1; y--) {
            multiplier = (y - y0) / (y1 - y0);
            c = c0 * (1 - multiplier) + c1 * multiplier;
            depth = z0 * (1 - multiplier) + z1 * multiplier;
            depthBufferCheck(x, y, depth, c);

            if (d < 0) {
                x++;
                d += inc;
            } else {
                d += dec;
            }

        }

    }
}

void Scene::triangleRasterization(Vec4 &vertex0 , Vec4 &vertex1, Vec4 &vertex2, Color &c0, Color &c1, Color &c2,Camera *cam, std::vector<std::vector<double> > &buffer){
    double x0 = vertex0.x, y0 = vertex0.y, z0 = vertex0.z;
    double x1 = vertex1.x, y1 = vertex1.y, z1 = vertex1.z;
    double x2 = vertex2.x, y2 = vertex2.y, z2 = vertex2.z;

    double d0_1 = x2*(y0-y1) + y2*(x1-x0) + x0*y1 - y0*x1;
    double d1_2 = x0*(y1-y2) + y0*(x2-x1) + x1*y2 - y1*x2;
    double d2_0 = x1*(y2-y0) + y1*(x0-x2) + x2*y0 - y2*x0;

    double x_min = min({x0, x1, x2});
    x_min = max(x_min, 0.0);

    double x_max = max({x0, x1, x2});
    x_max = min(x_max, double(cam->horRes-1));

    double y_min = min({y0, y1, y2});
    y_min = max(y_min, 0.0);

    double y_max = max({y0, y1, y2});
    y_max = min(y_max, double(cam->verRes-1));

    for (int y = y_min; y <= y_max ; ++y) {
        for (int x = x_min; x <= x_max ; ++x) {

            double a = ((x*(y1-y2) + y*(x2-x1) + x1*y2 - y1*x2)/d1_2);
            double b = ((x*(y2-y0) + y*(x0-x2) + x2*y0 - y2*x0)/d2_0);
            double c = ((x*(y0-y1) + y*(x1-x0) + x0*y1 - y0*x1)/d0_1);

            if (a >= 0 && b >= 0 && c >= 0){
                double depth = a*z0 + b*z1 + c*z2;
                if (buffer[x][y] >= depth){
                    buffer[x][y] = depth;
                    Color color(c0*a + c1*b + c2*c);
                    image[x][y] = color;
                }
            }
        }
    }
}

bool Scene::backfaceCheck(Vec4 vertex0, Vec4 vertex1, Vec4 vertex2){
    Vec3 e1_0 = Vec3((vertex1.x - vertex0.x), (vertex1.y - vertex0.y), (vertex1.z - vertex0.z), -1);
    Vec3 e2_0 = Vec3((vertex2.x - vertex0.x), (vertex2.y - vertex0.y), (vertex2.z - vertex0.z), -1);

    Vec3 normal = normalizeVec3(crossProductVec3(e1_0, e2_0));
    Vec3 temp = Vec3(vertex0.x, vertex0.y, vertex0.z, -1);
    return (dotProductVec3(normal,temp) < 0 );

}


/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera){
	// TODO: Implement this function
    /*
     * model transformations
     *
     * camera transformations
     *
     * clipping & culling (backface culling)
     *
     * projection transformation
     *
     * rasterization (also color calculations) (solid vs wireframe aaaa)
     *
     * */

    Matrix4 camera_transformed = cameraTransformation(camera);
    Matrix4 projected = projection(camera);
    Matrix4 viewport = viewportTransformation(camera);
    Matrix4 p_c= multiplyMatrixWithMatrix(projected, camera_transformed);

    for (Mesh *mesh : this->meshes) {
        Matrix4 model_transformed = modeling_transformation(mesh);
        Matrix4 M_p_c_m = multiplyMatrixWithMatrix(p_c, model_transformed);

        for (int i = 0; i < mesh->numberOfTriangles; ++i) {
            Triangle triangle = mesh->triangles[i];

            //triangle rasterization
            int index0 = triangle.vertexIds[0] - 1;
            int index1 = triangle.vertexIds[1] - 1;
            int index2 = triangle.vertexIds[2] - 1;

            Color color0 = *(this->colorsOfVertices[index0]);
            Color color1 = *(this->colorsOfVertices[index1]);
            Color color2 = *(this->colorsOfVertices[index2]);

            Vec3 *pointer0 = this->vertices[index0];
            Vec3 *pointer1 = this->vertices[index1];
            Vec3 *pointer2 = this->vertices[index2];

            Vec4 firstVertex = Vec4(pointer0->x, pointer0->y, pointer0->z, 1, pointer0->colorId);
            Vec4 secondVertex = Vec4(pointer1->x, pointer1->y, pointer1->z, 1, pointer1->colorId);
            Vec4 thirdVertex = Vec4(pointer2->x, pointer2->y, pointer2->z, 1, pointer2->colorId);

            firstVertex = Vec4(multiplyMatrixWithVec4(M_p_c_m, firstVertex));
            secondVertex = Vec4(multiplyMatrixWithVec4(M_p_c_m, secondVertex));
            thirdVertex = Vec4(multiplyMatrixWithVec4(M_p_c_m, thirdVertex));

            if (cullingEnabled) {
                bool backfacing = backfaceCheck(firstVertex, secondVertex, thirdVertex);

                if(backfacing){
                    continue;
                }
            }

            if (camera->projectionType == PERSPECTIVE_PROJECTION) {
                firstVertex = perspectiveDivide(firstVertex);
                secondVertex = perspectiveDivide(secondVertex);
                thirdVertex = perspectiveDivide(thirdVertex);
            }

            if (mesh->type == WIREFRAME_MESH) {
                //wireframe, do clipping

                Vec4 temp0 = firstVertex;
                Vec4 temp1 = secondVertex;
                Vec4 temp2 = thirdVertex;


                bool clipped_line0 = clipping(camera, firstVertex, secondVertex);
                bool clipped_line1 = clipping(camera, secondVertex, thirdVertex);
                bool clipped_line2 = clipping(camera, thirdVertex, firstVertex);

                firstVertex = multiplyMatrixWithVec4(viewport, firstVertex);
                secondVertex = multiplyMatrixWithVec4(viewport, secondVertex);
                thirdVertex = multiplyMatrixWithVec4(viewport, thirdVertex);

                temp0 = multiplyMatrixWithVec4(viewport, temp0);
                temp1 = multiplyMatrixWithVec4(viewport, temp1);
                temp2 = multiplyMatrixWithVec4(viewport, temp2);

                if (clipped_line0) {
                    // line rasterization firstVertex, secondVertex
                    lineRasterization(firstVertex, secondVertex, color0, color1, camera, depthBuffer);
                }

                if (clipped_line1) {
                    // line rasterization secondVertex, thirdVertex
                    lineRasterization(temp1, thirdVertex, color1, color2, camera, depthBuffer);

                }

                if (clipped_line2) {
                    // line rasterization thirdVertex, firstVertex
                    lineRasterization(temp2, temp0, color2, color0, camera, depthBuffer);

                }

            } else {
                //solid
                firstVertex = multiplyMatrixWithVec4(viewport, firstVertex);
                secondVertex = multiplyMatrixWithVec4(viewport, secondVertex);
                thirdVertex = multiplyMatrixWithVec4(viewport, thirdVertex);
                // rasterize triangle
                triangleRasterization(firstVertex, secondVertex, thirdVertex, color0, color1, color2, camera, depthBuffer);

            }
        }
    }
}
