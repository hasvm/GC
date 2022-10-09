#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include "Ponto.h"
#include "PontoText.h"
#include "tinyxml2.cpp"
#include "tinyxml2.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif
#include <IL/il.h> 

using namespace std;
using namespace tinyxml2;

class Transform {
public:
	float x, y, z;
public:
	Transform(float a, float b, float c){
		x = a;
		y = b;
		z = c;
	}
public:
	virtual void vf() {}
};

class Rotate : public Transform {
public:
	float angle;
	float time;
public:
	Rotate(float a, float b, float c, float ag, float t) : Transform(a, b, c) {		
		angle = ag;
		time = t;
	}
};

class Scale : public Transform {
public:
	Scale(float a, float b, float c) : Transform(a, b, c) {
	}
};

class Translate {
public:
	bool empty = true;
	float time;
	vector<Ponto> catmullPoints;
	float catmulls[50][3];
};

class Material {
public:
	float *diffuse = NULL, *specular = NULL, *emissive = NULL, *ambient = NULL;
};

class Light {
public:
	string type;
	float pos[4];
	Material color;
	float cons = NULL, quad = NULL, linear = NULL, spotCutOff = NULL, *spotDirection = NULL, spotExponent = NULL;
	unsigned int id;
};

class Group {
public:
	vector<vector<Ponto>> models;
	vector<int> numPontos; //vector com o numero de pontos de cada modelo(se for 1 o modelo é desenhado sem VBOs)
	vector<vector<Ponto>> normals;
	vector<int> numNormals;
	vector<vector<PontoText>> texturas;
	vector<GLuint> texture;
	vector<Transform*> transforms; //scales/rotates
	Translate translate; //translate
	vector<Ponto> colors;	
	vector<Light> lights;
	vector<Material> material;
	vector<Group> nodes;
};

Group group; //group das figuras e transformações a desenhar

unsigned int nodeN = 0;
//eixos
int xflag = 0;
int yflag = 0;
int zflag = 0;
//coordenadas para a camera 
float r = 20.0;
float px = 0, py = 0, pz = 0;
float dx = 0, dy = 0, dz = 0;
float gt = 0.0, clocks = 0.0;
float alpha = 0;
float beta = M_PI / 4;

GLenum DRAWING_MODE = GL_FILL;

Translate vectorToArray(Translate trans) {
	//transforma o vetor num array. 
	for (unsigned int it = 0; it < trans.catmullPoints.size(); it++) {
		Ponto p = trans.catmullPoints.at(it);
		trans.catmulls[it][0] = p.x;
		trans.catmulls[it][1] = p.y;
		trans.catmulls[it][2] = p.z;
	}
	return trans;
}

void multMatrixVector(float *m, float *v, float *res) {
	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k)
			res[j] += v[k] * m[j * 4 + k];
	}
}

void getCatmullRomPoint(float t, float *p0, float *p1, float *p2, float *p3, float *pos, float *deriv) {

	// catmull-rom matrix
	float m[4][4] = { {-0.5f,  1.5f, -1.5f,  0.5f},
						{ 1.0f, -2.5f,  2.0f, -0.5f},
						{-0.5f,  0.0f,  0.5f,  0.0f},
						{ 0.0f,  1.0f,  0.0f,  0.0f} };

	// Compute A = M * P
	float a[3][4];

	float px[4] = { p0[0], p1[0], p2[0], p3[0] };
	float py[4] = { p0[1], p1[1], p2[1], p3[1] };
	float pz[4] = { p0[2], p1[2], p2[2], p3[2] };
	
	multMatrixVector(*m, px, a[0]);
	multMatrixVector(*m, py, a[1]);
	multMatrixVector(*m, pz, a[2]);

	// Compute pos = T * A
	float tv[4] = { t*t*t, t*t, t, 1 };
	float tvd[4] = { 3 * t*t, 2 * t, 1, 0 };

	pos[0] = tv[0] * a[0][0] + tv[1] * a[0][1] + tv[2] * a[0][2] + tv[3] * a[0][3];
	pos[1] = tv[0] * a[1][0] + tv[1] * a[1][1] + tv[2] * a[1][2] + tv[3] * a[1][3];
	pos[2] = tv[0] * a[2][0] + tv[1] * a[2][1] + tv[2] * a[2][2] + tv[3] * a[2][3];

	// compute deriv = T' * A
	deriv[0] = tvd[0] * a[0][0] + tvd[1] * a[0][1] + tvd[2] * a[0][2] + tvd[3] * a[0][3];
	deriv[1] = tvd[0] * a[1][0] + tvd[1] * a[1][1] + tvd[2] * a[1][2] + tvd[3] * a[1][3];
	deriv[2] = tvd[0] * a[2][0] + tvd[1] * a[2][1] + tvd[2] * a[2][2] + tvd[3] * a[2][3];
}

void getGlobalCatmullRomPoint(float gt, float *pos, float *deriv, Translate trans) {
	int pointCount = trans.catmullPoints.size();
	float t = gt * pointCount; // this is the real global t
	int index = floor(t);  // which segment
	t = t - index; // where within  the segment.

	// indices store the points
	int indices[4];
	indices[0] = (index + pointCount - 1) % pointCount;
	indices[1] = (indices[0] + 1) % pointCount;
	indices[2] = (indices[1] + 1) % pointCount;
	indices[3] = (indices[2] + 1) % pointCount;
	getCatmullRomPoint(t, trans.catmulls[indices[0]], trans.catmulls[indices[1]], trans.catmulls[indices[2]], trans.catmulls[indices[3]], pos, deriv);
}

void renderCatmullRomCurve(Translate trans) {
	// desenhar a curva usando segmentos de reta - GL_LINE_LOOP
	float pos[3], deriv[3];
	glBegin(GL_LINE_LOOP);
	for (float gt = 0; gt <= 1; gt += 0.001) {
		getGlobalCatmullRomPoint(gt, pos, deriv, trans);
		glColor3f(0.15, 0.15, 0.15);
		glVertex3f(pos[0], pos[1], pos[2]);
	}
	glEnd();
}

void renderCatmullTranslate(Translate trans) {
	float pos[3], deriv[3];
	clocks = glutGet(GLUT_ELAPSED_TIME);
	gt = fmod(clocks, (float)(trans.time * 1000)) / (trans.time * 1000);
	getGlobalCatmullRomPoint(gt, pos, deriv, trans);
	glTranslatef(pos[0], pos[1], pos[2]);
}

void renderRotate(Rotate rot) {
	if (rot.angle != -1)
		glRotatef(rot.angle, rot.x, rot.y, rot.z);
	else if (rot.time != -1) {
		clocks = glutGet(GLUT_ELAPSED_TIME);
		float angle = 360 * (fmod(clocks, (float)(rot.time * 1000)) / (rot.time * 1000));
		glRotatef(angle, rot.x, rot.y, rot.z);
	}
}

void lights(Group g) {
	vector<Light>::iterator it;
	for (it = g.lights.begin(); it != g.lights.end(); it++) {
		glLightfv(GL_LIGHT0 + it->id, GL_POSITION, it->pos);

		if (it->color.diffuse)
			glLightfv(GL_LIGHT0 + it->id, GL_DIFFUSE, it->color.diffuse);
		if (it->color.specular)
			glLightfv(GL_LIGHT0 + it->id, GL_SPECULAR, it->color.specular);
		if (it->color.emissive)
			glLightfv(GL_LIGHT0 + it->id, GL_EMISSION, it->color.emissive);
		if (it->color.ambient)
			glLightfv(GL_LIGHT0 + it->id, GL_AMBIENT, it->color.ambient);

		if (strcmp(it->type.c_str(), "POINT") == 0) {
			if (it->cons)
				glLightf(GL_LIGHT0 + it->id, GL_CONSTANT_ATTENUATION, it->cons);
			if (it->linear)
				glLightf(GL_LIGHT0 + it->id, GL_LINEAR_ATTENUATION, it->linear);
			if (it->quad)
				glLightf(GL_LIGHT0 + it->id, GL_QUADRATIC_ATTENUATION, it->quad);
		}

		if (strcmp(it->type.c_str(), "SPOTLIGHT") == 0) {
			if (it->spotCutOff)
				glLightf(GL_LIGHT0 + it->id, GL_SPOT_CUTOFF, it->spotCutOff);
			if (it->spotExponent)
				glLightfv(GL_LIGHT0 + it->id, GL_SPOT_DIRECTION, it->spotDirection);
			if (it->spotExponent)
				glLightf(GL_LIGHT0 + it->id, GL_SPOT_EXPONENT, it->spotExponent);
		}
	}
}

int loadTexture(string s) {
	unsigned int t, tw, th;
	unsigned char *texData;
	unsigned int texID;

	ilInit();
	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_LOWER_LEFT);
	ilGenImages(1, &t);
	ilBindImage(t);
	ilLoadImage((ILstring)s.c_str());
	tw = ilGetInteger(IL_IMAGE_WIDTH);
	th = ilGetInteger(IL_IMAGE_HEIGHT);
	ilConvertImage(IL_RGBA, IL_UNSIGNED_BYTE);
	texData = ilGetData();

	glGenTextures(1, &texID);
	glBindTexture(GL_TEXTURE_2D, texID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

	gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, tw, th, GL_RGBA, GL_UNSIGNED_BYTE, texData);

	glBindTexture(GL_TEXTURE_2D, 0);
	return texID;
}

vector<Ponto> getFigure(string figure, int *i, vector<Ponto> *normals, vector<PontoText> *textures) {
	vector<Ponto> coords;
	string line;
	ifstream file;
	file.open(figure);
	// Verifica se é para desenhar com VBO ou triangles e envia valor pelo i para a getGroup.
	int numSpaces, nVertices, n;
	if (getline(file, line)) {
		file.close();
		numSpaces = count(line.begin(), line.end(), ' ');
		if (numSpaces == 0) {
			// Criar nova stream para poder voltar atrás e registar o valor da 1ª linha e as coordenadas.
			ifstream fileVBO;
			fileVBO.open(figure);
			fileVBO >> nVertices;
			*i = nVertices;
			/* Guarda os vértices */
			n = 0;
			//Coordenadas
			while (n < nVertices / 3) {
				Ponto coord;
				fileVBO >> coord.x >> coord.y >> coord.z;
				coords.push_back(coord);
				n++;
			}
			n = 0;
			//Normais
			while (n < nVertices / 3) {
				Ponto norm;
				fileVBO >> norm.x >> norm.y >> norm.z;
				(*normals).push_back(norm);
				n++;
			}
			//Pts de textura
			n = 0;
			while (n < nVertices / 3 && !fileVBO.eof()) {
				PontoText tex;
				fileVBO >> tex.x >> tex.y;
				(*textures).push_back(tex);
				n++;
			}
			fileVBO.close();
			return coords;
		}
		// Se tiver espaços na 1ª linha é pq são coordenadas, então não é para desenhar com VBO, logo i=1
		else *i = 1;
	}

	// Simplemente regista coordenadas da figura, nos restantes casos (quando é para desenhar com triângulos).
	ifstream fileGetLines;
	fileGetLines.open(figure);
	int x = 0;
	// Conta número de linhas para poder dividir por 3 e saber onde estão as coordenadas, as normais e os pts de textura.
	while (!fileGetLines.eof()) {
		Ponto coord;
		fileGetLines >> coord.x >> coord.y >> coord.z;
		x++;
	}
	fileGetLines.close();

	ifstream fileTriang;
	fileTriang.open(figure);
	n = 0;
	//Coordenadas
	while (n < x / 3) {
		Ponto coord;
		fileTriang >> coord.x >> coord.y >> coord.z;
		coords.push_back(coord);
		n++;
	}
	n = 0;
	//Normais
	while (n < x / 3) {
		Ponto norm;
		fileTriang >> norm.x >> norm.y >> norm.z;
		(*normals).push_back(norm);
		n++;
	}
	n = 0;
	//Pts de textura
	while (n < x / 3) {
		PontoText tex;
		fileTriang >> tex.x >> tex.y;
		(*textures).push_back(tex);
		n++;
	}
	fileTriang.close();
	return coords;
}

Group load(XMLElement* node) {
	string path = "../../3DFiles/"; //localização dos ficheiros 3D para desenhar
	XMLElement* gptr = (node->FirstChildElement()); //apontador para o filho do node group(group pointer)
	int i = 0;
	Group g;

	for (; gptr != nullptr; gptr = gptr->NextSiblingElement()) {
		string tag = gptr->Value();
		if (strcmp(tag.c_str(), "translate") == 0) {
			XMLElement* translateNode = gptr->FirstChildElement();
			g.translate.empty = false;
			g.translate.time = gptr->DoubleAttribute("time");
			for (; translateNode != nullptr; translateNode = translateNode->NextSiblingElement()) {
				g.translate.catmullPoints.push_back(Ponto(translateNode->DoubleAttribute("X"), translateNode->DoubleAttribute("Y"), translateNode->DoubleAttribute("Z")));
			}
			g.translate = vectorToArray(g.translate);
		}
		else if (strcmp(tag.c_str(), "rotate") == 0) {
			if (gptr->DoubleAttribute("angle")) 
				g.transforms.push_back(new Rotate(gptr->DoubleAttribute("axisX"), gptr->DoubleAttribute("axisY"), gptr->DoubleAttribute("axisZ"), gptr->DoubleAttribute("angle"), -1));
			else if (gptr->DoubleAttribute("time")) 
				g.transforms.push_back(new Rotate(gptr->DoubleAttribute("axisX"), gptr->DoubleAttribute("axisY"), gptr->DoubleAttribute("axisZ"), -1, gptr->DoubleAttribute("time")));
		}
		else if (strcmp(tag.c_str(), "scale") == 0) {
			g.transforms.push_back(new Scale(gptr->DoubleAttribute("X"), gptr->DoubleAttribute("Y"), gptr->DoubleAttribute("Z")));
		}
		else if (strcmp(tag.c_str(), "color") == 0) {
			g.colors.push_back(Ponto(gptr->DoubleAttribute("R"), gptr->DoubleAttribute("G"), gptr->DoubleAttribute("B")));
		}
		else if (strcmp(tag.c_str(), "models") == 0) {
			XMLElement* modelsNode = gptr->FirstChildElement();
			for (; modelsNode != NULL; modelsNode = modelsNode->NextSiblingElement()) {
				string figureName = path + modelsNode->Attribute("file");

				if (modelsNode->Attribute("texture")) {
					string textureName = path + modelsNode->Attribute("texture");
					g.texture.push_back(loadTexture(textureName));
				}
				if (modelsNode->Attribute("ambR")) {
					Material m; float* f;
					f = (float*)malloc(sizeof(float) * 4);
					f[0] = modelsNode->DoubleAttribute("ambR");
					f[1] = modelsNode->DoubleAttribute("ambG");
					f[2] = modelsNode->DoubleAttribute("ambB");
					f[3] = 1.0f;
					m.emissive = f;
					g.material.push_back(m);
				}
				if (modelsNode->Attribute("specR")) {
					Material m; float* f;
					f = (float*)malloc(sizeof(float) * 4);
					f[0] = modelsNode->DoubleAttribute("specR");
					f[1] = modelsNode->DoubleAttribute("specG");
					f[2] = modelsNode->DoubleAttribute("specB");
					f[3] = 1.0f;
					m.emissive = f;
					g.material.push_back(m);
				}

				if (modelsNode->Attribute("diffR")) {
					Material m; float* f;
					f = (float*)malloc(sizeof(float) * 4);
					f[0] = modelsNode->DoubleAttribute("diffR");
					f[1] = modelsNode->DoubleAttribute("diffG");
					f[2] = modelsNode->DoubleAttribute("diffB");
					f[3] = 1.0f;
					m.emissive = f;
					g.material.push_back(m);
				}

				if (modelsNode->Attribute("emiR")) {
					Material m; float* f;
					f = (float*)malloc(sizeof(float) * 4);
					f[0] = modelsNode->DoubleAttribute("emiR");
					f[1] = modelsNode->DoubleAttribute("emiG");
					f[2] = modelsNode->DoubleAttribute("emiB");
					f[3] = 1.0f;
					m.emissive = f;
					g.material.push_back(m);
				}
				int nVertex;
				vector<PontoText> textures;
				vector<Ponto> normals;
				vector<Ponto> coords = getFigure(figureName, &nVertex, &normals, &textures);
				g.models.push_back(coords);
				g.numPontos.push_back(nVertex);
				g.normals.push_back(normals);
				g.numNormals.push_back(normals.size());
				g.texturas.push_back(textures);
				i++;
			}
		}
		else if (strcmp(tag.c_str(), "lights") == 0) {
			XMLElement* lightsNode = gptr->FirstChildElement();
			glEnable(GL_LIGHTING);
			for (; lightsNode && nodeN < 8; lightsNode = lightsNode->NextSiblingElement(), nodeN++) {
				glEnable(GL_LIGHT0 + nodeN);
				Light l;
				l.id = nodeN;

				if (lightsNode->Attribute("type")) 
					l.type = lightsNode->Attribute("type");

				l.pos[0] = lightsNode->DoubleAttribute("X");
				l.pos[1] = lightsNode->DoubleAttribute("Y");
				l.pos[2] = lightsNode->DoubleAttribute("Z");
				l.pos[3] = 1.0f;

				if (strcmp(l.type.c_str(), "DIRECTIONAL") == 0) l.pos[3] = 0.0f;

				if (lightsNode->Attribute("diffR")) {
					l.color.diffuse = (float*)malloc(sizeof(float) * 4);
					l.color.diffuse[0] = lightsNode->DoubleAttribute("diffR");
					l.color.diffuse[1] = lightsNode->DoubleAttribute("diffG");
					l.color.diffuse[2] = lightsNode->DoubleAttribute("diffG");
					l.color.diffuse[3] = 1.0f;
				}
				if (lightsNode->Attribute("specR")) {
					l.color.specular = (float*)malloc(sizeof(float) * 4);
					l.color.specular[0] = lightsNode->DoubleAttribute("specR");
					l.color.specular[1] = lightsNode->DoubleAttribute("specG");
					l.color.specular[2] = lightsNode->DoubleAttribute("specB");
					l.color.specular[3] = 1.0f;
				}
				if (lightsNode->Attribute("emiR")) {
					l.color.emissive = (float*)malloc(sizeof(float) * 4);
					l.color.emissive[0] = lightsNode->DoubleAttribute("emiR");
					l.color.emissive[1] = lightsNode->DoubleAttribute("emiG");
					l.color.emissive[2] = lightsNode->DoubleAttribute("emiB");
					l.color.emissive[3] = 1.0f;
				}
				if (lightsNode->Attribute("ambR")) {
					l.color.ambient = (float*)malloc(sizeof(float) * 4);
					l.color.ambient[0] = lightsNode->DoubleAttribute("ambR");
					l.color.ambient[1] = lightsNode->DoubleAttribute("ambG");
					l.color.ambient[2] = lightsNode->DoubleAttribute("ambB");
					l.color.ambient[3] = 1.0f;
				}
				if (lightsNode->Attribute("const"))
					l.cons = lightsNode->DoubleAttribute("const");
				if (lightsNode->Attribute("linear"))
					l.linear = lightsNode->DoubleAttribute("linear");
				if (lightsNode->Attribute("quad"))
					l.quad = lightsNode->DoubleAttribute("quad");
				if (lightsNode->Attribute("spotDirX")) {
					l.spotDirection = (float*)malloc(sizeof(float) * 4);
					l.spotDirection[0] = lightsNode->DoubleAttribute("spotDirX");
					l.spotDirection[1] = lightsNode->DoubleAttribute("spotDirY");
					l.spotDirection[2] = lightsNode->DoubleAttribute("spotDirZ");
				}
				if (lightsNode->Attribute("spotCutOff"))
					l.spotCutOff = lightsNode->DoubleAttribute("spotCutOff");
				if (lightsNode->Attribute("spotExponent"))
					l.spotExponent = lightsNode->DoubleAttribute("spotExponent");

				g.lights.push_back(l);
			}
		}
		else if (strcmp(tag.c_str(), "group") == 0) {
			g.nodes.push_back(load(gptr));
		}
	}
	return g;
}

void drawFiguresVBOs(vector<Ponto> coords, vector<Ponto> normals, vector<PontoText> textures, int nPoints, int nNormals) {
	GLuint buffers[3];
	glGenBuffers(3, buffers);

	float *bufVertex = new float[nPoints + 5];
	float *bufNormal = new float[nNormals * 3];
	float *bufTextures = new float[nNormals * 3];

	vector<Ponto>::iterator it_coords;
	int it = 0;
	for (it_coords = coords.begin(); it_coords != coords.end(); it_coords++) {
		bufVertex[it++] = it_coords->x;
		bufVertex[it++] = it_coords->y;
		bufVertex[it++] = it_coords->z;
	}

	vector<Ponto>::iterator it_normals;
	int i = 0;
	for (it_normals = normals.begin(); it_normals != normals.end(); it_normals++) {
		bufNormal[i++] = it_normals->x;
		bufNormal[i++] = it_normals->y;
		bufNormal[i++] = it_normals->z;
	}

	vector<PontoText>::iterator it_textures;
	int n = 0;
	for (it_textures = textures.begin(); it_textures != textures.end(); it_textures++) {
		bufTextures[n++] = it_textures->x;
		bufTextures[n++] = it_textures->y;
	}

	glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
	glBufferData(GL_ARRAY_BUFFER, it * sizeof(float), bufVertex, GL_STATIC_DRAW);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
	glBufferData(GL_ARRAY_BUFFER, i * sizeof(float), bufNormal, GL_STATIC_DRAW);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
	glBufferData(GL_ARRAY_BUFFER, n * sizeof(float), bufTextures, GL_STATIC_DRAW);
	glTexCoordPointer(2, GL_FLOAT, 0, 0);

	glDrawArrays(GL_TRIANGLE_STRIP, 0, nPoints / 3);

	delete[] bufVertex;
	delete[] bufNormal;
	delete[] bufTextures;

	glDeleteBuffers(3, buffers);
}

void drawFiguresTriangles(vector<Ponto> figures, vector<Ponto> normais) {
	glBegin(GL_TRIANGLES);
	vector<Ponto>::iterator it_coords;
	vector<Ponto>::iterator it_normais;
	for (it_coords = figures.begin(), it_normais = normais.begin();
		it_coords != figures.end() && it_normais != normais.end();
		it_coords++, it_normais++) {
		glVertex3f(it_coords->x, it_coords->y, it_coords->z);
		glNormal3f(it_coords->x, it_coords->y, it_coords->z);
	}
	glEnd();
}

void display(Group g) {
	glPushMatrix();
	vector<Material>::iterator itmat = g.material.begin();
	vector<GLuint>::iterator ittexture = g.texture.begin();
	vector<vector<PontoText>>::iterator ittexturas = g.texturas.begin(); //array pontos textura
	vector<vector<Ponto>>::iterator itnormais = g.normals.begin(); //array pontos normais
	vector<Ponto>::iterator itcolor = g.colors.begin(); //iterador para percorrer o vector com as cores de cada modelo
	vector<int>::iterator itnumPontos = g.numPontos.begin();
	vector<int>::iterator itnumNormals = g.numNormals.begin();
	
	//Transformações
	if (!g.translate.empty) {
		renderCatmullRomCurve(g.translate);
		renderCatmullTranslate(g.translate);
		g.translate.catmullPoints.clear();
		g.translate.catmullPoints.shrink_to_fit();
	}
	for (vector<Transform*>::iterator it = g.transforms.begin(); it != g.transforms.end(); it++) {
		if (Rotate *r = dynamic_cast<Rotate*>(*it)) {
			renderRotate(*r);
		}
		if (Scale *s = dynamic_cast<Scale*>(*it)) {
			glScalef(s->x, s->y, s->z);
		}
	}
	//Modelos
	for (vector<vector<Ponto>>::iterator it = g.models.begin(); it != g.models.end(); it++) {
		glColor3f(itcolor->x, itcolor->y, itcolor->z);
		if (itmat != g.material.end()) {
			if (itmat->ambient) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, itmat->ambient); 
			if (itmat->diffuse) glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, itmat->diffuse);
			if (itmat->emissive) glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, itmat->emissive);
			if (itmat->specular) glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, itmat->specular);
		}
		if ((*itnumPontos) > 1) {
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, (*ittexture));
			drawFiguresVBOs(*it, *itnormais, *ittexturas, (*itnumPontos), (*itnumNormals));
			glBindTexture(GL_TEXTURE_2D, 0);
			glDisable(GL_TEXTURE_2D);
		}
		else if ((*itnumPontos) == 1) drawFiguresTriangles(*it, *itnormais);  //sem VBOs
		//itcolor++; itnumPontos++; itnumNormals++; itmat++; itnormais++; ittexturas++; ittexture++;
	}
	for (vector<Group>::iterator itg = g.nodes.begin(); itg != g.nodes.end(); itg++)
		display(*itg);
	lights(g);
	glPopMatrix();

	//reiniciar para valores default do glMaterial para nao aplicar aos futuros objetos
	float emission[4] = { 0,0,0,1 };
	float ambient[4] = { 0.2,0.2,0.2,1.0 };
	float specular[4] = { 0,0,0,1 };
	float diffuse[4] = { 0.8,0.8,0.8,1.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
}

void changeSize(int w, int h) {

	// Prevent xt divide by zero, when window is too short
	// (you cant make xt window with zero width).
	if (h == 0)
		h = 1;

	// compute window's aspect ratio 
	float ratio = w * 1.0 / h;

	// Set the projection matrix as current
	glMatrixMode(GL_PROJECTION);
	// Load Identity Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set perspective
	gluPerspective(45.0f, ratio, 1.0f, 1000.0f);

	// return to the model view matrix DRAWING_MODE
	glMatrixMode(GL_MODELVIEW);
}

void renderScene(void) {

	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPolygonMode(GL_FRONT_AND_BACK, DRAWING_MODE); //forma como as figuras vão ser desenhadas

	px = r * cosf(beta) * sinf(alpha);
	py = r * sinf(beta);
	pz = r * cosf(beta) * cosf(alpha);

	// set the camera
	glLoadIdentity();
	gluLookAt(px + dx, py + dy, pz + dz,
		dx, dy, dz,
		0.0f, 1.0f, 0.0f);
	
	//Para desenhar os eixos
	{
		//eixo x
		glColor3f(0, 0, 1);
		glBegin(GL_LINES);
		if (xflag) {
			glVertex3f(0, 0, 0);
			glVertex3f(200, 0, 0);
		}
		//eixo y
		if (yflag) {
			glVertex3f(0, 0, 0);
			glVertex3f(0, 200, 0);
		}
		//eixo z
		if (zflag) {
			glVertex3f(0, 0, 0);
			glVertex3f(0, 0, 200);
		}
		glEnd();
	}

	display(group); //desenha as figuras
	// End of frame
	glutSwapBuffers();
}

void keyboardEvents(unsigned char key, int x, int y) {
	switch (key) {
	//para alterar o drawing mode
	case '1':
		DRAWING_MODE = GL_FILL;
		glutPostRedisplay();
		break;
	case '2':
		DRAWING_MODE = GL_LINE;
		glutPostRedisplay();
		break;
	case '3':
		DRAWING_MODE = GL_POINT;
		glutPostRedisplay();
		break;
	//para desenhar/apagar os eixos
	case 'x':
		if (xflag == 1) 
			xflag = 0;
		else xflag = 1;
		glutPostRedisplay();
		break;
	case 'y':
		if (yflag == 1)
			yflag = 0;
		else yflag = 1;
		glutPostRedisplay();
		break;
	case 'z':
		if (zflag == 1)
			zflag = 0;
		else zflag = 1;
		glutPostRedisplay();
		break;
	//para mover a camera
	case 'w':
		dy -= 6.9;
		glutPostRedisplay();
		break;
	case 's':
		dy += 6.9;
		glutPostRedisplay();
		break;
	case 'a':
		dx += 6.9;
		glutPostRedisplay();
		break;
	case 'd':
		dx -= 6.9;
		glutPostRedisplay();
		break;
	case 'r':
		dz += 6.9;
		glutPostRedisplay();
		break;
	case 't':
		dz -= 6.9;
		glutPostRedisplay();
		break;
	case '+':
		r -= 6.9f;
		glutPostRedisplay();
		break;
	case '-':
		r += 6.9f;
		glutPostRedisplay();
		break;
	//para terminar o programa
	case 27:
		exit(0);
		glutPostRedisplay();
		break;
	default:
		break;
	}
}

void specialEvents(int key_code, int x, int y) {
	GLfloat colors[][3] = { { 0.0f, 0.0f, 0.0f},{1.0f, 1.0f, 1.0f },{0.0f, 1.0f, 0.0f },{0.0f, 0.0f, 1.0f },{0.0f, 1.0f, 1.0f },{1.0f, 0.0f, 1.0f }};
	//Black/White/Green/Blue/LightBlue/Purple
	static int back = 0;
	switch (key_code) {
	case GLUT_KEY_F1:
		back = (back + 1) % 6;
		glClearColor(colors[back][0], colors[back][1], colors[back][2], 1.0f);
		glutPostRedisplay();
		break;
	case GLUT_KEY_UP:
		if (beta + 0.1 > 1.5) beta = 1.5;
		else beta += 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_DOWN:
		if (beta - 0.1 < -1.5) beta = -1.5;
		else beta -= 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_LEFT:
		alpha += 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_RIGHT:
		alpha -= 0.1;
		glutPostRedisplay();
		break;
	default:
		break;
	}
}

void helpscreen() {
	printf("Prima 1/2/3 para alterar o Draw Mode das figuras entre fill/line/point respetivamente.\n");
	printf("Prima X,Y e Z para ativar e desativar os respetivos eixos.\n");
	printf("Prima F1 para alterar a background color.\n");
	printf("Prima W/S para alterar a coordenada Y do parametro lookAt da camera.\n");
	printf("Prima A/D para alterar a coordenada X do parametro lookAt da camera.\n");
	printf("Prima R/T para alterar a coordenada Z do parametro lookAt da camera.\n");
	printf("Prima as diferentes setas para alterar a posicao da camera.\n");
	printf("Prima +/- para alterar a distancia entre a camera e a figura.\n");
	printf("Prima Escape para terminar o programa.\n");
}

void printGroup(Group g) {
	
	printf("--------------------------GROUP------------------------\n");
	
	printf("Pontos:\n");
	vector<Material>::iterator itmat = g.material.begin();
	for (vector<vector<Ponto> >::iterator it = g.models.begin(); it != g.models.end(); it++) {
		printf("--------------------fig-------------\n");
		for (vector<Ponto>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
			printf("x:%f, y:%f, z:%f\n", it2->x, it2->y, it2->z);
			break;
		}
		itmat++;
	}
	for(vector<int>::iterator it = g.numPontos.begin(); it != g.numPontos.end(); it++)
		printf("pontos: %d\n", (*it));

	for (vector<int>::iterator it = g.numNormals.begin(); it != g.numNormals.end(); it++)
		printf("normais: %d\n", (*it));

	for (vector<GLuint>::iterator it = g.texture.begin(); it != g.texture.end(); it++)
		printf("ID: %d\n", *it);

	for (vector<Light>::iterator it = g.lights.begin(); it != g.lights.end();it++) {
		printf("light!!!!!!!!!!!!!!!!\n");
		printf("id:%d\n", it->id);
		printf("type: %s\n", it->type.c_str());
		printf("pos:%f\n", it->pos[3]);
		if(it->color.diffuse) printf("diffuse\n");
		if (it->color.specular) printf("specular\n");
		if (it->color.emissive) printf("emissive\n");
		if (it->color.ambient) printf("ambient\n");
		if (it->cons) printf("cons\n");
		if (it->linear) printf("linear\n");
		if (it->quad) printf("quad\n");
		if (it->spotCutOff) printf("spotcutoof\n");
		if (it->spotExponent) printf("spotexponen\n");
	}

	for(vector<Group>::iterator it3 = g.nodes.begin(); it3 != g.nodes.end(); it3++)
		printGroup(*it3);
}

void initGL() {
	// alguns settings para OpenGL
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);

	glClearColor(0, 0, 0, 0);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
}

int main(int argc, char **argv) {
	// init GLUT and the window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(200, 0);
	glutInitWindowSize(800, 800);
	glutCreateWindow("Engine");

	XMLDocument doc;
	XMLError fich = doc.LoadFile("../scene.xml"); //abre ficheiro XML com as figuras a desenhar
	if (fich != XML_SUCCESS) {
		printf("Erro no ficheiro XML\n");
		return -1;
	}
	XMLElement* ptr = doc.FirstChildElement("scene"); //proot é um apontador para o filho do node scene(<group>)
	if (ptr == nullptr) return -1;
	group = load(ptr);
	//printGroup(group);
	helpscreen();
	// Required callback registry 
	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);
	glutReshapeFunc(changeSize);

	#ifndef __APPLE__
		glewInit();
	#endif
	initGL();

	// Callback registration for keyboard processing
	glutKeyboardFunc(keyboardEvents);
	glutSpecialFunc(specialEvents);

	// enter GLUT's main cycle
	glutMainLoop();

	return 1;
}


