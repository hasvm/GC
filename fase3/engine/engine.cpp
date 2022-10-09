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
#include "tinyxml2.cpp"
#include "tinyxml2.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

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

class Group {
public:
	vector<vector<Ponto>> models;
	vector<Transform*> transforms; //scales/rotates
	Translate translate; //translate
	vector<Ponto> colors;	
	vector<int> numPontos; //vector com o numero de pontos de cada modelo(se for 0 o modelo é desenhado sem VBOs)
	vector<Group> nodes;
};

Group group; //group das figuras e transformações a desenhar

//eixos
int xflag = 1;
int yflag = 1;
int zflag = 1;
//coordenadas para a camera 
float r = 17;
float px, py, pz = 0;
float dx = 5, dy = 0, dz = -5;
float gt = 0.0, clocks = 0.0;
float alpha = M_PI / 6;
float beta = M_PI / 6;

GLenum DRAWING_MODE = GL_FILL;
GLuint buffers[1];

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

Group load(XMLElement* node) {
	string path = "../../3DFiles/"; //localização dos ficheiros 3D para desenhar
	XMLElement* gptr = (node->FirstChildElement()); //apontador para o filho do node group(group pointer)
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
			XMLElement* mptr = gptr->FirstChildElement(); //apontador para os filhos do models
			for (; mptr != nullptr; mptr = mptr->NextSiblingElement()) {
				int numSpaces, nVertices;
				string newFigure = path + mptr->Attribute("file");
				string line;
				ifstream file;
				file.open(newFigure);
				vector<Ponto> modelpoints;
				if (getline(file, line)) {
					file.close();
					numSpaces = count(line.begin(), line.end(), ' ');
					//VBO
					if (numSpaces == 0) {
						ifstream fileVBO;
						fileVBO.open(newFigure);
						fileVBO >> nVertices;
						/* Guarda os vértices */
						while (!fileVBO.eof()) {
							Ponto newC;
							fileVBO >> newC.x >> newC.y >> newC.z;
							modelpoints.push_back(newC);
						}
						fileVBO.close();
						g.numPontos.push_back(nVertices);
						g.models.push_back(modelpoints);
					}
					//Triangles
					else {
						ifstream fileTriang;
						fileTriang.open(newFigure);
						while (!fileTriang.eof()) { //percorre até chegar ao end of file
							Ponto ponto;
							fileTriang >> ponto.x >> ponto.y >> ponto.z;
							modelpoints.push_back(ponto); //guarda todos os pontos de uma figura num vector de pontos(modelpoints)
						}
						g.numPontos.push_back(0);
						g.models.push_back(modelpoints);
					}
				}
			}
		}
		else if (strcmp(tag.c_str(), "group") == 0) {
			g.nodes.push_back(load(gptr));
		}
	}
	return g;
}

void drawFiguresVBOs(vector<Ponto> figures, int nPoints) {
	glGenBuffers(1, buffers);
	glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Ponto) * figures.size(), &figures[0], GL_STATIC_DRAW);
	glVertexPointer(3, GL_FLOAT, 0, 0);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, figures.size());
	figures.clear();
	glDeleteBuffers(1, buffers);
}

void drawFiguresTriangles(vector<Ponto> figures) {
	vector<Ponto>::iterator it_fig;
	glBegin(GL_TRIANGLES);
	for (it_fig = figures.begin(); it_fig != figures.end(); it_fig++) 
		glVertex3f(it_fig->x, it_fig->y, it_fig->z);
	glEnd();
}

void display(Group g) {
	glPushMatrix();
	vector<Ponto>::iterator itcolor = g.colors.begin(); //iterador para percorrer o vector com as cores de cada modelo
	vector<int>::iterator itnum = g.numPontos.begin();

	//Transformações
	for (vector<Transform*>::iterator it = g.transforms.begin(); it != g.transforms.end(); it++) {
		if (Rotate *r = dynamic_cast<Rotate*>(*it)) {
			renderRotate(*r);
		}
		if (Scale *s = dynamic_cast<Scale*>(*it)) {
			glScalef(s->x, s->y, s->z);
		}
	}
	if (!g.translate.empty) {
		renderCatmullRomCurve(g.translate);
		renderCatmullTranslate(g.translate);
		g.translate.catmullPoints.clear();
		g.translate.catmullPoints.shrink_to_fit();
	}
	//Modelos
	for (vector<vector<Ponto> >::iterator it = g.models.begin(); it != g.models.end(); it++) {
		glColor3f(itcolor->x,itcolor->y,itcolor->z);
		if ((*itnum) > 0) drawFiguresVBOs(*it, *itnum);//com VBOs
		else drawFiguresTriangles(*it);                //sem VBOs
		itcolor++; //para percorrer todo o vector de cores
		itnum++;   //para percorrer o vector com o numero de vertices de cada modelo
	}
    
	//Nodos do group
	for(vector<Group>::iterator itg = g.nodes.begin(); itg != g.nodes.end(); itg++) 
		display(*itg);

	glPopMatrix();
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

	px = r * cos(beta) * sin(alpha);
	py = r * sin(beta);
	pz = r * cos(beta) * cos(alpha);

	// set the camera
	glLoadIdentity();
	gluLookAt(px, py, pz,
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

	glPolygonMode(GL_FRONT_AND_BACK,DRAWING_MODE); //forma como as figuras vão ser desenhadas
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
	/*
	printf("Pontos:\n");
	for (vector<vector<Ponto> >::iterator it = g.models.begin(); it != g.models.end(); it++)
		for (vector<Ponto>::iterator it2 = it->begin(); it2 != it->end(); it2++)
			printf("x:%f, y:%f, z:%f\n", it2->x, it2->y, it2->z);
	*/
	
	printf("Numero de pontos:\n");
	for (vector<int>::iterator it = g.numPontos.begin(); it != g.numPontos.end(); it++)
		printf("pontos: %d\n",*it);
	
	printf("Transformacoes:\n");
	for (vector<Transform*>::iterator it2 = g.transforms.begin(); it2 != g.transforms.end(); it2++) {
		if (Rotate *r = dynamic_cast<Rotate*>(*it2))
			printf("Rotate: time:%f, angle:%f, x:%f, y:%f, z:%f\n",r->time,r->angle,r->x,r->y,r->z);
		if (Scale *s = dynamic_cast<Scale*>(*it2))
			printf("Scale: x:%f, y:%f, z:%f\n", s->x, s->y, s->z);
	}
	if (!g.translate.empty) {
		printf("Translate: %f\n", g.translate.time);
		for (vector<Ponto>::iterator it = g.translate.catmullPoints.begin(); it != g.translate.catmullPoints.end(); it++)
			printf("time: %f, x: %f, y: %f, z: %f\n", g.translate.time, it->x, it->y, it->z);
	}

	for(vector<Group>::iterator it3 = g.nodes.begin(); it3 != g.nodes.end(); it3++)
		printGroup(*it3);
}

int main(int argc, char **argv) {

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

	// init GLUT and the window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(200,0);
	glutInitWindowSize(800, 800);
	glutCreateWindow("Engine");
	glEnableClientState(GL_VERTEX_ARRAY);
	helpscreen();

	#ifndef __APPLE__
		glewInit();
	#endif

	// Required callback registry 
	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);
	glutReshapeFunc(changeSize);

	// Callback registration for keyboard processing
	glutKeyboardFunc(keyboardEvents);
	glutSpecialFunc(specialEvents);

	//  OpenGL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	// enter GLUT's main cycle
	glutMainLoop();

	return 1;
}


