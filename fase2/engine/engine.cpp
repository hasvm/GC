#define GL_SILENCE_DEPRECATION

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
#include <GL/glut.h>
#endif

using namespace std;
using namespace tinyxml2;

class Transform {
public:
	float x, y, z;
	float angle;
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
public:
	Rotate(float a, float b, float c, float ag) : Transform(a, b, c) {		
		float angle = ag;
	}
};

class Translate : public Transform {
public:
	Translate(float a, float b, float c) : Transform(a, b, c) {
	}
};

class Scale : public Transform {
public:
	Scale(float a, float b, float c) : Transform(a, b, c) {
	}
};

class Group {
public:
	vector<vector<Ponto>> models;
	vector<vector<Ponto>> orbits;
	vector<Transform*> transforms;
	vector<Ponto> colors;
	vector<Group> nodes;
};

Group group; //group das figuras e transformações a desenhar

//eixos
int xflag = 1;
int yflag = 1;
int zflag = 1;
//coordenadas para a camera 
float r = 200;
float px, py, pz = 0;
float dx = 5, dy = 0, dz = -5;
float alpha = M_PI / 6;
float beta = M_PI / 6;

GLenum DRAWING_MODE = GL_FILL;

Group load(XMLElement* node) {
	string path = "../../3DFiles/"; //localização dos ficheiros 3D para desenhar
	XMLElement* gptr = (node->FirstChildElement()); //apontador para o filho do node group(group pointer)
	Group g;
//	int i = 1;
	for (; gptr != nullptr; gptr = gptr->NextSiblingElement()) {
		string tag = gptr->Value();
		if (strcmp(tag.c_str(), "translate") == 0) {
			g.transforms.push_back(new Translate(gptr->DoubleAttribute("X"), gptr->DoubleAttribute("Y"), gptr->DoubleAttribute("Z")));
		}
		else if (strcmp(tag.c_str(), "rotate") == 0) {
			printf("r\n");
			g.transforms.push_back(new Rotate(gptr->DoubleAttribute("X"), gptr->DoubleAttribute("Y"),
				gptr->DoubleAttribute("Z"), gptr->DoubleAttribute("angle")));
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
				string newFigure = path + mptr->Attribute("file");
				ifstream file;
				file.open(newFigure);
				vector<Ponto> modelpoints;
				while (!file.eof()) { //percorre até chegar ao end of file
					Ponto ponto;
					file >> ponto.x >> ponto.y >> ponto.z;
					modelpoints.push_back(ponto); //guarda todos os pontos de uma figura num vector de pontos(modelpoints)
				}
				if (strcmp(gptr->Attribute("type"), "figure") == 0) 
					g.models.push_back(modelpoints); //adicion	a o vector com os pontos de uma figura ao vector de figuras
				if (strcmp(gptr->Attribute("type"), "orbit") == 0) 
					g.orbits.push_back(modelpoints); //adiciona o vector com os pontos de uma figura ao vector de orbitas
			}
		}
		else if (strcmp(tag.c_str(), "group") == 0) {
			g.nodes.push_back(load(gptr));
		}
	}
	return g;
}

void display(Group g) {
	glPushMatrix();
	vector<Ponto>::iterator itcolor = g.colors.begin(); //iterador para percorrer o vector com as cores de cada modelo

	//Transformações
	for (vector<Transform*>::iterator it = g.transforms.begin(); it != g.transforms.end(); it++) {
		if (Rotate *r = dynamic_cast<Rotate*>(*it))
			glRotatef(r->angle, r->x, r->y, r->z);
		if (Translate *t = dynamic_cast<Translate*>(*it))
			glTranslatef(t->x, t->y, t->z);
		if (Scale *s = dynamic_cast<Scale*>(*it))
			glScalef(s->x, s->y, s->z);
	}

	//Modelos
	glBegin(GL_TRIANGLES);
	for (vector<vector<Ponto> >::iterator it = g.models.begin(); it != g.models.end(); it++) {
		glColor3f(itcolor->x,itcolor->y,itcolor->z);
		for (vector<Ponto>::iterator it2 = it->begin(); it2 != (it->end())-1; it2++) 
			glVertex3f(it2->x, it2->y, it2->z); // vertice(x,y,z)
		itcolor++; //para percorrer todo o vector de cores
	}
	glEnd();
 
	//Orbitas
	glColor3f(0.4f, 0.4f, 0.4f);
	for (vector<vector<Ponto> >::iterator it = g.orbits.begin(); it != g.orbits.end(); it++) {
		glBegin(GL_LINE_STRIP);
		for (vector<Ponto>::iterator it2 = it->begin(); it2 != (it->end())-1; it2++) {
			glVertex3f(it2->x, it2->y, it2->z); // vertice(x,y,z)
		}
		glEnd();
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
	printf("Pontos:\n");
	for (vector<vector<Ponto> >::iterator it = g.models.begin(); it != g.models.end(); it++)
		for (vector<Ponto>::iterator it2 = it->begin(); it2 != it->end(); it2++)
			printf("x:%f, y:%f, z:%f\n", it2->x, it2->y, it2->z);

	printf("Orbitas:\n");
	for (vector<vector<Ponto> >::iterator it = g.orbits.begin(); it != g.orbits.end(); it++)
		for (vector<Ponto>::iterator it2 = it->begin(); it2 != it->end(); it2++)
			printf("x:%f, y:%f, z:%f\n", it2->x, it2->y, it2->z);

	printf("Transformacoes:\n");
	for (vector<Transform*>::iterator it2 = g.transforms.begin(); it2 != g.transforms.end(); it2++) {
		if (Rotate *r = dynamic_cast<Rotate*>(*it2))
			printf("Rotate: angle:%f, x:%f, y:%f, z:%f",r->angle,r->x,r->y,r->z);
		if (Translate *t = dynamic_cast<Translate*>(*it2))
			printf("Translate: x:%f, y:%f, z:%f", t->x, t->y, t->z);
		if (Scale *s = dynamic_cast<Scale*>(*it2))
			printf("Scale: x:%f, y:%f, z:%f", s->x, s->y, s->z);
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
	helpscreen();

	// Required callback registry 
	glutDisplayFunc(renderScene);
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


