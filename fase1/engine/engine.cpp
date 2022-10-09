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

class Figura {
	public:
		vector<Ponto> fig; //vector com o vector dos pontos de uma figura
};

vector<Figura> figuras; //vector global das figuras que serão desenhadas

//eixos
int xflag = 1;
int yflag = 1;
int zflag = 1;
//coordenadas para a camera 
float r = 10;
float px, py, pz = 0;
float dx, dy, dz = 0;
float alpha = M_PI / 6;
float beta = M_PI / 6;

GLenum DRAWING_MODE = GL_FILL;

void load() {
	vector<string> figuresToLoad; //vector com os nomes das figuras presentes no ficheiro XML para serem desenhadas
	string path = "../../3DFiles/"; //localização dos ficheiros 3D para desenhar

	XMLDocument doc;
	XMLError load = doc.LoadFile("../scene.xml"); //abre ficheiro XML com as figuras a desenhar
	if (load != XML_SUCCESS) {
		printf("Erro no ficheiro XML\n");
		return;
	}
	
	XMLNode *pRoot = doc.FirstChildElement("scene");
	if (pRoot == nullptr) return;
	XMLElement *sceneFigures = pRoot->FirstChildElement("model");
	for (; sceneFigures != nullptr; sceneFigures = sceneFigures->NextSiblingElement("model")) {
		string newFigure = path + sceneFigures->Attribute("file");
		figuresToLoad.push_back(newFigure);
	}

	for (vector<string>::iterator it = figuresToLoad.begin(); it != figuresToLoad.end(); ++it) {
		Figura newFig;
		ifstream file;
		file.open(*it);
	
		while (!file.eof()) { //percorre até chegar ao end of file
			Ponto ponto;
			file >> ponto.x >> ponto.y >> ponto.z;
			newFig.fig.push_back(ponto); //coloca as coordenadas de um ponto no vector dos pontos da figura NewFig
		}
		
		figuras.push_back(newFig); //coloca no vector global uma Figura(newFig) para ser desenhada
	}
	
}

void display() {
	vector<Figura>::iterator it_figuras;
	int color = 0;
	int triangulo = 0;
	GLfloat colors[][3] = {{1.0f,0.0f,0.0f} ,{0.0f, 1.0f, 0.0f },{0.0f, 0.0f, 1.0f },{0.0f, 1.0f, 1.0f },
	{1.0f, 0.0f, 1.0f }, {1.1f,1.0f,0.0f } };

	for (it_figuras = figuras.begin(); it_figuras != figuras.end(); it_figuras++) {
		vector<Ponto> coord = it_figuras->fig; //coordenadas dos pontos de uma figura
		vector<Ponto>::iterator it_pontos; //iterador para percorrer os pontos de uma figura 

		glBegin(GL_TRIANGLES);
		//glColor3f(1.0, 0.0, 0.0);
		for (it_pontos = coord.begin(); it_pontos != coord.end(); it_pontos++,triangulo++) {
			color = rand() % 7;
			glColor3f(colors[color][0], colors[color][1], colors[color][2]); //pinta cada vértice de uma cor escolhida aleatoriamente do array colors
			glVertex3f(it_pontos->x, it_pontos->y, it_pontos->z); // vertice(x,y,z)
		}
		glEnd();
	}
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
		glBegin(GL_LINES);
		glColor3f(0.0, 0.0, 1.0);
		//eixo x
		if (xflag) {
			glVertex3f(0, 0, 0);
			glVertex3f(10, 0, 0);
		}
		//eixo y
		if (yflag) {
			glVertex3f(0, 0, 0);
			glVertex3f(0, 10, 0);
		}
		//eixo z
		if (zflag) {
			glVertex3f(0, 0, 0);
			glVertex3f(0, 0, 10);
		}
		glEnd();
	}

	glPolygonMode(GL_FRONT_AND_BACK,DRAWING_MODE); //forma como as figuras vão ser desenhadas
	load(); //carrega as figuras
	display(); //desenha as figuras
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
		dy -= 0.1;
		glutPostRedisplay();
		break;
	case 's':
		dy += 0.1;
		glutPostRedisplay();
		break;
	case 'a':
		dx += 0.1;
		glutPostRedisplay();
		break;
	case 'd':
		dx -= 0.1;
		glutPostRedisplay();
		break;
	case 'r':
		dz += 0.1;
		glutPostRedisplay();
		break;
	case 't':
		dz -= 0.1;
		glutPostRedisplay();
		break;
	case '+':
		r -= 0.2f;
		glutPostRedisplay();
		break;
	case '-':
		r += 0.2f;
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

int main(int argc, char **argv) {
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


