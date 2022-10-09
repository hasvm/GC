#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include "Ponto.h"
#include "PontoText.h"

using namespace std;
unsigned int *patches;
float *controlPoints;
int patchesNum;
int controlPointsNum;

void plane(float largura, float comprimento, string filename) {
	FILE *f;
	f = fopen(filename.c_str(), "w");

	if (f != NULL) {

		vector<Ponto> pontos;

		for (float i = -(largura / 2); i <= (largura / 2); i++) {
			for (float j = -(comprimento / 2); j <= (comprimento / 2); j++) {
				pontos.push_back(Ponto(j, 0.0, i));
				pontos.push_back(Ponto(j, 0.0, i + 1));
			}
		}

		fprintf(f, "%lu\n", pontos.size() * 3); //imprimir nº de pontos no ficheiro.
		for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
			fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
		}
	}

	fclose(f);
}

void box(float x, float y, float z, int divs, string filename) {
	FILE *f;
	f = fopen(filename.c_str(), "w");

	if (f != NULL) {
		vector<Ponto> pontos;
		float px, py, pz;

		// face frontal
		for (px = -x / 2; px < x / 2; px += x / divs) {
			for (py = y / 2; py > -y / 2; py -= y / divs) {
				pontos.push_back(Ponto(px, py, z / 2));
				pontos.push_back(Ponto(px, py - (y / divs), z / 2));
				pontos.push_back(Ponto(px + (x / divs), py, z / 2));

				pontos.push_back(Ponto(px, py - (y / divs), z / 2));
				pontos.push_back(Ponto(px + (x / divs), py - (y / divs), z / 2));
				pontos.push_back(Ponto(px + (x / divs), py, z / 2));
			}
		}

		//face traseira
		for (px = -x / 2; px < x / 2; px += x / divs) {
			for (float py = y / 2; py > -y / 2; py -= y / divs) {
				pontos.push_back(Ponto(px, py, -z / 2));
				pontos.push_back(Ponto(px + (x / divs), py, -z / 2));
				pontos.push_back(Ponto(px, py - (y / divs), -z / 2));

				pontos.push_back(Ponto(px + (x / divs), py, -z / 2));
				pontos.push_back(Ponto(px + (x / divs), py - (y / divs), -z / 2));
				pontos.push_back(Ponto(px, py - (y / divs), -z / 2));
			}
		}

		//face superior
		for (pz = z / 2; pz > -z / 2; pz -= z / divs) {
			for (px = -x / 2; px < x / 2; px += x / divs) {
				pontos.push_back(Ponto(px, y / 2, pz));
				pontos.push_back(Ponto(px + (x / divs), y / 2, pz));
				pontos.push_back(Ponto(px, y / 2, pz - (z / divs)));

				pontos.push_back(Ponto(px, y / 2, pz - (z / divs)));
				pontos.push_back(Ponto(px + (x / divs), y / 2, pz));
				pontos.push_back(Ponto(px + (x / divs), y / 2, pz - (z / divs)));
			}
		}

		//face inferior
		for (px = -x / 2; px < x / 2; px += x / divs) {
			for (float pz = -z / 2; pz < z / 2; pz += z / divs) {
				pontos.push_back(Ponto(px, -y / 2, pz));
				pontos.push_back(Ponto(px + (x / divs), -y / 2, pz));
				pontos.push_back(Ponto(px, -y / 2, pz + (z / divs)));

				pontos.push_back(Ponto(px + (x / divs), -y / 2, pz));
				pontos.push_back(Ponto(px + (x / divs), -y / 2, pz + (z / divs)));
				pontos.push_back(Ponto(px, -y / 2, pz + (z / divs)));
			}
		}
		//face lateral esquerda
		for (pz = z / 2; pz > -z / 2; pz -= z / divs) {
			for (py = y / 2; py > -y / 2; py -= y / divs) {
				pontos.push_back(Ponto(-x / 2, py, pz));
				pontos.push_back(Ponto(-x / 2, py, pz - (z / divs)));
				pontos.push_back(Ponto(-x / 2, py - (y / divs), pz));

				pontos.push_back(Ponto(-x / 2, py, pz - (z / divs)));
				pontos.push_back(Ponto(-x / 2, py - (y / divs), pz - (z / divs)));
				pontos.push_back(Ponto(-x / 2, py - (y / divs), pz));
			}
		}

		//face lateral direita
		for (py = y / 2; py > -y / 2; py -= y / divs) {
			for (pz = z / 2; pz > -z / 2; pz -= z / divs) {
				pontos.push_back(Ponto(x / 2, py, pz));
				pontos.push_back(Ponto(x / 2, py - (y / divs), pz));
				pontos.push_back(Ponto(x / 2, py, pz - (z / divs)));

				pontos.push_back(Ponto(x / 2, py, pz - (z / divs)));
				pontos.push_back(Ponto(x / 2, py - (y / divs), pz));
				pontos.push_back(Ponto(x / 2, py - (y / divs), pz - (z / divs)));
			}
		}

		for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
			fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
		}
	}
}

// A funcionar para VBOs e a gerar normais e pontos de textura
void sphere(float raio, int slices, int stacks, string filename) {
	FILE *f;
	f = fopen(filename.c_str(), "w");
	float x, y, z;

	if (f != NULL) {

		vector<Ponto> pontos;
		vector<Ponto> normais;
		vector<PontoText> texturas;
		float alpha = 2 * M_PI / slices;
		float beta = M_PI / stacks;
		float itStacks = 1.0f / ((float)stacks);
		float itSlices = 1.0f / ((float)slices);
		float cSlices = 0;
		float cStacks = 0;

		//metade superior da esfera
		for (int stack = stacks / 2; stack >= 0; stack--) {
			for (int slice = 0; slice <= slices; slice++) {
				pontos.push_back(Ponto(raio * cosf(stack*beta) * sinf(slice*alpha), raio * sinf(stack*beta), raio * cosf(stack*beta) * cosf(slice*alpha)));
				pontos.push_back(Ponto(raio * cosf((stack - 1)*beta) * sinf(slice*alpha), raio * sinf((stack - 1)*beta), raio * cosf((stack - 1)*beta) * cosf(slice*alpha)));

				normais.push_back(Ponto(cosf(stack*beta) * sinf(slice*alpha), sinf(stack*beta), cosf(stack*beta) * cosf(slice*alpha)));
				normais.push_back(Ponto(cosf((stack - 1)*beta) * sinf(slice*alpha), sinf((stack - 1)*beta), cosf((stack - 1)*beta) * cosf(slice*alpha)));

				texturas.push_back(PontoText(cSlices, cStacks));
				texturas.push_back(PontoText(cSlices, cStacks - itStacks));

				cSlices = cSlices + itSlices;
			}
			cStacks = cStacks - itStacks;
			cSlices = 0.0;
		}

		cSlices = 0;
		cStacks = 0.5;

		//metade inferior da esfera
		for (int stack = 0; stack <= stacks / 2; stack++) {
			for (int slice = 0; slice <= slices; slice++) {
				pontos.push_back(Ponto(raio * cosf(-stack * beta) * sinf(slice*alpha), raio * sinf(-stack * beta), raio * cosf(-stack * beta) * cosf(slice*alpha)));
				pontos.push_back(Ponto(raio * cosf(-(stack + 1)*beta) * sinf(slice*alpha), raio * sinf(-(stack + 1)*beta), raio * cosf(-(stack + 1)*beta) * cosf(slice*alpha)));

				normais.push_back(Ponto(cosf(-stack * beta) * sinf(slice*alpha), sinf(-stack * beta), cosf(-stack * beta) * cosf(slice*alpha)));
				normais.push_back(Ponto(cosf(-(stack + 1)*beta) * sinf(slice*alpha), sinf(-(stack + 1)*beta), cosf(-(stack + 1)*beta) * cosf(slice*alpha)));

				texturas.push_back(PontoText(cSlices, cStacks));
				texturas.push_back(PontoText(cSlices, cStacks - itStacks));

				cSlices = cSlices + itSlices;
			}
			cStacks = cStacks - itStacks;
			cSlices = 0.0;
		}
		//imprime pontos
		fprintf(f, "%lu\n", pontos.size() * 3);
		for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
			fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
		}
		for (unsigned int normal = 0; normal < normais.size(); normal++) {
			fprintf(f, "%f %f %f\n", normais[normal].getX(), normais[normal].getY(), normais[normal].getZ());
		}
		for (unsigned int textura = 0; textura < texturas.size(); textura++) {
			fprintf(f, "%f %f\n", texturas[textura].getX(), texturas[textura].getY());
		}
	}
}

void cone(float radius, float height, int slices, int stacks, string filename) {
	FILE *f;
	f = fopen(filename.c_str(), "w");

	if (f != NULL) {
		vector<Ponto> pontos;
		int slice, stack;

		float alpha = 2 * M_PI / slices;
		//		float beta = M_PI / stacks;

				// base
		for (slice = 0; slice < slices; slice++) {

			pontos.push_back(Ponto(radius * sin(alpha * slice), 0, radius * cos(alpha *slice)));
			pontos.push_back(Ponto(0, 0, 0));
			pontos.push_back(Ponto(radius * sin(alpha * (slice + 1)), 0, radius * cos(alpha *(slice + 1))));

			float novoRaio, alturaStack, nextStack, nextRaio;
			//lados do cone
			for (stack = 0; stack < stacks; stack++) {
				alturaStack = (height / stacks) * stack;
				novoRaio = (height - alturaStack) * radius / height;
				nextStack = (height / stacks) * (stack + 1);
				nextRaio = (height - nextStack) * radius / height;


				//triangulo cima
				pontos.push_back(Ponto(nextRaio * sin(alpha*slice), nextStack, nextRaio * cos(alpha*slice)));
				pontos.push_back(Ponto(novoRaio * sin(alpha*slice), alturaStack, novoRaio * cos(alpha*slice)));
				pontos.push_back(Ponto(nextRaio * sin(alpha*(slice + 1)), nextStack, nextRaio * cos(alpha*(slice + 1))));
				// traingulo baixo
				pontos.push_back(Ponto(novoRaio * sin(alpha*slice), alturaStack, novoRaio * cos(alpha*slice)));
				pontos.push_back(Ponto(novoRaio * sin(alpha*(slice + 1)), alturaStack, novoRaio * cos(alpha*(slice + 1))));
				pontos.push_back(Ponto(nextRaio * sin(alpha *(slice + 1)), nextStack, nextRaio * cos(alpha *(slice + 1))));

			}

			//ponta do cone
			pontos.push_back(Ponto(0.0, height, 0.0));
			pontos.push_back(Ponto(novoRaio * sin(alpha*slice), alturaStack, novoRaio * cos(alpha*slice)));
			pontos.push_back(Ponto(novoRaio * sin(alpha *(slice + 1)), alturaStack, novoRaio * cos(alpha *(slice + 1))));

		}

		for (unsigned int ponto = 0; ponto < pontos.size() - 1; ponto++) {
			fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
		}
	}
}

// A funcionar para VBOs e a gerar normais e pontos de textura
void anel(float raioint, float raioext, int slices, string filename) {
	FILE *f;
	f = fopen(filename.c_str(), "w");

	if (f != NULL) {
		// vetor onde os pontos do plano vão ser guardados
		vector<Ponto> pontos1, pontos2;
		vector<Ponto> normais;
		float x, z;
		// numero de divisões do circulo
		float alpha = 2 * M_PI / slices;

		// circunferência interior
		for (float div = 0; div <= 2 * M_PI + alpha; div += alpha) {
			x = raioint * sin(div);
			z = raioint * cos(div);
			pontos1.push_back(Ponto(x, 0, z));
		}
		// circunferência exterior
		for (float div = 0; div <= 2 * M_PI + alpha; div += alpha) {
			x = raioext * sin(div);
			z = raioext * cos(div);
			pontos2.push_back(Ponto(x, 0, z));

		}

		// printa o numeros de pontos
		fprintf(f, "%lu\n", (pontos1.size() + pontos2.size()) * 6);
		// PONTOS
		// parte de baixo do anel
		for (unsigned int i = 0; i < pontos1.size(); i++) {
			fprintf(f, "%f %f %f\n", pontos1[i].getX(), pontos1[i].getY(), pontos1[i].getZ());
			fprintf(f, "%f %f %f\n", pontos2[i].getX(), pontos2[i].getY(), pontos2[i].getZ());
		}
		//parte de cima do anel
		for (unsigned int i = 0; i < pontos1.size(); i++) {
			fprintf(f, "%f %f %f\n", pontos2[i].getX(), pontos2[i].getY(), pontos2[i].getZ());
			fprintf(f, "%f %f %f\n", pontos1[i].getX(), pontos1[i].getY(), pontos1[i].getZ());
		}
		// NORMAIS
		for (unsigned int i = 0; i < pontos1.size(); i++) {
			fprintf(f, "%f %f %f\n", pontos1[i].getX(), 1.0f, pontos1[i].getZ());
			fprintf(f, "%f %f %f\n", pontos2[i].getX(), 1.0f, pontos2[i].getZ());
		}
		for (unsigned int i = 0; i < pontos1.size(); i++) {
			fprintf(f, "%f %f %f\n", pontos1[i].getX(), -1.0f, pontos1[i].getZ());
			fprintf(f, "%f %f %f\n", pontos2[i].getX(), -1.0f, pontos2[i].getZ());
		}
		// PONTOS DE TEXTURA
		for (unsigned int i = 0; i < pontos1.size(); i++) {
			fprintf(f, "%f %f\n", 0.1f, 0.5f);
			fprintf(f, "%f %f\n", 0.9f, 0.5f);
		}
		for (unsigned int i = 0; i < pontos1.size(); i++) {
			fprintf(f, "%f %f\n", 0.1f, 0.5f);
			fprintf(f, "%f %f\n", 0.9f, 0.5f);
		}
	}
}

void circumference(float raio, int slices, string filename) {
	FILE *f;
	f = fopen(filename.c_str(), "w");

	if (f != NULL) {
		// vetor onde os pontos do plano vão ser guardados
		vector<Ponto> pontos;
		float x, z;
		// numero de divisões do circulo
		float alpha = 2 * M_PI / slices;

		// pontos da circunferência
		for (float div = 0; div <= 2 * M_PI + alpha; div += alpha) {
			x = raio * sin(div);
			z = raio * cos(div);
			pontos.push_back(Ponto(x, 0, z));
		}

		fprintf(f, "%lu\n", pontos.size() * 3);
		for (unsigned int i = 0; i < pontos.size(); i++) {
			fprintf(f, "%f %f %f\n", pontos[i].getX(), pontos[i].getY(), pontos[i].getZ());
		}
	}
}

void patchToArrays(string filename) {

	ifstream file(filename.c_str());
	string firstLine, line;

	// Colocar valor dos indices de cada patch no array.
	getline(file, firstLine);
	patchesNum = atoi(firstLine.c_str());
	patches = (unsigned int*)malloc(sizeof(unsigned int) * 16 * patchesNum);
	for (int p = 0; p < patchesNum; p++) {
		getline(file, line);
		istringstream indexes(line);
		string indexP;
		for (int i = 0; i < 16 && getline(indexes, indexP, ','); i++) { //getline(char, streamsize, delimitador)
			patches[p * 16 + i] = atoi(indexP.c_str());
		}
	}

	// Colocar valor dos pontos de controlo no array.
	getline(file, firstLine);
	controlPointsNum = atoi(firstLine.c_str());
	controlPoints = (float *)malloc(sizeof(float) * 3 * controlPointsNum);
	for (int cp = 0; cp < controlPointsNum; cp++) {
		getline(file, line);
		istringstream indexes(line);
		string indexCP;
		for (int i = 0; i < 3 && getline(indexes, indexCP, ','); i++) {
			controlPoints[cp * 3 + i] = (float)atof(indexCP.c_str());
		}
	}
}

Ponto pontoBezier(int p, float u, float v) {
	Ponto ponto = Ponto(0.0, 0.0, 0.0);

	// Polinomio de Bernstein
	float bernsteinU[4] = { powf(1 - u, 3), 3 * u * powf(1 - u, 2), 3 * powf(u, 2) * (1 - u), powf(u, 3) };
	float bernsteinV[4] = { powf(1 - v, 3), 3 * v * powf(1 - v, 2), 3 * powf(v, 2) * (1 - v), powf(v, 3) };
	//powf é a potencia (de floats) sendo que o 1º arg é a base e o 2º é o expoente.

	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < 4; i++) {

			int indexPatch = j * 4 + i;
			int indexCP = patches[p * 16 + indexPatch];
			ponto = Ponto(ponto.getX() + controlPoints[indexCP * 3 + 0] * bernsteinU[j] * bernsteinV[i],
				ponto.getY() + controlPoints[indexCP * 3 + 1] * bernsteinU[j] * bernsteinV[i],
				ponto.getZ() + controlPoints[indexCP * 3 + 2] * bernsteinU[j] * bernsteinV[i]);
		}
	}

	return ponto;
}

Ponto pontoBezierNormal(int p, float u, float v) {
	Ponto ponto = Ponto(0.0, 0.0, 0.0);

	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < 4; i++) {

			int indexPatch = j * 4 + i;
			int indexCP = patches[p * 16 + indexPatch];
			ponto = Ponto(ponto.getX() + controlPoints[indexCP * 3 + 0],
				ponto.getY() + controlPoints[indexCP * 3 + 1],
				ponto.getZ() + controlPoints[indexCP * 3 + 2]);
		}
	}
	return ponto;
}

// A funcionar para Triangles
void bezier(string patchFile, string filename, int tesselation) {
	patchToArrays(patchFile);
	vector<Ponto> pontos;
	vector<Ponto> normais;
	vector<PontoText> texturas;
	float inc = 1.0f / tesselation;

	for (int p = 0; p < patchesNum; p++) {

		for (int tv = 0; tv < tesselation; tv++) {
			float v = (float)tv / tesselation;

			for (int tu = 0; tu < tesselation; tu++) {
				float u = (float)tu / tesselation;

				// triângulo superior
				pontos.push_back(pontoBezier(p, (u + inc), (v + inc)));
				pontos.push_back(pontoBezier(p, u, v));
				pontos.push_back(pontoBezier(p, u, (v + inc)));

				// triângulo inferior
				pontos.push_back(pontoBezier(p, (u + inc), (v + inc)));
				pontos.push_back(pontoBezier(p, (u + inc), v));
				pontos.push_back(pontoBezier(p, u, v));

				normais.push_back(pontoBezierNormal(p, (u + inc), (v + inc)));
				normais.push_back(pontoBezierNormal(p, u, v));
				normais.push_back(pontoBezierNormal(p, (u + inc), v));

				normais.push_back(pontoBezierNormal(p, (u + inc), (v + inc)));
				normais.push_back(pontoBezierNormal(p, u, (v + inc)));
				normais.push_back(pontoBezierNormal(p, u, v));

				texturas.push_back(PontoText((u + inc), (v + inc)));
				texturas.push_back(PontoText(u, v));
				texturas.push_back(PontoText(u, (v + inc)));

				texturas.push_back(PontoText((u + inc), (v + inc)));
				texturas.push_back(PontoText((u + inc), v));
				texturas.push_back(PontoText(u, v));
			}
		}
	}

	FILE *f;
	f = fopen(filename.c_str(), "w");

	// Printar os pontos no ficheiro .3d
	for (int ponto = 0; ponto < pontos.size(); ponto++) {
		fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
	}
	for (int ponto = 0; ponto < normais.size(); ponto++) {
		fprintf(f, "%f %f %f\n", normais[ponto].getX(), normais[ponto].getY(), normais[ponto].getZ());
	}
	for (int ponto = 0; ponto < texturas.size(); ponto++) {
		fprintf(f, "%f %f\n", texturas[ponto].getX(), texturas[ponto].getY());
	}

	fclose(f);

	free(patches);
	free(controlPoints);
}

int main(int argc,const char* argv[]) {

	argc = 8;
	argv[1] = "anel"; argv[2] = "1.1";argv[3] = "1.3";argv[4] = "100";argv[5] = "anel.3d";

	/*
	argc = 5;
	argv[1] = "circumference"; argv[2] = "1";argv[3] = "50";argv[4] = "orbit.3d";
	
	argc = 5;
	argv[1] = "bezier"; argv[2] = "teapot.patch"; argv[3] = "teapot.3d"; argv[4] = "5";

    argc=8;
    argv[1]="anel"; argv[2]="1.1";argv[3]="1.3";argv[4]="100";argv[5]="anel.3d";
	
	//sphere - raidus - slices - stacks - filename
	argc = 6;
	argv[1] = "sphere"; argv[2] = "1";argv[3] = "50";argv[4] = "50";argv[5] = "sphere.3d";
	*/

	if (argc < 2) {
		printf("Parametros invalidos\n");
		return -1;
	}
	string tipo = argv[1];
	string path = "../../3DFiles/";

	if (tipo == "plane" && argc == 5) {
		float largura = stof(argv[2]);
		float comprimento = stof(argv[3]);
		string filename = argv[4];
		filename = path + filename;
		plane(largura, comprimento, filename);
		return 0;
	}
	
	if (tipo == "box" && argc == 7) {
		float x = stof(argv[2]);
		float y = stof(argv[3]);
		float z = stof(argv[4]);
		int divs = stoi(argv[5]);
		string filename = argv[6];
		filename = path + filename;
		box(x, y, z, divs, filename);
		return 0;
	}

	if (tipo == "sphere" && argc == 6) {
		float radius = stof(argv[2]);
		int slices = stoi(argv[3]);
		int stacks = stoi(argv[4]);
		string filename = argv[5];
		filename = path + filename;
		sphere(radius, slices, stacks, filename);
		return 0;
	}

	if (tipo == "cone" && argc == 7) {
		float radius = stof(argv[2]);
		float height = stof(argv[3]);
		int slices = stoi(argv[4]);
		int stacks = stoi(argv[5]);
		string filename = argv[6];
		filename = path + filename;
		cone(radius, height, slices, stacks, filename);
		return 0;
	}

	if (tipo == "circumference" && argc == 5) {
		float radius = stof(argv[2]);
		int slices = stoi(argv[3]);
		string filename = argv[4];
		filename = path + filename;
		circumference(radius, slices, filename);
		return 0;
	}
    
    if (tipo == "anel" && argc == 8) {
        float raioInt=stof(argv[2]);
        float raioExt=stof(argv[3]);
        int slices=stof(argv[4]);
        string filename=argv[5];
 
        filename = path + filename;
        anel(raioInt,raioExt, slices, filename);
        return 0;
    }

	if (tipo == "bezier" && argc == 5) {
		string patchFile = argv[2];
		string filename = argv[3];
		int tesselation = stoi(argv[4]);
		patchFile = path + patchFile;
		filename = path + filename;

		bezier(patchFile, filename, tesselation);

		return 0;
	}

	printf("Parametros Invalidos\n");
	return -1;
}
