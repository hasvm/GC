#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "Ponto.h"

using namespace std;

void plane(float l, string filename) {
	FILE *f;
	f = fopen(filename.c_str(), "w");

	if (f != NULL) {
		// vetor onde os pontos do plano vão ser guardados
		vector<Ponto> pontos;
		float lado = l / 2;

		//primeiro triangulo, virado para cima
		pontos.push_back(Ponto(lado, 0.0, lado));
		pontos.push_back(Ponto(lado, 0.0, -lado));
		pontos.push_back(Ponto(-lado, 0.0, lado));
		//segundo triangulo, virado para cima
		pontos.push_back(Ponto(-lado, 0.0, -lado));
		pontos.push_back(Ponto(-lado, 0.0, lado));
		pontos.push_back(Ponto(lado, 0.0, -lado));

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

void sphere(float radius, int slices, int stacks, string filename){
	FILE *f;
	f = fopen(filename.c_str(), "w");

	if (f != NULL) {
		std::vector<Ponto> pontos;

		float alpha = 2 * M_PI / slices;
		float beta = M_PI / stacks;

		for (int slice = 0; slice <= slices; slice++) {
			//topo da esfera
			float aux1 = (M_PI / 2) - beta;
			pontos.push_back(Ponto(0, radius, 0));
			pontos.push_back(Ponto(radius * cosf(aux1) * sinf(slice*alpha), radius * sinf(aux1), radius * cosf(aux1) * cosf(slice*alpha)));
			pontos.push_back(Ponto(radius * cosf(aux1) * sinf((slice + 1)*alpha), radius * sinf(aux1), radius * cosf(aux1) * cosf((slice + 1)*alpha)));
			//base da esfera
			float aux2 = (-M_PI / 2) + beta;
			pontos.push_back(Ponto(0, -radius, 0));
			pontos.push_back(Ponto(radius * cosf(aux2) * sinf((slice + 1)*alpha), radius * sinf(aux2), radius * cosf(aux2) * cosf((slice + 1)*alpha)));
			pontos.push_back(Ponto(radius * cosf(aux2) * sinf(slice*alpha), radius * sinf(aux2), radius * cosf(aux2) * cosf(slice*alpha)));

			//lados da esfera
			for (int stack = 0; stack < stacks / 2 - 1; stack++) {
				//metade superior
				//triangulo inferior
				pontos.push_back(Ponto(radius * cosf(stack*beta) * sinf(slice*alpha), radius * sinf(stack*beta), radius * cosf(stack*beta) * cosf(slice*alpha)));
				pontos.push_back(Ponto(radius * cosf(stack*beta) * sinf((slice + 1)*alpha), radius * sinf(stack*beta), radius * cosf(stack*beta) * cosf((slice + 1)*alpha)));
				pontos.push_back(Ponto(radius * cosf((stack + 1)*beta) * sinf((slice + 1)*alpha), radius * sinf((stack + 1)*beta), radius * cosf((stack + 1)*beta) * cosf((slice + 1)*alpha)));
				//triangulo superior
				pontos.push_back(Ponto(radius * cosf((stack + 1)*beta) * sinf(slice*alpha), radius * sinf((stack + 1)*beta), radius * cosf((stack + 1)*beta) * cosf(slice*alpha)));
				pontos.push_back(Ponto(radius * cosf(stack*beta) * sinf(slice*alpha), radius * sinf(stack*beta), radius * cosf(stack*beta) * cosf(slice*alpha)));
				pontos.push_back(Ponto(radius * cosf((stack + 1)*beta) * sinf((slice + 1)*alpha), radius * sinf((stack + 1)*beta), radius * cosf((stack + 1)*beta) * cosf((slice + 1)*alpha)));

				//metade inferior
				//triangulo inferior
				pontos.push_back(Ponto(radius * cosf(-(stack + 1)*beta) * sinf(slice*alpha), radius * sinf(-(stack + 1)*beta), radius * cosf(-(stack + 1)*beta) * cosf(slice*alpha)));
				pontos.push_back(Ponto(radius * cosf(-(stack + 1)*beta) * sinf((slice + 1)*alpha), radius * sinf(-(stack + 1)*beta), radius * cosf(-(stack + 1)*beta) * cosf((slice + 1)*alpha)));
				pontos.push_back(Ponto(radius * cosf(-stack * beta) * sinf((slice + 1)*alpha), radius * sinf(-stack * beta), radius * cosf(-stack * beta) * cosf((slice + 1)*alpha)));
				//triangulo superior
				pontos.push_back(Ponto(radius * cosf(-stack * beta) * sinf(slice*alpha), radius * sinf(-stack * beta), radius * cosf(-stack * beta) * cosf(slice*alpha)));
				pontos.push_back(Ponto(radius * cosf(-(stack + 1)*beta) * sinf(slice*alpha), radius * sinf(-(stack + 1)*beta), radius * cosf(-(stack + 1)*beta) * cosf(slice*alpha)));
				pontos.push_back(Ponto(radius * cosf(-stack * beta) * sinf((slice + 1)*alpha), radius * sinf(-stack * beta), radius * cosf(-stack * beta) * cosf((slice + 1)*alpha)));
			}
		}
		for (unsigned int ponto = 0; ponto < pontos.size() - 1; ponto++) {
			fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
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
		float beta = M_PI / stacks;

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

int main(int argc, char* argv[]) {
	
	/*argc = 4; //para testar
	argv[1] = "plane";
	argv[2] = "4";
	argv[3] = "plane.3d";*/

	argc = 7;
	argv[1] = "box"; argv[2] = "4.5";argv[3] = "4.5";argv[4] = "4.5";argv[5] = "2"; argv[6] = "box.3d";
	

	//argc = 6;
	//sphere - raidus - slices - stacks - filename
	//argv[1] = "sphere"; argv[2] = "1.8";argv[3] = "14";argv[4] = "14";argv[5] = "sphere.3d";
	

	/*argc = 7;
	//cone - radius - height - slices - stacks - filename
	argv[1] = "cone"; argv[2] = "1.5";argv[3] = "3";argv[4] = "12";argv[5] = "12"; argv[6] = "cone.3d";*/

	if (argc < 2) {
		printf("Parametros invalidos\n");
		return -1;
	}
	string tipo = argv[1];
	string path = "../../3DFiles/";

	if (tipo == "plane" && argc == 4) {
		float lado = stof(argv[2]);
		string filename = argv[3];
		filename = path + filename;
		plane(lado, filename);
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

	printf("Parametros Invalidos\n");
	return -1;
}