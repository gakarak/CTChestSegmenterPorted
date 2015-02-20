// ver 2.0

#include "stdafx.h"

Point3D ParseVertice(char* s);
Face3 ParseFace(char* s);
Edge ParseEdge(char* s);

/////////////////////////////////////////////////////////////////////////////////////////////////////

void ReadObj(char* filename, Point3D* vertices, Face3* facets, Edge* edges, int* verticesCount, int* facesCount, int* edgesCount){

	char buf[1000];
	FILE* f = fopen(filename, "rt");
	if (!f){
		sprintf(buf, "Failed to open file '%s' for reading", filename);
		throw buf;
	}

	(*verticesCount) = 0;
	(*facesCount) = 0;
	(*edgesCount) = 0;
	Point3D vert;
	Face3 fac;
	Edge ed;

	while (!feof(f)){
		fgets(buf, 1000, f);

		if (buf[0] == 'v'){
			vert = ParseVertice(buf);
			vertices[(*verticesCount)].x = vert.x;
			vertices[(*verticesCount)].y = vert.y;
			vertices[(*verticesCount)].z = vert.z;
			(*verticesCount)++;
		}

		if (buf[0] == 'f'){
			fac = ParseFace(buf);
			facets[(*facesCount)].verticeIndex1 = fac.verticeIndex1;
			facets[(*facesCount)].verticeIndex2 = fac.verticeIndex2;
			facets[(*facesCount)].verticeIndex3 = fac.verticeIndex3;
			(*facesCount)++;
		}

		if (buf[0] == 'e'){
			ed = ParseEdge(buf);
			edges[(*edgesCount)].verticeIndex1 = ed.verticeIndex1;
			edges[(*edgesCount)].verticeIndex2 = ed.verticeIndex2;
			(*edgesCount)++;
		}
	}

	fclose(f);
}

void WriteObj(char* filename, Point3D* vertices, Face3* facets, int verticesCount, int facetCount){
	int i, counter;

	FILE* f = fopen(filename, "wt");

	counter = 0;
	for (i = 0; i < verticesCount; i++){
		fprintf(f, "v %f %f %f\n", vertices[i].x, vertices[i].y, vertices[i].z);
		counter++;
	}
	fprintf(f, "# %i vertices\n\n", counter);

	counter = 0;
	fprintf(f, "g default\n");
	for (i = 0; i < facetCount; i++){
		if (facets[i].selected){
			fprintf(f, "f %i %i %i\n", facets[i].verticeIndex1 + 1, facets[i].verticeIndex2 + 1, facets[i].verticeIndex3 + 1);
			counter++;
		}
	}
	fprintf(f, "# %i faces\n\n", counter);

	fclose(f);
}

Point3D ParseVertice(char* s){
	Point3D vert;

	char buf[100];
	int j, i = 1;

	j = 0;
	while (s[i] == ' ')
		i++;
	while (s[i] != ' '){
		buf[j++] = s[i];
		i++;
	}
	buf[j] = 0;
	vert.x = (float)atof(buf);

	j = 0;
	while (s[i] == ' ')
		i++;
	while (s[i] != ' '){
		buf[j++] = s[i];
		i++;
	}
	buf[j] = 0;
	vert.y = (float)atof(buf);

	j = 0;
	while (s[i] == ' ')
		i++;
	while (s[i] != ' ' && i < strlen(s)){
		buf[j++] = s[i];
		i++;
	}
	buf[j] = 0;
	vert.z = (float)atof(buf);

	return vert;
}

Edge ParseEdge(char* s){
	Edge ed;

	char buf[10];
	int j, i = 1;

	j = 0;
	while (s[i] == ' ')
		i++;
	while (s[i] != ' '){
		buf[j++] = s[i];
		i++;
	}
	buf[j] = 0;
	ed.verticeIndex1 = atoi(buf) - 1;

	j = 0;
	while (s[i] == ' ')
		i++;
	while (s[i] != ' ' && i < strlen(s)){
		buf[j++] = s[i];
		i++;
	}
	buf[j] = 0;
	ed.verticeIndex2 = atoi(buf) - 1;

	return ed;
}

Face3 ParseFace(char* s){
	Face3 fac;

	char buf[10];
	int j, i = 1;

	j = 0;
	while (s[i] == ' ')
		i++;
	while (s[i] != ' '){
		buf[j++] = s[i];
		i++;
	}
	buf[j] = 0;
	fac.verticeIndex1 = atoi(buf) - 1;

	j = 0;
	while (s[i] == ' ')
		i++;
	while (s[i] != ' '){
		buf[j++] = s[i];
		i++;
	}
	buf[j] = 0;
	fac.verticeIndex2 = atoi(buf) - 1;

	j = 0;
	while (s[i] == ' ')
		i++;
	while (s[i] != ' ' && i < strlen(s)){
		buf[j++] = s[i];
		i++;
	}
	buf[j] = 0;
	fac.verticeIndex3 = atoi(buf) - 1;

	return fac;
}

Point3D Point3D_(float x, float y, float z){
	Point3D p;
	p.x = x;
	p.y = y;
	p.z = z;
	return p;
}

float Dot(Point3D p1, Point3D p2){
	return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

Point3D Cross(Point3D p1, Point3D p2){
	Point3D p;
	p.x = p1.y * p2.z - p2.y * p1.z;
	p.y = p1.z * p2.x - p2.z * p1.x;
	p.z = p1.x * p2.y - p2.x * p1.y;
	return p;
}

float MixProd(Point3D p1, Point3D p2, Point3D p3){
	return Dot(p1, Cross(p2, p3));
}