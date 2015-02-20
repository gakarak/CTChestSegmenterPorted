// ver 2.0

struct Point3D{
	float x;
	float y;
	float z;

	bool isBorder;
};

struct Face3{
	int verticeIndex1;
	int verticeIndex2;
	int verticeIndex3;

	float area;
	bool selected;
};

struct Edge{
	int verticeIndex1;
	int verticeIndex2;
};

void ReadObj(char* filename, struct Point3D* vertices, struct Face3* facets, struct Edge* edges, int* verticesCount, int* facesCount, int* edgesCount);

void WriteObj(char* filename, struct Point3D* vertices, struct Face3* facets, int verticesCount, int facesCount, int edgesCount);

Point3D Point3D_(float x, float y, float z);

float Dot(Point3D p1, Point3D p2);

Point3D Cross(Point3D p1, Point3D p2);

float MixProd(Point3D p1, Point3D p2, Point3D p3);