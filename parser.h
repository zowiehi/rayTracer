#ifndef PARSER_H
#define PARSER_H

typedef struct {
  int kind; // 0 = camera, 1 = sphere, 2 = plane
  union {
    struct {
      double width;
      double height;
    } camera;
    struct {
      double diffuse[3];
      double specular[3];
      double position[3];
      double radius;
      double reflect;
      double refract;
      double ior;
    } sphere;
    struct {
      double diffuse[3];
      double specular[3];
      double position[3];
      double normal[3];
      double reflect;
      double refract;
      double ior;
    } plane;
    struct {
      int type;
      double color[3];
      double position[3];
      double direction[3];
      double radial[3];
      double angular;
      double theta;
    } light;
  };
} Object;

double zerovec[3];

int getSize();
int next_c(FILE* json);
void expect_c(FILE* json, int d);
void skip_ws(FILE* json);
char* next_string(FILE* json);
double next_number(FILE* json);
double* next_vector(FILE* json);
Object** read_scene(char* filename);


#endif
