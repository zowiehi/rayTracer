/*
Created by Zowie Haugaard
10/20/16
Ray casting and lighting model
  This Ray caster creates a .ppm image from a .json style list of objects. This
  project has been modified from the origional ray caster to include support for
  a simple lighting model including spot and point lights

 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <ctype.h>
 #include <math.h>
 #include "parser.h"

 #define MAX_DEPTH 4

 double* trace(double* Ro, double* Rd, Object** objects, Object** lights, int depth);

//The shinyness component
 double ns = 20;

//clamp values between 0 and 1
double clamp(double v){
  if( v > 1) v = 1;
  else if(v < 0) v = 0;
  return v;
  }
//generate a vector of zero values
double* zvec(){
  double* vec = malloc(sizeof(double)*3);
  vec[0] = 0;
  vec[1] = 0;
  vec[2] = 0;
  return vec;
}

static inline double rad_to_deg(double rad) {
    return rad * (180.0 / M_PI);
}

//returnthe square of the supplied number
static inline double sqr(double v) {
   return v*v;
 }

//return the length of a vector
 static inline double len(double* v) {
   double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
   return len;
 }

//return the dot product of two double vectors
static inline double dot(double* v1, double* v2){
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

//return the cross product of two vectors
static inline double* cross(double* v1, double* v2){
  double* result = malloc(sizeof(double)*3);
  result[0] = v1[1]*v2[2] - v1[2]*v2[1];
  result[1] = v1[2]*v2[0] - v1[0]*v2[2];
  result[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return result;
}

//multply the coresponding values of two vectors by one another
static inline double* mult(double* v1, double* v2){
  double* result = malloc(sizeof(double)*3);
  result[0] = v1[0]*v2[0];
  result[1] = v1[1]*v2[1];
  result[2] = v1[2]*v2[2];
  return result;
}


//Subtract one vector from another
static inline double* sub(double* v1, double* v2){
  double* result = malloc(sizeof(double)*3);
  result[0] = v1[0] - v2[0];
  result[1] = v1[1] - v2[1];
  result[2] = v1[2] - v2[2];
  return result;
}

//add two vectors together
static inline double* add(double* v1, double* v2){
  double* result = malloc(sizeof(double)*3);
  result[0] = v1[0] + v2[0];
  result[1] = v1[1] + v2[1];
  result[2] = v1[2] + v2[2];
  return result;
}

//Scale a vector by some constant t
static inline double* scale(double t, double* v){
  double* result = malloc(sizeof(double)*3);
  result[0] = t * v[0];
  result[1] = t * v[1];
  result[2] = t * v[2];
  return result;
}

  //normalize a vector
 static inline double normalize(double* v) {
   double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
   v[0] /= len;
   v[1] /= len;
   v[2] /= len;
   return len;
 }

 //get the distance between two points, ie the length of the difference
 static inline double dist(double* v1, double* v2){
   return normalize(sub(v1, v2));
 }

 static inline double mix(double a, double b, double mix)
{
    return b * mix + a * (1 - mix);
}

//Calculate the radial attenuation in our lighting model
 double frad(double* rad, double d){

   double det = ((rad[0] * sqr(d)) + (rad[1] * d) + rad[2]);
   if(det == 0) return 0;
   else return 1 / det;
 }

 //Calculate the angular attenuation in our lighting model
 double fang(double ang, double theta, double* Rdn, double* dir){
   theta = rad_to_deg(theta);
   if(dir[0] == 0 && dir[1] == 0 && dir[2] == 0) return 1.0;
   double val = dot(dir, Rdn);

   if (val < cos(theta)) {
    return 0.0;
  }
   else return pow(val, ang);
 }

 //calculate the diffuse component of the lighting model
 double* dif(double* Kd, double* Il, double* N, double* L){
   double dots = dot(N,L);
   if(dots > 0)return scale(dots, mult(Kd, Il));
   else return zvec();
 }

 //calculate the specular component of the lighting model
 double* spec(double* Ks, double* Il, double* V, double* R, double* N, double* L){
   double dots = dot(V,R);
   double dots2 = dot(N,L);
   if (dots > 0 && dots2 > 0)return scale(pow(dots,ns), mult(Ks, Il));
   else return zvec();
 }

 //check for ray-plane intersections
 double plane_intersection(double* Ro, double* Rd,
          double* P, double* n){
            double a, b;

            a = dot(n, sub(P, Ro));
            b = dot(n, Rd);

            return a/b;
           }
  //check for ray-sphere intersections
 double sphere_intersection(double* Ro, double* Rd,
 			     double* C, double r){

             double a = (sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]));
             double b = (2 * (Ro[0] * Rd[0] - Rd[0] * C[0] + Ro[1] * Rd[1] - Rd[1] * C[1] + Ro[2] * Rd[2] - Rd[2] * C[2]));
             double c = sqr(Ro[0]) - 2*Ro[0]*C[0] + sqr(C[0]) + Ro[1] - 2*Ro[1]*C[1] + sqr(C[1]) + sqr(Ro[2]) - 2*Ro[2]*C[2] + sqr(C[2]) - sqr(r);

             //using the quadratic formula
             double det = sqr(b) - 4 * a * c;
             if (det < 0) return -1;

             det = sqrt(det);

             //only return the positive result
             double t0 = (-b - det) / (2*a);
             if (t0 > 0) return t0;

             double t1 = (-b + det) / (2*a);
             if (t1 > 0) return t1;

             return -1;

}

double* reflect(double* Ron, double* Rd, Object* best_obj, Object** objects, Object** lights, int depth){
  double* N;
  double* R;
  if(best_obj->kind == 1)	N = sub(Ron, best_obj->sphere.position);
  else if(best_obj->kind == 2) N = best_obj->plane.normal;
  normalize(N);
  R = add(scale(-2 * dot(Rd,N), N), Rd);
  normalize(R);
  if(best_obj->kind == 1) return scale(best_obj->sphere.reflect, trace(Ron, R, objects, lights, depth));
  if(best_obj->kind == 2) {
    return scale(best_obj->plane.reflect, trace(Ron, R, objects, lights, depth));
  }
}

double* refract(double* Ron, double* Rd, Object* best_obj, Object** objects, Object** lights, int depth){
  double* N;
  double R;
  double n;
  int in = 0;

  printf("Refract\n");
  if(best_obj->kind == 1)	N = sub(Ron, best_obj->sphere.position);
  else if(best_obj->kind == 2) N = best_obj->plane.normal;

  normalize(N);

  if(dot(Rd, N) > 0){
    N = sub(zvec(), N);
    in = 1;
  }

  if(best_obj->kind == 1)	n = 1/best_obj->sphere.ior;
  else if(best_obj->kind == 2) n = 1/best_obj->plane.ior;

  R = -dot(Rd,N);

  double c = sqrt( 1 - pow(n,2) * (1 - pow(R,2)));

  double* ref = add(scale(n, Rd), scale((n * R - c), N));
  normalize(ref);
  if(best_obj->kind == 1) return trace(sub(Ron, N), ref, objects, lights, depth);
  if(best_obj->kind == 2) {
    return trace(Ron, ref, objects, lights, depth);
  }
}




double* trace(double* Ro, double* Rd, Object** objects, Object** lights, int depth){

  double* curcolor = zvec();
  double* reflectColor = zvec();
  double* refractColor = zvec();
  double best_t = INFINITY;
  Object *best_obj = NULL;
  //Loop through objects and check for collisions with each
  for (int i=0; objects[i] != NULL; i += 1) {
    double t = 0;
    // run appropriate collision detection based on object type
    switch (objects[i]->kind) {
      case 0:
        break;
      case 1:
        t = sphere_intersection(Ro, Rd, objects[i]->sphere.position, objects[i]->sphere.radius);
        break;
      case 2:
        t = plane_intersection(Ro, Rd, objects[i]->plane.position, objects[i]->plane.normal);
        break;
      case 3:
        break;
      default:
        fprintf(stderr, "Error: unsupported object primative\n");
        exit(1);
    }

    //If we have a valid collision that is closer than the next best, set our
    //  best collision point and closest object
    if (t > 0 && t < best_t) {
      best_t = t;
      best_obj = objects[i];
    }
  }


  //If we didnt get a collision, go to the next loop and color the pixel black
  if (best_obj == NULL) {
    return curcolor;
  };
  //Get the diffuse color from the object
  double* difcol = malloc(sizeof(double)*3);
  if(best_obj->kind == 1) difcol = best_obj->sphere.diffuse;
  else if(best_obj->kind == 2) difcol = best_obj->plane.diffuse;

  //Get the specular color from the object
  double* speccol = malloc(sizeof(double)*3);
  if(best_obj->kind == 1) speccol = best_obj->sphere.diffuse;
  else if(best_obj->kind == 2) speccol = best_obj->plane.diffuse;


  //the point on the object where we got a collision
  double* Ron = add(scale(best_t, Rd), Ro);
  double* N;
  //get the normal of the object
  if(best_obj->kind == 1)	N = sub(Ron, best_obj->sphere.position);
  else if(best_obj->kind == 2)N = best_obj->plane.normal; // plane
  normalize(N);

  //loop through each of the lght objects
  for (int j=0; lights[j] != NULL; j+=1) {

    //the ray from that point to the current light
    double* Rdn = sub(lights[j]->light.position, Ron);
    //The distance to the light
    double dist = normalize(Rdn);
    //simple boolean to see if we are in shadow
    double closest_shadow_object = 0;
    //loop through Objects and see if any are between the collision point and the light
    for (int k=0; objects[k] != NULL; k+=1) {
      double t = 0;
      //we want to exclude the current object
      if (objects[k] == best_obj) continue;
      //detect collisions
      switch(objects[k]->kind) {
        case 0:
          break;
        case 1:
          t = sphere_intersection(Ron, Rdn, objects[k]->sphere.position, objects[k]->sphere.radius);
          break;
        case 2:
          t = plane_intersection(Ron, Rdn, objects[k]->plane.position, objects[k]->plane.normal);
          break;
        case 3:
          break;
      default:
          // ERROR
          break;
        }
        //the object was either behind the current object or behind the light
        if (t > dist || t <= 0) {
          continue;
        }
        //The object was between the collision point and the light
        else if (t < dist && t > 0) closest_shadow_object = 1;
      }
      //only do this if there is no shadow
      if (closest_shadow_object == 0) {
        double* L;
        double* R;
        double* V;
        double* Kd;
        double* Ks;

        L = Rdn; // light_position - Ron;

        R = sub(scale(2 * dot(L,N), N), L);
        normalize(R);
        V = scale(-1, Rd);
        //diffuse component
        Kd = dif(difcol, lights[j]->light.color, N, L); // uses object's diffuse color
        //specular component
        Ks = spec(speccol, lights[j]->light.color, V, R, N, L);; // uses object's specular color

        //our simple lighting model
        curcolor[0] += frad(lights[j]->light.radial, dist) * fang(lights[j]->light.angular, lights[j]->light.theta, L, lights[j]->light.direction) * (Kd[0] + Ks[0]);
        curcolor[1] += frad(lights[j]->light.radial, dist) * fang(lights[j]->light.angular, lights[j]->light.theta, L, lights[j]->light.direction) * (Kd[1] + Ks[1]);
        curcolor[2] += frad(lights[j]->light.radial, dist) * fang(lights[j]->light.angular, lights[j]->light.theta, L, lights[j]->light.direction) * (Kd[2] + Ks[2]);

      }
    }
    int inside = 0;
    if(depth < MAX_DEPTH){
      if(dot(Rd, N) < 0) {
        N = sub(zvec(), N);
        inside = 1;
      }
      double frat = -dot(Rd, N);
      double fresnel = mix(pow(1 - frat, 5), 1, 1);

      double refr;
      double refl;


      if(best_obj->kind == 1)	refr = best_obj->sphere.refract;
      else if(best_obj->kind == 2) refr = best_obj->plane.refract;

      if(best_obj->kind == 1)	refl = best_obj->sphere.reflect;
      else if(best_obj->kind == 2) refl = best_obj->plane.reflect;

      printf("depth %d\n", depth);
      depth += 1;
      if(refl > 0) reflectColor = reflect(Ron, Rd, best_obj, objects, lights, depth);
      if (refr > 0) refractColor = refract(Ron, Rd, best_obj, objects, lights, depth);

      curcolor[0] += ((fresnel * reflectColor[0]) + (refractColor[0] * (1 - fresnel) * refr)) * curcolor[0];
      curcolor[1] += ((fresnel * reflectColor[1]) + (refractColor[1] * (1 - fresnel) * refr)) * curcolor[1];
      curcolor[2] += ((fresnel * reflectColor[2]) + (refractColor[2] * (1 - fresnel) * refr)) * curcolor[2];
    }
    curcolor[0] = clamp(curcolor[0]);
    curcolor[1] = clamp(curcolor[1]);
    curcolor[2] = clamp(curcolor[2]);
    return curcolor;
}
 //represents a single pixel object
 typedef struct RGB {
   unsigned char r, g, b;
 }RGBPix;

 //this struct is used to store the entire image data, along with the width and height
 typedef struct Image {
   int width, height;
   RGBPix *data;
 }PPMImage;

 int main(int argc, char *argv[]) {

   //Make sure the right number of arguments was supplied
   if(argc < 5 || argc > 5) {
     fprintf(stderr, "usage: ./raycast width height object-file.json output-file.ppm \n");
     exit(1);
   }

   //check the object file
   char* inFile = argv[3];
   if(strstr(inFile, ".json") == NULL){
     fprintf(stderr, "Please provide a valid json file\n");
     exit(1);
   }

   //read in the scene file, terminate the array
   Object** objects = read_scene(inFile);
   objects[getSize()] = NULL;



   //create an array of just lights
   Object** lights = malloc(sizeof(Object*)*2);
   int Li = 0;
   //pull lights from our objects array and put them in the new one
   for (int i=0; objects[i] != NULL; i += 1){
     if(objects[i]->kind == 3){
       lights[Li] = malloc(sizeof(Object));
       lights[Li] = objects[i];
       Li += 1;
       lights = realloc(lights, sizeof(*lights)*(Li+2));
     }
   }
   lights[Li] = NULL;


   double *curcolor = malloc(sizeof(double)*3);
   static unsigned char backcolor[3] = {85, 85, 85};

   //get width and height from the input
   int width = atoi(argv[1]);
   int height = atoi(argv[2]);

   if(width < 1 || height < 1){
     printf("Please provide valid width and height values\n");
   }


   //create and allociate space for the ppmimage
   PPMImage image;
   image.width = width;
   image.height = height;
   image.data = malloc(sizeof(RGBPix)* width * height);

   //camera origin
   double cx = 0;
   double cy = 0;

   //view plane dimension variables
   float h;
   float w;

   //First, get the camera object. Then set the view plane dimensions based on
   // the camera's width and height
   for (int i=0; objects[i] != NULL; i += 1){
      if(objects[i]->kind == 0){
        if(objects[i]->camera.width > 0 && objects[i]->camera.height > 0){
          w = objects[i]->camera.width;
          h = objects[i]->camera.height;
        }

      }
   }

   //check the ppm image file
   char* outName = argv[4];
   if(strstr(outName, ".ppm") == NULL) {
     fprintf(stderr, "Please provide a valid .ppm file name to write to\n");
     exit(1);
   }

   //open the ppm image
   FILE *outFile = fopen(outName, "wb");

   //Print the header info in the image file
   (void) fprintf(outFile, "P6\n%d %d\n255\n", width, height);

   //pixel width and height
   double pixheight = h / height;
   double pixwidth = w / width;

   //loop through each pixel
   for (int y = 0; y <= height; y += 1) {
     for (int x = 0; x <= width; x += 1) {
       //Ray casting origin
       double Ro[3] = {0, 0, 0};
       //Ray from Ro through the center of the current pixel
       double Rd[3] = {
 	       cx - (w/2) + pixwidth * (x + 0.5),
 	       cy - (h/2) + pixheight * (y + 0.5),
 	       1
       };
       //normalize the ray
      (void) normalize(Rd);


      curcolor = trace(Ro, Rd, objects, lights, 0);

      int curPix = (height - y) * width + x;

      //24-bit colors from our double color values
     static unsigned char outcolor[3];

     //clamp colors between 0 and 1, then multiply by 255
     outcolor[0] = (char) (255 * clamp(curcolor[0]));
     outcolor[1] = (char) (255 * clamp(curcolor[1]));
     outcolor[2] = (char) (255 * clamp(curcolor[2]));

     //Color the pixel in the buffer
    image.data[curPix].r = outcolor[0];
    image.data[curPix].g = outcolor[1];
    image.data[curPix].b = outcolor[2];



   }

 }

 //Write our data array to the image, then close it out
 fwrite(image.data, 3 * image.width, image.height, outFile);
 //close our file
 (void) fclose(outFile);

 return 0;
 }
