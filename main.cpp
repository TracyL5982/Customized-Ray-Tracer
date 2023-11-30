#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>

#include "SETTINGS.h"

using namespace std;

int xRes = 800;
int yRes = 600;

//eye, lookAt, up
VEC3 eye(0.0, 0.0, 0.0);
VEC3 lookAt(0.0, 0.0, 1.0);
VEC3 up(0.0, 1.0, 0.0);

//sphere 0
VEC3 center_s_0(-3.5, 0, 10);
double radius_0 = 3;
VEC3 color_s_0(1, 0.25, 0.25);

//sphere 1
VEC3 center_s_1(3.5, 0, 10);
double radius_1 = 3;
VEC3 color_s_1(0.25, 0.25, 1);

//sphere 2
VEC3 center_s_2(0, -1000, 10);
double radius_2 = 997;
VEC3 color_s_2(0.5, 0.5, 0.5);

//light 0
VEC3 l_0_pos(10, 3, 5);
VEC3 l_0_color(1, 1, 1);

//light 1
VEC3 l_1_pos(-10, 3, 7.5);
VEC3 l_1_color(0.5, 0, 0);

//light 2
VEC3 l_2_pos(10, 10, 5);
VEC3 l_2_color(1, 1, 1);

//light 3
VEC3 l_3_pos(-10, 10, 7.5);
// VEC3 l_3_color(1, 1, 1);
VEC3 l_3_color(0.5, 0.25, 0.25);

//index of refraction
double n_air = 1.0;
double n_glass = 1.5;

//Phong exponent
double phong_ex = 10;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readPPM(const string& filename, int& xRes, int& yRes, float*& values)
{
  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "rb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for reading." << endl;
    cout << " Make sure you're not trying to read from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  // get the dimensions
  unsigned char newline;
  fscanf(fp, "P6\n%d %d\n255%c", &xRes, &yRes, &newline);
  if (newline != '\n') {
    cout << " The header of " << filename.c_str() << " may be improperly formatted." << endl;
    cout << " The program will continue, but you may want to check your input. " << endl;
  }
  int totalCells = xRes * yRes;

  // grab the pixel values
  unsigned char* pixels = new unsigned char[3 * totalCells];
  fread(pixels, 1, totalCells * 3, fp);

  // copy to a nicer data type
  values = new float[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = pixels[i];

  // clean up
  delete[] pixels;
  fclose(fp);
  cout << " Read in file " << filename.c_str() << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void writePPM(const string& filename, int& xRes, int& yRes, float*& values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

//////////////////////////////////////////////////////////////////////////////////
// generate ray
//////////////////////////////////////////////////////////////////////////////////
VEC3 generateRay (int i, int j){
  //calculate l, r, t, b first
  double fov = 65;
  double scale = tan (fov * 0.5 * M_PI/180);
  double n = 1.0;
  double t = scale * n;
  double aspect = 800.0f/600.0f;
  double r = scale * aspect * n;
  double l = -r;
  double b = -t;
  //set up the coordinates
  VEC3 g = lookAt - eye;
  VEC3 w = -g.normalized();
  VEC3 u = up.cross(w);
  u.normalize();
  VEC3 v = w.cross(u);
  v.normalize();
  VEC3 s = (l + (r - l) * (i + 0.5f) / xRes) * u + (b + (t - b) * (j + 0.5f) / yRes) * v - w*n;
  VEC3 d = s.normalized();
  d[0] = - d[0];
  return d;
}

//////////////////////////////////////////////////////////////////////////////////
// check sphere intersection
//////////////////////////////////////////////////////////////////////////////////
bool checkSphereIntersection(VEC3 center, VEC3 d, VEC3 eye, double radius, double*t){
  double a = d.dot(d);
  double b = 2*d.dot(eye - center);
  double c = (eye - center).dot(eye-center)-radius*radius;
  if(b*b - 4*a*c < 0) return false;
  double t1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
  double t2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
  if(t1 > 0 || t2 > 0){
    if(t1 > 0 && t2 > 0){
      //if both roots > 0, take the smaller one
      if(t1 < t2){
        *t = t1;
      }else{
        *t = t2;
      }
      //if one positive root, one negative root, take the positive one
    }else{
      if(t1 > 0){
        *t = t1;
      }else{
        *t = t2;
      }
    }
    return true;
  }else{
    return false;
  }
}

//////////////////////////////////////////////////////////////////////////////////
// calculate determinant
//////////////////////////////////////////////////////////////////////////////////
double getDeterminant(VEC3 a, VEC3 b, VEC3 c){
  return a[0] * (b[1] * c[2] - b[2] * c[1]) -
         a[1] * (b[0] * c[2] - b[2] * c[0]) +
         a[2] * (b[0] * c[1] - b[1] * c[0]);
}

VEC3 get_beta_gamma_t(VEC3 a_b, VEC3 a_c, VEC3 d, VEC3 a_e){
  double A = getDeterminant(a_b, a_c, d);
  if (abs(A) < std::numeric_limits<double>::epsilon()) {
    // Return a VEC3 with negative t to indicate no intersection
    return VEC3(-1, 0, 0);
  }
  double beta_top = getDeterminant(a_e, a_c, d);
  double gamma_top = getDeterminant(a_b, a_e, d);
  double t_top = -getDeterminant(a_b, a_c, a_e);
  double beta = beta_top/A;
  double gamma = gamma_top/A;
  double t = -t_top / A;
  return VEC3(beta, gamma, t);
}

//////////////////////////////////////////////////////////////////////////////////
// check triangle intersection
//////////////////////////////////////////////////////////////////////////////////
bool checkTriangleIntersection(vector<VEC3>& vertices, VEC3 d, VEC3 eye, double*t){
  VEC3 a = vertices[0];
  VEC3 b = vertices[1];
  VEC3 c = vertices[2];
  VEC3 a_b = b - a;
  VEC3 a_c = c - a;
  VEC3 a_e = eye - a;
  // Plane
  /*
  P = (x, y, z);
  N = a_b.cross(a_c)
  d = vertices[0]

  Ray
  P_0 = eye
  t = t
  V = d
  */
  // Step 1: find t
  VEC3 planeNormal = (a_b.cross(a_c)).normalized();
  float t1 = (vertices[0] - eye).dot(planeNormal) / (d.dot(planeNormal));
  if(t1 <= 0){
    return false;
  }

  VEC3 intersectionPoint = eye + d * t1;
  // Step 2: find whether intersectionPoint is in triangle
  VEC3 bgt = get_beta_gamma_t(vertices[0], vertices[1], vertices[2], intersectionPoint);
  double beta = bgt[0];
  double gamma = bgt[1];
  double alpha = bgt[2];
  if(alpha < 0) return false;
  if(gamma < 0 || gamma > 1) return false;
  if(beta < 0 || beta > 1-gamma) return false;
  *t = t1;
  return true;
}

//////////////////////////////////////////////////////////////////////////////////
// get shading
//////////////////////////////////////////////////////////////////////////////////
VEC3 getShading(VEC3 center, VEC3 d, VEC3 eye, double radius, VEC3 lightColor, VEC3 lightPos){
  VEC3 color(0.0, 0.0, 0.0);
  double t = 0;
  if(checkSphereIntersection(center, d, eye, radius, &t) == false){
    return color;
  }
  //calculate r (where the light hits the sphere based on t)
  VEC3 r = eye + t * d;
  VEC3 n = (r - center).normalized();
  VEC3 lightDir = (lightPos - r).normalized();
  double diff = max(0.0, n.dot(lightDir));
  color = lightColor * diff;
  return color;
}

VEC3 getTriangleShading(vector<VEC3> &vertices, VEC3 d, VEC3 eye, VEC3 lightColor, VEC3 lightPos){
  VEC3 color(0.0, 0.0, 0.0);
  double t = 0;
  if(checkTriangleIntersection(vertices, d, eye, &t) == false){
    return color;
  }
  //calculate r (where the light hits the sphere based on t)
  VEC3 r = eye + t * d;
  //calculate n
  VEC3 edge1 = vertices[1] - vertices[0];
  VEC3 edge2 = vertices[2] - vertices[0];
  VEC3 n = edge1.cross(edge2).normalized();
  VEC3 lightDir = (lightPos - r).normalized();
  double diff = max(0.0, n.dot(lightDir));
  color = lightColor * diff;
  return color;
}


//////////////////////////////////////////////////////////////////////////////////
// get reflection direction
//////////////////////////////////////////////////////////////////////////////////
VEC3 getReflection (VEC3 center, VEC3 r, VEC3 d){
  VEC3 n = (r - center).normalized();
  VEC3 reflectDir = d - 2 * (n.dot(d)) * n;
  return reflectDir.normalized();
}

VEC3 triangle_getReflection (vector<VEC3>& vertices, VEC3 d){
  VEC3 edge1 = vertices[1] - vertices[0];
  VEC3 edge2 = vertices[2] - vertices[0];
  VEC3 n = edge1.cross(edge2).normalized();
  VEC3 reflectDir = d - 2 * (n.dot(d)) * n;
  return reflectDir.normalized();
}

//////////////////////////////////////////////////////////////////////////////////
// get refraction direction
//////////////////////////////////////////////////////////////////////////////////
VEC3 getRefraction(VEC3 center, VEC3 r, VEC3 d, VEC3 n, double n_in, double n_out) {
  // Normalize vectors
  d = d.normalized();
  n = n.normalized();

  double epsilon = 1e-4;

  // Correct the normal direction
  if (d.dot(n) > 0) {
    // We are inside the sphere, so we need to invert the normal
    n = -n;
    // Also swap the refraction indices
    swap(n_in, n_out);
  }

  double cos_i = -d.dot(n);
  double sin_t2 = (n_in / n_out) * (n_in / n_out) * (1.0 - cos_i * cos_i);

  if (sin_t2 > 1.0) {
    // Total internal reflection
    cout << "total internal reflection!" << endl;
    return VEC3(); // Or handle as needed
  }

  double cos_t = sqrt(1.0 - sin_t2);
  VEC3 refractDir = (n_in / n_out) * d + (n_in / n_out * cos_i - cos_t) * n;

  return refractDir.normalized();
}

//////////////////////////////////////////////////////////////////////////////////
// generate specular highlight
//////////////////////////////////////////////////////////////////////////////////
VEC3 getHighlight (VEC3 center, VEC3 d, VEC3 eye, double radius, VEC3 lightColor, VEC3 lightPos){
  VEC3 color(0.0, 0.0, 0.0);
  double t = 0;
  if(checkSphereIntersection(center, d, eye, radius, &t) == false){
    return color;
  }
  //calculate r (where the light hits the sphere based on t)
  VEC3 r = eye + t * d;
  VEC3 n = (r - center).normalized();
  VEC3 lightDir = (lightPos - r).normalized();
  VEC3 reflectDir = (- lightDir + 2 * (n.dot(lightDir))*n).normalized();
  double spec = max(0.0, reflectDir.dot(-d));
  double coefficient = pow(spec, phong_ex);
  color = lightColor * coefficient;
  return color;
}

VEC3 getTriangleHighlight (vector<VEC3> vertices, VEC3 d, VEC3 eye, VEC3 lightColor, VEC3 lightPos){
  VEC3 color(0.0, 0.0, 0.0);
  double t = 0;
  if(checkTriangleIntersection(vertices, d, eye, &t) == false){
    return color;
  }
  //calculate r (where the light hits the sphere based on t)
  VEC3 r = eye + t * d;
  //calculate n
  VEC3 edge1 = vertices[1] - vertices[0];
  VEC3 edge2 = vertices[2] - vertices[0];
  VEC3 n = edge1.cross(edge2).normalized();
  VEC3 lightDir = (lightPos - r).normalized();
  VEC3 reflectDir = (- lightDir + 2 * (n.dot(lightDir))*n).normalized();
  double spec = max(0.0, reflectDir.dot(-d));
  double coefficient = pow(spec, phong_ex);
  color = lightColor * coefficient;
  return color;
}

//////////////////////////////////////////////////////////////////////////////////
// check if the point is in shadow
//////////////////////////////////////////////////////////////////////////////////
bool isInShadow(int i, VEC3 point, VEC3 lightPos,
                const vector<VEC3>& centers, const vector<double>& radii){
  double epsilon = 1e-4;
  VEC3 toLight = (lightPos - point).normalized();
  // time for the shadow ray to travel from point to light source
  double lightDistance = (lightPos - point).norm();
  // Check whether the light is occluded
  double t1 = 0;
  point += epsilon * toLight;
  for (size_t k = 0; k < centers.size(); k++) {
    if(k == i) continue;
    if (checkSphereIntersection(centers[k], toLight, point, radii[k], &t1)) {
      if(t1 > 0 && t1 < lightDistance){
        return true;
      }
    }
  }
  return false;
}

bool isInSphereShadow(bool isSphere, int i, VEC3 point, VEC3 lightPos,
                const vector<VEC3>& centers, const vector<double>& radii){
  double epsilon = 1e-4;
  VEC3 toLight = (lightPos - point).normalized();
  //time for the shadow ray to travel from point to light source
  double lightDistance = (lightPos - point).norm();
  // Check whether the light is occluded
  double t1 = 0;
  point += epsilon * toLight;
  for (size_t k = 0; k < centers.size(); k++) {
    if(isSphere){
       if(k == i) continue;
    }
    if (checkSphereIntersection(centers[k], toLight, point, radii[k], &t1)) {
      if(t1 > 0 && t1 < lightDistance){
        return true;
      }
    }
  }
  return false;
}

bool isInTriangleShadow (bool isTriangle, int i, VEC3 point, VEC3 lightPos, vector<vector<VEC3> >& triangles){
  double epsilon = 1e-2;
  VEC3 toLight = (lightPos - point).normalized();
  //time for the shadow ray to travel from point to light source
  double lightDistance = (lightPos - point).norm();
  // Check whether the light is occluded
  double t1 = 0;
  point += epsilon * toLight;

  for (size_t j = 0; j < triangles.size(); j++){
    if(isTriangle){
       if(j == i) continue;
    }
    if (checkTriangleIntersection(triangles[j], toLight, point, &t1)){
      if(t1 > 0 && t1 < lightDistance){
        return true;
      }
    }
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////////
// Calculate the Fresnel Effect
//////////////////////////////////////////////////////////////////////////////////
pair<double, double> calculateFresnel(VEC3 incident, VEC3 normal, double n1, double n2) {
  incident.normalized();
  normal.normalized();
  double cosi = incident.dot(normal);
  double etai = n1, etat = n2;
  if (cosi > 0) { swap(etai, etat); }
  double sint = etai / etat * sqrt(max(0.0, 1 - cosi * cosi));

  if (sint >= 1) {
    // Total internal reflection
    return make_pair(1.0, 0.0);
  } else {
    double cost = sqrt(max(0.0, 1 - sint * sint));
    cosi = fabs(cosi);
    double Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
    double Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
    double kr = (Rs * Rs + Rp * Rp) / 2;
    return make_pair(kr, 1 - kr);
  }
}

//////////////////////////////////////////////////////////////////////////////////
// Get Reflection Color
//////////////////////////////////////////////////////////////////////////////////
VEC3 getReflectionColor (VEC3 d, VEC3 point,
  vector<VEC3>& centers, vector<double>& radii, vector<VEC3>&colors,
  vector<VEC3>& lightPositions, vector<VEC3>& lightColors){

  VEC3 reflectionColor (0.0, 0.0, 0.0);
  VEC3 reflectedRay = getReflection(centers[1], point, d);
  double reflectionMinTime = INFINITY;
  int reflectionClosestSphereIndex = -1;
  VEC3 reflectionClosestSphereColor;

  for(size_t i = 0; i < centers.size(); i++){
    if(i == 1) continue;
    double t = 0.0;
    if(checkSphereIntersection(centers[i], reflectedRay, point, radii[i], &t)
      && t < reflectionMinTime){
        reflectionMinTime = t;
        reflectionClosestSphereIndex = i;
        reflectionClosestSphereColor = colors[i];
      }
  }

  if(reflectionClosestSphereIndex != -1){
    VEC3 reflectionPoint = point + reflectionMinTime * reflectedRay;
    for (size_t j = 0; j < lightPositions.size(); j++) {
      if (!isInShadow(reflectionClosestSphereIndex, reflectionPoint, lightPositions[j], centers, radii)) {
        reflectionColor += getShading(centers[reflectionClosestSphereIndex], reflectedRay, point, radii[reflectionClosestSphereIndex], lightColors[j], lightPositions[j])
              + getHighlight(centers[reflectionClosestSphereIndex], reflectedRay, point, radii[reflectionClosestSphereIndex], lightColors[j], lightPositions[j]);
      }
    }
  }

  reflectionColor[0] = reflectionColor[0] * reflectionClosestSphereColor[0];
  reflectionColor[1] = reflectionColor[1] * reflectionClosestSphereColor[1];
  reflectionColor[2] = reflectionColor[2] * reflectionClosestSphereColor[2];

  return reflectionColor;
}

VEC3 sphere_getReflectionColor_recursive (VEC3 d, VEC3 point,
  vector<VEC3>& centers, vector<double>& radii, vector<VEC3>&colors, vector<vector<VEC3> >& triangles,
  vector<VEC3>& lightPositions, vector<VEC3>& lightColors, VEC3& newEye, VEC3& newDir){

  VEC3 reflectionColor (0.0, 0.0, 0.0);
  VEC3 reflectedRay = getReflection(centers[1], point, d);
  double reflectionMinTime = INFINITY;
  int reflectionClosestSphereIndex = -1;
  int reflectionClosestTriangleIndex = -1;
  VEC3 reflectionClosestSphereColor (0.0, 0.0, 0.0);

  for(size_t i = 0; i < centers.size(); i++){
    if(i == 1) continue;
    double t = 0.0;
    if(checkSphereIntersection(centers[i], reflectedRay, point, radii[i], &t)
      && t < reflectionMinTime){
        reflectionMinTime = t;
        reflectionClosestSphereIndex = i;
        reflectionClosestSphereColor = colors[i];
      }
  }

  //check triangle intersection
  for (size_t i = 0; i < triangles.size(); i++) {
    double t = 0;
    if(checkTriangleIntersection(triangles[i], d, eye, &t) && t < reflectionMinTime){
      reflectionMinTime = t;
      reflectionClosestSphereIndex = -1;
      reflectionClosestTriangleIndex = i;
    }
  }

  if(reflectionClosestSphereIndex != -1){
    VEC3 reflectionPoint = point + reflectionMinTime * reflectedRay;
    for (size_t j = 0; j < lightPositions.size(); j++) {
      if (!isInShadow(reflectionClosestSphereIndex, reflectionPoint, lightPositions[j], centers, radii)) {
        reflectionColor += getShading(centers[reflectionClosestSphereIndex], reflectedRay, point, radii[reflectionClosestSphereIndex], lightColors[j], lightPositions[j])
              + getHighlight(centers[reflectionClosestSphereIndex], reflectedRay, point, radii[reflectionClosestSphereIndex], lightColors[j], lightPositions[j]);
      }
    }
  }

  if(reflectionClosestTriangleIndex != -1){
    newEye = point;
    newDir = reflectedRay.normalized();
  }

  reflectionColor[0] = reflectionColor[0] * reflectionClosestSphereColor[0];
  reflectionColor[1] = reflectionColor[1] * reflectionClosestSphereColor[1];
  reflectionColor[2] = reflectionColor[2] * reflectionClosestSphereColor[2];

  return reflectionColor;
}


VEC3 triangle_getReflectionColor_recursive (VEC3 d, VEC3 point,int index,
  vector<VEC3>& centers, vector<double>& radii, vector<VEC3>&colors, vector<vector<VEC3> > triangles,
  vector<VEC3>& lightPositions, vector<VEC3>& lightColors, VEC3& newEye, VEC3& newDir){

  VEC3 reflectionColor (0.0, 0.0, 0.0);
  VEC3 reflectedRay = triangle_getReflection(triangles[index], d);
  double reflectionMinTime = INFINITY;
  int reflectionClosestSphereIndex = -1;
  VEC3 reflectionClosestSphereColor (0.0, 0.0, 0.0);

  for(size_t i = 0; i < centers.size(); i++){
    double t = 0.0;
    if(checkSphereIntersection(centers[i], reflectedRay, point, radii[i], &t)
      && t < reflectionMinTime){
        reflectionMinTime = t;
        reflectionClosestSphereIndex = i;
        reflectionClosestSphereColor = colors[i];
      }
  }

  if(reflectionClosestSphereIndex != -1){
    VEC3 reflectionPoint = point + reflectionMinTime * reflectedRay;
    if(reflectionClosestSphereIndex == 1){
      //hit the reflective sphere
      newEye = point;
      newDir = reflectedRay.normalized();
    }else{
      for (size_t j = 0; j < lightPositions.size(); j++) {
        if (!isInShadow(reflectionClosestSphereIndex, reflectionPoint, lightPositions[j], centers, radii)) {
          reflectionColor += getShading(centers[reflectionClosestSphereIndex], reflectedRay, point, radii[reflectionClosestSphereIndex], lightColors[j], lightPositions[j])
                + getHighlight(centers[reflectionClosestSphereIndex], reflectedRay, point, radii[reflectionClosestSphereIndex], lightColors[j], lightPositions[j]);
        }
      }
    }
  }
  reflectionColor[0] = reflectionColor[0] * reflectionClosestSphereColor[0];
  reflectionColor[1] = reflectionColor[1] * reflectionClosestSphereColor[1];
  reflectionColor[2] = reflectionColor[2] * reflectionClosestSphereColor[2];

  return reflectionColor;
}
//////////////////////////////////////////////////////////////////////////////////
// Get Refraction Color
//////////////////////////////////////////////////////////////////////////////////
VEC3 getRefractionColor (VEC3 d, VEC3 point,
  vector<VEC3>& centers, vector<double>& radii, vector<VEC3>&colors,
  vector<VEC3>& lightPositions, vector<VEC3>& lightColors){

  VEC3 refractionColor (0.0, 0.0, 0.0);
  double refractionMinTime = INFINITY;
  int refractionClosestSphereIndex = -1;
  VEC3 refractionClosestSphereColor;
  VEC3 n_in = (point - centers[1]).normalized();
  VEC3 refractedRay = getRefraction(centers[1], point, d, n_in, n_air, n_glass);
  double epsilon = 1e-4;
  point = point + refractedRay * epsilon;
  double t_inside = 0.0;

  if(checkSphereIntersection(centers[1], refractedRay, point, radii[1], &t_inside)){
    VEC3 refractPoint = point + t_inside * refractedRay;
    //refraction 2: glass --> air
    VEC3 n_out = (refractPoint - centers[1]).normalized();
    VEC3 refractedRay_2 = getRefraction(centers[1], refractPoint, refractedRay, n_out, n_air, n_glass);
    double epsilon = 1e-4;
    refractPoint = refractPoint + refractedRay_2 * epsilon;
    for(size_t i = 0; i < centers.size(); i++){
      if(i == 1) continue;
      double t = 0.0;
      if(checkSphereIntersection(centers[i], refractedRay_2, refractPoint, radii[i], &t)
        && t < refractionMinTime){
          refractionMinTime = t;
          refractionClosestSphereIndex = i;
          refractionClosestSphereColor = colors[i];
        }
    }
    //if refractedRay_2 does intersect something
    if(refractionClosestSphereIndex != -1){
      VEC3 refractionPoint = refractPoint + refractionMinTime * refractedRay_2;
      refractionPoint = refractionPoint - refractedRay_2 * epsilon;
      for (size_t j = 0; j < lightPositions.size(); j++) {
        if (!isInShadow(refractionClosestSphereIndex, refractionPoint, lightPositions[j], centers, radii)) {
          refractionColor += getShading(centers[refractionClosestSphereIndex], refractedRay_2, refractionPoint, radii[refractionClosestSphereIndex], lightColors[j], lightPositions[j])
                + getHighlight(centers[refractionClosestSphereIndex], refractedRay_2, refractionPoint, radii[refractionClosestSphereIndex], lightColors[j], lightPositions[j]);
        }
      }
    }
  }

  refractionColor[0] = refractionColor[0] * refractionClosestSphereColor[0];
  refractionColor[1] = refractionColor[1] * refractionClosestSphereColor[1];
  refractionColor[2] = refractionColor[2] * refractionClosestSphereColor[2];

  return refractionColor;
}

VEC3 sphere_getRefractionColor_recursive (VEC3 d, VEC3 point,
  vector<VEC3>& centers, vector<double>& radii, vector<VEC3>&colors, vector<vector <VEC3> >&triangles,
  vector<VEC3>& lightPositions, vector<VEC3>& lightColors, VEC3& newEye, VEC3& newDir){

  VEC3 refractionColor (0.0, 0.0, 0.0);
  double refractionMinTime = INFINITY;
  int refractionClosestSphereIndex = -1;
  int refractionClosestTriangleIndex = -1;
  VEC3 refractionClosestSphereColor;

  //refraction 1: air --> glass
  VEC3 n_in = (point - centers[1]).normalized();
  VEC3 refractedRay = getRefraction(centers[1], point, d, n_in, n_air, n_glass);
  double epsilon = 1e-4;
  point = point + refractedRay * epsilon;
  double t_inside = 0.0;

  if(checkSphereIntersection(centers[1], refractedRay, point, radii[1], &t_inside)){
    VEC3 refractPoint = point + t_inside * refractedRay;
    //refraction 2: glass --> air
    VEC3 n_out = (refractPoint - centers[1]).normalized();
    VEC3 refractedRay_2 = getRefraction(centers[1], refractPoint, refractedRay, n_out, n_air, n_glass);
    double epsilon = 1e-4;
    refractPoint = refractPoint + refractedRay_2 * epsilon;

    //check sphere intersection
    for(size_t i = 0; i < centers.size(); i++){
      if(i == 1) continue;
      double t = 0.0;
      if(checkSphereIntersection(centers[i], refractedRay_2, refractPoint, radii[i], &t)
        && t < refractionMinTime){
          refractionMinTime = t;
          refractionClosestSphereIndex = i;
          refractionClosestSphereColor = colors[i];
        }
    }

    //check triangle intersection
    for (size_t i = 0; i < triangles.size(); i++) {
      double t = 0;
      if(checkTriangleIntersection(triangles[i], d, eye, &t) && t < refractionMinTime){
        refractionMinTime = t;
        refractionClosestSphereIndex = -1;
        refractionClosestTriangleIndex = i;
      }
    }
    //if refractedRay_2 does intersect something
    if(refractionClosestSphereIndex != -1){
      if(refractionClosestSphereIndex == 1){
        //enter recursion loop
        newEye = refractPoint;
        newDir = refractedRay_2.normalized();
        return refractionColor;
      }
      else{
        VEC3 refractionPoint = refractPoint + refractionMinTime * refractedRay_2;
        refractionPoint = refractionPoint - refractedRay_2 * epsilon;
        for (size_t j = 0; j < lightPositions.size(); j++) {
          if (!isInShadow(refractionClosestSphereIndex, refractionPoint, lightPositions[j], centers, radii)) {
            refractionColor += getShading(centers[refractionClosestSphereIndex], refractedRay_2, refractionPoint, radii[refractionClosestSphereIndex], lightColors[j], lightPositions[j])
                  + getHighlight(centers[refractionClosestSphereIndex], refractedRay_2, refractionPoint, radii[refractionClosestSphereIndex], lightColors[j], lightPositions[j]);
          }
        }
        refractionColor[0] = refractionColor[0] * refractionClosestSphereColor[0];
        refractionColor[1] = refractionColor[1] * refractionClosestSphereColor[1];
        refractionColor[2] = refractionColor[2] * refractionClosestSphereColor[2];

        return refractionColor;
      }
    }
    else if (refractionClosestTriangleIndex != -1){
      newEye = refractPoint;
      newDir = refractedRay_2.normalized();;
      return refractionColor;
    }
  }
  return refractionColor;
}

//////////////////////////////////////////////////////////////////////////////////
// 1x.ppm
//////////////////////////////////////////////////////////////////////////////////
void make_1x_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];
  for (int y = 0; y < yRes; y++){
    for (int x = 0; x < xRes; x++){
      VEC3 d = generateRay(x,y);
      int index = (yRes -1 - y) * xRes + x;
      float red = 1.0f * d[0];
      values[index*3] = red*255.0f;
      values[index*3+1] = 0.0;
      values[index*3+2] = 0.0;
    }
  }
  writePPM("1x.ppm", xRes, yRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// 1xabs.ppm
//////////////////////////////////////////////////////////////////////////////////
void make_1xabs_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];
  for (int y = 0; y < yRes; y++){
    for (int x = 0; x < xRes; x++){
      VEC3 d = generateRay(x,y);
      int index = (yRes -1 - y) * xRes + x;
      float red = 1.0f * abs(d[0]);
      values[index*3] = red*255.0f;
      values[index*3+1] = 0.0;
      values[index*3+2] = 0.0;
    }
  }
   writePPM("1xabs.ppm", xRes, yRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// 1y.ppm
//////////////////////////////////////////////////////////////////////////////////
void make_1y_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];
  for (int y = 0; y < yRes; y++){
    for (int x = 0; x < xRes; x++){
      VEC3 d = generateRay(x,y);
      int index = (yRes -1 - y) * xRes + x;
      float green = 1.0f * d[1];
      values[index*3] = 0.0;
      values[index*3+1] = green * 255.0f;
      values[index*3+2] = 0.0;
    }
  }
   writePPM("1y.ppm", xRes, yRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// 1yabs.ppm
//////////////////////////////////////////////////////////////////////////////////
void make_1yabs_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];
  for (int y = 0; y < yRes; y++){
    for (int x = 0; x < xRes; x++){
      VEC3 d = generateRay(x,y);
      int index = (yRes -1 - y) * xRes + x;
      float green = 1.0f * abs(d[1]);
      values[index*3] = 0.0;
      values[index*3+1] = green * 255.0f;
      values[index*3+2] = 0.0;
    }
  }
  writePPM("1yabs.ppm", xRes, yRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// 2.ppm
//////////////////////////////////////////////////////////////////////////////////
void make_2_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];
  for (int y = 0; y < yRes; y++){
    for (int x = 0; x < xRes; x++){
      VEC3 d = generateRay(x,y);
      int index = (yRes -1 - y) * xRes + x;
      double t2 = 0;
      if (checkSphereIntersection(center_s_0, d, eye, radius_0, &t2)){
        values[index*3] = color_s_0[0] * 255.0f;
        values[index*3+1] = color_s_0[1] * 255.0f;
        values[index*3+2] = color_s_0[2] * 255.0f;
      }else if (checkSphereIntersection(center_s_1, d, eye, radius_1, &t2)) {
        values[index * 3] = color_s_1[0] * 255.0f;
        values[index * 3 + 1] = color_s_1[1] * 255.0f;
        values[index * 3 + 2] = color_s_1[2] * 255.0f;
      } else if (checkSphereIntersection(center_s_2, d, eye, radius_2, &t2)) {
        values[index * 3] = color_s_2[0] * 255.0f;
        values[index * 3 + 1] = color_s_2[1] * 255.0f;
        values[index * 3 + 2] = color_s_2[2] * 255.0f;
      }
    }
  }
  writePPM("2.ppm", xRes, yRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// 3.ppm
// shading
//////////////////////////////////////////////////////////////////////////////////
void make_3_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  //diffuse shading
  for (int y = 0; y < yRes; y++){
      for (int x = 0; x < xRes; x++){
        VEC3 d = generateRay(x,y);
        int index = (yRes -1 - y) * xRes + x;
        double t3 = 0;
        VEC3 color(0.0, 0.0, 0.0);
        VEC3 sphereColor;
        if (checkSphereIntersection(center_s_0, d, eye, radius_0, &t3)){
          color = getShading(center_s_0, d, eye, radius_0, l_0_color, l_0_pos);
          sphereColor = color_s_0;
        }else if (checkSphereIntersection(center_s_1, d, eye, radius_1, &t3)) {
          color = getShading(center_s_1, d, eye, radius_1, l_0_color, l_0_pos);
          sphereColor = color_s_1;
        } else if (checkSphereIntersection(center_s_2, d, eye, radius_2, &t3)) {
          color = getShading(center_s_2, d, eye, radius_2, l_0_color, l_0_pos);
          sphereColor = color_s_2;
        }
        if(color != VEC3(0.0, 0.0, 0.0)){
          color[0] = color[0] * sphereColor[0];
          color[1] = color[1] * sphereColor[1];
          color[2] = color[2] * sphereColor[2];
          values[index * 3] = color[0] * 255.0f;
          values[index * 3 + 1] = color[1] * 255.0f;
          values[index * 3 + 2] = color[2] * 255.0f;
      }
    }
  }
  writePPM("3.ppm", xRes, yRes, values);
}
//////////////////////////////////////////////////////////////////////////////////
// 4.ppm
// multiple lights
//////////////////////////////////////////////////////////////////////////////////
void make_4_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;
  for (int y = 0; y < yRes; y++){
      for (int x = 0; x < xRes; x++){
        VEC3 d = generateRay(x,y);
        int index = (yRes -1 - y) * xRes + x;
        double t3 = 0;
        VEC3 color(0.0, 0.0, 0.0);
        VEC3 sphereColor;
        if (checkSphereIntersection(center_s_0, d, eye, radius_0, &t3)){
          color = getShading(center_s_0, d, eye, radius_0, l_0_color, l_0_pos)
                + getShading(center_s_0, d, eye, radius_0, l_1_color, l_1_pos);
          sphereColor = color_s_0;
        }else if (checkSphereIntersection(center_s_1, d, eye, radius_1, &t3)) {
          color = getShading(center_s_1, d, eye, radius_1, l_0_color, l_0_pos)
                + getShading(center_s_1, d, eye, radius_1, l_1_color, l_1_pos);
          sphereColor = color_s_1;
        } else if (checkSphereIntersection(center_s_2, d, eye, radius_2, &t3)) {
          color = getShading(center_s_2, d, eye, radius_2, l_0_color, l_0_pos)
                + getShading(center_s_2, d, eye, radius_2, l_1_color, l_1_pos);
          sphereColor = color_s_2;
        }
        if(color != VEC3(0.0, 0.0, 0.0)){
          color[0] = color[0] * sphereColor[0];
          color[1] = color[1] * sphereColor[1];
          color[2] = color[2] * sphereColor[2];
          values[index * 3] = color[0] * 255.0f;
          values[index * 3 + 1] = color[1] * 255.0f;
          values[index * 3 + 2] = color[2] * 255.0f;
      }
    }
  }
  writePPM("4.ppm", xRes, yRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// 5.ppm
// specular highlights
//////////////////////////////////////////////////////////////////////////////////
void make_5_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;
  for (int y = 0; y < yRes; y++){
      for (int x = 0; x < xRes; x++){
        VEC3 d = generateRay(x,y);
        int index = (yRes -1 - y) * xRes + x;
        double t3 = 0;
        VEC3 color(0.0, 0.0, 0.0);
        VEC3 sphereColor;
        if (checkSphereIntersection(center_s_0, d, eye, radius_0, &t3)){
          color = getShading(center_s_0, d, eye, radius_0, l_0_color, l_0_pos)
                + getShading(center_s_0, d, eye, radius_0, l_1_color, l_1_pos)
                + getHighlight(center_s_0, d, eye, radius_0, l_0_color, l_0_pos)
                + getHighlight(center_s_0, d, eye, radius_0, l_1_color, l_1_pos);

          sphereColor = color_s_0;
        }else if (checkSphereIntersection(center_s_1, d, eye, radius_1, &t3)) {
          color = getShading(center_s_1, d, eye, radius_1, l_0_color, l_0_pos)
                + getShading(center_s_1, d, eye, radius_1, l_1_color, l_1_pos)
                + getHighlight(center_s_1, d, eye, radius_1, l_0_color, l_0_pos)
                + getHighlight(center_s_1, d, eye, radius_1, l_1_color, l_1_pos);
          sphereColor = color_s_1;
        } else if (checkSphereIntersection(center_s_2, d, eye, radius_2, &t3)) {
          color = getShading(center_s_2, d, eye, radius_2, l_0_color, l_0_pos)
                + getShading(center_s_2, d, eye, radius_2, l_1_color, l_1_pos)
                + getHighlight(center_s_2, d, eye, radius_2, l_0_color, l_0_pos)
                + getHighlight(center_s_2, d, eye, radius_2, l_1_color, l_1_pos);
          sphereColor = color_s_2;
        }
        if(color != VEC3(0.0, 0.0, 0.0)){
          color[0] = color[0] * sphereColor[0];
          color[1] = color[1] * sphereColor[1];
          color[2] = color[2] * sphereColor[2];
          values[index * 3] = color[0] * 255.0f;
          values[index * 3 + 1] = color[1] * 255.0f;
          values[index * 3 + 2] = color[2] * 255.0f;
      }
    }
  }

  writePPM("5.ppm", xRes, yRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// 6.ppm
// Shadows
//////////////////////////////////////////////////////////////////////////////////
void make_6_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  // Instead of using initializer lists, we push back the elements one by one.
  vector<VEC3> centers;
  centers.push_back(center_s_2);
  centers.push_back(center_s_0);
  centers.push_back(center_s_1);

  vector<double> radii;
  radii.push_back(radius_2);
  radii.push_back(radius_0);
  radii.push_back(radius_1);

  vector<VEC3> lightPositions;
  lightPositions.push_back(l_0_pos);
  lightPositions.push_back(l_1_pos);

  vector<VEC3> lightColors;
  lightColors.push_back(l_0_color);
  lightColors.push_back(l_1_color);

  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      VEC3 d = generateRay(x, y);
      int index = (yRes - 1 - y) * xRes + x;

      VEC3 color(0.0, 0.0, 0.0);
      VEC3 closestSphereColor;
      double minT = INFINITY;  // Use this to track the closest intersection
      int closestSphereIndex = -1;

      // Check intersections for each sphere
      for (size_t i = 0; i < centers.size(); i++) {
        double t = 0;
        if (checkSphereIntersection(centers[i], d, eye, radii[i], &t) && t < minT) {
          minT = t;
          closestSphereIndex = i;
          closestSphereColor = (i == 0) ? color_s_2 : (i == 1) ? color_s_0 : color_s_1;
        }
      }

      // If we found the closest intersection
      if (closestSphereIndex != -1) {
        VEC3 point = eye + minT * d;  // Closest point of intersection

        // Check each light source for this closest intersection
        for (size_t j = 0; j < lightPositions.size(); j++) {
          if (!isInShadow(closestSphereIndex, point, lightPositions[j], centers, radii)) {
            color += getShading(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j])
                  + getHighlight(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j]);
          }
        }
        color[0] = color[0] * closestSphereColor[0];
        color[1] = color[1] * closestSphereColor[1];
        color[2] = color[2] * closestSphereColor[2];
        values[index * 3] = color[0] * 255.0f;
        values[index * 3 + 1] = color[1] * 255.0f;
        values[index * 3 + 2] = color[2] * 255.0f;
      }
    }
  }
   writePPM("6.ppm", xRes, yRes, values);
}

 //////////////////////////////////////////////////////////////////////////////////
// 7.ppm
// Shadows with wall of white spheres
// making sphere 0 mirror
//////////////////////////////////////////////////////////////////////////////////
void make_7_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];

  // Initialize all pixel colors to black
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  // Define the spheres
  vector<VEC3> centers;
  vector<double> radii;
  vector<VEC3> colors;
  centers.push_back(center_s_1);
  centers.push_back(center_s_0);
  centers.push_back(center_s_2);
  radii.push_back(radius_1);
  radii.push_back(radius_0);
  radii.push_back(radius_2);
  colors.push_back(color_s_1);
  colors.push_back(color_s_0);
  colors.push_back(color_s_2);

  // Define the light sources
  vector<VEC3> lightPositions;
  vector<VEC3> lightColors;
  lightPositions.push_back(l_2_pos);
  lightPositions.push_back(l_3_pos);
  lightColors.push_back(l_2_color);
  lightColors.push_back(l_3_color);

  // Define the wall of white spheres
  VEC3 startSphere(-20, -2, 20);
  double wallRadius = 1.0;
  VEC3 wallColor(1, 1, 1);
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 10; ++j) {
      centers.push_back(startSphere + VEC3(2 * wallRadius * i, 2 * wallRadius * j, 0));
      radii.push_back(wallRadius);
      colors.push_back(wallColor);
    }
  }

  // Render loop
  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      VEC3 d = generateRay(x, y);
      int index = (yRes - 1 - y) * xRes + x;

      VEC3 color(0.0, 0.0, 0.0);
      VEC3 closestSphereColor;
      double minTime = INFINITY;
      int closestSphereIndex = -1;

      // Check intersections for each sphere
      for (size_t i = 0; i < centers.size(); i++) {
        double t = 0;
        if (checkSphereIntersection(centers[i], d, eye, radii[i], &t) && t < minTime) {
          minTime = t;
          closestSphereIndex = i;
          closestSphereColor = colors[i];
        }
      }

      // If we found the closest intersection
      if (closestSphereIndex != -1) {
        VEC3 point = eye + minTime * d;

        if(closestSphereIndex == 1){
          VEC3 reflectionColor = getReflectionColor(d, point, centers, radii, colors, lightPositions, lightColors);
          values[index * 3] = reflectionColor[0] * 255.0f;
          values[index * 3 + 1] = reflectionColor[1] * 255.0f;
          values[index * 3 + 2] = reflectionColor[2] * 255.0f;
        }

        else{
          // Check each light source for this closest intersection
          for (size_t j = 0; j < lightPositions.size(); j++) {
            if (!isInShadow(closestSphereIndex, point, lightPositions[j], centers, radii)) {
              color += getShading(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j])
                    + getHighlight(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j]);
            }
          }
          color[0] = color[0] * closestSphereColor[0];
          color[1] = color[1] * closestSphereColor[1];
          color[2] = color[2] * closestSphereColor[2];
          values[index * 3] = color[0] * 255.0f;
          values[index * 3 + 1] = color[1] * 255.0f;
          values[index * 3 + 2] = color[2] * 255.0f;
        }
      }
    }
  }

  writePPM("7.ppm", xRes, yRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// 8.ppm
// Refraction
//////////////////////////////////////////////////////////////////////////////////
void make_8_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];

  // Initialize all pixel colors to black
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  // Define the spheres
  vector<VEC3> centers;
  vector<double> radii;
  vector<VEC3> colors;
  centers.push_back(center_s_1);
  centers.push_back(center_s_0);
  centers.push_back(center_s_2);
  radii.push_back(radius_1);
  radii.push_back(radius_0);
  radii.push_back(radius_2);
  colors.push_back(color_s_1);
  colors.push_back(color_s_0);
  colors.push_back(color_s_2);

  // Define the light sources
  vector<VEC3> lightPositions;
  vector<VEC3> lightColors;
  lightPositions.push_back(l_2_pos);
  lightPositions.push_back(l_3_pos);
  lightColors.push_back(l_2_color);
  lightColors.push_back(l_3_color);

  // Define the wall of white spheres
  VEC3 startSphere(-20, -2, 20);
  double wallRadius = 1.0;
  VEC3 wallColor(1, 1, 1);
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 10; ++j) {
      centers.push_back(startSphere + VEC3(2 * wallRadius * i, 2 * wallRadius * j, 0));
      radii.push_back(wallRadius);
      colors.push_back(wallColor);
    }
  }

  // Render loop
  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      VEC3 d = generateRay(x, y);
      int index = (yRes - 1 - y) * xRes + x;

      VEC3 color(0.0, 0.0, 0.0);
      VEC3 closestSphereColor;
      double minTime = INFINITY;
      int closestSphereIndex = -1;

      // Check intersections for each sphere
      for (size_t i = 0; i < centers.size(); i++) {
        double t = 0;
        if (checkSphereIntersection(centers[i], d, eye, radii[i], &t) && t < minTime) {
          minTime = t;
          closestSphereIndex = i;
          closestSphereColor = colors[i];
        }
      }

      // If we found the closest intersection
      if (closestSphereIndex != -1) {
        VEC3 point = eye + minTime * d;
        //if light hit the glass sphere
        if(closestSphereIndex == 1){
          VEC3 refractionColor = getRefractionColor(d, point, centers, radii, colors, lightPositions, lightColors);
          values[index * 3] = refractionColor[0] * 255.0f;
          values[index * 3 + 1] = refractionColor[1] * 255.0f;
          values[index * 3 + 2] = refractionColor[2] * 255.0f;
        }

        else{
          // Check each light source for this closest intersection
          for (size_t j = 0; j < lightPositions.size(); j++) {
            if (!isInShadow(closestSphereIndex, point, lightPositions[j], centers, radii)) {
              color += getShading(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j])
                    + getHighlight(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j]);
            }
          }
          color[0] = color[0] * closestSphereColor[0];
          color[1] = color[1] * closestSphereColor[1];
          color[2] = color[2] * closestSphereColor[2];
          values[index * 3] = color[0] * 255.0f;
          values[index * 3 + 1] = color[1] * 255.0f;
          values[index * 3 + 2] = color[2] * 255.0f;
        }
      }
    }
  }

  writePPM("8.ppm", xRes, yRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// 9.ppm
// Fresnel Effects
//////////////////////////////////////////////////////////////////////////////////
void make_9_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];

  // Initialize all pixel colors to black
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  // Define the spheres
  vector<VEC3> centers;
  vector<double> radii;
  vector<VEC3> colors;
  centers.push_back(center_s_1);
  centers.push_back(center_s_0);
  centers.push_back(center_s_2);
  radii.push_back(radius_1);
  radii.push_back(radius_0);
  radii.push_back(radius_2);
  colors.push_back(color_s_1);
  colors.push_back(color_s_0);
  colors.push_back(color_s_2);

  // Define the light sources
  vector<VEC3> lightPositions;
  vector<VEC3> lightColors;
  lightPositions.push_back(l_2_pos);
  lightPositions.push_back(l_3_pos);
  lightColors.push_back(l_2_color);
  lightColors.push_back(l_3_color);

  // Define the wall of white spheres
  VEC3 startSphere(-20, -2, 20);
  double wallRadius = 1.0;
  VEC3 wallColor(1, 1, 1);
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 10; ++j) {
      centers.push_back(startSphere + VEC3(2 * wallRadius * i, 2 * wallRadius * j, 0));
      radii.push_back(wallRadius);
      colors.push_back(wallColor);
    }
  }

  // Render loop
  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      VEC3 d = generateRay(x, y);
      int index = (yRes - 1 - y) * xRes + x;

      VEC3 color(0.0, 0.0, 0.0);
      VEC3 closestSphereColor;
      double minTime = INFINITY;
      int closestSphereIndex = -1;

      // Check intersections for each sphere
      for (size_t i = 0; i < centers.size(); i++) {
        double t = 0;
        if (checkSphereIntersection(centers[i], d, eye, radii[i], &t) && t < minTime) {
          minTime = t;
          closestSphereIndex = i;
          closestSphereColor = colors[i];
        }
      }

      // If we found the closest intersection
      if (closestSphereIndex != -1) {
        VEC3 point = eye + minTime * d;

        //if light hit the glass sphere
        if(closestSphereIndex == 1){

          VEC3 n_in = (point - centers[1]).normalized();
          pair<double, double> fresnel_k = calculateFresnel(d, n_in, n_air, n_glass);

          //REFRACTION:
          VEC3 refractionColor = getRefractionColor(d, point, centers, radii, colors, lightPositions, lightColors);

          //REFLECTION:
          VEC3 reflectionColor = getReflectionColor(d, point, centers, radii, colors, lightPositions, lightColors);

          values[index * 3] = (fresnel_k.first * reflectionColor[0] + fresnel_k.second * refractionColor[0]) * 255.0f;
          values[index * 3 + 1] = (fresnel_k.first * reflectionColor[1] + fresnel_k.second * refractionColor[1]) * 255.0f;
          values[index * 3 + 2] = (fresnel_k.first * reflectionColor[2] + fresnel_k.second * refractionColor[2]) * 255.0f;
        }

        else{
          // Check each light source for this closest intersection
          for (size_t j = 0; j < lightPositions.size(); j++) {
            if (!isInShadow(closestSphereIndex, point, lightPositions[j], centers, radii)) {
              color += getShading(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j])
                    + getHighlight(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j]);
            }
          }
          color[0] = color[0] * closestSphereColor[0];
          color[1] = color[1] * closestSphereColor[1];
          color[2] = color[2] * closestSphereColor[2];
          values[index * 3] = color[0] * 255.0f;
          values[index * 3 + 1] = color[1] * 255.0f;
          values[index * 3 + 2] = color[2] * 255.0f;
        }
      }
    }
  }
  writePPM("9.ppm", xRes, yRes, values);
}


//////////////////////////////////////////////////////////////////////////////////
// 10.ppm
// Ray-triangle intersection
//////////////////////////////////////////////////////////////////////////////////
void make_10_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];

  // Initialize all pixel colors to black
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  // Define the spheres
  vector<VEC3> centers;
  vector<double> radii;
  vector<VEC3> colors;

  centers.push_back(center_s_2);
  centers.push_back(center_s_0);
  radii.push_back(radius_2);
  radii.push_back(radius_0);
  colors.push_back(color_s_2);
  colors.push_back(color_s_0);

  // Define the light sources
  vector<VEC3> lightPositions;
  vector<VEC3> lightColors;
  lightPositions.push_back(l_2_pos);
  lightPositions.push_back(l_3_pos);
  lightColors.push_back(l_2_color);
  lightColors.push_back(l_3_color);

  // Define the wall of white spheres
  VEC3 startSphere(-20, -2, 20);
  double wallRadius = 1.0;
  VEC3 wallColor(1, 1, 1);
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 10; ++j) {
      centers.push_back(startSphere + VEC3(2 * wallRadius * i, 2 * wallRadius * j, 0));
      radii.push_back(wallRadius);
      colors.push_back(wallColor);
    }
  }

  // Define the triangles
  VEC3 vertex_a(0.5, 3, 10);
  VEC3 vertex_b(0.5, -3, 10);
  VEC3 vertex_c(6.5, 3, 10);
  VEC3 vertex_d(6.5, -3, 10);

  // Translation to origin
  VEC3 center(3.5, 0, 10);
  vertex_a -= center;
  vertex_b -= center;
  vertex_c -= center;
  vertex_d -= center;

  // Rotate 45 degrees around the y axis
  MATRIX3 R_y_45;
  R_y_45.setZero();
  R_y_45(0, 0) = cos(M_PI_4); R_y_45(0, 2) = sin(M_PI_4);
  R_y_45(1, 1) = 1;
  R_y_45(2, 0) = -sin(M_PI_4); R_y_45(2, 2) = cos(M_PI_4);

  // Apply rotation
  vertex_a = R_y_45 * vertex_a;
  vertex_b = R_y_45 * vertex_b;
  vertex_c = R_y_45 * vertex_c;
  vertex_d = R_y_45 * vertex_d;

  // Translate back
  vertex_a += center;
  vertex_b += center;
  vertex_c += center;
  vertex_d += center;

  //triangle abc
  vector<VEC3> tri_abc;
  tri_abc.push_back(vertex_a);
  tri_abc.push_back(vertex_c);
  tri_abc.push_back(vertex_b);
  VEC3 color_tri_abc (0.0, 1.0, 0.0);

  //triangle bcd
  vector<VEC3> tri_bcd;
  tri_bcd.push_back(vertex_b);
  tri_bcd.push_back(vertex_c);
  tri_bcd.push_back(vertex_d);
  VEC3 color_tri_bcd (1.0, 0.0, 0.0);

  //vector of triangles
  vector<vector <VEC3> > triangles;
  vector<VEC3> triangle_colors;
  triangles.push_back(tri_abc);
  triangles.push_back(tri_bcd);
  triangle_colors.push_back (color_tri_abc);
  triangle_colors.push_back (color_tri_bcd);

  // Render loop
  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      VEC3 d = generateRay(x, y);
      int index = (yRes - 1 - y) * xRes + x;

      VEC3 color(0.0, 0.0, 0.0);
      VEC3 closestShapeColor;
      double minTime = INFINITY;
      int closestSphereIndex = -1;
      int closestTriangleIndex = -1;

      // Check intersections for each sphere
      for (size_t i = 0; i < centers.size(); i++) {
        double t = 0;
        if (checkSphereIntersection(centers[i], d, eye, radii[i], &t) && t < minTime) {
          minTime = t;
          closestSphereIndex = i;
          closestShapeColor = colors[i];
        }
      }

      // Check intersections for each triangle
      for (size_t i = 0; i < triangles.size(); i++) {
        double t = 0;
        if(checkTriangleIntersection(triangles[i], d, eye, &t) && t < minTime){
          minTime = t;
          closestSphereIndex = -1;
          closestTriangleIndex = i;
          closestShapeColor = triangle_colors[i];
        }
      }

      // If we found the closest sphere intersection
      if (closestSphereIndex != -1) {
        double epsilon = 1e-4;
        VEC3 point = eye + (minTime - epsilon) * d;

        //if light hit the glass sphere
        if(closestSphereIndex == 1){

          VEC3 n_in = (point - centers[1]).normalized();
          pair<double, double> fresnel_k = calculateFresnel(d, n_in, n_air, n_glass);

          //REFRACTION:
          VEC3 refractionColor = getRefractionColor(d, point, centers, radii, colors, lightPositions, lightColors);

          //REFLECTION:
          VEC3 reflectionColor = getReflectionColor(d, point, centers, radii, colors, lightPositions, lightColors);

          values[index * 3] = (fresnel_k.first * reflectionColor[0] + fresnel_k.second * refractionColor[0]) * 255.0f;
          values[index * 3 + 1] = (fresnel_k.first * reflectionColor[1] + fresnel_k.second * refractionColor[1]) * 255.0f;
          values[index * 3 + 2] = (fresnel_k.first * reflectionColor[2] + fresnel_k.second * refractionColor[2]) * 255.0f;
        }

        else{
          // Check each light source for this closest intersection
          for (size_t j = 0; j < lightPositions.size(); j++) {
            if (!isInSphereShadow(true, closestSphereIndex, point, lightPositions[j], centers, radii)
              && !isInTriangleShadow(false, closestSphereIndex, point, lightPositions[j], triangles)) {
              color += getShading(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j])
                    + getHighlight(centers[closestSphereIndex], d, eye, radii[closestSphereIndex], lightColors[j], lightPositions[j]);
            }
          }
          color[0] = color[0] * closestShapeColor[0];
          color[1] = color[1] * closestShapeColor[1];
          color[2] = color[2] * closestShapeColor[2];
          values[index * 3] = color[0] * 255.0f;
          values[index * 3 + 1] = color[1] * 255.0f;
          values[index * 3 + 2] = color[2] * 255.0f;
        }
      }

      //if we find the closest triangle intersection
      if(closestTriangleIndex != -1){
        VEC3 point = eye + minTime * d;
        for (size_t j = 0; j < lightPositions.size(); j++) {
            if (!isInSphereShadow(false, closestSphereIndex, point, lightPositions[j], centers, radii)
              && !isInTriangleShadow(false, closestSphereIndex, point, lightPositions[j], triangles)) {
              // color[0] = 1;
              // color[1] = 1;
              // color[2] = 1;
              color += getTriangleShading(triangles[closestTriangleIndex], d, eye, lightColors[j], lightPositions[j])
                    + getTriangleHighlight(triangles[closestTriangleIndex], d, eye, lightColors[j], lightPositions[j]);
            }
          }
          color[0] = color[0] * closestShapeColor[0];
          color[1] = color[1] * closestShapeColor[1];
          color[2] = color[2] * closestShapeColor[2];
          values[index * 3] = color[0] * 255.0f;
          values[index * 3 + 1] = color[1] * 255.0f;
          values[index * 3 + 2] = color[2] * 255.0f;
      }

    }
  }
  writePPM("10.ppm", xRes, yRes, values);
}

VEC3 recursion (VEC3 origin, VEC3 d, int& count, int index,
                vector<VEC3>& centers, vector<double>& radii, vector<VEC3>& colors,
                vector<vector<VEC3> > &triangles,
                vector<VEC3>& lightPositions, vector<VEC3>& lightColors){

  VEC3 color(0.0, 0.0, 0.0);
  VEC3 closestShapeColor;
  double minTime = INFINITY;
  int closestSphereIndex = -1;
  int closestTriangleIndex = -1;

  // Check intersections for each sphere
  for (size_t i = 0; i < centers.size(); i++) {
    double t = 0;
    if (checkSphereIntersection(centers[i], d, origin, radii[i], &t) && t < minTime) {
      minTime = t;
      closestSphereIndex = i;
      closestShapeColor = colors[i];
    }
  }

  // Check intersections for each triangle
  for (size_t i = 0; i < triangles.size(); i++) {
    double t = 0;
    if(checkTriangleIntersection(triangles[i], d, origin, &t) && t < minTime){
      minTime = t;
      closestSphereIndex = -1;
      closestTriangleIndex = i;
    }
  }

  // If we found the closest sphere intersection
  if (closestSphereIndex != -1) {
    VEC3 point = origin + minTime * d;

    //if light hit the glass sphere
    if(closestSphereIndex == 1){

      VEC3 n_in = (point - centers[1]).normalized();
      pair<double, double> fresnel_k = calculateFresnel(d, n_in, n_air, n_glass);

      VEC3 zero (0.0, 0.0, 0.0);

      //REFRACTION:
      VEC3 newEye1;
      VEC3 newDir1;
      VEC3 refractionColor = sphere_getRefractionColor_recursive(d, point, centers, radii, colors, triangles, lightPositions, lightColors, newEye1, newDir1);
      //REFLECTION:
      VEC3 newEye2;
      VEC3 newDir2;
      VEC3 reflectionColor = getReflectionColor(d, point, centers, radii, colors, lightPositions, lightColors);

      if(count == 9){
        color[0] = fresnel_k.first * reflectionColor[0] + fresnel_k.second * refractionColor[0];
        color[1] = fresnel_k.first * reflectionColor[1] + fresnel_k.second * refractionColor[1];
        color[2] = fresnel_k.first * reflectionColor[2] + fresnel_k.second * refractionColor[2];
        return color;
      }

      if(newEye1 != zero){
        count += 1;
        refractionColor = recursion(newEye1, newDir1, count, index, centers, radii, colors,
                  triangles, lightPositions, lightColors);
        }
      // VEC3 refractionColor = getRefractionColor(d, point, centers, radii, colors, lightPositions, lightColors);

      if(newEye1 != zero){
        count += 1;
        reflectionColor = recursion(newEye2, newDir2, count, index, centers, radii, colors,
                  triangles, lightPositions, lightColors);
        }

      color[0] = fresnel_k.first * reflectionColor[0] + fresnel_k.second * refractionColor[0];
      color[1] = fresnel_k.first * reflectionColor[1] + fresnel_k.second * refractionColor[1];
      color[2] = fresnel_k.first * reflectionColor[2] + fresnel_k.second * refractionColor[2];
      return color;
    }
    else{
      // Check each light source for this closest intersection
      for (size_t j = 0; j < lightPositions.size(); j++) {
        if (!isInSphereShadow(true, closestSphereIndex, point, lightPositions[j], centers, radii)
          && !isInTriangleShadow(false, closestTriangleIndex, point, lightPositions[j], triangles)) {
          color += getShading(centers[closestSphereIndex], d, origin, radii[closestSphereIndex], lightColors[j], lightPositions[j])
                + getHighlight(centers[closestSphereIndex], d, origin, radii[closestSphereIndex], lightColors[j], lightPositions[j]);
        }
      }
      color[0] = color[0] * closestShapeColor[0];
      color[1] = color[1] * closestShapeColor[1];
      color[2] = color[2] * closestShapeColor[2];
      return color;
    }

    if(closestTriangleIndex != -1){
      VEC3 newEye;
      VEC3 newDir;
      VEC3 triangle_reflection_color = triangle_getReflectionColor_recursive(d, origin, closestTriangleIndex, centers, radii, colors, triangles, lightPositions, lightColors, newEye, newDir);
      if(count == 9){
        return triangle_reflection_color;
      }
      count = count + 1;
      color = recursion(newEye, newDir, count, index, centers, radii, colors, triangles, lightPositions, lightColors);
      return color;
    }
  }

  //if we find the closest triangle intersection
  if(closestTriangleIndex != -1){
    VEC3 point = eye + minTime * d;
    for (size_t j = 0; j < lightPositions.size(); j++) {
        if (!isInSphereShadow(false, closestSphereIndex, point, lightPositions[j], centers, radii)
          && !isInTriangleShadow(true, closestSphereIndex, point, lightPositions[j], triangles)) {

          color += getTriangleShading(triangles[closestTriangleIndex], d, origin, lightColors[j], lightPositions[j])
                + getTriangleHighlight(triangles[closestTriangleIndex], d, origin, lightColors[j], lightPositions[j]);
        }
      }
  }

  return color;
}
//////////////////////////////////////////////////////////////////////////////////
// 11.ppm
// reflection recursion
//////////////////////////////////////////////////////////////////////////////////
void make_11_ppm(){
  int totalCells = xRes * yRes;
  float* values = new float[3 * totalCells];

  // Initialize all pixel colors to black
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  // Define the spheres
  vector<VEC3> centers;
  vector<double> radii;
  vector<VEC3> colors;

  centers.push_back(center_s_2);
  centers.push_back(center_s_0);
  radii.push_back(radius_2);
  radii.push_back(radius_0);
  colors.push_back(color_s_2);
  colors.push_back(color_s_0);

  // Define the light sources
  vector<VEC3> lightPositions;
  vector<VEC3> lightColors;
  lightPositions.push_back(l_2_pos);
  lightPositions.push_back(l_3_pos);
  lightColors.push_back(l_2_color);
  lightColors.push_back(l_3_color);

  // Define the wall of white spheres
  VEC3 startSphere(-20, -2, 20);
  double wallRadius = 1.0;
  VEC3 wallColor(1, 1, 1);
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 10; ++j) {
      centers.push_back(startSphere + VEC3(2 * wallRadius * i, 2 * wallRadius * j, 0));
      radii.push_back(wallRadius);
      colors.push_back(wallColor);
    }
  }

  // Define the triangles
  VEC3 vertex_a(0.5, 3, 10);
  VEC3 vertex_b(0.5, -3, 10);
  VEC3 vertex_c(6.5, 3, 10);
  VEC3 vertex_d(6.5, -3, 10);

  // Translation to origin
  VEC3 center(3.5, 0, 10);
  vertex_a -= center;
  vertex_b -= center;
  vertex_c -= center;
  vertex_d -= center;

  // Rotate 45 degrees around the y axis
  MATRIX3 R_y_45;
  R_y_45.setZero();
  R_y_45(0, 0) = cos(M_PI_4); R_y_45(0, 2) = sin(M_PI_4);
  R_y_45(1, 1) = 1;
  R_y_45(2, 0) = -sin(M_PI_4); R_y_45(2, 2) = cos(M_PI_4);

  // Apply rotation
  vertex_a = R_y_45 * vertex_a;
  vertex_b = R_y_45 * vertex_b;
  vertex_c = R_y_45 * vertex_c;
  vertex_d = R_y_45 * vertex_d;

  // Translate back
  vertex_a += center;
  vertex_b += center;
  vertex_c += center;
  vertex_d += center;

  //triangle abc
  vector<VEC3> tri_abc;
  tri_abc.push_back(vertex_a);
  tri_abc.push_back(vertex_c);
  tri_abc.push_back(vertex_b);
  VEC3 color_tri_abc (0.0, 1.0, 0.0);

  //triangle bcd
  vector<VEC3> tri_bcd;
  tri_bcd.push_back(vertex_b);
  tri_bcd.push_back(vertex_c);
  tri_bcd.push_back(vertex_d);
  VEC3 color_tri_bcd (1.0, 0.0, 0.0);

  //vector of triangles
  vector<vector <VEC3> > triangles;
  triangles.push_back(tri_abc);
  triangles.push_back(tri_bcd);

  // Render loop
  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      VEC3 d = generateRay(x, y);
      VEC3 origin = eye;
      int index = (yRes - 1 - y) * xRes + x;
      VEC3 color;
      int count = 0;
      color = recursion(origin, d, count, index, centers, radii, colors, triangles, lightPositions, lightColors);
      values[index * 3] = color[0] * 255.0f;
      values[index * 3 + 1] = color[1] * 255.0f;
      values[index * 3 + 2] = color[2] * 255.0f;
      }
    }

  writePPM("11.ppm", xRes, yRes, values);
}
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  make_1x_ppm();
  make_1xabs_ppm();
  make_1y_ppm();
  make_1yabs_ppm();
  make_2_ppm();
  make_3_ppm();
  make_4_ppm();
  make_5_ppm();
  make_6_ppm();
  make_7_ppm();
  make_8_ppm();
  make_9_ppm();
  make_10_ppm();
  make_11_ppm();
  return 0;
}
