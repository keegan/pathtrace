#include <iostream>
#include <windows.h> // for bitmap headers
#include <atomic>
#include <vector>
#include <thread>
#include <array>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>


const size_t c_imageWidth = 512;
const size_t c_imageHeight = 512;
// sampling parameters
const size_t c_samplesPerPixel = 1000;
const size_t c_numBounces = 5;
const float c_rayBounceEpsilon = 0.001f;

//=================================================================================
// Derived values
const size_t c_numPixels = c_imageWidth * c_imageHeight;
const float c_aspectRatio = float(c_imageWidth) / float(c_imageHeight);
//const float c_cameraHorizFOV = c_cameraVerticalFOV * c_aspectRatio;
//const float c_windowTop = tan(c_cameraVerticalFOV / 2.0f) * c_nearPlaneDistance;
//const float c_windowRight = tan(c_cameraHorizFOV / 2.0f) * c_nearPlaneDistance;
//const TVector3 c_cameraFwd = Normalize(c_cameraLookAt - c_cameraPos);
//const TVector3 c_cameraRight = Cross({ 0.0f, 1.0f, 0.0f }, c_cameraFwd);
// TVector3 c_cameraUp = Cross(c_cameraFwd, c_cameraRight);

//=================================================================================

typedef std::array<uint8_t, 3> TPixelBGRU8;
typedef std::array<float, 3> TPixelRGBF32;

std::vector<TPixelRGBF32> g_pixels;



struct Vector {
    double x,y,z;
    Vector(double x = 0, double y = 0, double z = 0){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Vector operator+(const Vector &other){
        return Vector(x + other.x, y + other.y, z + other.z);
    }
    Vector operator-(const Vector &other){
        return Vector(x - other.x, y - other.y, z - other.z);
    }
    Vector operator*(double factor){
        return Vector(x * factor, y * factor, z * factor);
    }
    double dot(const Vector &other){
        return (x * other.x) + (y * other.y) + (z * other.z);
    }
    Vector norm(){
        return *this = *this * (1.0 / sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z)));
    }
};

Vector randHemisphere(const Vector &n);

const Vector cam_pos = Vector(0, 0, -10);
const Vector cam_dir = Vector(0, 0, 1);

struct Ray {
    Vector start;
    Vector dir;
    Ray(Vector start = Vector(), Vector end = Vector()){
        this->start = start;
        this->dir = (end - start).norm();
    }
};

struct Color {
    double r, g, b;
    Color(double r = 0, double g = 0, double b = 0){
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color operator+(const Color &other){
        return Color(r + other.r, g + other.g, b + other.b);
    }
    Color operator*(double factor){
        return Color(r * factor, g * factor, b * factor);
    }
    Color operator*(const Color &other){
        return Color(r * other.r, g * other.g, b * other.b);
    }
};



class Thing {
public:
    Vector pos;
    Vector* shine_direction;
    Color em; // emmisivity (abbreviated bc i cant spell it)
    Color c; // object's actual (reflectance) color
    virtual float intersect(const Ray &r) = 0;
    //virtual Color colorAt(const Ray &r, float t) = 0;
    virtual Color emission(const Vector &point, Vector &incoming) = 0;
    virtual Vector normal(const Vector &point) = 0;
    virtual Vector randReflect(const Vector &point, Vector &incoming) = 0;
};

class Sphere: public Thing {
public:
    double r, r2;
    bool mirror;
    Sphere(Vector p, Color e, Color c, double radius, bool mirror=false, Vector* shine=NULL){
        this->pos = p;
        this->em = e;
        this->c = c;
        this->r = radius;
        this->r2 = r*r;
        this->mirror = mirror;
        this->shine_direction=shine;
    }
    float intersect(const Ray &r){
        Vector L = this->pos - r.start; // draw a vector from Ray start to sphere center
        float t_CA = L.dot(r.dir);


        float d2 = L.dot(L) - t_CA * t_CA;
        //if (d2 > r2) { return 0; } // 0 means no intersect (or object behind camera)
        float t_HC = sqrt(r2 - d2);
        float t0 = t_CA - t_HC;
        float t1 = t_CA + t_HC;

        if(t0 > t1){ // we want the closer intersection
            std::swap(t0, t1);
        }

        if (t0 < 0){ //if inside the sphere, we want the face in front of us
            t0 = t1;
            if(t0 < 0){
                return 0; // if both are negative, the whole sphere is behind us
            }
        }

        return t0; // return the time to point of intersection
    }
    virtual Color emission(const Vector &point, Vector &incoming){
        if(this->shine_direction == NULL){
            return this->em;
        }
        Vector copy = point;
        float cos2 = (copy - this-> pos).dot(*this->shine_direction);
        float cos = incoming.norm().dot(*this->shine_direction);
        if(cos2 < 0){ return Color(); }
        return this->em * cos2;
    }
    Vector normal(const Vector &point){
        return (this->pos - point).norm() * -1;
    }
    Vector randReflect(const Vector &point, Vector &incoming){
        Vector normal = this->normal(point);
        if(this->mirror){
            return incoming - (normal * (incoming.dot(normal) * 2));
        }
        return randHemisphere(normal);
    }

   /* Color colorAt(const Ray &r, float t){
        //At each stage of the recursion, we return: EmittedLight + 2 * RecursiveLight * Dot(Normal, RandomHemisphereAngle) * SurfaceDiffuseColor.
        Color result;
        result = result + this->em;


        return result;
    }*/
};

std::vector<Thing*> things;


inline float Clamp(float v, float min, float max)
{
    if (v < min)
        return min;
    else if (v > max)
        return max;
    else
        return v;
}

Vector randHemisphere(const Vector &n){
    float x = (float(rand()) / (RAND_MAX / 2)) - 1;
    float y = (float(rand()) / (RAND_MAX / 2)) - 1;
    float z = (float(rand()) / (RAND_MAX / 2)) - 1;

    Vector dir = Vector(x, y, z).norm();

    if(dir.dot(n) < 0){
        dir = dir * -1;
    }

    return dir;
}

Color tracePath(Ray r, int depth){
    //At each stage of the recursion, we return: EmittedLight + 2 * RecursiveLight * Dot(Normal, RandomHemisphereAngle) * SurfaceDiffuseColor.


    Color result = {0, 0, 0};

    if(depth > c_numBounces){
        return result; //recursion max
    }

    Thing* closest_thing = NULL;
    float min_time = 1e25; // a big magic number. eww TODO: not this
    for(Thing* t : things){
        float time = t->intersect(r);
        //std::cout << time;
        if(time > 0){ // we have a hit
            if(time < min_time){
                min_time = time;
                closest_thing = t;
            }
        }
    }

    if(closest_thing != NULL){
        //Color c = closest_thing->em;
        Ray reflect;
        reflect.start = r.start + (r.dir * min_time); // point we reflect from
        Vector thing_normal = closest_thing->normal(reflect.start);

        reflect.dir = closest_thing->randReflect(reflect.start, r.dir);

        reflect.start = reflect.start + (reflect.dir * 0.01); // epsilon for floating point errors (stay outside object)

        Color reflected_color = tracePath(reflect, depth + 1);
        float cos_theta = reflect.dir.dot(thing_normal);

        result = closest_thing->emission(reflect.start, reflect.dir);
        result = result + (closest_thing->c * reflected_color * (2 * cos_theta));

    }
    return result;
}

/*
double dist(vector a, vector b){
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;

    return sqrt(dx*dx + dy*dy + dz*dz);
}
 */

bool SaveImage (const char* fileName)
{
    // allocate memory for our bitmap BGR U8 image

    std::vector<TPixelBGRU8> outPixels;
    outPixels.resize(c_numPixels);

    //int* outPixels = new int(c_numPixels);

    // convert from RGB F32 to BGR U8
    for (size_t i = 0; i < c_numPixels; ++i)
    {
        const TPixelRGBF32& src = g_pixels[i];
        TPixelBGRU8& dest = outPixels[i];

        // apply gamma correction
        TPixelRGBF32 correctedPixel;
        correctedPixel[0] = powf(src[0], 1.0f / 2.2f);
        correctedPixel[1] = powf(src[1], 1.0f / 2.2f);
        correctedPixel[2] = powf(src[2], 1.0f / 2.2f);

        // clamp and convert
        dest[0] = uint8_t(Clamp(correctedPixel[2] * 255.0f, 0.0f, 255.0f));
        dest[1] = uint8_t(Clamp(correctedPixel[1] * 255.0f, 0.0f, 255.0f));
        dest[2] = uint8_t(Clamp(correctedPixel[0] * 255.0f, 0.0f, 255.0f));
    }

    // write the bitmap

    // open the file if we can
    FILE *file;
    fopen_s(&file, fileName, "wb");
    if (!file)
        return false;

    // make the header info
    BITMAPFILEHEADER header;
    BITMAPINFOHEADER infoHeader;

    header.bfType = 0x4D42;
    header.bfReserved1 = 0;
    header.bfReserved2 = 0;
    header.bfOffBits = 54;

    infoHeader.biSize = 40;
    infoHeader.biWidth = (LONG)c_imageWidth;
    infoHeader.biHeight = (LONG)c_imageHeight;
    infoHeader.biPlanes = 1;
    infoHeader.biBitCount = 24;
    infoHeader.biCompression = 0;
    infoHeader.biSizeImage = (DWORD)c_numPixels*3;
    infoHeader.biXPelsPerMeter = 0;
    infoHeader.biYPelsPerMeter = 0;
    infoHeader.biClrUsed = 0;
    infoHeader.biClrImportant = 0;

    header.bfSize = infoHeader.biSizeImage + header.bfOffBits;

    // write the data and close the file
    fwrite(&header, sizeof(header), 1, file);
    fwrite(&infoHeader, sizeof(infoHeader), 1, file);
    fwrite(&outPixels[0], infoHeader.biSizeImage, 1, file);
    fclose(file);
    return true;
}

TPixelRGBF32 renderPixel(float x, float y){
    Color acc;

    for(int i = 0; i < c_samplesPerPixel; i ++) {
        float jitter_x = (float(rand()) / (RAND_MAX / 2)) - 1.0;
        float jitter_y = (float(rand()) / (RAND_MAX / 2)) - 1.0;
        jitter_x *= 0.5;
        jitter_y *= 0.5;
        float real_x = -5 + ((x + jitter_x) * 10 / c_imageWidth);
        float real_y = -5 + ((y + jitter_y) * 10 / c_imageHeight);
        Vector pixel_pos = Vector(real_x, real_y, 0);
        Color c = tracePath(Ray(cam_pos, pixel_pos), 0);

        acc = acc + (c * (1.0 / c_samplesPerPixel));
    }

    return {float(acc.r), float(acc.g), float(acc.b)};

    return {0, 0, 0};
}

int main() {
    g_pixels.resize(c_numPixels);
    things.push_back(new Sphere(Vector(1, 4, 7), Color(1, 0.5, 0), Color(0.3, 0.3, 0.1), 1, false, new Vector(1, 0, 0)));
    things.push_back(new Sphere(Vector(0.1, 4, 7), Color(), Color(1, 1, 1), sqrt(2), true, new Vector(1, 0, 0)));
    things.push_back(new Sphere(Vector(5, 0, 7), Color(0, .1, 1), Color(1, 1, 1), 2));
    things.push_back(new Sphere(Vector(0, 0, 5), Color(1, .3, .3), Color(0.5,1,1), 3));
    things.push_back(new Sphere(Vector(-5, -3, 7), Color(0, 0, .1), Color(1, 1, 1), 2));
    things.push_back(new Sphere(Vector(0, -2003, 0), Color(.2, .2, .2), Color(1, 1, 1), 2000));
    things.push_back(new Sphere(Vector(-2005, 0, 0), Color(), Color(.8, 1, .8), 2000));
    things.push_back(new Sphere(Vector(2005, 0, 0), Color(), Color(0, 1, 1), 2000, true));
    things.push_back(new Sphere(Vector(0, 2009, 0), Color(), Color(0, 1, 1), 2000));
    things.push_back(new Sphere(Vector(0, 0, 2010), Color(), Color(1, 1, 1), 2000));


    for(int i = 0; i < c_numPixels; i ++){
        g_pixels[i] = {1, 1, 0};
    }


    for(int x = 0; x < c_imageWidth; x ++){
        for(int y = 0; y < c_imageHeight; y ++) {
            g_pixels[(y * c_imageWidth) + x] = renderPixel(x, y);
        }
        std::cout << x;
    }

    std::cout << "Hello, World!" << std::endl;

    std::cout << SaveImage("out.bmp");

    return 0;
}