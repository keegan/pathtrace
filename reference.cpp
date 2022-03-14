//
// Torbert, 8 February 2016
// K. Lanzillotta 2019
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//
#define M 600 //640
#define N 400 //480
#define AMBIENT 0.4

#define SQUARE 4

#define NONE -1

#define G 0.3
#define dt 0.0333333333333

#define TIME 5
//

typedef struct vector {
    double x;
    double y;
    double z;
} vector;


typedef enum Color {RED=0, BLUE=1, GREEN=2, FLOOR=3, BACK=5} Color;

typedef struct sphere {
    double r;
    vector v;
    enum Color c;
} sphere;

vector e = { 0.50 , 0.50 , -1.00 } ; // the eye
vector g = { 0.00 , 2.25 , -0.50 } ; // the light

const int c_vals[6][3] = {
        {211, 79, 23},      //RED
        {24, 102, 219},      //BLUE
        {124, 186, 67},      //GREEN
        {102, 127, 65},   //FLOOR0
        {122, 183, 31},    //FLOOR2
        {178, 247, 243}         //BACK
};

double mix[4];

void init(){
    mix[RED] = -1;
    mix[BLUE] = 0.3;
    mix[GREEN] = 0.001;
    mix[FLOOR] = 0.6;
}

sphere arr[4] = {
        {
                0.50,
                {0.00, 0.75, 1.25},
                RED
        },
        {
                3.00,
                {0.50, 1.50, 5.00},
                BLUE
        },
        {
                0.25,
                {1.00, 1.25, 0.50},
                GREEN
        },
        {
                20000.25,
                {0.50, -20000.00, 0.50},
                FLOOR
        }
};


void diff( vector* t , vector u , vector v ) // t = u - v
{
    t->x = u.x - v.x ;
    t->y = u.y - v.y ;
    t->z = u.z - v.z ;
}

double dot(vector* a, vector* b){
    double out = a->x * b->x;
    out += a->y * b->y;
    out += a->z * b->z;

    return out;
}

double dist(vector a, vector b){
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;

    return sqrt(dx*dx + dy*dy + dz*dz);
}

int* findcolor(int x, int y);
double intersect(vector* ray, sphere* dest, vector start);
int* castray(vector* ray, vector start);

int main(int argc, char** argv)
{
    init();
    //int* rgb[N][M]; // red-green-blue for each pixel
    double v = 0;
    int counter = 0;
    for(double t = 0; t < TIME; t += dt){
        arr[GREEN].v.y += v * dt;
        v -= G * dt;
        double d = dist(arr[GREEN].v, arr[FLOOR].v);
        if(d <= arr[FLOOR].r + arr[GREEN].r){
            v *= -1;
        }

        int*** rgb = malloc(N * sizeof(int*));
        for(int i = 0; i < M; i ++){
            rgb[i] = malloc(M * sizeof(int*));
        }
        //
        int y , x ;
        //
        FILE* fout ;
        //
        for( y = 0 ; y < N ; y++ )
        {
            for( x = 0 ; x < M ; x++)
            {
                rgb[y][x] = findcolor(x, y);
            }
        }
        //
        //
        //
        int length = snprintf(NULL, 0,"frames/f%04d.ppm", counter);
        char filename[length + 1];
        sprintf(filename, "frames/f%04d.ppm", counter);
        fout = fopen(filename , "w" ) ;
        //
        fprintf( fout , "P3\n" ) ;
        fprintf( fout , "%d %d\n" , M , N ) ;
        fprintf( fout , "255\n" ) ;
        //
        for( y = 0 ; y < N ; y++ )
        {
            for( x = 0 ; x < M ; x++)
            {
                fprintf( fout , "%d %d %d\n" ,
                         rgb[y][x][0] , rgb[y][x][1] , rgb[y][x][2] ) ;
            }
        }
        fclose( fout ) ;
        counter ++;
    }
    //
    return 0 ;
}

void normalize(vector* v){
    double mag = v->x * v->x;
    mag += v->y * v->y;
    mag += v->z * v->z;
    mag = sqrt(mag);

    v->x /= mag;
    v->y /= mag;
    v->z /= mag;
}

int* findcolor(int intx, int inty){
    double x, y;
    int r=0, g=0, b=0;
    int* out = malloc(3 * sizeof(int));
    for(int i = 0; i < SQUARE; i ++){
        for(int j = 0; j < SQUARE; j ++){
            x = intx + ((double)i / SQUARE) + (0.5 / SQUARE);
            x /= N;
            x -= 1.0/6;

            y = inty + ((double)j / SQUARE) + (0.5 / SQUARE);
            y /= N;
            y = 1.0 - y;

            vector* ray = malloc(sizeof(vector));
            vector point = {x, y, 0};

            diff(ray, point, e);

            normalize(ray);

            int* res = castray(ray, e);
            r += res[0];
            g += res[1];
            b += res[2];
        }
    }
    out[0] = r / (SQUARE * SQUARE);
    out[1] = g / (SQUARE * SQUARE);
    out[2] = b / (SQUARE * SQUARE);
    return out;
}
int* castray(vector* ray, vector start){
    int pos = 0;

    int min_color = BACK;
    int* out;
    double t = 1.2e+25;
    vector normal = {0, 0, 0};
    vector p = {0, 0, 0};
    for(int i = 0; i < 4; i ++){
        double newt = intersect(ray, &arr[i], start);
        if(newt < t && newt > 0){
            t = newt;
            min_color = arr[i].c;
        }
    }
    if(min_color != BACK){
        sphere s = arr[min_color];

        p.x = t * ray->x + start.x;
        p.y = t * ray->y + start.y;
        p.z = t * ray->z + start.z;


        normal.x = (p.x - s.v.x) / s.r;
        normal.y = (p.y - s.v.y) / s.r;
        normal.z = (p.z - s.v.z) / s.r;

        p.x += normal.x * 0.001;
        p.y += normal.y * 0.001;
        p.z += normal.z * 0.001;

        if(min_color == FLOOR){
            pos = abs(
                    ((int)( p.x * 10 + 10000)) % 2
                    - ((int)( p.z * 10 + 10000)) % 2
            );
        }

        vector* tolight = malloc(sizeof(vector));
        diff(tolight, g, p);
        normalize(tolight);

        for(int i = 0; i < 4; i ++){
            double sun = intersect(tolight, &arr[i], p);
            if(min_color != i && sun > 0){
                goto SHADOW;
            }
        }


        double brightness = dot(&normal, tolight);
        if(brightness < 0){
            brightness = 0;
        }
        brightness = brightness * (1.0 - AMBIENT);
        brightness = brightness + AMBIENT;
        out = malloc(3* sizeof(int));
        out[0] = c_vals[min_color + pos][0] * brightness;
        out[1] = c_vals[min_color + pos][1] * brightness;
        out[2] = c_vals[min_color + pos][2] * brightness;

        goto EXIT;
    }

    return c_vals[min_color];

    SHADOW:
    out = malloc(3* sizeof(int));
    out[0] = c_vals[min_color + pos][0] * AMBIENT;
    out[1] = c_vals[min_color + pos][1] * AMBIENT;
    out[2] = c_vals[min_color + pos][2] * AMBIENT;

    EXIT:

    if(mix[min_color] > NONE){

        vector* new = malloc(sizeof(vector));
        new->x = ray->x - 2 * (normal.x * dot(ray, &normal));
        new->y = ray->y - 2 * (normal.y * dot(ray, &normal));
        new->z = ray->z - 2 * (normal.z * dot(ray, &normal));

        normalize(new);

        int* reflect = castray(new, p);

        if(reflect == c_vals[BACK]){
            goto GETOUT;
        }


        out[0] *= mix[min_color];
        out[1] *= mix[min_color];
        out[2] *= mix[min_color];

        out[0] += (1-mix[min_color]) * reflect[0];
        out[1] += (1-mix[min_color]) * reflect[1];
        out[2] += (1-mix[min_color]) * reflect[2];
    }

    GETOUT:
    return out;
}

double intersect(vector* ray, sphere* dest, vector start){
    double dx, dy, dz;
    dx = start.x - dest->v.x;
    dy = start.y - dest->v.y;
    dz = start.z - dest->v.z;

    double a = (ray->x)*(ray->x) + (ray->y)*(ray->y) + (ray->z)*(ray->z);
    double b = 2 * (ray->x*dx + ray->y*dy + ray->z*dz);
    double c = dx*dx + dy*dy + dz*dz - (dest->r * dest->r);

    double desc = (b*b) - (4 * a * c);
    if(desc <= 0){return NONE;}

    double t1 = (     sqrt(desc) - b) / (2 * a);
    double t2 = (-1 * sqrt(desc) - b) / (2 * a);

    if(t1 < 0){return t2;}
    if(t2 < 0){return t1;}
    if(t1 < t2){return t1;}
    return t2;
}
//
// end of file
//
