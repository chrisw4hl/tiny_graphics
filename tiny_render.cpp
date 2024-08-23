#include "./tinyrenderer/tgaimage.cpp"
#include "./tinyrenderer/model.cpp"
#include "./tinyrenderer/geometry.cpp"
#include "tinyrenderer/tgaimage.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <ostream>
#include <utility>
#include <vector>

int width = 1000;
int height = 1000;
const TGAColor white = TGAColor({255, 255, 255, 255});
const TGAColor red   = TGAColor({0, 0,   255,   255});
const TGAColor blue   = TGAColor({255, 0,   0,   255});
const TGAColor green   = TGAColor({0, 255,   0,   255});

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color){
        
        std::vector<int> points_x;
        bool steep = false;
        if (std::abs(x0-x1)<std::abs(y0-y1)) { 
                std::swap(x0, y0); 
                std::swap(x1, y1); 
                steep = true; 
        }
        if (x1<x0){
                std::swap(x0,x1);
                std::swap(y0,y1);
        }

        int x = x0;
        int y = y0;
        bool y_change = false;
        float error = 0;
        float derror = std::abs((y1-y0)/float(x1-x0));
        while(x<x1){
                if(steep){
                        image.set(y,x,color);
                }else{
                        image.set(x,y, color);
                }
                x++;
                error += derror;
                if(error > .5){
                        y+=(y1>y0?1:-1);
                        error-=1.;
                }
                
        }
}

vec3 barycentric(double bary_coor[2][2], double r_diff[2]){
        double lambda_1 = bary_coor[0][0]*r_diff[0] + bary_coor[0][1]*r_diff[1];
        double lambda_2 = bary_coor[1][0]*r_diff[0] + bary_coor[1][1]*r_diff[1];
        double lambda_3 = 1. - lambda_1 - lambda_2;

        vec3 coord = vec3({lambda_1, lambda_2, lambda_3});
        return coord;
}

bool is_inside(vec3 coord){
        if ((coord.x < 0) || (coord.y < 0) || (coord.z < 0)){
                return false;
        }
        return true;
}

float viewport[4][4] = {{float(width/2.),0.,0.,0.},{0.,float(height/2.),0.,0.},{0.,0.,10.,0.},{0.,0.,0.,0.}};

std::vector<vec3> transform(std::vector<vec3> points, int y_angle, std::vector<std::vector<float>> viewport){
        std::vector<vec3>* result = new std::vector<vec3>;




        return *result;
}


void triangle(vec3 v0, vec3 v1, vec3 v2, TGAImage &image, TGAColor color){

        std::vector<vec3> vec{v0, v1, v2};
        std::sort(vec.begin(),vec.end(), [](vec3 a, vec3 b){
                                        return a.y < b.y;
                                });

        line(vec[0].x, vec[0].y, vec[1].x, vec[1].y, image, green);
        line(vec[0].x, vec[0].y, vec[2].x, vec[2].y, image, red);
        line(vec[1].x, vec[1].y, vec[2].x, vec[2].y, image, red);
        float low_slope = (vec[1].x-vec[0].x)/(vec[1].y-vec[0].y);
        float high_slope = (vec[2].x-vec[0].x)/(vec[2].y-vec[0].y);
        int low = vec[0].y;
        int mid = vec[1].y;
        int high = vec[2].y;

        for(int i = 1; i<= mid-low; i++){
                
                int x_low = low_slope*i + vec[0].x;
                int x_high = high_slope*i + vec[0].x;
                line(x_low+1, i+low, x_high+1, i+low, image, white);

        }
        

        float top_slope = (vec[2].x-vec[1].x)/(vec[2].y-vec[1].y);
        for(int i = mid; i<= high; i++){
                
                int x_low = high_slope*(i-mid) + vec[0].x + high_slope*(mid-low);
                int x_high = top_slope*(i-mid) + vec[1].x;
                line(x_low+1, i, x_high+1, i, image, white);

        }
        
        std::cout<<"done"<<std::endl;
}

void triangle2(vec3 v0, vec3 v1, vec3 v2, TGAImage &image, TGAColor color, float* zbuffer){ 
        std::vector<vec3> vec_vert_sort{v0, v1, v2};
        std::vector<vec3> vec_horiz_sort{v0, v1, v2};
        std::sort(vec_vert_sort.begin(),vec_vert_sort.end(), [](vec3 a, vec3 b){
                                                                   return a.y < b.y;
                                                                });
        std::sort(vec_horiz_sort.begin(),vec_horiz_sort.end(), [](vec3 a, vec3 b){
                                                                return a.x < b.x;
                                                                });

        double bary_convert[2][2] = { {v0.x - v2.x, v1.x - v2.x}, {v0.y - v2.y, v1.y - v2.y}};
        double det = bary_convert[0][0]*bary_convert[1][1] - bary_convert[1][0]*bary_convert[0][1];

        if (det == 0){
                std::cout<<"error, singular matrix"<<std::endl;
        }
        double bary_convert_inv[2][2] = {{1/det*bary_convert[1][1], -1/det*bary_convert[0][1]}, 
                {-1/det*bary_convert[1][0], 1/det*bary_convert[0][0]}};

        for(int i = vec_horiz_sort[0].x; i< vec_horiz_sort[2].x; i++){
                for(int j = vec_vert_sort[0].y; j< vec_vert_sort[2].y; j++){
                        double r_diff[2] = {double(i) - v2.x, double(j)- v2.y};
                        vec3 bary_coord = barycentric(bary_convert_inv, r_diff);
                        bool test = is_inside(bary_coord);
                        if (test){
                                image.set(i,j, color);
                        }
                }
        }

       
}
TGAColor interp_texture(vec3 bary_coor, std::vector<vec2> texture_coor, Model &model){
        int x = int(bary_coor.x*texture_coor[0].x*width+bary_coor.y*texture_coor[1].x*width+bary_coor.z*texture_coor[2].x*width);
        int y = int(bary_coor.x*texture_coor[0].y*height+bary_coor.y*texture_coor[1].y*height+bary_coor.z*texture_coor[2].y*height);
        //std::cout<<texture_coor[0].x<<","<<texture_coor[1].x<<std::endl;
        return model.diffuse().get(x,y);
}

void triangle3(vec3 v0, vec3 v1, vec3 v2, TGAImage &image, TGAColor color, float *zbuffer, std::vector<vec2> texture_coor, Model &model){ 
        double bary_convert[2][2] = { {v0.x - v2.x, v1.x - v2.x}, {v0.y - v2.y, v1.y - v2.y}};
        
        double det = bary_convert[0][0]*bary_convert[1][1] - bary_convert[1][0]*bary_convert[0][1];

        if (det == 0){
                std::cout<<"error, singular matrix"<<std::endl;
        }
        double bary_convert_inv[2][2] = {{1/det*bary_convert[1][1], -1/det*bary_convert[0][1]}, 
                {-1/det*bary_convert[1][0], 1/det*bary_convert[0][0]}};

        for(int i = std::min({v0.x,v1.x,v2.x}); i< std::max({v0.x,v1.x,v2.x}); i++){
                for(int j = std::min({v0.y,v1.y,v2.y}); j< std::max({v0.y,v1.y,v2.y}); j++){
                        double r_diff[2] = {double(i) - v2.x, double(j)- v2.y};
                        vec3 bary_coor = barycentric(bary_convert_inv, r_diff);
                        bool test = is_inside(bary_coor);
                        if(!test){
                                continue;
                        }
                        TGAColor texture_color = interp_texture(bary_coor, texture_coor, model);
                        float z= v0.z*bary_coor[0]+v1.z*bary_coor[1]+v2.z*bary_coor[2];
                        if ((zbuffer[i+width*j]<z)){
                                zbuffer[i+ width*j] = z;
                                image.set(i,j, texture_color);
                        }

                }
        }

       
}

vec3 normalize(vec3 v0, vec3 v1){
        vec3 ret = vec3({(v1.z*v0.y-v0.z*v1.y), (v0.z*v1.x-v0.x*v1.z), (v1.y*v0.x-v1.x*v0.y)});
        float len = std::sqrt(ret.x*ret.x+ret.y*ret.y+ret.z*ret.z);
        return vec3({ret.x/len,ret.y/len,ret.z/len});
}

int main(int argc, char** argv){
        std::cout<<"here"<<std::endl;
        int camera_distance = width;
        TGAImage image(width, height, TGAImage::RGB);
        float* zbuffer = new float[width*height];
        Model face_model("./tinyrenderer/obj/african_head/african_head.obj");
        Model* model = &face_model;
        for (int i=0; i<model->nfaces(); i++) { 
                vec3 v0 = model->vert(i,0); 
                vec3 v1 = model->vert(i,1); 
                vec3 v2 = model->vert(i,2);
                vec3 v0_screen = vec3({(v0.x+1.)*width/2, (v0.y+1.)*height/2, v0.z});  
                vec3 v1_screen = vec3({(v1.x+1.)*width/2, (v1.y+1.)*height/2, v1.z}); 
                vec3 v2_screen = vec3({(v2.x+1.)*width/2, (v2.y+1.)*height/2, v2.z});
                float intensity = normalize(v2-v0, v1-v0)*vec3({0.,0.,-1.})*255;
                uint8_t res = uint8_t(intensity);
                TGAColor shaded_white = TGAColor({res,res,res,255});
                if(intensity>0){
                        vec2 test = model->uv(1,0);
                        std::vector<vec2> texture_coor = {model->uv(i,0), model->uv(i,1), model->uv(i,2)};
                        triangle3(v0_screen, v1_screen, v2_screen, image, shaded_white, zbuffer, texture_coor, face_model);                          
                }
        }
        image.write_tga_file("output_1.tga");
        return 0;
}
