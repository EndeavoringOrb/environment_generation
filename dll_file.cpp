#include "PerlinNoise.hpp"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <conio.h>
#include <vector>

using namespace std;

extern "C"{
    __declspec( dllexport ) bool isinCircle(const vector<float>& point, const vector<float>& center, const float radius) {
        return (pow(pow((point[0]-center[0]),2)+pow(point[1]-center[1],2)+pow(point[2]-center[2],2),0.333333333333333333333333333) <= radius) ? true : false ;
    }
}

extern "C"{
    __declspec( dllexport ) vector<vector<float>> getPoints(vector<float>& camera_position, const int size = 2, const int substeps = 32, const siv::PerlinNoise::seed_type seed = 123456u, const float squish_num = 1, const float noise_amp = 1, const float scale = 1, const float view_radius = 1) {
        
        // const siv::PerlinNoise::seed_type seed = 123456u;

        const siv::PerlinNoise perlin{ seed };
        
        vector<vector<float>> points;
        vector<vector<float>> output_points;
        vector<float> point;
        
        for (int z = -size/2*substeps+camera_position[0]*substeps; z < size*substeps-size/2*substeps+camera_position[2]*substeps; ++z)
        {
            for (int x = -size/2*substeps+camera_position[0]*substeps; x < size*substeps-size/2*substeps+camera_position[0]*substeps; ++x)
            {
                const float noise = perlin.normalizedOctave2D((x*squish_num/substeps),(z*squish_num/substeps),4);

                points.push_back({x/static_cast<float>(substeps), noise*noise_amp, z/static_cast<float>(substeps)});
            }
        }
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            if (isinCircle(points[i], camera_position, view_radius))
            {
                output_points.push_back(points[i]);
            }
        }

        return output_points;
    }
}