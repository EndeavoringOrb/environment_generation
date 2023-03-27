#include "PerlinNoise.hpp"
#include <SFML/Graphics.hpp>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <conio.h>
#include <vector>

using namespace std;

sf::Vector3f rotatePoint3D(sf::Vector3f& p, double angleX, double angleY, double angleZ) {
    // point should be formatted: x, y, z
    // angles should be in radians
    
    // Calculate sin and cos values of the angles
    double sinX = std::sin(angleX);
    double cosX = std::cos(angleX);
    double sinY = std::sin(angleY);
    double cosY = std::cos(angleY);
    double sinZ = std::sin(angleZ);
    double cosZ = std::cos(angleZ);

    // Rotation matrix
    vector<vector<double>> R = {
        { cosY*cosZ, cosX*sinZ + sinX*sinY*cosZ, sinX*sinZ - cosX*sinY*cosZ },
        { -cosY*sinZ, cosX*cosZ - sinX*sinY*sinZ, sinX*cosZ + cosX*sinY*sinZ },
        { sinY, -sinX*cosY, cosX*cosY }
    };

    // Perform the rotation
    vector<double> result(3);
    for (int i = 0; i < 3; i++) {
        result[i] = p.x*R[i][0] + p.y*R[i][1] + p.z*R[i][2];
    }

    // Return the rotated point
    return sf::Vector3f(result[0],result[1],result[2]);
}

bool isInFieldOfView(double hfov_deg, double vfov_deg, double cam_x, double cam_y, double cam_z, double point_x, double point_y, double point_z) {
    // Convert field of view angles to radians
    double hfov_rad = hfov_deg * M_PI / 180.0;
    double vfov_rad = vfov_deg * M_PI / 180.0;
    
    // Calculate vector from camera to point
    double dx = point_x - cam_x;
    double dy = point_y - cam_y;
    double dz = point_z - cam_z;
    
    // Calculate horizontal and vertical angles from camera to point
    double yaw = atan2(dy, dx);
    double pitch = asin(dz / sqrt(dx*dx + dy*dy + dz*dz));
    
    // Check if point is within field of view
    if (abs(yaw) <= hfov_rad / 2 && abs(pitch) <= vfov_rad / 2) {
        return true;
    } else {
        return false;
        cout << "false" << endl;
    }
}

int main()
{
	const siv::PerlinNoise::seed_type seed = 123456u;

	const siv::PerlinNoise perlin{ seed };

    cout << "Enter size of plane. Must be at least 2." << endl;
    int size;
    cin >> size;

    cout << "Enter number of substeps. Enter 1 if you want no substeps." << endl;
    int substeps;
    cin >> substeps;

    cout << "Enter noise amplifying constant." << endl;
    int noise_amp;
    cin >> noise_amp;

    //camera stuff
    int hFOV = 90; //115
    int vFOV = 90; //180
    int camera_x = 0;
    int camera_y = 0;
    int camera_z = 0;
    vector<sf::Vector2f> inview_points;


    const int x_iters = size;
    //const int y_iters = 5;
    const int z_iters = size;

    const int window_x = 800;
    const int window_y = 600;

    double scale = ((window_y>window_x) ? window_x: window_y)/size;
    cout << "scale: " << scale << endl;

    int offset_x = (size%2==0 ? window_x/2 : window_x/2-scale/2);//-scale*size/2;
    int offset_y = window_y/2;//-scale*size/2;

    const double squish_num = 1;

    const double pi = 3.141592653589793238;

    const double rotation_amount = 3 * pi / 180;

    // Create an empty vector of 3D points
    vector<sf::Vector3f> points;

    cout << "generating points" << endl;

	for (int z = -z_iters/2*substeps; z < z_iters*substeps-z_iters/2*substeps; ++z)
	{
		for (int x = -x_iters/2*substeps; x < x_iters*substeps-x_iters/2*substeps; ++x)
		{
            //cout << (x*squish_num) << " " << (z*squish_num) << endl;
            const double noise = perlin.normalizedOctave2D((x*squish_num/substeps),(z*squish_num/substeps),4);
            //const double noise = perlin.octave2D_01((x * 0.01), (y * 0.01), 4);

            points.push_back(sf::Vector3f(x/static_cast<double>(substeps), noise*noise_amp, z/static_cast<double>(substeps)));
            //cout << noise << '\t' << endl;
		}
        //cout << '\n' << endl;
	}

    //scale = (x_iters>z_iters) ? scale/x_iters : scale/z_iters;

    cout << "Finished generating " << points.size() << " points.\nPress any key to continue..." << endl;
    cin.get();

    // Create the SFML window
    sf::RenderWindow window(sf::VideoMode(800, 600), "SFML window");

    // Create a vertex array to hold the points
    sf::VertexArray vertices(sf::Points, points.size());
    for (std::size_t i = 0; i < points.size(); ++i)
    {

        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
    }

    // Start the main loop
    while (window.isOpen())
    {
        // Handle events
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed)
            {
                //cout << event.key.code << endl;
                if (event.key.code == sf::Keyboard::Equal)
                {
                    scale *= 1.5;
                    cout << "scale: " << scale << endl;
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == 56)
                {
                    scale /= 1.5;
                    cout << "scale: " << scale << endl;
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::Up)
                {
                    offset_y += 5;
                    cout << "offset y: " << offset_y << endl;
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::Down)
                {
                    offset_y -= 5;
                    cout << "offset y: " << offset_y << endl;
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::Left)
                {
                    offset_x += 5;
                    cout << "offset x: " << offset_x << endl;
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::Right)
                {
                    offset_x -= 5;
                    cout << "offset x: " << offset_x << endl;
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::A)
                {
                    // Up arrow key pressed
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        points[i] = rotatePoint3D(points[i],rotation_amount,0,0);
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::Q)
                {
                    // Down arrow key pressed
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        points[i] = rotatePoint3D(points[i],-rotation_amount,0,0);
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::S)
                {
                    // Left arrow key pressed
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        points[i] = rotatePoint3D(points[i],0,rotation_amount,0);
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::W)
                {
                    // Right arrow key pressed
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        points[i] = rotatePoint3D(points[i],0,-rotation_amount,0);
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::D)
                {
                    // Left arrow key pressed
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        points[i] = rotatePoint3D(points[i],0,0,rotation_amount);
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
                else if (event.key.code == sf::Keyboard::E)
                {
                    // Right arrow key pressed
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        points[i] = rotatePoint3D(points[i],0,0,-rotation_amount);
                        vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
                    }
                }
            }
        }

        // Clear the window
        window.clear();

        // Draw the vertices
        inview_points.clear();
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            if (isInFieldOfView(hFOV,vFOV,camera_x,camera_y,camera_z,points[i].x,points[i].y,points[i].z))
            {
                inview_points.push_back(sf::Vertex(sf::Vector2f(points[i].x, points[i].y)));
            }
            //vertices[i] = sf::Vertex(sf::Vector2f(points[i].x*scale+offset_x, points[i].y*scale+offset_y));
        }
        window.draw(inview_points);

        // Display the window
        window.display();
    }

    return 0;
}