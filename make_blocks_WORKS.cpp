#include "PerlinNoise.hpp"
#include <SFML/Graphics.hpp>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <conio.h>
#include <vector>

using namespace std;

sf::Vector3f rotatePoint3D(sf::Vector3f& point, sf::Vector3f& center, float angleX, float angleY, float angleZ) {
    // point should be formatted: x, y, z
    // angles should be in radians
    sf::Vector3f p = point - center;

    
    // Calculate sin and cos values of the angles
    float sinX = std::sin(angleX);
    float cosX = std::cos(angleX);
    float sinY = std::sin(angleY);
    float cosY = std::cos(angleY);
    float sinZ = std::sin(angleZ);
    float cosZ = std::cos(angleZ);
 
    // Rotation matrix
    vector<vector<float>> R = {
        { cosY*cosZ, cosX*sinZ + sinX*sinY*cosZ, sinX*sinZ - cosX*sinY*cosZ },
        { -cosY*sinZ, cosX*cosZ - sinX*sinY*sinZ, sinX*cosZ + cosX*sinY*sinZ },
        { sinY, -sinX*cosY, cosX*cosY }
    };

    // Perform the rotation
    vector<float> result(3);
    for (int i = 0; i < 3; i++) {
        result[i] = p.x*R[i][0] + p.y*R[i][1] + p.z*R[i][2];
    }

    // Return the rotated point
    return sf::Vector3f(result[0],result[1],result[2])+center;
}

// Converts degrees to radians
float toRadians(float degrees) {
    return degrees * M_PI / 180.0;
}

// Helper cross function for getViewMatrix
sf::Vector3f cross(const sf::Vector3f& v1, const sf::Vector3f& v2) {
    return sf::Vector3f(
        v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x
    );
}

// Helper dot function for getViewMatrix
float dot(const sf::Vector3f& v1, const sf::Vector3f& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// Calculates the view matrix
void getViewMatrix(const sf::Vector3f& cameraPos, const sf::Vector3f& cameraDir, float viewMatrix[16]) {
    sf::Vector3f up = sf::Vector3f(0.0f, 1.0f, 0.0f);
    sf::Vector3f zAxis = sf::Vector3f(cameraDir.x, cameraDir.y, cameraDir.z);
    float mag = sqrt(zAxis.x*zAxis.x + zAxis.y*zAxis.y + zAxis.z*zAxis.z);
    if (mag > 0.0f) {
        zAxis = sf::Vector3f(zAxis.x / mag, zAxis.y / mag, zAxis.z / mag);
        sf::Vector3f xAxis = cross(up, zAxis);
        sf::Vector3f yAxis = cross(zAxis, xAxis);
        viewMatrix[0] = xAxis.x; viewMatrix[4] = xAxis.y; viewMatrix[8]  = xAxis.z; viewMatrix[12] = -dot(xAxis, cameraPos);
        viewMatrix[1] = yAxis.x; viewMatrix[5] = yAxis.y; viewMatrix[9]  = yAxis.z; viewMatrix[13] = -dot(yAxis, cameraPos);
        viewMatrix[2] = zAxis.x; viewMatrix[6] = zAxis.y; viewMatrix[10] = zAxis.z; viewMatrix[14] = -dot(zAxis, cameraPos);
        viewMatrix[3] = 0.0;     viewMatrix[7] = 0.0;     viewMatrix[11] = 0.0;     viewMatrix[15] = 1.0;
    }
    else {
        // zAxis has zero magnitude, set xAxis and yAxis to arbitrary non-parallel vectors
        sf::Vector3f xAxis = sf::Vector3f(1.0f, 0.0f, 0.0f);
        if (zAxis.x == 0.0f && zAxis.z == 0.0f) {
            xAxis = sf::Vector3f(0.0f, 0.0f, 1.0f);
        }
        sf::Vector3f yAxis = cross(zAxis, xAxis);
        viewMatrix[0] = xAxis.x; viewMatrix[4] = xAxis.y; viewMatrix[8]  = xAxis.z; viewMatrix[12] = -dot(xAxis, cameraPos);
        viewMatrix[1] = yAxis.x; viewMatrix[5] = yAxis.y; viewMatrix[9]  = yAxis.z; viewMatrix[13] = -dot(yAxis, cameraPos);
        viewMatrix[2] = zAxis.x; viewMatrix[6] = zAxis.y; viewMatrix[10] = zAxis.z; viewMatrix[14] = -dot(zAxis, cameraPos);
        viewMatrix[3] = 0.0;     viewMatrix[7] = 0.0;     viewMatrix[11] = 0.0;     viewMatrix[15] = 1.0;
    }
}

// Calculates the projection matrix
void getProjectionMatrix(float fov, float aspectRatio, float nearPlane, float farPlane, float projectionMatrix[16]) {
    float tanHalfFov = tan(toRadians(fov/2.0));
    projectionMatrix[0] = 1.0 / (aspectRatio * tanHalfFov); projectionMatrix[4] = 0.0;              projectionMatrix[8]  = 0.0;                                       projectionMatrix[12] = 0.0;
    projectionMatrix[1] = 0.0;                                projectionMatrix[5] = 1.0 / tanHalfFov; projectionMatrix[9]  = 0.0;                                       projectionMatrix[13] = 0.0;
    projectionMatrix[2] = 0.0;                                projectionMatrix[6] = 0.0;              projectionMatrix[10] = -(farPlane + nearPlane) / (farPlane - nearPlane); projectionMatrix[14] = -2.0 * farPlane * nearPlane / (farPlane - nearPlane);
    projectionMatrix[3] = 0.0;                                projectionMatrix[7] = 0.0;              projectionMatrix[11] = -1.0;                                      projectionMatrix[15] = 0.0;
}

// Transforms a 3D point by a matrix
sf::Vector3f transformPoint(const sf::Vector3f& point, const float matrix[16]) {
    float x = point.x * matrix[0] + point.y * matrix[4] + point.z * matrix[8]  + matrix[12];
    float y = point.x * matrix[1] + point.y * matrix[5] + point.z * matrix[9]  + matrix[13];
    float z = point.x * matrix[2] + point.y * matrix[6] + point.z * matrix[10] + matrix[14];
    return sf::Vector3f(x, y, z);
}

// Projects a 3D point to 2D screen coordinates
sf::Vector2f projectPoint(const sf::Vector3f& point, const float viewMatrix[16], const float projectionMatrix[16], int screenWidth, int screenHeight) {
    sf::Vector3f viewPos = transformPoint(point, viewMatrix);
    sf::Vector3f projPos = transformPoint(viewPos, projectionMatrix);
    float x = screenWidth * (projPos.x + 1.0f) / 2.0f;
    float y = screenHeight * (1.0f - projPos.y) / 2.0f;
    return sf::Vector2f(x, y);
}

// Renders a vector of 3D points as seen from a camera
vector<sf::Vector2f> renderPoints(const vector<sf::Vector3f>& points, const sf::Vector3f& cameraPos, const sf::Vector3f& cameraDir, float fov, int screenWidth, int screenHeight) {
    float viewMatrix[16], projectionMatrix[16];
    getViewMatrix(cameraPos, cameraDir, viewMatrix);
    getProjectionMatrix(fov, (float)screenWidth / (float)screenHeight, 0.00001f, 10.0f, projectionMatrix);
    vector<sf::Vector2f> screenPoints;
    for (int i = 0; i < points.size(); i++) {
        screenPoints.push_back(projectPoint(points[i], viewMatrix, projectionMatrix, screenWidth, screenHeight));
    }
    return screenPoints;
}


vector<sf::Vector3f> getNewPoints(const float& x, const float& z, const float& x_steps, const float& z_steps, const siv::PerlinNoise& perlin, const float& squish_num, const int& substeps, const int& noise_amp) {
    vector<sf::Vector3f> new_points;
    
    if (x_steps >= 0 and z_steps >= 0)
    {
        for (int i = z; i < z+z_steps; ++i)
        {
            for (int j = x; j < x+x_steps; ++j)
            {
                //std::cout << (x*squish_num) << " " << (z*squish_num) << endl;
                const float noise = perlin.normalizedOctave2D((j*squish_num/substeps),(i*squish_num/substeps),4);
                //const float noise = perlin.octave2D_01((x * 0.01), (y * 0.01), 4);

                new_points.push_back(sf::Vector3f(x/static_cast<float>(substeps), noise*noise_amp, z/static_cast<float>(substeps)));
                //std::cout << noise << '\t' << endl;
            }
            //std::cout << '\n' << endl;
        }
    }
    else if (x_steps >= 0 and z_steps <= 0)
    {
        for (int i = z; i > z+z_steps; --i)
        {
            for (int j = x; j < x+x_steps; ++j)
            {
                //std::cout << (x*squish_num) << " " << (z*squish_num) << endl;
                const float noise = perlin.normalizedOctave2D((j*squish_num/substeps),(i*squish_num/substeps),4);
                //const float noise = perlin.octave2D_01((x * 0.01), (y * 0.01), 4);

                new_points.push_back(sf::Vector3f(x/static_cast<float>(substeps), noise*noise_amp, z/static_cast<float>(substeps)));
                //std::cout << noise << '\t' << endl;
            }
            //std::cout << '\n' << endl;
        }
    }
    else if (x_steps <= 0 and z_steps >= 0)
    {
        for (int i = z; i < z+z_steps; ++i)
        {
            for (int j = x; j > x+x_steps; --j)
            {
                //std::cout << (x*squish_num) << " " << (z*squish_num) << endl;
                const float noise = perlin.normalizedOctave2D((j*squish_num/substeps),(i*squish_num/substeps),4);
                //const float noise = perlin.octave2D_01((x * 0.01), (y * 0.01), 4);

                new_points.push_back(sf::Vector3f(x/static_cast<float>(substeps), noise*noise_amp, z/static_cast<float>(substeps)));
                //std::cout << noise << '\t' << endl;
            }
            //std::cout << '\n' << endl;
        }
    }
    else if (x_steps <= 0 and z_steps <= 0)
    {
        for (int i = z; i > z+z_steps; --i)
        {
            for (int j = x; j > x+x_steps; --j)
            {
                //std::cout << (x*squish_num) << " " << (z*squish_num) << endl;
                const float noise = perlin.normalizedOctave2D((j*squish_num/substeps),(i*squish_num/substeps),4);
                //const float noise = perlin.octave2D_01((x * 0.01), (y * 0.01), 4);

                new_points.push_back(sf::Vector3f(x/static_cast<float>(substeps), noise*noise_amp, z/static_cast<float>(substeps)));
                //std::cout << noise << '\t' << endl;
            }
            //std::cout << '\n' << endl;
        }
    }

    return new_points;
}

bool isinCircle(const sf::Vector3f& point, const sf::Vector3f& center, const float radius) {
    return (pow(pow((point.x-center.x),2)+pow(point.y-center.y,2)+pow(point.z-center.z,2),0.333333333333333333333333333) <= radius) ? true : false ;
}

vector<sf::Vector3f> getPoints(sf::Vector3f& camera_position, const int size = 2, const int substeps = 32, const siv::PerlinNoise::seed_type seed = 123456u, const float squish_num = 1, const float noise_amp = 1, const float scale = 1, const float view_radius = 1) {
    
    // const siv::PerlinNoise::seed_type seed = 123456u;

	const siv::PerlinNoise perlin{ seed };
    
    vector<sf::Vector3f> points;
    vector<sf::Vector3f> output_points;
    sf::Vector3f point;
    
    for (int z = -size/2*substeps+camera_position.z*substeps; z < size*substeps-size/2*substeps+camera_position.z*substeps; ++z)
    {
        for (int x = -size/2*substeps+camera_position.x*substeps; x < size*substeps-size/2*substeps+camera_position.x*substeps; ++x)
        {
            const float noise = perlin.normalizedOctave2D((x*squish_num/substeps),(z*squish_num/substeps),4);

            points.push_back(sf::Vector3f(x/static_cast<float>(substeps), noise*noise_amp, z/static_cast<float>(substeps)));
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


sf::Color getColor(const float& Y_POS) {
    return sf::Color(
        255*((Y_POS+1)/2), // red
        255*(1-(Y_POS+1)/2), // green
        0, // blue
        255 // alpha
        );
}


void addPoint(const int& i, sf::VertexArray& vertices, const vector<sf::Vector3f>& points, const sf::Vector3f& point, const sf::Vector3f& camera_position, const float& view_radius, const float& scale, const int& offset_x, const int& offset_y) {
    if (isinCircle(point, camera_position, view_radius))
    {
        vertices.append(sf::Vertex(sf::Vector2f(point.x*scale+offset_x-camera_position.x*scale, point.y*scale+offset_y),getColor(points[i].y)));
        //vertices[-1].color = sf::Color(255, 0, 0, 255); //*points[i].y
    }
}

sf::Vector2f calculate_screen_coordinate(float h_fov_deg, float v_fov_deg, const sf::Vector3f& camera_pos, const sf::Vector3f& camera_dir, const sf::Vector3f& point_pos, const sf::Vector2f& screen_size) {
    // Convert fov from degrees to radians
    float h_fov = std::tan(h_fov_deg / 2 * M_PI / 180) * 2;
    float v_fov = std::tan(v_fov_deg / 2 * M_PI / 180) * 2;

    // Calculate the vector from the camera to the point
    sf::Vector3f camera_to_point_vector = point_pos - camera_pos;

    // Calculate the magnitude of the camera to point vector
    float magnitude = std::sqrt(camera_to_point_vector.x * camera_to_point_vector.x + camera_to_point_vector.y * camera_to_point_vector.y + camera_to_point_vector.z * camera_to_point_vector.z);

    // Calculate the angle between the camera to point vector and the camera's forward direction
    sf::Vector3f forward_vector = -camera_dir;
    float cos_theta = (camera_to_point_vector.x * forward_vector.x + camera_to_point_vector.y * forward_vector.y + camera_to_point_vector.z * forward_vector.z) / magnitude;
    float theta = std::acos(cos_theta);

    // Calculate the screen coordinate of the point
    float tan_h_fov = std::tan(h_fov / 2);
    float tan_v_fov = std::tan(v_fov / 2);
    float x = (camera_to_point_vector.x / magnitude) / (tan_h_fov * screen_size.x / 2);
    float y = (camera_to_point_vector.y / magnitude) / (tan_v_fov * screen_size.y / 2);

    // Return the screen coordinate as an sf::Vector2f
    return sf::Vector2f((x + 1) * screen_size.x / 2, (1 - y) * screen_size.y / 2);
}

int main()
{
	const siv::PerlinNoise::seed_type seed = 123456u;

	const siv::PerlinNoise perlin{ seed };

    std::cout << "Enter size of plane." << endl;
    int size;
    std::cin >> size;

    std::cout << "Enter number of substeps. Enter 1 if you want no substeps." << endl;
    int substeps;
    std::cin >> substeps;

    std::cout << "Enter noise amplifying constant." << endl;
    int noise_amp;
    std::cin >> noise_amp;


    const int x_iters = size;
    const int z_iters = size;

    const int window_x = 800;
    const int window_y = 600;

    //camera stuff

    sf::Vector3f camera_position = sf::Vector3f(0.0f,0.0f,0.0f);
    const float view_radius = 1;
    cout << "Radius: " << view_radius << endl;

    const float scaling_constant = 1.1;
    float scale = ((window_y>window_x) ? window_x: window_y)/size;
    float old_scale;

    int offset_x = (size%2==0 ? window_x/2 : window_x/2-scale/2);//-scale*size/2;
    int offset_y = window_y/2;//-scale*size/2;

    const float squish_num = 1;

    const float pi = 3.141592653589793238;

    const float rotation_amount = 3 * pi / 180;
    const float move_amount = 0.01;

    sf::Vector3f points_rotation = sf::Vector3f(0.0f,0.0f,0.0f);

    sf::Vector3f point;

    // Create an empty vector of 3D points
    vector<sf::Vector3f> points;

    int max_z = z_iters*substeps-z_iters/2*substeps-1;
    int min_z = -z_iters/2*substeps;
    int max_x = x_iters*substeps-x_iters/2*substeps-1;
    int min_x = -x_iters/2*substeps;


	for (int z = -z_iters/2*substeps; z < z_iters*substeps-z_iters/2*substeps; ++z)
	{
		for (int x = -x_iters/2*substeps; x < x_iters*substeps-x_iters/2*substeps; ++x)
		{
            //std::cout << (x*squish_num) << " " << (z*squish_num) << endl;
            const float noise = perlin.normalizedOctave2D((x*squish_num/substeps),(z*squish_num/substeps),4);
            //const float noise = perlin.octave2D_01((x * 0.01), (y * 0.01), 4);

            points.push_back(sf::Vector3f(x/static_cast<float>(substeps), noise*noise_amp, z/static_cast<float>(substeps)));
            //std::cout << noise << '\t' << endl;
		}
        //std::cout << '\n' << endl;
	}

    //scale = (x_iters>z_iters) ? scale/x_iters : scale/z_iters;

    std::cout << "Finished generating " << points.size() << " points.\nPress any key to continue..." << endl;
    std::cin.get();

    // Create the SFML window
    sf::RenderWindow window(sf::VideoMode(window_x, window_y), "SFML window");

    // Create a vertex array to hold the points
    //vector<sf::Vector2f> screenPoints = renderPoints(points,camera_position,camera_direction,FOV,window_x,window_y);

    sf::VertexArray vertices(sf::Points);
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        if (isinCircle(points[i], camera_position, view_radius))
        {
            vertices.append(sf::Vertex(sf::Vector2f(point.x*scale+offset_x-camera_position.x*scale, point.y*scale+offset_y),getColor(points[i].y)));
        }
    }
    cout << "Vertex amount: " << vertices.getVertexCount() << endl;
    cout << "Scale: " << scale << "\r";

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
                if (event.key.code == sf::Keyboard::Num1)
                {
                    camera_position = sf::Vector3f(camera_position.x+move_amount,camera_position.y,camera_position.z);
                    points.clear();
                    for (int z = -z_iters/2*substeps+camera_position.z*substeps; z < z_iters*substeps-z_iters/2*substeps+camera_position.z*substeps; ++z)
                    {
                        for (int x = -x_iters/2*substeps+camera_position.x*substeps; x < x_iters*substeps-x_iters/2*substeps+camera_position.x*substeps; ++x)
                        {
                            const float noise = perlin.normalizedOctave2D((x*squish_num/substeps),(z*squish_num/substeps),4);

                            points.push_back(sf::Vector3f(x/static_cast<float>(substeps), noise*noise_amp, z/static_cast<float>(substeps)));
                        }
                    }
                    vertices.clear();
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        point = rotatePoint3D(points[i],camera_position,points_rotation.x,points_rotation.y,points_rotation.z);
                        addPoint(i, vertices, points, point, camera_position, view_radius, scale, offset_x, offset_y);
                    }
                }
                else if (event.key.code == sf::Keyboard::Num2)
                {
                    camera_position = sf::Vector3f(camera_position.x-move_amount,camera_position.y,camera_position.z);
                    points.clear();
                    for (int z = -z_iters/2*substeps+camera_position.z*substeps; z < z_iters*substeps-z_iters/2*substeps+camera_position.z*substeps; ++z)
                    {
                        for (int x = -x_iters/2*substeps+camera_position.x*substeps; x < x_iters*substeps-x_iters/2*substeps+camera_position.x*substeps; ++x)
                        {
                            const float noise = perlin.normalizedOctave2D((x*squish_num/substeps),(z*squish_num/substeps),4);

                            points.push_back(sf::Vector3f(x/static_cast<float>(substeps), noise*noise_amp, z/static_cast<float>(substeps)));
                        }
                    }
                    vertices.clear();
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        point = rotatePoint3D(points[i],camera_position,points_rotation.x,points_rotation.y,points_rotation.z);
                        addPoint(i, vertices, points, point, camera_position, view_radius, scale, offset_x, offset_y);
                    }
                }
                if (event.key.code == sf::Keyboard::Equal)
                {
                    old_scale = scale;
                    scale *= scaling_constant;
                    for (std::size_t i = 0; i < vertices.getVertexCount(); ++i)
                    {
                        vertices[i] = sf::Vertex(sf::Vector2f(((vertices[i].position.x-offset_x)/old_scale)*scale+offset_x, ((vertices[i].position.y-offset_y)/old_scale)*scale+offset_y), getColor(points[i].y));
                    }
                    cout << "Scale: " << scale << "\r";
                }
                else if (event.key.code == 56)
                {
                    old_scale = scale;
                    scale /= scaling_constant;
                    for (std::size_t i = 0; i < vertices.getVertexCount(); ++i)
                    {
                        vertices[i] = sf::Vertex(sf::Vector2f(((vertices[i].position.x-offset_x)/old_scale)*scale+offset_x, ((vertices[i].position.y-offset_y)/old_scale)*scale+offset_y), getColor(points[i].y));
                    }
                    cout << "Scale: " << scale << "\r";
                }
                else if (event.key.code == sf::Keyboard::A)
                {
                    vertices.clear();
                    points_rotation = sf::Vector3f(points_rotation.x+rotation_amount, points_rotation.y, points_rotation.z);
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        point = rotatePoint3D(points[i],camera_position,points_rotation.x,points_rotation.y,points_rotation.z);
                        addPoint(i, vertices, points, point, camera_position, view_radius, scale, offset_x, offset_y);
                    }
                }
                else if (event.key.code == sf::Keyboard::Q)
                {
                    vertices.clear();
                    points_rotation = sf::Vector3f(points_rotation.x-rotation_amount, points_rotation.y, points_rotation.z);
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        point = rotatePoint3D(points[i],camera_position,points_rotation.x,points_rotation.y,points_rotation.z);
                        addPoint(i, vertices, points, point, camera_position, view_radius, scale, offset_x, offset_y);
                    }
                }
                else if (event.key.code == sf::Keyboard::S)
                {
                    vertices.clear();
                    points_rotation = sf::Vector3f(points_rotation.x, points_rotation.y+rotation_amount, points_rotation.z);
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        point = rotatePoint3D(points[i],camera_position,points_rotation.x,points_rotation.y,points_rotation.z);
                        addPoint(i, vertices, points, point, camera_position, view_radius, scale, offset_x, offset_y);
                    }
                }
                else if (event.key.code == sf::Keyboard::W)
                {
                    vertices.clear();
                    points_rotation = sf::Vector3f(points_rotation.x, points_rotation.y-rotation_amount, points_rotation.z);
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        point = rotatePoint3D(points[i],camera_position,points_rotation.x,points_rotation.y,points_rotation.z);
                        if (isinCircle(point, camera_position, view_radius))
                        addPoint(i, vertices, points, point, camera_position, view_radius, scale, offset_x, offset_y);
                    }
                }
                else if (event.key.code == sf::Keyboard::D)
                {
                    vertices.clear();
                    points_rotation = sf::Vector3f(points_rotation.x, points_rotation.y, points_rotation.z+rotation_amount);
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        point = rotatePoint3D(points[i],camera_position,points_rotation.x,points_rotation.y,points_rotation.z);
                        addPoint(i, vertices, points, point, camera_position, view_radius, scale, offset_x, offset_y);
                    }
                }
                else if (event.key.code == sf::Keyboard::E)
                {
                    vertices.clear();
                    points_rotation = sf::Vector3f(points_rotation.x, points_rotation.y, points_rotation.z-rotation_amount);
                    for (std::size_t i = 0; i < points.size(); ++i)
                    {
                        point = rotatePoint3D(points[i],camera_position,points_rotation.x,points_rotation.y,points_rotation.z);
                        addPoint(i, vertices, points, point, camera_position, view_radius, scale, offset_x, offset_y);
                    }
                }
            }
        }

        // Clear the window
        window.clear();

        // Draw the vertices
        window.draw(vertices);

        // Display the window
        window.display();
    }

    return 0;
}