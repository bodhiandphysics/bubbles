#ifndef MARCHINGCUBES_HPP
#define MARCHINGCUBES_HPP
#include <assert.h>
#include <atomic>
#include <vector>
#include <unordered_map>
#include <optional>
#include <glm/glm.hpp>
#include "table.h"
#include <omp.h>
#include <random>

bool is_present(int color, int* colors, int clength); 

struct coord {

	long x;
	long y;
	long z;

	coord(long x, long y, long z): x(x), y(y), z(z) {}
	bool operator==(const coord &other) const {

		return (x == other.x && y == other.y && z == other.z);
	}
};

template<>
struct std::hash<coord> {
	const long MINDEX = 100000; 
	size_t operator()(const coord& c) const {
		return std::hash<long>{}((c.x%MINDEX) + (c.y%MINDEX)*MINDEX + (c.z%MINDEX)*MINDEX*MINDEX);
	}
};

struct fcoord {

	float x;
	float y;
	float z;
	int color;

	fcoord(float x, float y, float , int color): x(x), y(y), z(z), color(color) {}
	bool operator==(const fcoord& other) const {

		return (x == other.x && y == other.y && z == other.z && color == other.color);
	}
};

template<>
struct std::hash<fcoord> {
	const long MINDEX = 100000; 
	size_t operator()(const fcoord& f) const {
		return std::hash<float>{}((fmod(f.x, MINDEX) + fmod(f.y, MINDEX)*MINDEX + fmod(f.z, MINDEX)*MINDEX*MINDEX)* f.color);
	}
};


struct Pot{

	int color;
	int new_color;
	bool on_border = false;
	bool already_draw = false;
	long x;
	long y; 
	long z;
	
};




template<typename T>
struct Grid {

	size_t x_max;
	size_t y_max;
	size_t z_max;
	T* data;


	Grid(size_t x, size_t y, size_t z): x_max(x), y_max(y), z_max(z) {

		data = (T*) malloc((x_max * y_max * z_max) * sizeof(T));
	}

	~Grid() {free(data);}

	T* operator()(long x, long y, long z) {

		if (x > x_max || y > y_max || z > z_max)
			return nullptr;

		long index = x + x_max*y + y_max*x_max*z;
		if (index >= x_max * y_max * z_max)
			return nullptr;
		else {
			return &data[index];
		}

	}
};

struct Simulation {

	long x_max, y_max, z_max;


	const int FREESPACE = 0;
	const int BOUNDARY = INT_MAX;

	float ALPHA = 1;
	float BETA = 1;
	Grid<Pot> grid;
	float scale = 1.f;
	std::vector<int> colors;
	std::vector<glm::vec3> color_map;
	int num_colors;

	std::vector<float> P0V0;
	std::vector<float> pressures;

	coord* border;
	std::atomic<int> border_length;
	std::atomic<int> bborder_length;
	coord* bborder;

	float surface_area = 0;

	// random stuff

	std::vector<std::mt19937> randoms;

	std::unordered_map<fcoord, long> point_map;
	std::vector<glm::vec3> points;
	int num_points = 0;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec3> cols;
	std::vector<GLint> indexes; 

	void flip(long x, long y, long z);
	void flip_cells();
	void rebuild();
	void confirm_flip();
	float get_de();
	void run_sim_iteration();

	u_char get_index_flipped(int color, long x, long y, long z, int fcolor, long fx, long fy, long fz, bool);
	u_char get_index(int color, long x, long y, long z, bool old);
	float get_surface_area(long x, long y, long z, bool old = true, 
							      bool flipped = false, int fcolor = 0, 
								  int fx = 0, int fy = 0, int fz = 0);
	bool on_border(long x, long y, long z);
	bool on_border(coord at);
	bool next_to_border(long x, long y, long z);
	bool next_to_border(coord at);

	glm::vec3 get_point(char edge, glm::vec3 cube_pos);
	float roll01();

	Simulation(size_t x_max, 
					   size_t y_max, 
					   size_t z_max,
					   float scale,  
					   std::vector<glm::vec3> color_map, 
					   std::vector<float> the_pressures,
					   std::vector<fcoord> bubbles,
					   std::vector<float> radii);
	~Simulation();
	void render();


private:
	std::vector<std::atomic<int> > volumes;


};



#endif