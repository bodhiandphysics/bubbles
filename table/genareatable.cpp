#include <glm/glm.hpp>
#include <iostream>
#include "table.h"
#include <assert.h>


glm::vec3 get_point(char edge) {

	assert(edge < 12 && edge >= 0);

	switch(edge) {

	case 0: 
		return glm::vec3{.5f,0.f,0.f};
	case 1: 
		return glm::vec3{1.f,0.f,.5f};
	case 2: 
		return glm::vec3{.5f,0.f,1.f};
	case 3: 
		return glm::vec3{0.f,0.f,0.5f};
	case 4: 
		return glm::vec3{.5f,1.f,0.f};
	case 5: 
		return glm::vec3{1.f,1.f,.5f};
	case 6: 
		return glm::vec3{.5f,1.f,1.f};
	case 7: 
		return glm::vec3{0.f,1.f,0.5f};	
	case 8: 
		return glm::vec3{0.f,.5f,0.f};
	case 9: 
		return glm::vec3{1.f,.5f,0.f};
	case 10: 
		return glm::vec3{1.f,0.5f,1.f};
	case 11: 
		return glm::vec3{0.f,0.5f,1.f};
	}

	return glm::vec3(0);  // should never happen
}
int main(int argc, char const *argv[])
{
		
	float areas[256];

	for (int i = 0; i < 256; i++) {

		float area = 0;

		for (int j = 0; j < 16; j+= 3) {

			if (TRIANGLE_TABLE[i][j] == -1)
				break;

			glm::vec3 p1 = get_point(TRIANGLE_TABLE[i][j]);
			glm::vec3 p2 = get_point(TRIANGLE_TABLE[i][j+1]);
			glm::vec3 p3 = get_point(TRIANGLE_TABLE[i][j+2]);

			glm::vec3 s1 = p2 - p1;
			glm::vec3 s2 = p3 - p1;

			area += glm::length(glm::cross(s1, s2)) / 2.0;
			assert(area >= 0);
		}

		areas[i] = area;
	}



	std::cout << "{\n";
	for (int i = 0; i < 256; i++)
		std::cout << areas[i] << ",\n";

	std::cout <<"}\n"; 
	return 0;
}