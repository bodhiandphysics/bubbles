#include "marchingcubes.h"
#include <unordered_map>


glm::vec3 Simulation::get_point(char edge, glm::vec3 cube_pos) {

	assert(edge < 12 && edge >= 0);

	switch(edge) {

	case 0: 
		return cube_pos + (scale)*glm::vec3{.5f,0.f,0.f};
	case 1: 
		return cube_pos + (scale)*glm::vec3{1.f,0.f,.5f};
	case 2: 
		return cube_pos + (scale)*glm::vec3{.5f,0.f,1.f};
	case 3: 
		return cube_pos + (scale)*glm::vec3{0.f,0.f,0.5f};
	case 4: 
		return cube_pos + (scale)*glm::vec3{.5f,1.f,0.f};
	case 5: 
		return cube_pos + (scale)*glm::vec3{1.f,1.f,.5f};
	case 6: 
		return cube_pos + (scale)*glm::vec3{.5f,1.f,1.f};
	case 7: 
		return cube_pos + (scale)*glm::vec3{0.f,1.f,0.5f};	
	case 8: 
		return cube_pos + (scale)*glm::vec3{0.f,.5f,0.f};
	case 9: 
		return cube_pos + (scale)*glm::vec3{1.f,.5f,0.f};
	case 10: 
		return cube_pos + (scale)*glm::vec3{1.f,0.5f,1.f};
	case 11: 
		return cube_pos + (scale)*glm::vec3{0.f,0.5f,1.f};
	}

	return glm::vec3(0);  // should never happen
}

void Simulation::render() {

	points.clear();
	normals.clear();
	cols.clear();

	std::unordered_map<fcoord, size_t> point_map;
	std::vector<float> angles; 


	for (int i = 0; i < border_length; i++) {

		auto cell = grid(border[i].x, border[i].y, border[i].z);

		if (!cell) continue;

		int local_colors[8];
		int num_local_colors = 0;

		for (int i = 0; i <=1; i++) {
			for (int j = 0; j <= 1; j++) {
				for (int k = 0; k <= 1; k++) {
					auto other_cell = grid(cell->x + i, cell->y + j, cell->z + k);
					if (other_cell && !is_present(other_cell->color, local_colors, num_local_colors) && other_cell->color != BOUNDARY)
						local_colors[num_local_colors++] = other_cell->color;
				}
			}
		}


		for (int colorp = 0; colorp < num_local_colors; colorp++) {

			int color = local_colors[colorp];
			if (color == FREESPACE) continue;

			glm::vec3 c = color_map[color];
			const char* table_p = TRIANGLE_TABLE[get_index(color, cell->x,cell->y,cell->z, false)];

			glm::vec3 pos{cell->x*scale, cell->y*scale, cell->z*scale};
			while (*table_p != -1) {

				glm::vec3 p1 = get_point(*table_p, pos);
				glm::vec3 p2 = get_point(*(table_p+1), pos);
				glm::vec3 p3 = get_point(*(table_p+2), pos);

				glm::vec3 normal = glm::normalize(glm::cross(p2-p1, p3-p1));

				float theta1 = acos(glm::dot(p2 - p1, p3 - p1)/(glm::length(p2 - p1)*glm::length(p3 - p1)));
				float theta2 = acos(glm::dot(p3 - p2, p1 - p2)/(glm::length(p3 - p2)*glm::length(p1 - p2)));
				float theta3 = acos(glm::dot(p2 - p3, p1 - p3)/(glm::length(p2 - p3)*glm::length(p1 - p3)));

				auto fcoord1 = fcoord(p1.x, p1.y, p1.z, color);
				auto fcoord2 = fcoord(p2.x, p2.y, p2.z, color);
				auto fcoord3 = fcoord(p3.x, p3.y, p3.z, color);

				if (theta1 < 0) theta1 += M_PI;
				if (theta2 < 0) theta2 += M_PI;
				if (theta3 < 0) theta3 += M_PI;

				fcoord const* coords[3] = {&fcoord1, &fcoord2, &fcoord3};
				glm::vec3* ps[3] = {&p1, &p2, &p3};
				float* angs[3] {&theta1, &theta2, &theta3};

				for (int i = 0; i < 3; i++) {
					if (point_map.find(*coords[i]) == point_map.end()) {

					points.push_back(*ps[i]);
					cols.push_back(c);
					normals.push_back(normal * *angs[i]);
					angles.push_back(*angs[i]);

					point_map[*coords[i]] = points.size() - 1;
					indexes.push_back(points.size() - 1); 
					}
					else {
						size_t index = point_map[*coords[i]];
						normals[index] += normal * *angs[i];
						angles[index] += *angs[i];
						indexes.push_back(GLint(index));
					}

				}

				table_p += 3;
			}

		}

		
	}

	for (int i = 0; i < normals.size(); i++) {

		normals[i] /= angles[i];
		normals[i] = glm::normalize(normals[i]);
	}

	
}