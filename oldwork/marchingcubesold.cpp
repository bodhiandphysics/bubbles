
#include "marchingcubes.h"

VoxelGrid::VoxelGrid(size_t x_max, size_t y_max, size_t z_max): grid(x_max, y_max, z_max) {

	for (int i = 0; i < x_max * y_max * z_max; i++)
		grid.data[i].color = 0;

	colors.push_back(1);
	color_map.push_back(glm::vec3{0.f,0.f,0.f});
	color_map.push_back(glm::vec3{1.f, 0.f, 0.f});
}

VoxelGrid::VoxelGrid(float x_max, float y_max, float z_max, float ascale): VoxelGrid(long(x_max / ascale), long(y_max / ascale), long(z_max / ascale)) {
	scale = ascale;
}	

u_char VoxelGrid::get_index(uint color, long x, long y, long z) {

	u_char retval = 0;

	auto cell = grid(x,y,z);
	if (cell && cell->color == color) retval |= 0b00000001;
	
	cell = grid(x + 1,y,z);
	if (cell && cell->color == color) retval |= 0b00000010;
	
	cell = grid(x + 1,y,z + 1);	
	if (cell && cell->color == color) retval |= 0b00000100;
	
	cell = grid(x ,y,z + 1);
	if (cell && cell->color == color) retval |= 0b00001000;
	
	cell = grid(x,y + 1,z);
	if (cell && cell->color == color) retval |= 0b00010000;
	
	cell = grid(x + 1,y + 1,z);
	if (cell && cell->color == color) retval |= 0b00100000;
	
	cell = grid(x + 1,y + 1,z + 1);
	if (cell && cell->color == color) retval |= 0b01000000;
	
	cell = grid(x,y + 1,z + 1);
	if (cell && cell->color == color) retval |= 0b10000000;

	return retval;

}
u_char VoxelGrid::get_index_flipped(int color, long x, long y, long z, int fcolor, long fx, long fy, long fz, bool old) {

	u_char retval = 0;

	int check_color = 0; 

	auto cell = grid(x,y,z);
	if (cell) {
		if (coord(x, y, z) == coord(fx, fy, z))
			check_color = fcolor;
		else check_color = old? cell->old_color: cell->new_color;
	} else check_color = BOUNDARY_COLOR;
	if (check_color == color) retval |= 0b00000001;
		
	cell = grid(x + 1,y,z);
	if (cell) {
		if (coord(x + 1, y, z) == coord(fx, fy, z))
			check_color = fcolor;
		else check_color = old? cell->old_color: cell->new_color;
		} else check_color = BOUNDARY_COLOR;
		if (check_color == color) retval |= 0b00000010;

	
	cell = grid(x + 1,y,z + 1);	
	if (cell) {
		if (coord(x + 1, y, z + 1) == coord(fx, fy, z))
			check_color = fcolor;
		else check_color = old? cell->old_color: cell->new_color;
	} else check_color = BOUNDARY_COLOR;
	if (check_color == color) retval |= 0b00000100;
	
	
	cell = grid(x ,y,z + 1);
	if (cell) {
		if (coord(x, y, z + 1) == coord(fx, fy, z))
			check_color = fcolor;
		else check_color = old? cell->old_color: cell->new_color;
	} else check_color = BOUNDARY_COLOR;
	if (check_color == color) retval |= 0b00001000;
	
	
	
	cell = grid(x,y + 1,z);
	if (cell) {
		if (coord(x, y + 1, z) == coord(fx, fy, z))
			check_color = fcolor;
		else check_color = old? cell->old_color: cell->new_color;
	} else check_color = BOUNDARY_COLOR;
	if (check_color == color) retval |= 0b00010000;
	
	
	
	cell = grid(x + 1,y + 1,z);
	if (cell) {
		if (coord(x + 1, y + 1, z) == coord(fx, fy, z))
			check_color = fcolor;
		else check_color = old? cell->old_color: cell->new_color;
	} else check_color = BOUNDARY_COLOR;
	if (check_color == color) retval |= 0b00100000;
	
	
	cell = grid(x + 1,y + 1,z + 1);
	if (cell) {
		if (coord(x + 1, y + 1, z + 1) == coord(fx, fy, z))
			check_color = fcolor;
		else check_color = old? cell->old_color: cell->new_color;
	} else check_color = BOUNDARY_COLOR;
	if (check_color == color) retval |= 0b01000000;
	
	
	cell = grid(x,y + 1,z + 1);
	if (cell) {
		if (coord(x, y + 1, z + 1) == coord(fx, fy, z))
			check_color = fcolor;
		else check_color = old? cell->old_color: cell->new_color;
	} else check_color = BOUNDARY_COLOR;
	if (check_color == color) retval |= 0b10000000;
	

	return retval;

}

glm::vec3 VoxelGrid::get_point(char edge, glm::vec3 cube_pos) {

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

void VoxelGrid::get_triangles(long x, long y, long z, int color) {

	const glm::vec3 c = color_map[color];
	const char* table_p = TRIANGLE_TABLE[get_index(color, x,y,z)];

	glm::vec3 pos{x*scale, y*scale, z*scale};
	while (*table_p != -1) {

		glm::vec3 p1 = get_point(*table_p, pos);
		glm::vec3 p2 = get_point(*(table_p+1), pos);
		glm::vec3 p3 = get_point(*(table_p+2), pos);

		glm::vec3 normal = glm::normalize(glm::cross(p2-p1, p3-p1));

		points.push_back(p1);
		normals.push_back(normal);
		cols.push_back(c);	
		points.push_back(p2);
		normals.push_back(normal);
		cols.push_back(c);	
		points.push_back(p3);
		normals.push_back(normal);
		cols.push_back(c);

		table_p += 3;
	}
}

float VoxelGrid::get_surface_area(long x, long y, long z, int fcolor, int fx, int fy, int fz, bool old = true) {	

	float sa_return = 0; 

	for (int i = 0; i < num_colors; i++) {

		int color = colors[i];

		if (color == 0 || color == BOUNDARY_COLOR) continue; 

		const char* table_p = TRIANGLE_TABLE[get_index_flipped(color, x,y,z, fcolor, fx, fy, fz, bool old)];

		glm::vec3 pos{x*scale, y*scale, z*scale};
		while (*table_p != -1) {

			glm::vec3 p1 = get_point(*table_p, pos);
			glm::vec3 p2 = get_point(*(table_p+1), pos);
			glm::vec3 p3 = get_point(*(table_p+2), pos);

			float sa = glm::length(glm::cross(p2 - p1, p3 - p1)) / 2; 
			sa_return += sa;	
			table_p += 3;
		}

	}

	return sa_return;

}


bool VoxelGrid::on_border(long x, long y, long z) {

	for (long k = z-1; k <= z+1; k++)
		for (long j = y - 1; j <= y+1; j++)
			for (long i = x-1; i <= x+1; i++) {

				if (i == x && j == y && k = z) continue;
				if (!grid(i, j, k)) continue;

				if ((grid(i, j, k)).value.color != grid(x, y, z).value.color)
					return true;

			}

	return false;
}

bool VoxelGrid::on_border(coord at) {return on_border(at.x, at.y, at.z);}

bool VoxelGrid::next_to_border(long x, long y, long z) {

	for (long k = z-1; k <= z+1; k++)
		for (long j = y - 1; j <= y+1; j++)
			for (long i = x-1; i <= x+1; i++) {

				if (i == x && j == y && k = z) continue;
				if (!grid(i, j, k)) continue;

				if ((grid(i, j, k)).value.on_border);
					return true;

			}
			
	return false;
}

bool VoxelGrid::next_to_border(coord at) {return next_to_border(at.x, at.y, at.z);}

void VoxelGrid::rebuild() {

	std::atomic<int> delta_p[NUM_COLORS];

	std::atomic<int> new_border_length = 0;
	for (auto& i: delta_p) i = 0;

	coord* new_border = (coord*) malloc(x_max * y_max * z_max*sizeof(coord));
	coord* new_bborder = (coord*) malloc(x_max * y_max * z_max*sizeof(coord));
	std::atomic<bool>* place_bb = (std::atomic<bool>*) calloc(x_max * y_max * z_max*sizeof(coord));

	int bl = border_length;

	#pragma omp parralel

	{

		int thread_num = omp_get_thread_num();
		int nthreads = omp_gem_num_threads();

		int delta_pt[NUM_COLORS]
		for (int i = thread_num; i < bl; i+= nthreads) {
			if (on_border(border[i])) 
				new_border[new_border_length++] = border[i];

			auto cell = grid(border_i);
			if (cell->color != cell->new_color) {
				delta_pt[cell->color]--;
				delta_pt[cell->new_color]++;
				cell->color = cell->new_color;
				cell->old_color = cell->new_color;
			}
			
		}

		int bbl = bborder_length;

		for (int i = thread_num i < bbl; nthreads++) {
			if (on_border(bborder[i]))
				new_border[new_bborder_length++] = bborder[i];

			coord current = bborder[i];
		
			for (int x = current.x - 1; x <= current.x + 1; x++)
				for (int y = current.y - 1; y <= current.y + 1; y++)
					for (int z = current.z - 1; z <= current.z + 1; z++) {

					if (grid(x, y, z) && next_to_border(x, y, z) && !(place_bb[x + x_max*y + x_max*y_max*z].exchange(true))){ 
						bborder[new_bborder_length++] = coord(x, y, z);
					}
			}
		}

		for (auto& i: delta_pt)
			delta_p += delta_pt; 
	}

	//synch

	free(border);
	free(bborder);
	free(place_bb);

	border = new_border;
	bborder = new_bborder;

	border_length = new_border_length;
	bborder_length = new_bborder_length;

}

bool is_present(int color, int* colors, int clength) {

	int* p = colors;
	while(*p++ != -1) 
		if (p->first == color) return true;
	return false;
}

void flip(long x, long y, long z) {

	int colors[NUM_COLORS + 1];
	float probs[NUM_COLORS + 1];
	float sa[NUM_COLORS + 1]
	int num_colors = 0;
	colors[0]->first = -1;
	auto cell = grid(x, y, z);
	if (!cell) return;
	int old_color = cell->color;
	for (long k = z-1; k <= z+1; k++)
		for (long j = y - 1; j <= y+1; j++)
			for (long i = x-1; i <= x+1; i++) {

				auto other = grid(i, j, k);
				if (other && !is_present(other->color, colors, NUM_COLORS + 1) && other->color !=  BOUNDARY_COLOR) {
					colors[num_colors++] = color;
				}
	}

	int p = 0;

	float sum_likelihood = 0;

	float new_sa = 0;
	while(colors[p] != -1) {
		int new_color = colors[p];
		new_sa += get_surface_area(x,     y,     z,     colors, num_colors, new_color, x, y, z);
		new_sa += get_surface_area(x,     y,     z - 1, colors, num_colors, new_color, x, y, z);
		new_sa += get_surface_area(x,     y - 1, z,     colors, num_colors, new_color, x, y, z);
		new_sa += get_surface_area(x,     y - 1, z - 1, colors, num_colors, new_color, x, y, z);
		new_sa += get_surface_area(x - 1, y,     z,     colors, num_colors, new_color, x, y, z);
		new_sa += get_surface_area(x - 1, y,     z - 1, colors, num_colors, new_color, x, y, z);
		new_sa += get_surface_area(x - 1, y - 1, z,     colors, num_colors, new_color, x, y, z);
		new_sa += get_surface_area(x - 1, y - 1, z - 1, colors, num_colors, new_color, x, y, z);

		sa[p] = new_sa;

		float delta_sa = new_sa - cell->old_sa;

		float dp = pressures[old_color] - pressures[new_color];

		float delta_e_local = ALPHA*dp + BETA* delta_sa;

		float likelihood = exp(-delta_e_local);

		probs[p] = likelihood;
		sum_likelihood += likelihood;
		color_p++;
	}

	p = 0;

	while(colors[p] != -1) {

		probs[p] /= sum_likelihood;
		p++;
	}

	float roll = roll01();

	float aggprob = 0;

	p = 0;

	while(colors[p] != -1) {

		agprob += probs[p];
		if (roll < aggprob) {

			if (cell->new_color != colors[p]) {
				cell->new_color = colors[p];
				cell->old_sa = sa[p];
			}
			break;
		}
	}
}

void confirm_flip() {

	#pragma omp parralel
	{

		int tn = omp_get_thread_num();
		int nt = omp_gem_num_threads();
		for (int i = tn; i < border_length; i += nt;)
			grid(cell.x, cell.y, cell.z)->old_color = grid(cell.x, cell.y, cell.z)->new_color;
    }	
}

float get_de() {

	std::atomic<float> return_eng = 0;

	#pragma omp parralel

	{

		float return_eng_t = 0;

		int tn = omp_get_thread_num();
		int nt = omp_gem_num_threads();

		for (int i = tn; i < border_length; i += nt) {

			auto cell = grid(border[i].x, border[i].y, border[i].z);
			float sa_o = get_surface_area(border[i].x, border[i].y, border[i].z, true); 
			float sa_new = get_surface_area(border[i].x, border[i].y, border[i].z, false);

			float dsa = sa_new - sa_o;

			float dp = [cell->new_color] - pressures[cell->old_color];

			return_eng_t += alpha*dp + beta*dsa;

		}

		for (int i = tn; i < border_length; i += nt) {

			auto cell = grid(bborder[i].x, bborder[i].y, bborder[i].z);
			float sa_o = get_surface_area(bborder[i].x, bborder[i].y, bborder[i].z, true); 
			float sa_new = get_surface_area(bborder[i].x, bborder[i].y, bborder[i].z, false);

			float dsa = sa_new - sa_o;

			float dp = pressures[cell->new_color] - pressures[cell->old_color];

			return_eng_t += alpha*dp + beta*dsa;

		}

		return_eng += return_eng_t;

	}

	return return_eng;
}

void VoxelGrid::flip_cells() {

	#pragma omp parralel

	{
		int tn = omp_get_thread_num();
		int nt = omp_gem_num_threads();

		for (int i = tn; i < border_length; i += nt)
			flip(border[i].x, border[i].y, border[i].z);
	}
}

void VoxelGrid::run_sim_iteration() {

	float old_de = get_de();

	for (int i = 0; i < 30; i ++) {

		flip_cells();

		float new_de = get_de();

		float prob = std::max(1.f, exp(-(new_de - old_de)));

		if (roll01() < prob) {
			confirm_flip();
			old_de = new_de;
		}
	}

	rebuild();
}



void VoxelGrid::render(int color) {

	points.clear();
	normals.clear();
	colors.clear();
	surface_area = 0;

	for (auto pot: border) {

		auto index = pot - grid.data; 
		long i = index % grid.x_max;
		long j = (index - i) % grid.y_max;
		long k = index - i - j;

		get_triangles(i, j, k , color);
	}
}