#include "marchingcubes.h"
#include "table.h"

// returns the marching cubes index of a set of 8 grid cells starting at xyz
// based on its colors... allows applying a flipped cell in a different color


u_char Simulation::get_index_flipped(int color, long x, long y, long z, int fcolor, long fx, long fy, long fz, bool old) {
	
	u_char retval = 0;

	int check_color = 0; 

	for (int i = 0; i <= 1; i++)
		for (int j = 0; j <= 1; j++)
			for (int k = 0; k <= 1; k++) {
				auto cell = grid(x + i,y + j,z + k);
				if (cell) { 
					if (coord(x + i, y + j, z + k) == coord(fx, fy, z))
						check_color = fcolor;
					else check_color = old? cell->color: cell->new_color;
				} else check_color = BOUNDARY;
				if (color == check_color) retval |= MASK_TABLE[i][j][k];
	}
	return retval;

}

u_char Simulation::get_index(int color, long x, long y, long z, bool old) {

	u_char retval = 0;
	int check_color = 0;

	for (int i = 0; i <= 1; i++)
		for (int j = 0; j <= 1; j++)
			for (int k = 0; k <= 1; k++) {
				auto cell = grid(x + i,y + j,z + k);
				if (cell)
					check_color = old? cell->color: cell->new_color;
				else check_color = BOUNDARY;
				if (color == check_color) retval |= MASK_TABLE[i][j][k];
	}
	return retval;
}


float Simulation::get_surface_area(long x, long y, long z, bool old, 
							      bool flipped, int fcolor, 
								  int fx, int fy, int fz) {	

	int local_colors[8];
	int num_local_colors = 0;
	float sa_return = 0.f;


	for (int i = 0; i <= 1; i++)
		for (int j = 0; j <= 1; j++)
			for (int k = 0; k <=1; k++) {

				auto cell = grid(x + i, y + j, z + k);
				if (cell && old && is_present(cell->color, local_colors, num_local_colors))
					colors[num_local_colors++] == cell->color;
				if (cell && !old && is_present(cell->new_color, local_colors, num_local_colors))
					colors[num_local_colors++] == cell->new_color;
				
	}

	for (int i = 0; i < num_local_colors; i++) {

		int color = local_colors[i];

		if (color == 0 || color == BOUNDARY) continue; 

		if (flipped)
			sa_return += scale * scale * AREA_TABLE[get_index_flipped(color, x,y,z, fcolor, fx, fy, fz, old)]/2;
		else
			sa_return += scale * scale * AREA_TABLE[get_index(color, x, y, z, old)]/2;

	}

	return sa_return;

}

bool Simulation::on_border(long x, long y, long z) {

	int check_color = BOUNDARY;
	if (auto cell = grid(x, y, z))
		check_color = cell->color;
	for (long k = z-1; k <= z+1; k++)
		for (long j = y - 1; j <= y+1; j++)
			for (long i = x-1; i <= x+1; i++) {

				if (i == x && j == y && k == z) continue;
				if (!grid(i, j, k)) continue;

				if ((grid(i, j, k))->color != check_color)
					return true;

			}

	return false;
}

bool Simulation::on_border(coord at) {return on_border(at.x, at.y, at.z);}

bool Simulation::next_to_border(long x, long y, long z) {

	for (long k = 0; k <= 1; k++)
		for (long j = 0; j <= 1; j++)
			for (long i = 0; i <= x; i++) {

				if (i == 0 && j == 0 && k == 0 ) continue;
				if (!grid(i, j, k)) continue;

				if ((grid(x + i, y + j, z + k))->on_border);
					return true;

			}
			
	return false;
}

bool Simulation::next_to_border(coord at) {return next_to_border(at.x, at.y, at.z);}


