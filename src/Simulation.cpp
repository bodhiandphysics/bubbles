
#include "marchingcubes.h"


bool is_present(int color, int* colors, int clength) {

	for (int i = 0; i < clength; i++) {
		if (colors[i] == color) return true;
	}
	return false;
}

Simulation::Simulation(size_t x_max, 
					   size_t y_max, 
					   size_t z_max,
					   float scale,  
					   std::vector<glm::vec3> color_map, 
					   std::vector<float> the_pressures,
					   std::vector<fcoord> bubbles,
					   std::vector<float> radii):  x_max(x_max),
												   y_max(y_max),
												   z_max(z_max),
												   border_length(0),
												   bborder_length(0),
												   grid(x_max, y_max, z_max), 
												   scale(scale), color_map(color_map), 
												   num_colors(color_map.size()) {




	int num_threads = std::min(omp_get_num_threads(), 10);
	std::random_device r;

	for (int i = 0; i < num_threads; i++) {
		randoms.emplace_back(r());
	}

	volumes = (std::atomic<int>*) malloc(num_colors * sizeof(std::atomic<int>));
	for(int i = 0; i < num_colors; i++)
		volumes[i] = 0;

	int free_volume = 0;
	for (long i = 0; i < x_max; i++) {
		for (long j = 0; j < y_max; j++) {
			for (long k = 0; k < z_max; k++) {

				auto cell = grid(i, j, k);
				cell->on_border = false;
				cell->x = i;
				cell->y = j;
				cell->z = k;
				if (i == 0 ||
					i == x_max - 1 ||
					j == 0 ||
					j == y_max - 1 ||
					k == 0 ||
					k == z_max - 1) {
					cell->color = BOUNDARY;
					cell->new_color = BOUNDARY;
					cell->old_color = BOUNDARY;
				}
				else {
					cell->color = FREESPACE;
					cell->new_color = FREESPACE;
					cell->old_color = FREESPACE;
					free_volume++;
				} 
			}
		}	
	}
	volumes[FREESPACE] = free_volume;		


	for (int bnum = 0; bnum < num_colors - 1; bnum++) {

		for (long i = 0; i < x_max; i++) {
			for (long j = 0; j < y_max; j++) {
				for (long k = 0; k < z_max; k++) {

					auto cell = grid(i, j, k);

					float xcoord = i*scale;
					float ycoord = j*scale;
					float zcoord = k*scale;
					auto bubble_coord = bubbles[bnum];
					float x0 = bubble_coord.x;
					float y0 = bubble_coord.y;
					float z0 = bubble_coord.z;
					float r2 = pow(radii[bnum], 2);

					if (pow(xcoord - x0, 2) + pow(ycoord - y0, 2) + pow(zcoord - z0, 2) < r2) {

						assert(cell->color == FREESPACE || cell->color == BOUNDARY);
						if (cell->color != BOUNDARY) {
							volumes[cell->color]--;
							cell->color = bnum + 1; // first color is always freespace
							cell->new_color = bnum + 1;
							cell->old_color = bnum + 1;
							volumes[bnum + 1]++;
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < num_colors; i++) {

		P0V0.push_back(the_pressures[i] * volumes[i]);
		pressures.push_back(the_pressures[i]);
	}

	std::vector<coord> first_border;
	for (long i = 0; i < x_max; i++) {
		for (long j = 0; j < y_max; j++) {
			for (long k = 0; k < z_max; k++) {

				auto cell = grid(i, j, k);
				if (on_border(i, j, k)) {
					cell->on_border = true;
					first_border.push_back(coord(i, j, k));
					border_length++; 
				}
			}
		}	
	}

	border = (coord*) malloc(first_border.size() * sizeof(coord));
	for (int i = 0; i < first_border.size(); i++)
		border[i] = first_border[i];

	std::vector<coord> first_bborder;
	for (long i = 0; i < x_max; i++)
			for (long j = 0; j < y_max; j++)
				for (long k = 0; k < z_max; k++) {

					auto cell = grid(i, j, k);
					if (!cell->on_border && next_to_border(i, j, k)) {
						first_bborder.push_back(coord(i,j,k));
						bborder_length++; 
					}
	}

	bborder = (coord*) malloc(first_bborder.size() * sizeof(coord));
	for (int i = 0; i < first_border.size(); i++)
		bborder[i] = first_bborder[i];

}

float Simulation::roll01() {

	int tn = omp_get_thread_num();

	auto& engine = randoms[tn % randoms.size()];

	std::uniform_real_distribution<float> distr(0.f, 1.f);

	return distr(engine);

}

void Simulation::rebuild() {

	std::atomic<int> new_border_length = 0;
	std::atomic<int> new_bborder_length = 0;

	coord* new_border = (coord*) malloc(x_max * y_max * z_max*sizeof(coord));
	coord* new_bborder = (coord*) malloc(x_max * y_max * z_max*sizeof(coord));
	std::atomic<bool>* already_placed_bb = (std::atomic<bool>*) calloc(x_max * y_max * z_max, sizeof(std::atomic<bool>));
	std::atomic<bool>* already_placed_b = (std::atomic<bool>*) calloc(x_max * y_max * z_max, sizeof(std::atomic<bool>));

	int bl = border_length;

	#pragma omp parralel

	{
		int nthreads = omp_get_num_threads();
		int thread_num = omp_get_thread_num();
		int* dvolumes = (int*) calloc(num_colors, num_colors * sizeof(int));
		int span = border_length / nthreads;
		int start = thread_num * span;
		int end = thread_num == nthreads - 1? border_length.load(): (thread_num+1) * span;
		for (int i = start; i < end; i++) {

			coord current = border[i];
			auto* cell = grid(current.x, current.y, current.z);
			cell->on_border = false;
			for (int x = current.x - 1; x <= current.x + 1; x++) {
				for (int y = current.y - 1; y <= current.y + 1; y++){
					for (int z = current.z; z <= current.z + 1; z++) {

						auto new_border_cell = grid(x, y, z);
						if (!already_placed_b[x + y*x_max + z*x_max*y_max].load(std::memory_order_relaxed) &&
						    !already_placed_b[x + y*x_max + z*x_max*y_max].exchange(true, std::memory_order_acquire) &&
						    on_border(x, y, z, false)) {
							new_border[new_border_length++] = coord(x, y, z);
						} 
					}
				}
			}

			
			if (cell && cell->old_color != cell->new_color) {
				dvolumes[cell->old_color]--;
				dvolumes[cell->new_color]++;
				cell->old_color = cell->new_color;
				cell->color = cell->new_color;
			}
			
		}

		for (int i = 0; i < num_colors; i++)
			volumes[i] += dvolumes[i];
		free(dvolumes);	
	}

	for (int i = 0; i < new_border_length; i++) {

		if (auto cell = grid(new_border[i].x, new_border[i].y, new_border[i].z)) cell->on_border = true;
	}

	//synch

	#pragma omp parrallel

	{	
		int thread_num = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int span = floor(new_border_length/nthreads);
		int start = thread_num*(span - 1);
		int end = thread_num == nthreads - 1? new_border_length.load(): thread_num*span;

		for (int i = start; i < end; i++) {

			coord current = new_border[i];
			for (int x = current.x - 1; x <= current.x; x++) {
				for (int y = current.y - 1; y <= current.y; y++) {
					for (int z = current.z - 1; z <= current.z; z++) {

						if (grid(x, y, z) && 
							!grid(x, y, z)->on_border && 
							!already_placed_bb[x + y*x_max + z*x_max*y_max].load(std::memory_order_relaxed) &&
							!already_placed_bb[x + x_max*y + x_max*y_max*z].exchange(true, std::memory_order_acquire)) { 
							bborder[new_bborder_length++] = coord(x, y, z);
						}
					}
				}
			}
		}
	}

		
	

	//synch
	
	for (int i = 0; i < num_colors; i++)
			pressures[i] = P0V0[i]/volumes[i]; 

	free(border);
	free(bborder);
	free(already_placed_bb);
	free(already_placed_b);

	border = new_border;
	bborder = new_bborder;

	border_length = new_border_length.load();
	bborder_length = new_bborder_length.load();

}

void Simulation::flip(long x, long y, long z) {

	int local_colors[27];
	float probs[27];
	float sa[27];
	int num_local_colors = 0;
	auto cell = grid(x, y, z);
	if (!cell) return;
	int old_color = cell->color;
	if (old_color == BOUNDARY) return;
	for (long k = z-1; k <= z+1; k++)
		for (long j = y - 1; j <= y+1; j++)
			for (long i = x-1; i <= x+1; i++) {

				auto other = grid(i, j, k);
				if (other && other->color != BOUNDARY && !is_present(other->color, local_colors, num_local_colors)) {
					local_colors[num_local_colors++] = other->color;
				}
	}

	float sum_likelihood = 0;

	float delta_sa = 0;
	for (int c = 0; c < num_local_colors; c++) {
		int new_color = local_colors[c];
		for (int i = 0; i <= 1; i++) {
			for (int j = 0; j <= 1; j++) {
				for (int k = 0; k <= 1; k ++) {
					delta_sa += get_surface_area(x - i, y - j, z - k, false, true, new_color, x, y, z);
				}
			}	
		}

		float dp = scale*scale*scale*(pressures[old_color] - pressures[new_color]);

		float delta_e_local = ALPHA*dp + BETA*pow(delta_sa, .5);

		float likelihood = exp(-delta_e_local);

		probs[c] = likelihood;
		sum_likelihood += likelihood;
	}

	for (int c = 0; c < num_local_colors; c++)
		probs[c] /= sum_likelihood;

	float roll = roll01();

	float agprob = 0;

	for (int c = 0; c < num_local_colors; c++) {

		agprob += probs[c];
		if (roll < agprob) {

			if (cell->new_color != local_colors[c]) {
				cell->new_color = local_colors[c];
			}
			break;
		}
	}
}



float Simulation::get_de() {

	std::atomic<float> return_eng = 0;

	#pragma omp parralel

	{

		float return_eng_t = 0;

		int tn = omp_get_thread_num();
		int nt = omp_get_num_threads();

		for (int i = tn; i < border_length; i += nt) {

			auto cell = grid(border[i].x, border[i].y, border[i].z);
			float sa_o = get_surface_area(border[i].x, border[i].y, border[i].z, true); 
			float sa_new = get_surface_area(border[i].x, border[i].y, border[i].z, false);

			float dsa = sa_new - sa_o;

			float dp = scale*scale*scale*(pressures[cell->color] - pressures[cell->new_color]);

			return_eng_t += ALPHA*dp + BETA*pow(dsa, .5);

		}
		for (int i = tn; i < bborder_length; i += nt) {

			auto cell = grid(bborder[i].x, bborder[i].y, bborder[i].z);
			float sa_o = get_surface_area(bborder[i].x, bborder[i].y, bborder[i].z, true); 
			float sa_new = get_surface_area(bborder[i].x, bborder[i].y, bborder[i].z, false);

			float dsa = sa_new - sa_o;

			return_eng_t +=  BETA*dsa;

		}

		return_eng = return_eng.load() + return_eng_t;

	}

	return return_eng;
}

void Simulation::flip_cells() {

	#pragma omp parralel

	{
		int tn = omp_get_thread_num();
		int nt = omp_get_num_threads();

		for (int i = tn; i < border_length; i += nt) {
				flip(border[i].x, border[i].y, border[i].z);
		}
	}
}

void Simulation::confirm_flip() {

	#pragma omp parralel
	{

		int tn = omp_get_thread_num();
		int nt = omp_get_num_threads();
		for (int i = tn; i < border_length; i += nt) {
			if(auto* cell = grid(border[i].x, border[i].y, border[i].z)) 
				cell->color = cell->new_color;
		}
    }	
}

void Simulation::run_sim_iteration() {

	float old_de = 10000000;

	for (int i = 0; i < 30; i++) {

		flip_cells();

		float new_de = get_de();

		float prob = std::min(1.0, exp(-(new_de - old_de)));

		if (roll01() < prob) {
			confirm_flip();
			old_de = new_de;
		}
	}

	rebuild();
}
