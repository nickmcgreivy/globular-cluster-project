#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
 
//double G = 6.67259 * pow(10,-11);
double G = 1;
using namespace std;

int n = 20;
vector<double > masses(n);
vector<vector<double> > positions(n);
vector<vector<double> > past_positions(n);


vector<double> r1(3);
vector<double> r2(3);
vector<double> r3(3);
double t = 0;
double h = 0.05;




//	Force function, takes in two vectors and returns the force the second
//	exerts on the first as a vector.

vector<double> force(vector<double> &r1, vector<double> &r2, float mass1, float mass2) {

	vector<double> f(3);

	double constant = G*mass1*mass2 / pow( pow((r2[0] - r1[0]),2) + pow((r2[1] - r1[1]), 2) + pow((r2[2] - r1[2]), 2) , 3/2 );

	for (int i = 0; i < 3; i++){

	f[i] = constant * (r2[i] - r1[i]);
	
	}
	
	return f;
}



// Takes a vector containing all objects in the system, and calculates net force on all of them.

vector<vector<double> > net_force(vector<vector<double> > &all_objs) {

	vector<vector<double> > f(n);

	

	
	for (int j = 0; j < n; j++) {
		
		vector<double> force_total(3);
		vector<double> temp_force(3);

		for (int x = 0; x < n; x++) {
			
			if (not (j == x)) {

				temp_force = force(all_objs[j], all_objs[x], masses[j], masses[x]);

				for (int i = 0; i < 3; i ++) {
					force_total[i] = force_total[i] + temp_force[i];
				}

			}}
			f[j] = force_total;
		}


	return f;
}


void verlet_step(vector<vector<double> > &r, vector<vector<double> > &r_prev, double &t, double &h) {

	vector<vector<double> > next_r(n);
	vector<vector<double> > all_forces(n);

	all_forces = net_force(r);

	for (int i = 0; i < n; i++) {
		
		vector<double> pos(3);
		next_r[i] = pos;
		
		for (int j = 0; j < 3; j++){
			
			next_r[i][j] = 2*r[i][j] - r_prev[i][j] + ((pow(h,2) / masses[i]) * all_forces[i][j]);


		
		}
		
	}

	past_positions = positions;
	positions = next_r;

}



void init() {

	for (int i = 0; i < n; i++) {

	
		vector<double> r(3);
		r[0] = 10*i;
		r[1] = -10*i;
		r[2] = 0;

		vector<double> r_prev(3);
		r_prev[0] = 10*i;
		r_prev[1] = -10*i + 0.05;
		r_prev[2] = 0;

		positions[i] = r;
		past_positions[i] = r_prev;
		masses[i] = 1;

	}
		/*/
		r[0] = 5;
		r[1] = -5;
		r[2] = 0;


		r_prev[0] = 5;
		r_prev[1] = -5 - 0.025;
		r_prev[2] = 0;

		positions[1] = r;
		past_positions[1] = r_prev;
		masses[1] = 2;



		r[0] = -5;
		r[1] = -5;
		r[2] = 0;


		r_prev[0] = -5 - 0.05;
		r_prev[1] = -5;
		r_prev[2] = 0;

		positions[2] = r;
		past_positions[2] = r_prev;
		masses[2] = 1;


		r[0] = -5;
		r[1] = 5;
		r[2] = 0;


		r_prev[0] = -5;
		r_prev[1] = 5 + 0.05;
		r_prev[2] = 0;

		positions[3] = r;
		past_positions[3] = r_prev;
		masses[3] = 1;
		
		/*/
	//}

}




int main() {

	init();

	ofstream myFile;
	myFile.open("test_sim_values.csv");


	for (int l = 0; l < 8000; l++) {

		string s = "";

		for (int i = 0; i < n; i++) {
			

			for (int j = 0; j < 3; j++) {
				s += to_string(positions[i][j]) + ","; 


			}

		}

		myFile << s << endl;

		verlet_step(positions, past_positions, t, h);
		t += h;


		}
	
}




