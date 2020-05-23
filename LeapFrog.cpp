#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
using namespace std;


// Declares universal time step "h" and universal time "t"

double h = 0.005;
double t = 0;


//Overload the + , - , * operators and create a to_string function for <double> vectors
vector<double> operator+(const vector<double> &v1,const vector<double> &v2) {

	int size = v1.size();
	vector<double> v3(size);

	for (int i = 0; i < size; i ++) {
		v3[i] = v1[i] + v2[i];
	}

	return v3;
}

vector<double> operator-(const vector<double> &v1,const vector<double> &v2) {

	int size = v1.size();
	vector<double> v3(size);

	for (int i = 0; i < size; i ++) {
		v3[i] = v1[i] - v2[i];
	}

	return v3;	
}

vector<double> operator*(const double num, const vector<double> &v1) {

	int size = v1.size();
	vector<double> v3(size);

	for (int i = 0; i < size; i ++) {
		v3[i] = v1[i] * num;
	}

	return v3;
}

vector<vector<double> > operator*(const double num, const vector<vector<double> > &v1) {

	int size = v1.size();
	vector<vector<double> > v3(size);

	for (int i = 0; i < size; i ++) {
		v3[i] = num* v1[i];
	}

	return v3;
}

double magnitude(const vector<double> &v) {

	int size = v.size();
	double sum = 0;


	for (int i = 0; i < size; i ++) {

		sum += pow( v[i] , 2 );

	}

	return sqrt( sum );

}

// Makes it easier to print vectors as comma separated values  
string to_string(const vector<double> &v) {

	int size = v.size(); 
	string s = "";

	s = s + to_string(v[0]);

	for (int i = 1; i < size; i ++) {

		s = s + "," + to_string(v[i]);		

	}

	return s;
}




// Single particle. Keeps track of position, momentum, mass, and any other attributes 

struct Particle {

	vector<double> position;
	vector<double> curr_momentum;
	vector<double> past_momentum;
	double mass;

	Particle(double x,double y,double z, double vx, double vy, double vz, double mass_) {

		position.push_back(x);
		position.push_back(y);
		position.push_back(z);



		curr_momentum.push_back(mass_*vx);
		curr_momentum.push_back(mass_*vy);
		curr_momentum.push_back(mass_*vz);

		past_momentum.push_back(mass_*vx);
		past_momentum.push_back(mass_*vy);
		past_momentum.push_back(mass_*vz);

		mass = mass_;

	}

	void change_pos(const vector<double> &delta_pos) {
	
		position = position + delta_pos;
	
	}

	void change_mom(const vector<double> &delta_mom) {
		
		past_momentum = curr_momentum;
		curr_momentum = curr_momentum + delta_mom;
	
	}

	vector<double> get_mom() {

		vector<double> mom;
		mom = (.5) * ( past_momentum + curr_momentum );
		return mom;
	
	}

	string output_data() {

		string s = "";

		s = s + to_string(position) + ",";

		s = s + to_string(get_mom());

		return s;

	}

};

struct System {

	vector<Particle> sys;

	System() {
	}

	void add_particle(Particle &p) {

		sys.push_back(p);

	}

	Particle& operator[](int index) {

		return sys[index];

	}

	int size() {
		return sys.size();
	}

	string output_data() {

		string s = "";

		for (int i = 0; i < sys.size(); i++) {

			s = s + sys[i].output_data() + "  ,  ";


		}

		return s;

	}

};


// Force functions can be made more general, however I assume it is a force due to particle interactions
// Force ON p1 DUE TO p2, represented as a vector

//Gravitational Constant

double G = 1;

vector<double> force(const Particle &p1, const Particle &p2) {

	vector<double> f(3); 

	// This is specifically gravity 
	// ( G * m * m / |r|^3 ) * r

	f = ( G*p1.mass*p2.mass / ( pow( magnitude(p2.position - p1.position) , 3 ) ) )  *  (p2.position - p1.position);

	return f;

}

// Returns a vector, holding vector<double>, representing the force on each particle in the system.

vector<vector<double> > net_force(System sys) { 
	
	int size = sys.size();

	vector<vector<double> > net_f(size);
	
	for(int i = 0; i < size; i ++) {

		vector<double> net_force_on_i(3);
		
		for(int j = 0; j < size; j ++) {

			if (not (i == j)) {

				net_force_on_i = net_force_on_i + force(sys[i], sys[j]);

			}
		}

		net_f[i] = net_force_on_i;

	}

	return net_f;

}

// The driver for the leap-frogging algorithm
// For the first time step in leap-frogging, all velocities need to be advanced by h/2

bool already_advanced = false;

void leap_frog(System &sys) {

	int size = sys.size();

	vector<vector<double> > f;
	f = net_force(sys);

	if( not already_advanced ) {

		for (int i = 0; i < size; i++)
			
			sys[i].change_mom( (h/2) * f[i] );

		already_advanced = true;

	}

	else {

		for (int i = 0; i < size; i++ ) {

			sys[i].change_pos( (h / sys[i].mass) * sys[i].curr_momentum );

			sys[i].change_mom( h * f[i] );
		
		}

		t += h;
	}
}


int main() {

	ofstream myFile;
	myFile.open("test_sim_values.csv");

	System sys;

	Particle* x;
	for (int i = 0; i < 360; i = i + 20) {

		float rads1 = i * (3.14 / 180);

		x = new Particle(20*cos(rads1),20*sin(rads1),0,0,0,0,10);
		sys.add_particle(*x);

	}

	for (int i = 0; i < 4000; i++) {

		myFile << sys.output_data() + "time:" + to_string(t) << endl;
		cout << "Simulated: " + to_string(i) << endl;
		leap_frog(sys);

	}

	return 0;
}

