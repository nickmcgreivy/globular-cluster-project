#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>
#include "NewFiles/integrators.cpp"
#include "NewFiles/forces.cpp"

int main()
{

	std::ofstream myFile;
	myFile.open("NewFiles/Output/Data.csv");

	Eigen::Matrix<double,2,3> positions;
	Eigen::Matrix<double,2,3> velocities;
	Eigen::VectorXf masses(2);

	positions <<	-3,0,0,
					3,0,0;


	velocities <<	0,-.3,0,
					0,.3,0;


	masses << 1,1;

	double num_particles = positions.rows();

	ODE_Integrator::Forces::Gravity::set_masses(masses);

	std::cout << positions << std::endl;
	
	double t = 0;
	double dt = 0.5;
	
	ODE_Integrator::Integrators::ArrayPair output;

	for (int i = 0; i < num_particles; i++)
	{
		myFile << positions(i,0) << "," << positions(i,1) << "," << positions(i,2) << ",";
	}

	myFile << "time:" << "," << t << std::endl;


	while (t < 500)
	{

		t += dt;
		output = ODE_Integrator::Integrate(positions, velocities, dt, ODE_Integrator::Forces::gravity, "rk");
		
		positions = output[0];
		velocities = output[1];

		for (int i = 0; i < num_particles; i++)
		{
			myFile << positions(i,0) << "," << positions(i,1) << "," << positions(i,2) << ",";
		}

		myFile << "time:" << "," << t << std::endl;

		//std::cout << velocities << std:: endl << "time: " << t << std::endl;
	}

	std::cout << ODE_Integrator::Integrators::RungaKutta::c_m << std::endl;

	return 0;

}