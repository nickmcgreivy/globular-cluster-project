#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>

namespace ODE_Integrator
{

	//Namespace that includes some example forces
	namespace Forces
	{
		namespace Gravity
		{
			void set_masses(Eigen::VectorXf &masses_);
		
			Eigen::VectorXd calc_gravity(int index_of_1, const Eigen::VectorXd pos1, const Eigen::VectorXd vel1, int index_of_2, const Eigen::VectorXd pos2, const Eigen::VectorXd vel2);

		}

		
		Eigen::MatrixXd gravity(const Eigen::MatrixXd &positions, const Eigen::MatrixXd &velocities, double t);
	
	}

}
