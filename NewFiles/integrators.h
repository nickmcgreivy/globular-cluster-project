#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>


namespace ODE_Integrator
{
	
	double t = 0;
	Eigen::MatrixXd (*calc_accelerations)(const Eigen::MatrixXd &positions, const Eigen::MatrixXd &velocities, double t);


	//Namespace that includes all of the integration methods
	namespace Integrators
	{	
		
		//Used to return the position and velocity arrays from the main integrator function
		struct ArrayPair 
		{
			
			Eigen::MatrixXd array1;
			Eigen::MatrixXd array2;

			ArrayPair(Eigen::MatrixXd &arr1, Eigen::MatrixXd &arr2) 
			{
				array1 = arr1;
				array2 = arr2;
			}
			

			ArrayPair()
			{

			}

			Eigen::MatrixXd &operator[](int i)
			{
				assert(i == 0 or i == 1);
				if (i == 0)
				{
					return array1;
				}
				if (i == 1)
				{
					return array2;
				}
			}
		
		};	

		void leap_frog(Eigen::MatrixXd &positions, Eigen::MatrixXd &velocities, double dt);

	}


	ODE_Integrator::Integrators::ArrayPair Integrate(Eigen::MatrixXd positions, Eigen::MatrixXd velocities, double dt, Eigen::MatrixXd (*calc_acc)(const Eigen::MatrixXd &positions, const Eigen::MatrixXd &velocities, double t), std::string method);	


}