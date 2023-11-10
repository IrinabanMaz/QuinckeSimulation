#include"ElectroHydroSim.h"
#include "fluidtests.h"

int main()
{
	testFluid1();
	testFluid2();
	testFluid3();
	testFluid4();
	testFluid5();

	std::fstream f("q-data.txt");

	std::vector<double> epsminus = { 1.2 , 0.8};
	std::vector<double> sigmaminus = { 0.8 , 1.2 };

	std::vector<int> p = { 6, 8 , 10 , 12, 14 , 16 };
	std::vector<double> maxqs = { 0.0 , 0.0 ,0.0 ,0.0 ,0.0, 0.0 , 0.0 };
	

	LaplaceParticleSystem q;
	for (int i = 0; i < 6; i++)
	{
		ElectroHydroSim sim(RectCoord(0.0, 0.0, 5.0), 1.0, epsminus, 1.0, sigmaminus, p[i]);

		std::vector<SphereCoordSystem> spheres;

		spheres.push_back(RectCoord(0.0, 0.0, 1.1));
		spheres.push_back(RectCoord(0.0, 0.0, -1.1));
		//spheres.push_back(RectCoord(0.0, 3.0, 0.0));
		//spheres.push_back(RectCoord(0.0, -3.0, 0.0));

		sim.initialize(spheres);

		LaplaceParticleSystem lastq;
	
		double dqdt = 0.0;
		double lastdqdt = 100;

		for (int j = 0; j < 100; j++)
		{
			

			
		
			sim.solveElectric();
			sim.timeStepOhmicForward(0.01);

			std::cout << " Polarization at timestep " << j << " for particle 1:\n {";
			for (int j = 0; j < sim.q.particles[0].consts.NUMGLNODES; j++)
				std::cout << sim.q.particles[0].data[j][0] << ", ";
			std::cout << "} \n\n";

			std::cout << " Polarization at timestep " << j << " for particle 2:\n {";
			for (int j = 0; j < sim.q.particles[1].consts.NUMGLNODES; j++)
				std::cout << sim.q.particles[1].data[j][0] << ", ";

			std::cout << "} \n\n";

			
			/*
			std::cout << " Polarization at timestep " << i << " for particle 3:\n {";
			for (int j = 0; j < sim.q.particles[2].consts.NUMTRAPNODES; j++)
				std::cout << sim.q.particles[2].data[sim.q.particles[2].consts.NUMGLNODES / 4][j] << ", ";

			std::cout << "} \n\n";

			std::cout << " Polarization at timestep " << i << " for particle 4:\n {";
			for (int j = 0; j < sim.q.particles[3].consts.NUMTRAPNODES; j++)
				std::cout << sim.q.particles[3].data[sim.q.particles[3].consts.NUMGLNODES / 4][j] << ", ";

			std::cout << "} \n\n";
			*/
		


			std::cout << "Change in q:" << (dqdt = sim.OhmicCurrent.norm()) << "\n";

			if (dqdt > lastdqdt)
				break;
		
			lastdqdt = dqdt;
			q = lastq = sim.q;
		}
		maxqs[i] = lastq.particles[0].data[0][0];
	}
	
	std::cout << "Convergence of max of q: {";
	for (int i = 0; i < 5; i++)
		std::cout << abs(maxqs[i] - maxqs[5]) << " ";

	std::cout << "}\n";
	
	SphericalVectorField identity(&SphereToRect);
	SphereData particle1points = discretize(&identity, 
		                                    &q.particles[0].system, 
		                                     q.particles[0].consts);
	
	
	f << "Particle 1 points: \n" << particle1points;

	f << " Particle 1 charge data: \n";
	f << q.particles[0].data;

	SphereData particle2points = discretize(&identity, 
		                                    &q.particles[1].system, 
		                                    q.particles[1].consts);
	f << "Particle 2 points: \n" << particle2points;



	f << " Particle 2 charge data: \n";
	f << q.particles[1].data;

	f.close();


	return 0;
}

