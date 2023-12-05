#include"ElectroHydroSim.h"
#include "fluidtests.h"

int main()
{
/*
	std::fstream outputdata("output_data_in_plane8.txt");
	std::array<double, 14> rs = { 3.0 , 2.9 , 2.8 , 2.7 , 2.6 , 2.5 , 2.4 , 2.3, 2.25 , 2.2,2.15,2.1,2.05 , 2.01 };
	std::array<double, 14> MUFs = { 1.53416156 , 1.55461831 ,  1.57664602, 1.60039150, 1.62599574 , 1.65356957 ,
								   1.68314045,1.71452988,1.73072508, 1.74703222, 1.76311004, 1.77826951, 1.79070892,1.79223228 };
	std::array<double, 14> MOmegaFs = { 0.140187608,0.149081378,0.158630515,0.168801210,0.179481457,0.190411103,
										0.201041187 ,0.210224103 ,0.213577673 ,0.215412007 , 0.214756269 , 0.209731929 , 0.195478184 , 0.159607490 };

	std::time_t t0 = time(0);
	for (int i = 0; i < 6; i++)
	{
		testMobilityWilson3SphereInPlane(i, 8, outputdata);
		outputdata << " total time elapsed: " << time(0) - t0 << "\n";

		std::cout << "s = " << rs[i] << "\n";
		std::cout << " total time elapsed: " << time(0) - t0 << "\n";
	}

	outputdata.close();

	//testMobilityWilson3SpherePerp(2.05, 1.79070892, 0.195478184,8);
	 //SolveMobilitySingleSphere(RectCoord(0.0, 0.0, 6.0 * PI), RectCoord(), 8);
	 */

	//testFluid1();
	//testFluid2();
	//testFluid3();
	//testFluid4();

	std::fstream f("q-data.txt");

	std::vector<double> epsminus = { 1.2 , 0.8};
	std::vector<double> sigmaminus = { 0.8 , 1.2 };

	std::vector<int> p = {8,10,12};
	std::vector<double> maxqs;
	

	LaplaceParticleSystem q;
	for (int i = 0; i < p.size(); i++)
	{
		ElectroHydroSim sim(RectCoord(0.0 , 0.0, 4.0),
			                1.0, 
			                epsminus, 
			                1.0, 
			                sigmaminus,
			                p[i]);

	
		std::vector<SphereCoordSystem> spheres;

		spheres.push_back(RectCoord(4.1, 0.0, 0.0));
		spheres.push_back(RectCoord(0.0, 0.0, 0.0));
		//spheres.push_back(RectCoord(0.0, 3.0, 0.0));
		//spheres.push_back(RectCoord(0.0, -3.0, 0.0));

		sim.initialize(spheres);
		LaplaceParticleSystem lastq;
	
		double dqdt1 = 0.0;
		double dqdt2 = 0.0;
		double lastdqdt1 = 100;
		double lastdqdt2 = 100;

		RectCoord particle1force, particle2force, particle1torque, particle2torque;

		for (int j = 0; j < 100; j++)
		{
			

			
			sim.solveElectric();
			sim.computeElectricForces();
			sim.solveFluid();
			sim.timeStepForward(0.01);
			

			size_t nodes = sim.q.particles[0].consts.NUMGLNODES;
			std::cout << " Polarization at timestep " << j << " for particle 1:\n";
		
				std::cout << sim.q.particles[0].data << "\n\n";
			

			nodes = sim.q.particles[1].consts.NUMGLNODES;
			std::cout << " Polarization at timestep " << j << " for particle 2:\n {";
							std::cout << sim.q.particles[1].data << "\n\n";


			
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
		


			std::cout << "Change in q:" << (dqdt1 = sim.OhmicCurrent.norm()) << "\n";

			if (dqdt1 > lastdqdt1)
				break;
		
			lastdqdt1 = dqdt1;
			q = lastq = sim.q;

			particle1force = Integrate(&sim.fluidSolution,
				&sim.spheres[0],
				sim.electricTraction.particles[0].consts);

			particle2force = Integrate(&sim.fluidSolution,
				&sim.spheres[1],
				sim.electricTraction.particles[1].consts);

			SphereData fluidData1 = discretize(&sim.fluidSolution,
				&spheres[0],
				sim.fluidSolution.particles[0].consts);
			particle1torque = rotationCoefficient(fluidData1,
				spheres[0],
				sim.electricTraction.particles[0].consts);

			SphereData fluidData2 = discretize(&sim.fluidSolution,
												&spheres[1],
												sim.fluidSolution.particles[1].consts);
			particle1torque = rotationCoefficient(fluidData2,
													spheres[1],
				sim.electricTraction.particles[1].consts);

			 particle2torque = rotationCoefficient(sim.electricTraction.particles[1].data,
				spheres[1],
				sim.electricTraction.particles[1].consts);

			std::cout << "Net velocity on particles : " << particle1force << "\n" << particle2force << "\n\n";
			std::cout << "Net angular velocity on particles : " << particle1torque << " \n" << particle2torque << "\n\n";

			std::cout << " Center of particles" << sim.spheres[0].center << "\n"
				<< sim.spheres[1].center << "\n\n";
			
		}
		maxqs.push_back( lastq.particles[0].data[0][0]);
		
	}
	
	std::cout << "Convergence of max of q: {";
	for (int i = 0; i < maxqs.size() - 1; i++)
		std::cout << abs(maxqs[i] - maxqs.front()) << " ";

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

