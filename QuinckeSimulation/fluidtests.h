#pragma once
#include "C:\Users\David\source\repos\IrinabanMaz\SphericalHarmonics\MobilityProblem.h"

void testFluid1()
{
	std::fstream f("testfluid1.txt");

	const double PI = 4.0 * atan(1.0);

	std::vector<RectCoord> forces;
	std::vector<RectCoord> torques;
	std::vector<SphereCoordSystem> spheres;

	RectCoord F1(-6.0 * PI, 0.0, 0.0);
	forces.push_back(F1);
	torques.push_back(RectCoord());

	SphereCoordSystem sphere1(RectCoord(2.1, 0.0, 0.0));
	spheres.push_back(sphere1);

	RectCoord F2(0.0 , 0.0 , 0.0);
	forces.push_back(F2);
	torques.push_back(RectCoord());

	SphereCoordSystem sphere2;
	spheres.push_back(sphere2);

	StokesParticleSystem solution = SolveMobilityManySphere(spheres,
		forces,
		torques,
		16,
		StokesParticleSystem(),
		30,
		1e-5);

	RectCoord U1 = Integrate(&solution, &sphere1, solution.particles[0].consts) / (4.0 * PI);
	RectCoord U2 = Integrate(&solution, &sphere2, solution.particles[1].consts) / (4.0 * PI);

	SphereData flow1 = discretize(&solution, &sphere1, solution.particles[0].consts);
	SphereData flow2 = discretize(&solution, &sphere2, solution.particles[1].consts);
	RectCoord Omega1 = rotationCoefficient(flow1, sphere1, solution.particles[0].consts) * 3.0 / 8.0 / PI;
	RectCoord Omega2 = rotationCoefficient(flow2, sphere2, solution.particles[1].consts) * 3.0 / 8.0 / PI;

	SphericalVectorField identity(&SphereToRect);
	SphereData particle1points = discretize(&identity, &sphere1, solution.particles[0].consts);
	f << "Particle 1 points: \n" << particle1points;


	f << "Particle 1 velocity: " << U1 << "\n";
	f << "Particle 1 angular velocity: " << Omega1 << "\n\n";

	f << " Particle 1 velocity data: \n";
	f << discretize(&solution, &sphere1, solution.particles[0].consts);

	ScalSphereData forcedensity(solution.particles[0].consts.NUMTRAPNODES, 
		                        solution.particles[0].consts.NUMGLNODES);
	for (int i = 0; i < forcedensity.size(); i++)
		for (int j = 0; j < forcedensity[0].size(); j++)
			forcedensity[i][j] = norm(solution.particles[0].data[i][j]);
	f << "Particle 1 force data: \n";
	f << forcedensity;

	SphereData particle2points = discretize(&identity, &sphere2, solution.particles[1].consts);
	f << "Particle 2 points: \n" << particle2points;


	f << "Particle 2 velocity: " << U2 << "\n";
	f << "Particle 2 angular velocity: " << Omega2 << "\n\n";

	f << " Particle 2 velocity data: \n";
	f << discretize(&solution, &sphere2, solution.particles[1].consts);

	for (int i = 0; i < forcedensity.size(); i++)
		for (int j = 0; j < forcedensity[0].size(); j++)
			forcedensity[i][j] = norm(solution.particles[1].data[i][j]);
	f << "Particle 2 force data: \n";
	f << forcedensity;


	f.close();
}

void testFluid2()
{
	std::fstream f("testfluid2.txt");

	const double PI = 4.0 * atan(1.0);

	std::vector<RectCoord> forces;
	std::vector<RectCoord> torques;
	std::vector<SphereCoordSystem> spheres;

	RectCoord F1(0.0, 6.0 * PI, 0.0);
	forces.push_back(F1);
	torques.push_back(RectCoord());

	SphereCoordSystem sphere1(RectCoord(2.1, 0.0, 0.0));
	spheres.push_back(sphere1);

	RectCoord F2(0.0 , 0.0 , 0.0);
	forces.push_back(F2);
	torques.push_back(RectCoord());

	SphereCoordSystem sphere2;
	spheres.push_back(sphere2);

	StokesParticleSystem solution = SolveMobilityManySphere(spheres,
		forces,
		torques,
		16,
		StokesParticleSystem(),
		30,
		1e-5);

	RectCoord U1 = Integrate(&solution, &sphere1, solution.particles[0].consts) / (4.0 * PI);
	RectCoord U2 = Integrate(&solution, &sphere2, solution.particles[1].consts) / (4.0 * PI);

	SphereData flow1 = discretize(&solution, &sphere1, solution.particles[0].consts);
	SphereData flow2 = discretize(&solution, &sphere2, solution.particles[1].consts);
	RectCoord Omega1 = rotationCoefficient(flow1, sphere1, solution.particles[0].consts) * 3.0 / 8.0 / PI;
	RectCoord Omega2 = rotationCoefficient(flow2, sphere2, solution.particles[1].consts) * 3.0 / 8.0 / PI;

	SphericalVectorField identity(&SphereToRect);
	SphereData particle1points = discretize(&identity, &sphere1, solution.particles[0].consts);
	f << "Particle 1 points: \n" << particle1points;


	f << "Particle 1 velocity: " << U1 << "\n";
	f << "Particle 1 angular velocity: " << Omega1 << "\n\n";

	f << " Particle 1 velocity data: \n";
	f << discretize(&solution, &sphere1, solution.particles[0].consts);

	ScalSphereData forcedensity(solution.particles[0].consts.NUMTRAPNODES,
		solution.particles[0].consts.NUMGLNODES);
	for (int i = 0; i < forcedensity.size(); i++)
		for (int j = 0; j < forcedensity[0].size(); j++)
			forcedensity[i][j] = norm(solution.particles[0].data[i][j]);
	f << "Particle 1 force data: \n";
	f << forcedensity;

	SphereData particle2points = discretize(&identity, &sphere2, solution.particles[1].consts);
	f << "Particle 2 points: \n" << particle2points;


	f << "Particle 2 velocity: " << U2 << "\n";
	f << "Particle 2 angular velocity: " << Omega2 << "\n\n";

	f << " Particle 2 velocity data: \n";
	f << discretize(&solution, &sphere2, solution.particles[1].consts);

	for (int i = 0; i < forcedensity.size(); i++)
		for (int j = 0; j < forcedensity[0].size(); j++)
			forcedensity[i][j] = norm(solution.particles[1].data[i][j]);
	f << "Particle 2 force data: \n";
	f << forcedensity;


	f.close();
}



void testFluid3()
{
	std::fstream f("testfluid3.txt");

	const double PI = 4.0 * atan(1.0);

	std::vector<RectCoord> forces;
	std::vector<RectCoord> torques;
	std::vector<SphereCoordSystem> spheres;

	RectCoord F1;
	forces.push_back(F1);
	RectCoord T1(0.0, 8.0 * PI, 0.0);
	torques.push_back(T1);

	SphereCoordSystem sphere1(RectCoord(2.1, 0.0, 0.0));
	spheres.push_back(sphere1);

	RectCoord F2;
	forces.push_back(F2);
	RectCoord T2(0.0, 0.0, 0.0);
	torques.push_back(T2);

	SphereCoordSystem sphere2;
	spheres.push_back(sphere2);

	StokesParticleSystem solution = SolveMobilityManySphere(spheres,
		forces,
		torques,
		16,
		StokesParticleSystem(),
		30,
		1e-5);

	RectCoord U1 = Integrate(&solution, &sphere1, solution.particles[0].consts) / (4.0 * PI);
	RectCoord U2 = Integrate(&solution, &sphere2, solution.particles[1].consts) / (4.0 * PI);


	SphereData flow1 = discretize(&solution, &sphere1, solution.particles[0].consts);
	RectCoord Omega1 = rotationCoefficient(flow1, sphere1, solution.particles[0].consts) *(3.0 / 8.0 / PI);

	SphereData flow2 = discretize(&solution, &sphere2, solution.particles[1].consts);
	RectCoord Omega2 = rotationCoefficient(flow2, sphere2, solution.particles[1].consts) * (3.0 / 8.0 / PI);


	SphericalVectorField identity(&SphereToRect);
	SphereData particle1points = discretize(&identity, &sphere1, solution.particles[0].consts);
	f << "Particle 1 points: \n" << particle1points;


	f << "Particle 1 velocity: " << U1 << "\n";
	f << "Particle 1 angular velocity: " << Omega1 << "\n\n";

	f << " Particle 1 velocity data: \n";
	f << discretize(&solution, &sphere1, solution.particles[0].consts);

	ScalSphereData forcedensity(solution.particles[0].consts.NUMTRAPNODES,
		solution.particles[0].consts.NUMGLNODES);
	for (int i = 0; i < forcedensity.size(); i++)
		for (int j = 0; j < forcedensity[0].size(); j++)
			forcedensity[i][j] = norm(solution.particles[0].data[i][j]);
	f << "Particle 1 force data: \n";
	f << forcedensity;

	SphereData particle2points = discretize(&identity, &sphere2, solution.particles[1].consts);
	f << "Particle 2 points: \n" << particle2points;


	f << "Particle 2 velocity: " << U2 << "\n";
	f << "Particle 2 angular velocity: " << Omega2 << "\n\n";

	f << " Particle 2 velocity data: \n";
	f << discretize(&solution, &sphere2, solution.particles[1].consts);

	for (int i = 0; i < forcedensity.size(); i++)
		for (int j = 0; j < forcedensity[0].size(); j++)
			forcedensity[i][j] = norm(solution.particles[1].data[i][j]);
	f << "Particle 2 force data: \n";
	f << forcedensity;


	f.close();
}


void testFluid4()
{
	std::fstream f("testfluid4.txt");

	const double PI = 4.0 * atan(1.0);

	std::vector<RectCoord> forces;
	std::vector<RectCoord> torques;
	std::vector<SphereCoordSystem> spheres;

	RectCoord F1;
	forces.push_back(F1);
	RectCoord T1(8.0 * PI, 0.0, 0.0);
	torques.push_back(T1);

	SphereCoordSystem sphere1(RectCoord(2.1, 0.0, 0.0));
	spheres.push_back(sphere1);

	RectCoord F2;
	forces.push_back(F2);
	RectCoord T2(0.0, 0.0, 0.0);
	torques.push_back(T2);

	SphereCoordSystem sphere2;
	spheres.push_back(sphere2);

	StokesParticleSystem solution = SolveMobilityManySphere(spheres,
		forces,
		torques,
		16,
		StokesParticleSystem(),
		30,
		1e-5);

	RectCoord U1 = Integrate(&solution, &sphere1, solution.particles[0].consts) / (4.0 * PI);
	RectCoord U2 = Integrate(&solution, &sphere2, solution.particles[1].consts) / (4.0 * PI);


	SphereData flow1 = discretize(&solution, &sphere1, solution.particles[0].consts);
	RectCoord Omega1 = rotationCoefficient(flow1, sphere1, solution.particles[0].consts) * 3.0 / 8.0 / PI;

	SphereData flow2 = discretize(&solution, &sphere2, solution.particles[1].consts);
	RectCoord Omega2 = rotationCoefficient(flow2, sphere2, solution.particles[1].consts) * 3.0 / 8.0 / PI;

	SphericalVectorField identity(&SphereToRect);
	SphereData particle1points = discretize(&identity, &sphere1, solution.particles[0].consts);
	f << "Particle 1 points: \n" << particle1points;


	f << "Particle 1 velocity: " << U1 << "\n";
	f << "Particle 1 angular velocity: " << Omega1 << "\n\n";

	f << " Particle 1 velocity data: \n";
	f << discretize(&solution, &sphere1, solution.particles[0].consts);

	ScalSphereData forcedensity(solution.particles[0].consts.NUMTRAPNODES,
		solution.particles[0].consts.NUMGLNODES);
	for (int i = 0; i < forcedensity.size(); i++)
		for (int j = 0; j < forcedensity[0].size(); j++)
			forcedensity[i][j] = norm(solution.particles[0].data[i][j]);
	f << "Particle 1 force data: \n";
	f << forcedensity;

	SphereData particle2points = discretize(&identity, &sphere2, solution.particles[1].consts);
	f << "Particle 2 points: \n" << particle2points;


	f << "Particle 2 velocity: " << U2 << "\n";
	f << "Particle 2 angular velocity: " << Omega2 << "\n\n";

	f << " Particle 2 velocity data: \n";
	f << discretize(&solution, &sphere2, solution.particles[1].consts);

	for (int i = 0; i < forcedensity.size(); i++)
		for (int j = 0; j < forcedensity[0].size(); j++)
			forcedensity[i][j] = norm(solution.particles[1].data[i][j]);
	f << "Particle 2 force data: \n";
	f << forcedensity;


	f.close();
}

void testMobilityWilson3SpherePerp(double distance, double trueMUF, double trueMOmegaF, int numterms, std::fstream& f)
{
	const double PI = 4.0 * std::atan(1.0);

	RectCoord north(0.0, 1.0, 0.0);
	RectCoord southwest(std::cos(7.0 / 6.0 * PI), std::sin(7.0 / 6.0 * PI), 0.0); //points on an equilateral trangle.
	RectCoord southeast(std::cos(-PI / 6.0), std::sin(-PI / 6.0), 0.0);

	double normalizer = norm(north - southeast); // get the side length

	std::vector<SphereCoordSystem> centers;
	centers.push_back(distance / normalizer * north);
	centers.push_back(distance / normalizer * southwest); // create centers with given side length.
	centers.push_back(distance / normalizer * southeast);

	//centers.push_back(RectCoord(0.0, 0.0, 10.0));

	std::vector<RectCoord> forces;
	forces.push_back(RectCoord(0.0, 0.0, 6.0 * PI)); //force on first particle, pointing to the origin.
	forces.push_back(RectCoord(0.0, 0.0, 6.0 * PI));                    //force on second particle = 0.
	forces.push_back(RectCoord(0.0, 0.0, 6.0 * PI));                    //force on third particle = 0.
	//forces.push_back(RectCoord());

	std::vector<RectCoord> torques;
	torques.push_back(RectCoord());
	torques.push_back(RectCoord());                  //zero torque on all particles.
	torques.push_back(RectCoord());
	//torques.push_back(RectCoord());

	time_t time = std::time(0);
	time_t time0 = time;

	StokesParticleSystem ps(centers, 1);
	//system = std::move(SolveMobilityManySphere(centers, forces, torques, 1,system, 30,1e-5));
   //std::cout << "Time for first solve (N = 1): " << std::time(0) - time << "\n";

	double U3 = 1.0, U3prev = 0.0;

	QuadConstants consts(2 * numterms + 1, 2 * numterms + 1);

	ps = StokesParticleSystem(centers, numterms);
	time = std::time(0);
	ps = SolveMobilityManySphere(centers, forces, torques, numterms, ps, 30, 1e-5);
	f << "Time for solve (N = " << numterms << "): " << std::time(0) - time << "\n";

	RectCoord particle2velocity = Integrate(&ps, &ps.particles[1].system, consts) / (4.0 * PI);
	RectCoord particle1rotation = rotationCoefficient(discretize(&ps, &ps.particles[0].system, consts), ps.particles[0].system, consts) * 0.375 / PI;
	RectCoord particle2rotation = rotationCoefficient(discretize(&ps, &ps.particles[1].system, consts), ps.particles[1].system, consts) * 0.375 / PI;
	RectCoord particle3rotation = rotationCoefficient(discretize(&ps, &ps.particles[2].system, consts), ps.particles[2].system, consts) * 0.375 / PI;

	f << "MUF = " << std::setprecision(16) << Integrate(&ps, &ps.particles[0].system, consts) / (4.0 * PI) << "\n";
	f << "MUF error = " << (Integrate(&ps, &ps.particles[0].system, consts) / (4.0 * PI) - RectCoord(0.0, 0.0, trueMUF)) / trueMUF << "\n";
	f << "MOmegaF = " << norm(particle2rotation) << "\n";
	f << "MOmegaF error = " << (norm(particle2rotation) - trueMOmegaF) / trueMOmegaF << "\n";




	//std::cout << "Force on particle 1 after solving: " << Integrate(&system.traction, 1.0, system.particles[0].center) <<"\n";
	//std::cout << "Force on particle 2 after solving: " << Integrate(&system.traction,  1.0, system.particles[1].center) <<"\n";
	//std::cout << "Force on particle 3 after solving: " << Integrate(&system.traction, 1.0, system.particles[2].center) <<"\n";


}

std::array<double, 14> torquers = { 3.0 , 2.9 , 2.8 , 2.7 , 2.6 , 2.5 , 2.4 , 2.3, 2.25 , 2.2,2.15,2.1,2.05 , 2.01 };
std::array<double, 14> torqueMUTs = { 0.140187608, 0.149081378 ,  0.158630515 ,  0.168801210 , 0.179481457 ,  0.190411103 ,  0.201041187 ,  0.210224103 ,
								 0.213577673 , 0.215412007 ,  0.214756269 , 0.209731929 , 0.195478184 ,  0.159607490 };
std::array<double, 14> torqueMOmegaTs = { 0.771951396 , 0.772261325 , 0.771969753 ,  0.770692987 ,  0.767809702 ,  0.762285310 , 0.752325593 ,  0.734617762 ,
									 0.720980855 ,  0.702342997 , 0.676197489 , 0.637667933 , 0.573859459 ,  0.455596646 };


void testMobilityWilson3SpherePerpTorque(int index, int numterms, std::fstream& f)
{

	std::array<double, 14> rs = { 3.0 , 2.9 , 2.8 , 2.7 , 2.6 , 2.5 , 2.4 , 2.3, 2.25 , 2.2,2.15,2.1,2.05 , 2.01 };
	std::array<double, 14> MUTs = { 0.140187608, 0.149081378 ,  0.158630515 ,  0.168801210 , 0.179481457 ,  0.190411103 ,  0.201041187 ,  0.210224103 ,
									 0.213577673 , 0.215412007 ,  0.214756269 , 0.209731929 , 0.195478184 ,  0.159607490 };
	std::array<double, 14> MOmegaTs = { 0.771951396 , 0.772261325 , 0.771969753 ,  0.770692987 ,  0.767809702 ,  0.762285310 , 0.752325593 ,  0.734617762 ,
										 0.720980855 ,  0.702342997 , 0.676197489 , 0.637667933 , 0.573859459 ,  0.455596646 };
	const double PI = 4.0 * std::atan(1.0);

	RectCoord north(0.0, 1.0, 0.0);
	RectCoord southwest(std::cos(7.0 / 6.0 * PI), std::sin(7.0 / 6.0 * PI), 0.0); //points on an equilateral trangle.
	RectCoord southeast(std::cos(-PI / 6.0), std::sin(-PI / 6.0), 0.0);

	double normalizer = norm(north - southeast); // get the side length

	std::vector<SphereCoordSystem> centers;
	centers.push_back(rs[index] / normalizer * north);
	centers.push_back(rs[index] / normalizer * southwest); // create centers with given side length.
	centers.push_back(rs[index] / normalizer * southeast);

	//centers.push_back(RectCoord(0.0, 0.0, 10.0));

	std::vector<RectCoord> forces;
	forces.push_back(RectCoord()); //force on first particle, pointing to the origin.
	forces.push_back(RectCoord());                    //force on second particle = 0.
	forces.push_back(RectCoord());                    //force on third particle = 0.
	//forces.push_back(RectCoord());

	std::vector<RectCoord> torques;
	torques.push_back(RectCoord(-6.0 * PI, 0.0, 0.0));
	torques.push_back(RectCoord(-6.0 * PI * std::sin(7.0 / 6.0 * PI), 6.0 * PI * std::cos(7.0 / 6.0 * PI), 0.0));                  //zero torque on all particles.
	torques.push_back(RectCoord(-6.0 * PI * std::sin(-PI / 6.0), 6.0 * PI * std::cos(-PI / 6.0), 0.0));
	//torques.push_back(RectCoord());

	time_t time = std::time(0);
	time_t time0 = time;

	StokesParticleSystem ps(centers, 1);
	//system = std::move(SolveMobilityManySphere(centers, forces, torques, 1,system, 30,1e-5));
   //std::cout << "Time for first solve (N = 1): " << std::time(0) - time << "\n";

	double U3 = 1.0, U3prev = 0.0;

	QuadConstants consts(2 * numterms + 1, numterms + 1);

	ps = StokesParticleSystem(centers, numterms);
	time = std::time(0);
	ps = SolveMobilityManySphere(centers, forces, torques, numterms, ps, 40, 1e-5);
	f << "Time for solve (N = " << numterms << "): " << std::time(0) - time << "\n";

	RectCoord particle2velocity = Integrate(&ps, &ps.particles[1].system, consts) / (4.0 * PI);
	RectCoord particle1rotation = rotationCoefficient(discretize(&ps, &ps.particles[0].system, consts), ps.particles[0].system, consts) * 0.375 / PI;
	RectCoord particle2rotation = rotationCoefficient(discretize(&ps, &ps.particles[1].system, consts), ps.particles[1].system,  consts) * 0.375 / PI;
	RectCoord particle3rotation = rotationCoefficient(discretize(&ps, &ps.particles[2].system, consts), ps.particles[2].system,  consts) * 0.375 / PI;
	f << "Error in MUT = " << (Integrate(&ps, &ps.particles[0].system, consts) / (4.0 * PI) - RectCoord(0.0, 0.0, MUTs[index])) / MUTs[index] << "\n";
	f << "Error in MOmegaT = " << (norm(particle1rotation) - MOmegaTs[index]) /MOmegaTs[index] << "\n";




	//std::cout << "Force on particle 1 after solving: " << Integrate(&system.traction, 1.0, system.particles[0].center) <<"\n";
	//std::cout << "Force on particle 2 after solving: " << Integrate(&system.traction,  1.0, system.particles[1].center) <<"\n";
	//std::cout << "Force on particle 3 after solving: " << Integrate(&system.traction, 1.0, system.particles[2].center) <<"\n";


}

void testMobilityWilson3SphereInPlane(int index, int numterms, std::fstream& f)
{

	
	std::array<double, 6> rs = { 6.00 , 4.00 , 3.00 , 2.50 , 2.10 , 2.01 };
	std::array<double, 6> U1s = { 0.99581 , 0.97964 , 0.93905 , 0.87765 , 0.73857 , 0.65528 };
	std::array<double, 6> U2s = { 0.21586 , 0.31859 , 0.41694 , 0.49545 , 0.59718 , 0.63461 };
	std::array<double, 6> U3s = { 0.05078 , 0.06925 , 0.07824 , 0.07393 , 0.03517 , 0.00498 };
	std::array<double, 6> Omegas = { 0.010159 , 0.021634 , 0.035022 , 0.045466 , 0.052035 , 0.037336 };
	
	const double PI = 4.0 * std::atan(1.0);

	RectCoord north(0.0, 1.0, 0.0);
	RectCoord southwest(std::cos(7.0 / 6.0 * PI), std::sin(7.0 / 6.0 * PI), 0.0); //points on an equilateral trangle.
	RectCoord southeast(std::cos(-PI / 6.0), std::sin(-PI / 6.0), 0.0);

	double normalizer = norm(north - southeast); // get the side length

	std::vector<SphereCoordSystem> centers;
	centers.push_back(rs[index] / normalizer * north);
	centers.push_back(rs[index] / normalizer * southwest); // create centers with given side length.
	centers.push_back(rs[index] / normalizer * southeast);

	//centers.push_back(RectCoord(0.0, 0.0, 10.0));

	std::vector<RectCoord> forces;
	forces.push_back(RectCoord(0.0, -6.0 * PI, 0.0)); //force on first particle, pointing to the origin.
	forces.push_back(RectCoord());                    //force on second particle = 0.
	forces.push_back(RectCoord());                    //force on third particle = 0.
	//forces.push_back(RectCoord());

	std::vector<RectCoord> torques;
	torques.push_back(RectCoord());
	torques.push_back(RectCoord());                  //zero torque on all particles.
	torques.push_back(RectCoord());
	//torques.push_back(RectCoord());

	time_t time = std::time(0);
	time_t time0 = time;

	StokesParticleSystem ps(centers, 1);
	//system = std::move(SolveMobilityManySphere(centers, forces, torques, 1,system, 30,1e-5));
   //std::cout << "Time for first solve (N = 1): " << std::time(0) - time << "\n";

	double U3 = 1.0, U3prev = 0.0;

	QuadConstants consts(2 * numterms + 1, numterms + 1);

	ps = StokesParticleSystem(centers, numterms);
	time = std::time(0);
	ps = SolveMobilityManySphere(centers, forces, torques, numterms, ps, 40, 1e-5);
	f << "s = " << rs[index];
	f << "Time for solve (N = " << numterms << "): " << std::time(0) - time << "\n";

	RectCoord particle1velocity = Integrate(&ps, &ps.particles[0].system, consts) / (4.0 * PI);
	RectCoord particle2velocity = Integrate(&ps, &ps.particles[1].system, consts) / (4.0 * PI);
	RectCoord particle1rotation = rotationCoefficient(discretize(&ps, &ps.particles[0].system, consts), ps.particles[0].system, consts) * 0.375 / PI;
	RectCoord particle2rotation = rotationCoefficient(discretize(&ps, &ps.particles[1].system, consts), ps.particles[1].system, consts) * 0.375 / PI;
	RectCoord particle3rotation = rotationCoefficient(discretize(&ps, &ps.particles[2].system, consts), ps.particles[2].system, consts) * 0.375 / PI;

	f << "U1 = " << -particle1velocity.y << "\n";
	f << "U1 error = " << std::abs(-particle1velocity.y - U1s[index]) / U1s[index] << "\n";
	f << "U2 = " << -particle2velocity.y << "\n";
	f << "U2 error = " << std::abs(-particle2velocity.y - U2s[index]) / U2s[index] << "\n";
	f << "U3 = " << -particle2velocity.x << "\n";
	f << "U3 error = " << std::abs(-particle2velocity.x - U3s[index]) / U3s[index] << "\n";
	f << "Omega = " << norm(particle2rotation) << "\n";
	f << "Omega error = " << std::abs(norm(particle2rotation) - Omegas[index]) / Omegas[index] << "\n";



	//std::cout << "Force on particle 1 after solving: " << Integrate(&system.traction, 1.0, system.particles[0].center) <<"\n";
	//std::cout << "Force on particle 2 after solving: " << Integrate(&system.traction,  1.0, system.particles[1].center) <<"\n";
	//std::cout << "Force on particle 3 after solving: " << Integrate(&system.traction, 1.0, system.particles[2].center) <<"\n";


}

void testMobilityWilson3SphereLinear(int index, int numterms, std::fstream& f)
{
	const double PI = 4.0 * std::atan(1.0);

	std::array<double, 5> rs = { 4.0 , 3.0 , 2.4 , 2.1 , 2.05 };
	std::array<double, 5> D0s = { 0.7194 , 0.6522 , 0.5946 , 0.5610 , 0.5561 };
	std::array<double, 5> D1s = { 0.7755 , 0.7177 , 0.6663, 0.6337 , 0.6277 };
	std::array<double, 5> Omegas = { 0.0586 , 0.1036, 0.1578 , 0.1919 , 0.1934 };

	RectCoord middle(0.0, 0.0, 0.0);
	RectCoord right(rs[index], 0.0, 0.0); //points on a line.
	RectCoord left(-rs[index], 0.0, 0.0);


	std::vector<SphereCoordSystem> centers;
	centers.push_back(middle);
	centers.push_back(right); // create centers with given side length.
	centers.push_back(left);

	//centers.push_back(RectCoord(0.0, 0.0, 10.0));

	std::vector<RectCoord> forces;
	forces.push_back(RectCoord(0.0, 6.0 * PI, 0.0)); //force on first particle, pointing to the origin.
	forces.push_back(RectCoord(0.0, 6.0 * PI, 0.0));                    //force on second particle = 0.
	forces.push_back(RectCoord(0.0, 6.0 * PI, 0.0));                    //force on third particle = 0.
	//forces.push_back(RectCoord());

	std::vector<RectCoord> torques;
	torques.push_back(RectCoord());
	torques.push_back(RectCoord());                  //zero torque on all particles.
	torques.push_back(RectCoord());
	//torques.push_back(RectCoord());

	time_t time = std::time(0);
	time_t time0 = time;

	StokesParticleSystem ps(centers, 1);
	//system = std::move(SolveMobilityManySphere(centers, forces, torques, 1,system, 30,1e-5));
   //std::cout << "Time for first solve (N = 1): " << std::time(0) - time << "\n";

	double U3 = 1.0, U3prev = 0.0;

	QuadConstants consts(2 * numterms + 1, numterms + 1);

	ps = StokesParticleSystem(centers, numterms);
	time = std::time(0);
	ps = SolveMobilityManySphere(centers, forces, torques, numterms, ps, 30, 1e-5);
	f << "Time for solve (p = " << numterms << " , s = " << rs[index] << "): " << std::time(0) - time << "\n";

	RectCoord particle1velocity = Integrate(&ps, &ps.particles[0].system, consts) / (4.0 * PI);
	RectCoord particle2velocity = Integrate(&ps, &ps.particles[1].system, consts) / (4.0 * PI);
	RectCoord particle2rotation = rotationCoefficient(discretize(&ps, &ps.particles[1].system, consts), ps.particles[1].system, consts) * 0.375 / PI;

	f << "D0 = " << particle1velocity.y << "\n";
	f << "D0 error = " << std::abs(1.0 / particle1velocity.y - D0s[index]) / D0s[index] << "\n";
	f << "D1 = " << particle2velocity.y << "\n";
	f << "D1 error = " << std::abs(1.0 / particle2velocity.y - D1s[index]) / D1s[index] << "\n";
	f << "Omega = " << norm(particle2rotation) << "\n";
	f << "Omega error = " << std::abs(norm(particle2rotation) - Omegas[index]) / Omegas[index] << "\n";

	//std::cout << "Force on particle 1 after solving: " << Integrate(&system.traction, 1.0, system.particles[0].center) <<"\n";
	//std::cout << "Force on particle 2 after solving: " << Integrate(&system.traction,  1.0, system.particles[1].center) <<"\n";
	//std::cout << "Force on particle 3 after solving: " << Integrate(&system.traction, 1.0, system.particles[2].center) <<"\n";


}