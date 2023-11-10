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

	RectCoord F2(6.0 * PI , 0.0 , 0.0);
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

	ScalSphereData forcedensity(33, 33);
	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
			forcedensity[i][j] = norm(solution.particles[0].data[i][j]);
	f << "Particle 1 force data: \n";
	f << forcedensity;

	SphereData particle2points = discretize(&identity, &sphere2, solution.particles[1].consts);
	f << "Particle 2 points: \n" << particle2points;


	f << "Particle 2 velocity: " << U2 << "\n";
	f << "Particle 2 angular velocity: " << Omega2 << "\n\n";

	f << " Particle 2 velocity data: \n";
	f << discretize(&solution, &sphere2, solution.particles[1].consts);

	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
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

	SphereCoordSystem sphere1(RectCoord(-2.1, 0.0, 0.0));
	spheres.push_back(sphere1);

	RectCoord F2(0.0 , -6.0 * PI , 0.0);
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

	ScalSphereData forcedensity(33, 33);
	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
			forcedensity[i][j] = norm(solution.particles[0].data[i][j]);
	f << "Particle 1 force data: \n";
	f << forcedensity;

	SphereData particle2points = discretize(&identity, &sphere2, solution.particles[1].consts);
	f << "Particle 1 points: \n" << particle2points;


	f << "Particle 2 velocity: " << U2 << "\n";
	f << "Particle 2 angular velocity: " << Omega2 << "\n\n";

	f << " Particle 2 velocity data: \n";
	f << discretize(&solution, &sphere2, solution.particles[1].consts);

	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
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

	SphereCoordSystem sphere1(RectCoord(-2.1, 0.0, 0.0));
	spheres.push_back(sphere1);

	RectCoord F2;
	forces.push_back(F2);
	RectCoord T2(0.0, 8.0 * PI, 0.0);
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

	ScalSphereData forcedensity(33, 33);
	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
			forcedensity[i][j] = norm(solution.particles[0].data[i][j]);
	f << "Particle 1 force data: \n";
	f << forcedensity;

	SphereData particle2points = discretize(&identity, &sphere2, solution.particles[1].consts);
	f << "Particle 2 points: \n" << particle2points;


	f << "Particle 2 velocity: " << U2 << "\n";
	f << "Particle 2 angular velocity: " << Omega2 << "\n\n";

	f << " Particle 2 velocity data: \n";
	f << discretize(&solution, &sphere2, solution.particles[1].consts);

	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
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
	RectCoord T1(0.0 ,8.0 * PI, 0.0);
	torques.push_back(T1);

	SphereCoordSystem sphere1(RectCoord(-2.1, 0.0, 0.0), RectCoord(1.0, 0.0, 0.0));
	spheres.push_back(sphere1);

	RectCoord F2;
	forces.push_back(F2);
	RectCoord T2(0.0, -8.0 * PI, 0.0);
	torques.push_back(T2);

	SphereCoordSystem sphere2(RectCoord(), RectCoord(-1.0, 0.0, 0.0));
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


	f << "Particle 1 velocity: " << U1 << "\n";
	f << "Particle 1 angular velocity: " << Omega1 << "\n\n";

	SphericalVectorField identity(&SphereToRect);
	SphereData particle1points = discretize(&identity, &sphere1, solution.particles[0].consts);
	f << "Particle 1 points: \n" << particle1points;


	f << "Particle 1 velocity: " << U1 << "\n";
	f << "Particle 1 angular velocity: " << Omega1 << "\n\n";

	f << " Particle 1 velocity data: \n";
	f << discretize(&solution, &sphere1, solution.particles[0].consts);

	ScalSphereData forcedensity(33, 33);
	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
			forcedensity[i][j] = norm(solution.particles[0].data[i][j]);
	f << "Particle 1 force data: \n";
	f << forcedensity;

	SphereData particle2points = discretize(&identity, &sphere2, solution.particles[1].consts);
	f << "Particle 2 points: \n" << particle2points;


	f << "Particle 2 velocity: " << U2 << "\n";
	f << "Particle 2 angular velocity: " << Omega2 << "\n\n";

	f << " Particle 2 velocity data: \n";
	f << discretize(&solution, &sphere2, solution.particles[1].consts);

	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
			forcedensity[i][j] = norm(solution.particles[1].data[i][j]);
	f << "Particle 2 force data: \n";
	f << forcedensity;


	f.close();
}

void testFluid5()
{
	std::fstream f("testfluid5.txt");

	const double PI = 4.0 * atan(1.0);

	std::vector<RectCoord> forces;
	std::vector<RectCoord> torques;
	std::vector<SphereCoordSystem> spheres;

	RectCoord F1;
	forces.push_back(F1);
	RectCoord T1(8.0 * PI, 0.0, 0.0);
	torques.push_back(T1);

	SphereCoordSystem sphere1(RectCoord(-2.1, 0.0, 0.0), RectCoord(0.0, 0.0, -1.0));
	spheres.push_back(sphere1);

	RectCoord F2;
	forces.push_back(F2);
	RectCoord T2(0.0, 8.0 * PI, 0.0);
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

	ScalSphereData forcedensity(33, 33);
	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
			forcedensity[i][j] = norm(solution.particles[0].data[i][j]);
	f << "Particle 1 force data: \n";
	f << forcedensity;

	SphereData particle2points = discretize(&identity, &sphere2, solution.particles[1].consts);
	f << "Particle 2 points: \n" << particle2points;


	f << "Particle 2 velocity: " << U2 << "\n";
	f << "Particle 2 angular velocity: " << Omega2 << "\n\n";

	f << " Particle 2 velocity data: \n";
	f << discretize(&solution, &sphere2, solution.particles[1].consts);

	for (int i = 0; i < 33; i++)
		for (int j = 0; j < 33; j++)
			forcedensity[i][j] = norm(solution.particles[1].data[i][j]);
	f << "Particle 2 force data: \n";
	f << forcedensity;

	f.close();
}
