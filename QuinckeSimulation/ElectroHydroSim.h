#pragma once
#include "C:\Users\David\source\repos\ElectrostaticsQuincke\ElectrostaticsQuincke\LaplaceManyParticles.h"
#include "C:\Users\David\source\repos\IrinabanMaz\SphericalHarmonics\MobilityProblem.h"

class ElectroHydroSim
{
public:

	LaplaceParticleSystem electricSolution;
	LaplaceParticleSystem q;
	LaplaceParticleSystem OhmicCurrent;
	LaplaceParticleSystem forceDensity;

	double fluidConductivity;
	std::vector<double> particleConductivities;
	double fluidPermitivity;
	std::vector<double> particlePermitivities;
	RectCoord appliedField;


	std::vector<RectCoord> forces;
	std::vector<RectCoord> torques;


	int numSeriesTerms;

	std::vector<SphereCoordSystem> spheres;

	ElectroHydroSim(RectCoord E0, double epsplus, std::vector<double> epsminus, double sigmaplus, std::vector<double>  sigmaminus, int nst) :
		fluidPermitivity(epsplus),
		particlePermitivities(epsminus),
		fluidConductivity(sigmaplus),
		particleConductivities(sigmaminus),
		appliedField(E0),
		numSeriesTerms(nst) {}

	void initialize(std::vector<SphereCoordSystem> cs)
	{
		spheres = cs;
		q = LaplaceParticleSystem(spheres, numSeriesTerms);
		electricSolution = q;
	}

	void solveElectric()
	{
		int size = spheres.size();


		
		LPSQuinckeOperator L(fluidPermitivity);

		LaplaceParticleSystem appliedFieldTerm(spheres, numSeriesTerms);
		appliedFieldTerm.setPermitivities(particlePermitivities);

		LaplaceParticleSystem rhs;

		for (int i = 0; i < appliedFieldTerm.particles.size(); i++)
		{
			double lambda = (particlePermitivities[i] - fluidPermitivity) / (particlePermitivities[i] + fluidPermitivity);
			NormalField field(appliedField);
			ScalSphereData data = discretize(&field,
				&appliedFieldTerm.particles[i].system,
				appliedFieldTerm.particles[i].consts);

			appliedFieldTerm.particles[i] = LaplaceParticle(data, numSeriesTerms, appliedFieldTerm.particles[i].system) * lambda;
			
			rhs.push_back(q.particles[i] * (particlePermitivities[i] + fluidPermitivity) 
				          + appliedFieldTerm.particles[i]);

			rhs.particles[i].setPermitivity(particlePermitivities[i]);
		}

		std::cout << "Norm of rhs in GMRES solve: " << rhs.norm() << "\n";
		LPSIdentityPreconditioner I;

		time_t t0 = time(0);
		GMRES(&L, &electricSolution, &rhs, &I, 30, 1, 1e-5);
		std::cout << "Time to solve: " << time(0) - t0 << std::endl;

		//std::cout << Integrate(&BIEsolution, 1.0, BIEsolution.particles[0].center) * 0.25 / PI << std::endl;
		//std::cout << Integrate(&rh, 1.0, rh.particles[0].center) * 0.25 / PI << std::endl;

	}

	void updateOhmicCurrent()
	{

		std::vector<double> permjumpinv;
		permjumpinv.resize(particlePermitivities.size());
		for (int i = 0; i < electricSolution.particles.size(); i++)
			permjumpinv[i] = 1.0 / (fluidPermitivity - particlePermitivities[i]);
		LaplaceParticleSystem EnPlus = (q - electricSolution * particlePermitivities)
			* permjumpinv;

		LaplaceParticleSystem EnMinus = (q - electricSolution * fluidPermitivity)
			* permjumpinv;



	    OhmicCurrent = EnPlus * fluidConductivity - EnMinus * particleConductivities;
	}

	void timeStepOhmicForward(double dt)
	{
		updateOhmicCurrent();

		double timescale = particlePermitivities[0] / particleConductivities[0];

		q += OhmicCurrent * (-timescale * dt);
	}

	void timeStepOhmicBackwardFP(double dt , double tol)
	{
		
		updateOhmicCurrent();

		double timescale = particlePermitivities[0] / particleConductivities[0];

		LaplaceParticleSystem qstart = q;
		LaplaceParticleSystem qlast = q + OhmicCurrent * (-timescale * dt);

		int i = 1;
		double err = (qlast - q).norm();
		std::cout << "Using implicit Euler with fixed point iteration: \n";
			while (err > tol)
			{
				qlast = q;
				solveElectric();

				updateOhmicCurrent();

				q = qstart + OhmicCurrent * (-timescale * dt);

				err = OhmicCurrent.norm();

				std::cout << "Error in fixed point method: " << err << "\n\n";
			}

	
	}


};
