#pragma once
#include "ElectrostaticsQuincke/ElectrostaticsQuincke/LaplaceManyParticles.h"
#include "StokesSolver/MobilityProblem.h"


LaplaceParticleSystem operator *(const StokesParticleSystem& Ps, const StokesParticleSystem& Qs)
{
	LaplaceParticleSystem temp(Ps.getN());

	for (int i = 0; i < temp.particles.size(); i++)
	{
		ScalSphereData vals = Ps.particles[i].data * Qs.particles[i].data;

		temp.particles[i] = LaplaceParticle(vals,
							Ps.particles[i].numSeriesTerms,
							Ps.particles[i].system);
	}
	return temp;
}

StokesParticleSystem operator *(const LaplaceParticleSystem& Ps, const StokesParticleSystem& Qs)
{
	StokesParticleSystem temp(Ps.getN());

	for (int i = 0; i < temp.particles.size(); i++)
	{
		SphereData vals = Ps.particles[i].data * Qs.particles[i].data;

		temp.particles[i] = StokesParticle(vals,
			Ps.particles[i].system,
			Ps.particles[i].numSeriesTerms);
	}
	return temp;
}


class ConstantField : public SphericalVectorField {
private:
	RectCoord direction;
public:
	ConstantField(const RectCoord& d) : direction(d)
	{}

	RectCoord operator()(const SphereCoord& x) const
	{
		return direction;
	}
};

class ElectroHydroSim
{
public:

	LaplaceParticleSystem electricSolution;
	LaplaceParticleSystem q;
	
	LaplaceParticleSystem OhmicCurrent;

	double fluidConductivity;
	std::vector<double> particleConductivities;
	double fluidPermitivity;
	std::vector<double> particlePermitivities;
	RectCoord appliedField;


	std::vector<RectCoord> forces;
	std::vector<RectCoord> torques;


	StokesParticleSystem electricTraction;
	StokesParticleSystem fluidSolution;
	LaplaceParticleSystem flowAdvection;

	std::vector<RectCoord> Us;
	std::vector<RectCoord> Omegas;

	int numSeriesTerms;

	std::vector<SphereCoordSystem> spheres;

	ElectroHydroSim(RectCoord E0, double epsplus, std::vector<double> epsminus, double sigmaplus, std::vector<double>  sigmaminus, int nst) :
		fluidPermitivity(epsplus),
		particlePermitivities(epsminus),
		fluidConductivity(sigmaplus),
		particleConductivities(sigmaminus),
		appliedField(E0),
		numSeriesTerms(nst){}

	void initialize(std::vector<SphereCoordSystem> cs)
	{
		spheres = cs;
		Us.resize(cs.size());
		Omegas.resize(cs.size());
		q = LaplaceParticleSystem(spheres, numSeriesTerms);
		flowAdvection = LaplaceParticleSystem(spheres, numSeriesTerms);
		electricTraction = StokesParticleSystem(spheres, numSeriesTerms);
		fluidSolution = StokesParticleSystem(spheres, numSeriesTerms);
		q.setPermitivities(particlePermitivities);
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

	void updateAdvection()
	{
		double PI = 4.0 * atan(1.0);
		for (int i = 0; i < flowAdvection.getN(); i++)
		{
			

			ForceBalance rigidFlow(4.0 * PI * Us[i], PI / 0.375 * Omegas[i]);

			SphereData rigidFlowDiscrete = discretize(&rigidFlow,
				&spheres[i],
				fluidSolution.particles[i].consts);
			/*
			SphereData flow = discretize(&fluidSolution,
				&spheres[i],
				fluidSolution.particles[i].consts);
				*/

			SurfaceGrad gradq(&q.particles[i].fourierdata);

			SphereData gradqDiscrete = discretize(&gradq,
				&spheres[i],
				q.particles[i].consts);

			ScalSphereData advectionData = rigidFlowDiscrete * gradqDiscrete;

			flowAdvection.particles[i] = LaplaceParticle(advectionData, numSeriesTerms , spheres[i]);
		}
	}

	void timeStepOhmicForward(double dt)
	{
		updateOhmicCurrent();

		double timescale = particlePermitivities[0] / particleConductivities[0];

		q += OhmicCurrent * (-timescale * dt);
	}

	void timeStepOhmicAB(double dt)
	{
		static LaplaceParticleSystem om2 = LaplaceParticleSystem(spheres , numSeriesTerms);

		updateOhmicCurrent();
		if (numSeriesTerms != om2.particles[0].numSeriesTerms)
			om2 = OhmicCurrent;
		double timescale = particlePermitivities[0] / particleConductivities[0];
		
		q += -OhmicCurrent * 1.5 * dt * timescale + om2 * 0.5 * timescale * dt;
		om2 = OhmicCurrent;
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

	void timeStepForward(double dt)
	{
		updateOhmicCurrent();
		updateAdvection();

		double timescale = particlePermitivities[0] / particleConductivities[0];

		q += -(OhmicCurrent + flowAdvection) * timescale * dt;

		for (int i = 0; i < spheres.size(); i++)
		{
			spheres[i] = SphereCoordSystem(spheres[i].center + Us[i] * timescale * dt);
		}
	}

	void computeElectricForces()
	{
		StokesParticleSystem E(spheres, numSeriesTerms);
		StokesParticleSystem n(spheres, numSeriesTerms);
		StokesParticleSystem E0(spheres , numSeriesTerms);
		SphericalVectorField er(&e_r);

		for (int i = 0; i < E.particles.size(); i++)
		{
			time_t t0 = time(0);
			E.particles[i] = StokesParticle(discretize(&electricSolution.E, 
				                        &spheres[i],
				                        electricSolution.particles[i].consts),
				                        spheres[i] , numSeriesTerms);
			
			n.particles[i] = StokesParticle(discretize(&er,
											&spheres[i],
											electricSolution.particles[i].consts),
											spheres[i],
											numSeriesTerms);

			ConstantField e0(appliedField);

			E0.particles[i] = StokesParticle(discretize(&e0,
				&spheres[i],
				electricSolution.particles[i].consts),
				spheres[i],
				numSeriesTerms);


			std::cout << "Field on particle " << i + 1 << " computed. Time elapsed: " << time(0) - t0 << "\n";
		}

		StokesParticleSystem externalTraction = ((E0 * n) * E0 - (E0 * E0) * n * 0.5) * fluidPermitivity;
		electricTraction = ((E *n)* E  - (E * E) * n * 0.5) * fluidPermitivity + externalTraction;

		for (int i = 0; i < electricTraction.particles.size(); i++)
		{
			RectCoord F = IntegrateDiscrete(electricTraction.particles[i].data,
				&electricTraction.particles[i].system,
				electricTraction.particles[i].consts);
			RectCoord T = rotationCoefficient(electricTraction.particles[i].data,
				electricTraction.particles[i].system,
				electricTraction.particles[i].consts);

			std::cout << "Net Force on particle " << i + 1 << ": " << F << "\n";
			std::cout << "Net Torque on particle " << i + 1 << ": " << T << "\n";

		}
	}

	void solveFluid()
	{
		PSLHSOperator L;


		const double PI = 4.0 * atan(1.0);
		

		fluidSolution = fluidSolution -electricTraction;
		PSRHSOperator R;
		StokesParticleSystem rhs = - (R * electricTraction);

		SPSIdentityPreconditioner I;

		GMRES(&L, &fluidSolution, &rhs, &I, 30, 1, 1e-5);



		StokesParticleSystem soln = fluidSolution + electricTraction;

		for (int i = 0; i < fluidSolution.particles.size(); i++)
		{
			SphereData flowdata = discretize(&fluidSolution,
				&spheres[i],
				fluidSolution.particles[i].consts);

			 Us[i] = IntegrateDiscrete(flowdata,
				&spheres[i],
				fluidSolution.particles[i].consts) / 4.0 * PI;

			 Omegas[i] = rotationCoefficient(flowdata,
				spheres[i],
				fluidSolution.particles[i].consts) * 0.375 / PI;
		}
	}
	

};
