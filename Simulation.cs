using Resonance;
using System.Diagnostics;
using System.Numerics;
using System.Xml.Linq;


namespace Resonance
{
    public abstract class Simulation
    {
        protected Star MainStar;
        protected Planet[] Planets;
        protected Particle[] Particles;

        public Simulation(Star mainStar, Planet[] planets, int particlesCount)
        {
            MainStar = mainStar;
            Planets = planets;
            CreateParticles(particlesCount);
        }

        public Simulation(Star mainStar, Planet planet, int particlesCount)
        {
            MainStar = mainStar;
            Planets = new Planet[1] { planet };
            CreateParticles(particlesCount);
        }

        public Simulation(Star mainStar, int particlesCount)
        {
            MainStar = mainStar;
            Planets = new Planet[0];
            CreateParticles(particlesCount);
        }

        public int PlanetsCount { get { return Planets.Length; } }

        public int ParticlesCount { get { return Particles.Length; } }

        public void CreateParticles(int count)
        {
            if (count == 0)
            {
                Particles = new Particle[0];
                return;
            }

            BasicKeplerianNormal normal;
            switch (count)
            {
                case 0:
                    normal = new(new(Constants.AU, 0, 0), new(Constants.AU / 2, 1, Math.PI / 2));
                    break;
                case 1:
                    BasicKeplerian planet = Planets[0].ToBasicKeplerian(MainStar);
                    normal = new(new(planet.SemiMajorAxis, 0, 0),
                                 new(planet.SemiMajorAxis / 2, 1, Math.PI / 2));
                    break;
                default:
                    normal = PlanetsKeplerianDistribution();
                    break;
            }

            Random random = new();

            Particles = new Particle[count];
            for (int i = 0; i < count; i++)
                Particles[i] = Particle.FromKeplerian(MainStar, random.NextKeplerian(normal, MainStar.Mass));
        }

        private BasicKeplerianNormal PlanetsKeplerianDistribution()
        {
            double length = PlanetsCount;
            BasicKeplerian planetKeplerian;

            BasicKeplerianNormal normal = new();

            foreach (Planet planet in Planets)
            {
                planetKeplerian = planet.ToBasicKeplerian(MainStar);

                normal.Mu.SemiMajorAxis += planetKeplerian.SemiMajorAxis;
                normal.Mu.Eccentricity += planetKeplerian.Eccentricity;
                normal.Mu.Inclination += planetKeplerian.Inclination;
            }

            normal.Mu.SemiMajorAxis /= length;
            normal.Mu.Eccentricity /= length;
            normal.Mu.Inclination /= length;

            foreach (Planet planet in Planets)
            {
                planetKeplerian = planet.ToBasicKeplerian(MainStar);

                normal.Sigma.SemiMajorAxis += Math.Pow(planetKeplerian.SemiMajorAxis - normal.Mu.SemiMajorAxis, 2);
                normal.Sigma.Eccentricity += Math.Pow(planetKeplerian.Eccentricity - normal.Mu.Eccentricity, 2);
                normal.Sigma.Inclination += Math.Pow(planetKeplerian.Inclination - normal.Mu.Inclination, 2);
            }

            normal.Sigma.SemiMajorAxis = Math.Sqrt(normal.Sigma.SemiMajorAxis);
            normal.Sigma.Eccentricity = Math.Sqrt(normal.Sigma.Eccentricity);
            normal.Sigma.Inclination = Math.Sqrt(normal.Sigma.Inclination);

            return normal;
        }

        private (Vector3D star, Vector3D[] planets) AccelerationStarAndPlanets()
        {
            Vector3D mainStarAcc = Vector3D.Zero;
            Vector3D[] planetsAcc = new Vector3D[PlanetsCount];

            for (int i = 0; i < PlanetsCount; i++)
            {
                Vector3D delta = Planets[i].Pos - MainStar.Pos;
                Vector3D aoverm = Constants.G / Math.Pow(delta.Magnitude, 3) * delta;
                mainStarAcc += Planets[i].Mass * aoverm;
                planetsAcc[i] = -MainStar.Mass * aoverm;
            }

            for (int i = 0; i < PlanetsCount; i++)
            {
                for (int j = i + 1; j < PlanetsCount; j++)
                {
                    Vector3D delta = Planets[j].Pos - Planets[i].Pos;
                    Vector3D aoverm = Constants.G / Math.Pow(delta.Magnitude, 3) * delta;
                    planetsAcc[i] += Planets[j].Mass * aoverm;
                    planetsAcc[j] -= Planets[i].Mass * aoverm;
                }
            }

            return (mainStarAcc, planetsAcc);
        }

        private Vector3D[] AccelerationParticles()
        {
            Vector3D[] accelerations = new Vector3D[ParticlesCount];

            for (int i = 0; i < ParticlesCount; i++)
            {
                Vector3D delta = MainStar.Pos - Particles[i].Pos;
                accelerations[i] = MainStar.Mass * Constants.G / Math.Pow(delta.Magnitude, 3) * delta;
                foreach (Planet planet in Planets)
                {
                    delta = planet.Pos - Particles[i].Pos;
                    accelerations[i] += planet.Mass * Constants.G / Math.Pow(delta.Magnitude, 3) * delta;
                }
            }

            return accelerations;
        }

        private void UpdatePositions(double dt, Vector3D mainStarAcc, Vector3D[] planetsAcc, Vector3D[] particlesAcc)
        {
            MainStar.Pos += dt * (MainStar.Vel + dt * mainStarAcc / 2);
            for (int i = 0; i < PlanetsCount; i++)
                Planets[i].Pos += dt * (Planets[i].Vel + dt * planetsAcc[i] / 2);
            for (int i = 0; i < ParticlesCount; i++)
                Particles[i].Pos += dt * (Particles[i].Vel + dt * particlesAcc[i] / 2);
        }

        private void UpdateVelocities(double dt, Vector3D mainStarAcc, Vector3D[] planetsAcc, Vector3D[] particlesAcc,
                                                 Vector3D newMainStarAcc, Vector3D[] newPlanetsAcc, Vector3D[] newParticlesAcc)
        {
            MainStar.Vel += dt * (mainStarAcc + newMainStarAcc) / 2;
            for (int i = 0; i < PlanetsCount; i++)
                Planets[i].Vel += dt * (planetsAcc[i] + newPlanetsAcc[i]) / 2;
            for (int i = 0; i < ParticlesCount; i++)
                Particles[i].Vel += dt * (particlesAcc[i] + newParticlesAcc[i]) / 2;
        }

        public void NextStep(double dt)  // https://gamedev.stackexchange.com/questions/15708/how-can-i-implement-gravity
        {
            var (mainStarAcc, planetsAcc) = AccelerationStarAndPlanets();
            Vector3D[] particlesAcc = AccelerationParticles();

            UpdatePositions(dt, mainStarAcc, planetsAcc, particlesAcc);

            var (newMainStarAcc, newPlanetsAcc) = AccelerationStarAndPlanets();
            Vector3D[] newParticlesAcc = AccelerationParticles();

            UpdateVelocities(dt, mainStarAcc, planetsAcc, particlesAcc, newMainStarAcc, newPlanetsAcc, newParticlesAcc);

            MainStar.Pos += MainStar.Vel * dt;
            MainStar.Vel += mainStarAcc * dt;
        }

        public abstract void Run(int steps, double dt, int substeps);
    }

    public class KeplerianSimulation : Simulation
    {
        public KeplerianSimulation(Star mainStar, int particlesCount = 0) : base(mainStar, particlesCount) { }

        public KeplerianSimulation(Star mainStar, Planet planet, int particlesCount = 0) : base(mainStar, planet, particlesCount) { }
        
        public KeplerianSimulation(Star mainStar, Planet[] planets, int particlesCount = 0) : base(mainStar, planets, particlesCount) { }

        public override void Run(int steps = 1000, double dt = 86400, int substeps = 1)
        {
            Stopwatch timer = new();
            timer.Start();

            for (int step = 0; step < steps; step++)
            {
                for (int i = 0; i < substeps; i++)
                    NextStep(dt);

                for (int i = 0; i < PlanetsCount; i++)
                {
                    BasicKeplerian keplerian = Keplerian.BasicFromCartesian(MainStar.Mass, Planets[i].Pos - MainStar.Pos, Planets[i].Vel - MainStar.Vel);
                }

                for (int i = 0; i < ParticlesCount; i++)
                {
                    BasicKeplerian keplerian = Keplerian.BasicFromCartesian(MainStar.Mass, Particles[i].Pos - MainStar.Pos, Particles[i].Vel - MainStar.Vel);
                }

                double remainingSeconds = ((double)steps / (step + 1) - 1) * timer.ElapsedMilliseconds / 1000;

                Console.Write($"\r{step + 1}/{steps}, " + Math.Round(remainingSeconds, 2).ToString() + " seconds remaining");
            }
              
            Console.WriteLine($"\r{steps}/{steps}, finished in " + Math.Round((double)timer.ElapsedMilliseconds / 1000, 2).ToString() + " seconds");
        }
    }
}
