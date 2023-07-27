namespace Resonance
{
    public class Body3D
    {
        public Vector3D Pos;
        public Vector3D Vel;

        public Body3D()
        {
            Pos = Vector3D.Zero;
            Vel = Vector3D.Zero;
        }

        public Body3D(Vector3D pos, Vector3D vel)
        {
            Pos = pos;
            Vel = vel;
        }
    }

    public class Star : Body3D
    {
        public double Mass;

        public Star(double mass) : base()
        {
            Mass = mass;
        }
    }

    public class Particle : Body3D
    {
        public Particle(Vector3D pos, Vector3D vel) : base(pos, vel) { }

        private double EccentricAnormalyNumeric(double meanArnormaly, double eccentricity, int iterations)
        {
            double eccentricAnormaly = meanArnormaly;
            for (int i = 0; i < iterations; i++)
            {
                eccentricAnormaly = eccentricAnormaly - (eccentricAnormaly - eccentricity * Math.Sin(eccentricAnormaly) - meanArnormaly) / (1 - eccentricity * Math.Cos(eccentricAnormaly));
            }
            return eccentricAnormaly;
        }

        public static Particle FromKeplerian(Star mainStar, Keplerian keplerian)
        {
            Body3D cartesian = keplerian.ToCartesian();
            return new(mainStar.Pos + cartesian.Pos, mainStar.Vel + cartesian.Vel);
        }

        public BasicKeplerian ToBasicKeplerian(Star mainStar)
        {
            double mu = Constants.G * mainStar.Mass;
            Vector3D posRelative = Pos - mainStar.Pos;
            Vector3D velRelative = Vel - mainStar.Vel;

            Vector3D angMomentum = Vector3D.Cross(posRelative, velRelative);

            Vector3D eccentricityVector = ((velRelative.Magnitude * velRelative.Magnitude - mu / posRelative.Magnitude) * posRelative - Vector3D.Dot(posRelative, velRelative) * velRelative) / mu;
            double eccentricity = eccentricityVector.Magnitude;
            double energy = velRelative.Magnitude * velRelative.Magnitude / 2 - mu / posRelative.Magnitude;

            double semiMajorAxis;
            if (Math.Abs(eccentricity - 1.0) > double.Epsilon)
                semiMajorAxis = -mu / (2 * energy);
            else
                semiMajorAxis = double.PositiveInfinity;

            double inclination = Math.Acos(angMomentum.z / angMomentum.Magnitude);

            return new(semiMajorAxis, eccentricity, inclination);
        }

        // https://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto
        public Keplerian ToKeplerian(Star mainStar)
        {
            double mu = Constants.G * mainStar.Mass;
            Vector3D posRelative = Pos - mainStar.Pos;
            Vector3D velRelative = Vel - mainStar.Vel;

            Vector3D angMomentum = Vector3D.Cross(posRelative, velRelative);
            Vector3D nodeVector = Vector3D.Cross(new Vector3D(0, 0, 1), velRelative);

            Vector3D eccentricityVector = ((velRelative.Magnitude * velRelative.Magnitude - mu / posRelative.Magnitude) * posRelative - Vector3D.Dot(posRelative, velRelative) * velRelative) / mu;
            double eccentricity = eccentricityVector.Magnitude;
            double energy = velRelative.Magnitude * velRelative.Magnitude / 2 - mu / posRelative.Magnitude;

            double semiMajorAxis, p;
            if (Math.Abs(eccentricity - 1.0) > double.Epsilon)
            {
                semiMajorAxis = -mu / (2 * energy);
                p = semiMajorAxis * (1 - eccentricity * eccentricity);
            }
            else
            {
                p = angMomentum.Magnitude * angMomentum.Magnitude / mu;
                semiMajorAxis = double.PositiveInfinity;
            }

            double inclination = Math.Acos(angMomentum.z / angMomentum.Magnitude);

            double longitudeAscending = Math.Acos(nodeVector.x / nodeVector.Magnitude);
            if (nodeVector.y < 0) longitudeAscending = 2 * Math.PI - longitudeAscending;

            double argumentPeriapsis = Math.Acos(Vector3D.Dot(nodeVector, eccentricityVector) / (nodeVector.Magnitude * eccentricity));
            if (eccentricityVector.z < 0) argumentPeriapsis = 2 * Math.PI - argumentPeriapsis;

            double trueAnormaly = Math.Acos(Vector3D.Dot(eccentricityVector, posRelative) / (eccentricity * posRelative.Magnitude));
            if (Vector3D.Dot(Pos, velRelative) < 0) trueAnormaly = 2 * Math.PI - trueAnormaly;

            return new(mainStar.Mass, semiMajorAxis, eccentricity, inclination, longitudeAscending, argumentPeriapsis, trueAnormaly);
        }
    }

    public class Planet : Particle
    {
        public double Mass;

        public Planet(double mass, Vector3D pos, Vector3D vel) : base(pos, vel)
        {
            Mass = mass;
        }

        public static Planet FromKeplerian(Keplerian keplerian,
                                               double planetMass)
        {
            Body3D cartesian = keplerian.ToCartesian();
            return new(planetMass, cartesian.Pos, cartesian.Vel);
        }
    }
}
