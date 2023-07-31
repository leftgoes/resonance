namespace Resonance
{
    public static class Elements
    {
        public const int SemiMajorAxis = 0;
        public const int Eccentricity = 1;
        public const int Inclination = 2;
    }

    public class BasicKeplerian
    {
        public double SemiMajorAxis;
        public double Eccentricity;
        public double Inclination;

        public BasicKeplerian(double semiMajorAxis, double eccentricity, double inclination)
        {
            SemiMajorAxis = semiMajorAxis;
            Eccentricity = eccentricity;
            Inclination = inclination;
        }

        public static BasicKeplerian Zero { get { return new(0.0, 0.0, 0.0); } }

        public bool IsEscaping { get { return Eccentricity >= 1.0; } }

        public static BasicKeplerian FromCartesian(double starMass, Vector3D pos, Vector3D vel)
        {
            double mu = Constants.G * starMass;
            Vector3D angMomentum = Vector3D.Cross(pos, vel);

            double eccentricity = (((vel.Magnitude * vel.Magnitude - mu / pos.Magnitude) * pos - Vector3D.Dot(pos, vel) * vel) / mu).Magnitude;

            double semiMajorAxis;
            if (Math.Abs(eccentricity - 1.0) > double.Epsilon)
                semiMajorAxis = -mu / (vel.Magnitude * vel.Magnitude - 2 * mu / pos.Magnitude);
            else
                semiMajorAxis = double.PositiveInfinity;

            double inclination = Math.Acos(angMomentum.z / angMomentum.Magnitude);

            return new(semiMajorAxis, eccentricity, inclination);
        }
    }
    
    public class BasicKeplerianNormal
    {
        public BasicKeplerian Mu;
        public BasicKeplerian Sigma;

        public BasicKeplerianNormal(BasicKeplerian mu, BasicKeplerian sigma)
        {
            Mu = mu;
            Sigma = sigma;
        }

        public BasicKeplerianNormal()
        {
            Mu = BasicKeplerian.Zero;
            Sigma = BasicKeplerian.Zero;
        }
    }

    public class Keplerian : BasicKeplerian
    {
        public double StarMass;
        public double LongitudeAscending;
        public double ArgumentPeriapsis;
        public double TrueAnomaly;

        public Keplerian(double starMass, double semiMajorAxis, double eccentricity, double inclination,
                         double longitudeAscending, double argumentPeriapsis, double trueAnomaly) : base(semiMajorAxis, eccentricity, inclination)
        {
            StarMass = starMass;
            LongitudeAscending = longitudeAscending;
            ArgumentPeriapsis = argumentPeriapsis;
            TrueAnomaly = trueAnomaly;
        }

        // https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
        public Body3D ToCartesian()
        {
            double mu = Constants.G * StarMass;

            double eccentricAnomaly = Math.Atan(Math.Sqrt(1 - Eccentricity * Eccentricity) * Math.Sin(TrueAnomaly) / (Eccentricity + Math.Cos(TrueAnomaly)));
            double distance = SemiMajorAxis * (1 - Eccentricity * Math.Cos(eccentricAnomaly));

            Vector3D oPos = distance * new Vector3D(Math.Cos(TrueAnomaly), Math.Sin(TrueAnomaly), 0);
            Vector3D oVel = Math.Sqrt(mu * SemiMajorAxis) / distance * new Vector3D(-Math.Sin(eccentricAnomaly),
                                                                                     Math.Sqrt(1 - Eccentricity * Eccentricity) * Math.Cos(eccentricAnomaly),
                                                                                     0);

            Vector3D pos = new(oPos.x * (Math.Cos(ArgumentPeriapsis) * Math.Cos(LongitudeAscending) - Math.Sin(ArgumentPeriapsis) * Math.Cos(Inclination) * Math.Sin(LongitudeAscending)) - oPos.y * (Math.Sin(ArgumentPeriapsis) * Math.Cos(LongitudeAscending) + Math.Cos(ArgumentPeriapsis) * Math.Cos(Inclination) * Math.Sin(LongitudeAscending)),
                               oPos.x * (Math.Cos(ArgumentPeriapsis) * Math.Sin(LongitudeAscending) + Math.Sin(ArgumentPeriapsis) * Math.Cos(Inclination) * Math.Cos(LongitudeAscending)) + oPos.y * (Math.Cos(ArgumentPeriapsis) * Math.Cos(Inclination) * Math.Cos(LongitudeAscending) - Math.Sin(ArgumentPeriapsis) * Math.Sin(LongitudeAscending)),
                               oPos.x * Math.Sin(ArgumentPeriapsis) * Math.Sin(Inclination) + oPos.y * Math.Cos(ArgumentPeriapsis) * Math.Sin(Inclination));
            Vector3D vel = new(oVel.x * (Math.Cos(ArgumentPeriapsis) * Math.Cos(LongitudeAscending) - Math.Sin(ArgumentPeriapsis) * Math.Cos(Inclination) * Math.Sin(LongitudeAscending)) - oVel.y * (Math.Sin(ArgumentPeriapsis) * Math.Cos(LongitudeAscending) + Math.Cos(ArgumentPeriapsis) * Math.Cos(Inclination) * Math.Sin(LongitudeAscending)),
                               oVel.x * (Math.Cos(ArgumentPeriapsis) * Math.Sin(LongitudeAscending) + Math.Sin(ArgumentPeriapsis) * Math.Cos(Inclination) * Math.Cos(LongitudeAscending)) + oVel.y * (Math.Cos(ArgumentPeriapsis) * Math.Cos(Inclination) * Math.Cos(LongitudeAscending) - Math.Sin(ArgumentPeriapsis) * Math.Sin(LongitudeAscending)),
                               oVel.x * Math.Sin(ArgumentPeriapsis) * Math.Sin(Inclination) + oVel.y * Math.Cos(ArgumentPeriapsis) * Math.Sin(Inclination));

            return new(pos, vel);
        }

        public new static Keplerian FromCartesian(double starMass, Vector3D pos, Vector3D vel)
        {
            double mu = Constants.G * starMass;
            Vector3D angMomentum = Vector3D.Cross(pos, vel);
            Vector3D nodeVector = Vector3D.Cross(Vector3D.ZUnit, vel);

            Vector3D eccentricityVector = ((vel.Magnitude * vel.Magnitude - mu / pos.Magnitude) * pos - Vector3D.Dot(pos, vel) * vel) / mu;
            double eccentricity = eccentricityVector.Magnitude;
            double energy = vel.Magnitude * vel.Magnitude / 2 - mu / pos.Magnitude;

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

            double trueAnormaly = Math.Acos(Vector3D.Dot(eccentricityVector, pos) / (eccentricity * pos.Magnitude));
            if (Vector3D.Dot(pos, vel) < 0) trueAnormaly = 2 * Math.PI - trueAnormaly;

            return new(starMass, semiMajorAxis, eccentricity, inclination, longitudeAscending, argumentPeriapsis, trueAnormaly);
        }

        public new static Keplerian Zero(double starMass) { return new(starMass, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); }
    }
}
