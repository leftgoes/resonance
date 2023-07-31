using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Resonance
{
    public class DoubleRange
    {
        public double Min;
        public double Max;

        public DoubleRange(double min, double max)
        {
            Min = min;
            Max = max;
        }

        public DoubleRange()
        {
            Min = 0.0;
            Max = 1.0;
        }

        public override string ToString()
        {
            return $"[{Min}, {Max}]";
        }
    }


    public abstract class SimulationData
    {
        public readonly int StepsCount;
        public readonly int PlanetCount;
        public readonly int ParticlesCount;

        public SimulationData(int stepsCount, int planetCount, int particlesCount)
        {
            StepsCount = stepsCount;
            PlanetCount = planetCount;
            ParticlesCount = particlesCount;
        }

        protected IEnumerable<(int step, int particleIndex)> ParticlesIndices()
        {
            for (int step = 0; step < StepsCount; step++) {
                for (int particleIndex = 0; particleIndex < ParticlesCount; particleIndex++) {
                    yield return (step, particleIndex);
                }
            }
        }

        protected abstract IEnumerable<(int step, int particleIndex)> NonNaNParticlesIndices();
    }

    public class BasicKeplerianData : SimulationData
    {
        double[,,] PlanetData;
        double[,,] ParticleData;

        public BasicKeplerianData(int stepsCount, int planetCount, int particlesCount) : base(stepsCount, planetCount, particlesCount)
        {
            PlanetData = new double[stepsCount, planetCount, 3];
            ParticleData = new double[stepsCount, particlesCount, 3];
        }

        protected override IEnumerable<(int step, int particleIndex)> NonNaNParticlesIndices()
        {
            for (int step = 0; step < StepsCount; step++)
            {
                for (int particleIndex = 0; particleIndex < ParticlesCount; particleIndex++)
                {
                    if (double.IsNaN(ParticleData[step, particleIndex, 0]))
                        continue;
                    yield return (step, particleIndex);
                }
            }
        }

        public void AddPlanetDatapoint(int step, int planetIndex, BasicKeplerian keplerian)
        {
            PlanetData[step, planetIndex, Elements.SemiMajorAxis] = keplerian.SemiMajorAxis;
            PlanetData[step, planetIndex, Elements.Eccentricity] = keplerian.Eccentricity;
            PlanetData[step, planetIndex, Elements.Inclination] = keplerian.Inclination;
        }

        public void AddParticleDatapoint(int step, int particleIndex, BasicKeplerian keplerian)
        {
            if (keplerian.IsEscaping)
            {
                ParticleData[step, particleIndex, Elements.SemiMajorAxis] = double.NaN;
                ParticleData[step, particleIndex, Elements.Eccentricity] = double.NaN;
                ParticleData[step, particleIndex, Elements.Inclination] = double.NaN;
            }
            else
            {
                ParticleData[step, particleIndex, Elements.SemiMajorAxis] = keplerian.SemiMajorAxis;
                ParticleData[step, particleIndex, Elements.Eccentricity] = keplerian.Eccentricity;
                ParticleData[step, particleIndex, Elements.Inclination] = keplerian.Inclination;
            }
        }

        private double LinMap(double x, double f1, double f2, double t1, double t2)
        {
            return (t1 - t2) * (x - f1) / (f2 - f1) + t1;
        }

        private double LinMap(double x, DoubleRange from, DoubleRange to)
        {
            return (to.Min - to.Max) * (x - from.Min) / (from.Max - from.Min) + to.Min;
        }

        private (DoubleRange x, DoubleRange y) ParticlesXYRanges(int xIndex, int yIndex, double percentile)
        {
            int length = NonNaNParticlesIndices().Count();
            double[] xFlattened = new double[length];
            double[] yFlattened = new double[length];

            int flattenedIndex = 0;
            foreach ((int step, int particleIndex) in NonNaNParticlesIndices()) {
                xFlattened[flattenedIndex] = ParticleData[step, particleIndex, xIndex];
                yFlattened[flattenedIndex] = ParticleData[step, particleIndex, yIndex];
                flattenedIndex++;
            }

            Array.Sort(xFlattened); Array.Sort(yFlattened);


            DoubleRange xRange = new(xFlattened[(int)(percentile / 100 * length)],
                                     xFlattened[(int)((1.0 - percentile / 100) * length)]);
            DoubleRange yRange = new(yFlattened[(int)(percentile / 100 * length)],
                                     yFlattened[(int)((1.0 - percentile / 100) * length)]);

            return (xRange, yRange);
        }

        public Videoframes ParticlesToVideoframes(int width, int height, int xIndex, int yIndex, double percentile)
        {
            Videoframes frames = new(width, height, StepsCount);

            var (xRange, yRange) = ParticlesXYRanges(xIndex, yIndex, percentile);

            Func<double, double> x2i = (double x) => LinMap(x, xRange, new(0.0, frames.Width));
            Func<double, double> y2j = (double y) => LinMap(y, yRange, new(frames.Height, 0.0));

            foreach ((int step, int particleIndex) in ParticlesIndices())
            {
                double x = ParticleData[step, particleIndex, xIndex];
                double y = ParticleData[step, particleIndex, yIndex];
                double i = x2i(x);
                double j = y2j(y);
                frames.DrawPoint(step, i, j);
            }
                

            return frames;
        }
    }
}
