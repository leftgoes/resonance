﻿namespace Resonance
{
    public class Vector3D
    {
        public double x, y, z;

        public Vector3D(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public override string ToString()
        {
            return $"[{x}, {y}, {z}]";
        }

        public static Vector3D operator +(Vector3D u, Vector3D v)
        {
            return new(u.x + v.x, u.y + v.y, u.z + v.z);
        }

        public static Vector3D operator -(Vector3D u, Vector3D v)
        {
            return new(u.x - v.x, u.y - v.y, u.z - v.z);
        }

        public static Vector3D operator *(double s, Vector3D v)
        {
            return new(s * v.x, s * v.y, s * v.z);
        }

        public static Vector3D operator *(Vector3D v, double s)
        {
            return new(s * v.x, s * v.y, s * v.z);
        }

        public static Vector3D operator /(Vector3D v, double s)
        {
            return new(v.x / s, v.y / s, v.z / s);
        }

        public static Vector3D Cross(Vector3D u, Vector3D v)
        {
            return new(u.y * v.z - u.z * v.y,
                       u.z * v.x - u.x * v.z,
                       u.x * v.y - u.y * v.x);
        }

        public static double Dot(Vector3D u, Vector3D v)
        {
            return u.x * v.x + u.y * v.y + u.z * v.z;
        }


        public static Vector3D Zero { get { return new(0, 0, 0); } }

        public static Vector3D XUnit { get { return new(1, 0, 0); } }

        public static Vector3D YUnit { get { return new(0, 1, 0); } }

        public static Vector3D ZUnit { get { return new(0, 0, 1); } }

        public double Magnitude { get { return Math.Sqrt(x * x + y * y + z * z); } }

        public Vector3D Normalized { get { return this / Magnitude; } }
    }
}
