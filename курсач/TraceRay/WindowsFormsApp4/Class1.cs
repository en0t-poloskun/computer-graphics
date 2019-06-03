using System;
using WindowsFormsApp4;

namespace WindowsFormsApp4
{
    internal class Sphere
    {
        public float radius;
        public float[] center, color;

        public Sphere(float[] cent, float r, float[] col)
        {
            center = cent;
            radius = r;
            color = col;
        }
    }
}
