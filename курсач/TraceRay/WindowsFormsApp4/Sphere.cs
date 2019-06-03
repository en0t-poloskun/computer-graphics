namespace WindowsFormsApp4
{
    public class Figure
    {
        public float specular, reflective, transparent;
        public float[] color;
        public int flag;

        public Figure(float[] c, float s, float r, float t, int f)
        {
            specular = s;
            reflective = r;
            color = c;
            transparent = t;
            flag = f;
        }
    }

    public class Sphere : Figure
    {
        public float radius;
        public float[] center;

        public Sphere(float[] cent, float rad, float[] c, float s, float r, float t, int f) : base(c, s, r, t, f)
        {
            center = cent;
            radius = rad;
        }
    }

    public class Triangle : Figure
    {
        public float[] point1, point2, point3;

        public Triangle(float[] p1, float[] p2, float[] p3, float[] c, float s, float r, float t, int f) : base(c, s, r, t, f)
        {
            point1 = p1;
            point2 = p2;
            point3 = p3;
        }
    }

    public class Cilinder : Figure
    {
        public float[] center;
        public float A, B;

        public Cilinder(float[] cent, float a, float b, float[] c, float s, float r, float t, int f) : base(c, s, r, t, f)
        {
            center = cent;
            A = a;
            B = b;
        }
    }

    public class Circle : Figure
    {
        public float[] center;
        public float R;
        public float[] point2;
        public float[] point3;

        public Circle(float[] cent, float rad, float[] c, float s, float r, float t, int f) : base(c, s, r, t, f)
        {
            center = cent;
            R = rad;
            point2 = new float[3];
            point3 = new float[3];

            point2[0] = center[0]+R;
            point2[1] = center[1];
            point2[2] = center[2];

            point3[0] = center[0];
            point3[1] = center[1];
            point3[2] = center[2]+R;
          
        }
    }


    public class Light
    {
        public float ltype, intensity;
        public float[] position;

        public Light(float t, float i, float[] p)
        {
            ltype = t;
            intensity = i;
            position = p;
        }
    }
}