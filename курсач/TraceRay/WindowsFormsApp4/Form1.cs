using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using WindowsFormsApp4;
using System.Diagnostics;


namespace WindowsFormsApp4
{
    public partial class Form1 : Form
    {
        public const float viewport_size = 1;
        public const float projection_plane_z = 1;
        public const float recursion_depth = 3;
        static Sphere s1 = new Sphere(new float[] { 0, 0, -2 }, (float)0.5, new float[] { 157, 161, 170 }, 500, (float)0.5, (float)1, 0);
        static Sphere s2 = new Sphere(new float[] { -2, 0, 4 }, 1, new float[] { 0, 255, 0 }, 10, (float)0.4, (float)0.1, 0);
        static Sphere s3 = new Sphere(new float[] { -3, 0, 3 }, 1, new float[] { 0, 0, 255 }, 500, (float)0.3, (float)0.1, 0);
        static Sphere s4 = new Sphere(new float[] { 0, -5001, 0 }, 5000, new float[] { 255, 255, 0 }, 1000, (float)0.5, (float)0, 0);
        static Light l1 = new Light(0, (float)0.2, null);
        static Light l2 = new Light(1, (float)0.6, new float[] { 0, 0, -5});
        static Light l3 = new Light(2, (float)0.2, new float[] { 1, 4, 4 });
        static Triangle t1 = new Triangle(new float[] { 0, -1, 2 }, new float[] { 1, -1, 3 }, new float[] { 0, 1, 2 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t2 = new Triangle(new float[] { 1, 1, 3 }, new float[] { 1, -1, 3 }, new float[] { 0, 1, 2 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t3 = new Triangle(new float[] { 0, -1, 2 }, new float[] { -1, -1, 3 }, new float[] { 0, 1, 2 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t4 = new Triangle(new float[] { -1, 1, 3 }, new float[] { -1, -1, 3 }, new float[] { 0, 1, 2 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t5 = new Triangle(new float[] { 0, -1, 4 }, new float[] { 1, -1, 3 }, new float[] { 0, 1, 4 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t6 = new Triangle(new float[] { 1, 1, 3 }, new float[] { 1, -1, 3 }, new float[] { 0, 1, 4 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t7 = new Triangle(new float[] { -0, -1, 4 }, new float[] { -1, -1, 3 }, new float[] { 0, 1, 4 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t8 = new Triangle(new float[] { -1, 1, 3 }, new float[] { -1, -1, 3 }, new float[] { 0, 1, 4 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t9 = new Triangle(new float[] { 0, 1, 2 }, new float[] { 1, 1, 3 }, new float[] { 0, 1, 4 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t10 = new Triangle(new float[] { 0, 1, 2 }, new float[] { -1, 1, 3 }, new float[] { 0, 1, 4 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t11 = new Triangle(new float[] { 0, -1, 2 }, new float[] { 1, -1, 3 }, new float[] { 0, -1, 4 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Triangle t12 = new Triangle(new float[] { 0, -1, 2 }, new float[] { -1, -1, 3 }, new float[] { 0, -1, 4 }, new float[] { 255, 0, 0 }, 500, (float)0.2, (float)0, 0);
        static Cilinder c1 = new Cilinder(new float[] { 3, 0, 3 }, 1, 1, new float[] { 0, 255, 0 }, 10, (float)0.4, (float)0, 0);
        //static Circle c11 = new Circle(new float[] { 2, 3, 4 }, 1, new float[] { 0, 255, 0 }, 10, (float)0.4);
        //static Circle c12 = new Circle(new float[] { 2, -1, 4 }, 1, new float[] { 0, 255, 0 }, 10, (float)0.4);
        public float[] background_color = { 0, 0, 0 };
        Sphere[] spheres = { s1, s3, s4};
        Light[] lights = { l1, l2, l3 };
        Triangle[] triangles = { t2, t1, t3, t4, t5, t6, t8, t7, t10, t9 };
        private Cilinder[] cilinders = { c1};
        Circle[] circles = { };
        public Form1()
        {
            InitializeComponent();
            float[] camera_position = { 0, 0, -5 };
            float[,] camera_rotation = new float[,] { { 1, 0, 0},
                                                        {0, 1, 0},
                                                        {0, 0, 1}};
            Draw(camera_position, camera_rotation);
        }


        private void Draw(float[] camera_position, float[,] camera_rotation)
        {

            Bitmap bmp = new Bitmap(pictureBox1.Width, pictureBox1.Height);
            Graphics graph = Graphics.FromImage(bmp);
            graph.TranslateTransform(pictureBox1.Width / 2, pictureBox1.Height / 2);
            graph.FillRectangle(new SolidBrush(Color.Black), new Rectangle(0, 0, 1, 1));
            float[] direction = new float[3];
            float[] color = new float[3];
            //float[] camera_position = { -2, 0, -5};
            //float[,] camera_rotation = new float[,] { { 1, 0, 0},
            //                                            {0, 1, 0},
            //                                            {0, 0, 1}};

            //Stopwatch stopwatch = new Stopwatch();
            //stopwatch.Restart();
            //for (int k = 1; k <= 10; k++)
            //{ 
            for (int x = -pictureBox1.Width / 2; x < pictureBox1.Width / 2; x++)
            {

                for (int y = -pictureBox1.Height; y < pictureBox1.Height / 2; y++)
                {

                    direction = CanvasToViewport(x, y, direction);
                    direction = MultiplyMV(camera_rotation, direction);
                    color = TraceRay(camera_position, direction, 1, float.MaxValue, recursion_depth);
                    color = Clamp(color, color);
                    if (x > -pictureBox1.Width / 2 && x < pictureBox1.Width / 2 && y > -pictureBox1.Height / 2 && y < pictureBox1.Height / 2)
                    {
                        graph.FillRectangle(new SolidBrush(Color.FromArgb((int)color[0], (int)color[1], (int)color[2])), new Rectangle(x, y, 1, 1));
                    }
                }
            }
            //}
            //stopwatch.Stop();

            pictureBox1.Image = bmp;
            bmp.RotateFlip(RotateFlipType.Rotate180FlipX);
        }
        public float[] Clamp(float[] vec, float[] res)
        {
            res[0] = Math.Min(255, Math.Max(0, vec[0]));
            res[1] = Math.Min(255, Math.Max(0, vec[1]));
            res[2] = Math.Min(255, Math.Max(0, vec[2]));
            return res;
        }

        public float[] MultiplyMV(float[,] mat, float[] vec)
        {
            float[] result = new float[3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    result[i] += vec[j] * mat[i, j];
                }
            }
            return result;
        }

        public float DotProduct(float[] v1, float[] v2)
        {
            return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        }
        public float[] Subtract(float[] v1, float[] v2, float[] res)
        {
            res[0] = v1[0] - v2[0];
            res[1] = v1[1] - v2[1];
            res[2] = v1[2] - v2[2];
            return res;
        }

        public float[] Multiply(float k, float[] vec, float[] res)
        {
            res[0] = k * vec[0];
            res[1] = k * vec[1];
            res[2] = k * vec[2];
            return res;
        }
        public float[] VecMultiply(float[] AB, float[] AC, float[] n)
        {
            n[0] = AB[1] * AC[2] - AC[1] * AB[2];
            n[1] = AB[2] * AC[0] - AC[2] * AB[0];
            n[2] = AB[0] * AC[1] - AC[0] * AB[1];
            return n;
        }

        public float[] Add(float[] v1, float[] v2, float[] res)
        {
            res[0] = v1[0] + v2[0];
            res[1] = v1[1] + v2[1];
            res[2] = v1[2] + v2[2];
            return res;
        }

        public float Length(float[] vec)
        {
            return (float)Math.Sqrt(DotProduct(vec, vec));

        }
        public float[] ReflectRay(float[] v1, float[] v2, float[] res)
        {
            float[] r1 = new float[3];
            r1 = Multiply(2 * DotProduct(v1, v2), v2, r1);
            res = Subtract(r1, v1, res);
            return res;
        }

        public float[] IntersectRaySphere(float[] origin, float[] direction, Sphere sphere, float[] res)
        {
            float[] oc = new float[3];
            oc = Subtract(origin, sphere.center, oc);
            float k1 = DotProduct(direction, direction);
            float k2 = 2 * DotProduct(oc, direction);
            float k3 = DotProduct(oc, oc) - sphere.radius * sphere.radius;

            float discriminant = k2 * k2 - 4 * k1 * k3;
            if (discriminant < 0)
            {
                res[0] = float.MaxValue;
                res[1] = float.MaxValue;
                return res;
            }

            float t1 = (-k2 + (float)Math.Sqrt(discriminant)) / (2 * k1);
            float t2 = (-k2 - (float)Math.Sqrt(discriminant)) / (2 * k1);
            res[0] = t1;
            res[1] = t2;
            return res;
        }

        public float[] IntersectRayCilinder(float[] origin, float[] direction, Cilinder cilinder, float[] res)
        {
            float gran1 = cilinder.center[1] + cilinder.B;
            float gran2 = cilinder.center[1] - cilinder.B;
            float a = direction[0] * direction[0] + direction[2] * direction[2];
            float b = 2 * ((origin[0] - cilinder.center[0]) * direction[0] + (origin[2] - cilinder.center[2]) * direction[2]);
            float c = (origin[0] - cilinder.center[0]) * (origin[0] - cilinder.center[0]) + (origin[2] - cilinder.center[2]) * (origin[2] - cilinder.center[2]) - cilinder.A * cilinder.A;
            float D = b * b - 4 * a * c;
            if (D < 0)
            {
                res[0] = float.MaxValue;
                res[1] = float.MaxValue;
                return res;
            }

            float t1 = (-b + (float)Math.Sqrt(D)) / (2 * a);
            float t2 = (-b - (float)Math.Sqrt(D)) / (2 * a);
            float M1 = origin[1] + direction[1] * t1;
            float M2 = origin[1] + direction[1] * t2;

            if (((M1 < gran1) && (M1 > gran2)) || ((M2 < gran1) && (M2 > gran2)))
            {
                res[0] = t1;
                res[1] = t2;
                return res;

            }
            else
            {
                res[0] = float.MaxValue;
                res[1] = float.MaxValue;
                return res;
            }
        }

        public float IntersectRayTriangle(float[] origin, float[] direction, Triangle triangle, float[] n, float[] AB)
        {

            float D;
            D = -(n[0] * (triangle.point1[0]) + n[1] * (triangle.point1[1]) + n[2] * (triangle.point1[2]));
            float t;
            if (n[0] * direction[0] + n[1] * direction[1] + n[2] * direction[2] == 0)
            {
                t = float.MaxValue;
                return t;
            }
            t = -(n[0] * origin[0] + n[1] * origin[1] + n[2] * origin[2] + D) / (n[0] * direction[0] + n[1] * direction[1] + n[2] * direction[2]);
            float[] M = new float[3];
            M[0] = origin[0] + direction[0] * t;
            M[1] = origin[1] + direction[1] * t;
            M[2] = origin[2] + direction[2] * t;
            float[] na = new float[3];
            float[] AM = new float[3];
            AM = Subtract(M, triangle.point1, AM);
            na = VecMultiply(AB, AM, n);
            float[] nb = new float[3];
            float[] BC = new float[3];
            BC = Subtract(triangle.point3, triangle.point2, BC);
            float[] BM = new float[3];
            BM = Subtract(M, triangle.point2, BM);
            nb = VecMultiply(BC, BM, nb);
            float[] nc = new float[3];
            float[] CA = new float[3];
            CA = Subtract(triangle.point1, triangle.point3, CA);
            float[] CM = new float[3];
            CM = Subtract(M, triangle.point3, CM);
            nc = VecMultiply(CA, CM, nc);
            if ((na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2] > 0) && (na[0] * nc[0] + na[1] * nc[1] + na[2] * nc[2] >= 0))
            {
                return t;
            }
            else
            {
                t = float.MaxValue;
                return t;
            }

        }

        public float IntersectRayCircle(float[] origin, float[] direction, Circle circle, float[] n)
        {

            float D;
            D = -(n[0] * (circle.center[0]) + n[1] * (circle.center[1]) + n[2] * circle.center[2]);
            float t;
            if (n[0] * direction[0] + n[1] * direction[1] + n[2] * direction[2] == 0)
            {
                t = float.MaxValue;
                return t;
            }
            t = -(n[0] * origin[0] + n[1] * origin[1] + n[2] * origin[2] + D) / (n[0] * direction[0] + n[1] * direction[1] + n[2] * direction[2]);
            float[] M = new float[3];
            M[0] = origin[0] + direction[0] * t;
            M[2] = origin[2] + direction[2] * t;

            if ((M[0] - circle.center[0]) * (M[0] - circle.center[0]) + (M[2] - circle.center[2]) * (M[2] - circle.center[2]) <= circle.R * circle.R)
            {
                return t;
            }

            else
            {
                t = float.MaxValue;
                return t;
            }


        }

        public float ComputeLighting(float[] point, float[] normal, float[] view, float specular)
        {
            float Intensity = 0;
            float lenght_n = Length(normal);
            float lenght_v = Length(view);
            for (int i = 0; i < lights.Length; i++)
            {
                Light light = lights[i];
                if (light.ltype == 0)
                {
                    Intensity += light.intensity;
                }
                else
                {
                    float[] vec_l = new float[3];
                    float t_max;
                    if (light.ltype == 1)
                    {
                        vec_l = Subtract(light.position, point, vec_l);
                        t_max = 1;
                    }
                    else
                    {
                        vec_l = light.position;
                        t_max = float.MaxValue;

                    }
                    var blocker = ClosestIntersection(point, vec_l, (float)0.001, t_max);
                    if (blocker.Item1 != null && blocker.Item2 != 0)
                    {
                        continue;
                    }
                    float n_dot_l;
                    n_dot_l = DotProduct(normal, vec_l);
                    if (n_dot_l > 0)
                    {
                        Intensity += light.intensity * n_dot_l / (lenght_n * Length(vec_l));
                    }

                    if (specular != -1)
                    {
                        float[] vec_r = new float[3];
                        float[] res = new float[3];
                        vec_r = ReflectRay(vec_l, normal, vec_r);
                        float r_dot_v = DotProduct(vec_r, view);
                        if (r_dot_v > 0)
                        {
                            Intensity += light.intensity * (float)Math.Pow(r_dot_v / (Length(vec_r) * lenght_v), specular);

                        }
                    }

                }

            }
            return Intensity;
        }

        public Tuple<Figure, float, float, float[]> ClosestIntersection(float[] origin, float[] direction, float min_t, float max_t)
        {
            float sec_point = 0;
            float[] n = new float[3];
            int[] flag = new int[2];
            float[] point = new float[3];
            float[] res = new float[3];
            float[] ress = new float[3];
            float[] normal = new float[3];
            float[] ts = new float[2];
            float[] tc = new float[2];
            float tcc;
            float tt;
            float closest_t = float.MaxValue;
            Figure closest_figure = null;
            for (int i = 0; i < spheres.Length; i++)
            {
                if (spheres[i].flag == 0)
                {
                    ts = IntersectRaySphere(origin, direction, spheres[i], ts);

                    if (ts[0] < closest_t && min_t < ts[0] && ts[0] < max_t)
                    {

                        closest_t = ts[0];
                        closest_figure = spheres[i];
                        flag[0] = 1;
                        flag[1] = i;
                        sec_point = ts[1];
                    }
                    if (ts[1] < closest_t && min_t < ts[1] && ts[1] < max_t)
                    {
                        closest_t = ts[1];
                        closest_figure = spheres[i];
                        flag[0] = 1;
                        flag[1] = i;
                        sec_point = ts[0];
                    }
                }
            }

            for (int i = 0; i < triangles.Length; i++)
            {
                if (triangles[i].flag == 0)
                {
                    float[] AB = new float[3];
                    AB = Subtract(triangles[i].point2, triangles[i].point1, AB);
                    float[] AC = new float[3];
                    AC = Subtract(triangles[i].point3, triangles[i].point1, AC);

                    n = VecMultiply(AB, AC, n);
                    


                    tt = IntersectRayTriangle(origin, direction, triangles[i], n, AB);
                    //normal = Multiply(1 / Length(n), n, n);

                    if (tt < closest_t && min_t < tt && tt < max_t)
                    {
                        closest_t = tt;
                        closest_figure = triangles[i];
                        flag[0] = 0;
                        if (i % 2 != 0)
                        {
                            flag[1] = i - 1;
                        }
                        else
                        {
                            flag[1] = i;
                        }

                        //normal = Multiply(1 / Length(n), n, n);

                    }
                }
            }
            for (int i = 0; i < cilinders.Length; i++)
            {
                if (cilinders[i].flag == 0)
                {
                    tc = IntersectRayCilinder(origin, direction, cilinders[i], tc);
                    if (tc[0] < closest_t && min_t < tc[0] && tc[0] < max_t)
                    {

                        closest_t = tc[0];
                        closest_figure = cilinders[i];
                        flag[0] = 2;
                        flag[1] = i;
                    }
                    if (tc[1] < closest_t && min_t < tc[1] && tc[1] < max_t)
                    {

                        closest_t = tc[1];
                        closest_figure = cilinders[i];
                        flag[0] = 2;
                        flag[1] = i;
                    }
                }
            }




            if (closest_figure != null)
            {
                if (flag[0] == 1)
                {

                    point = Add(origin, Multiply(closest_t, direction, res), point);
                    normal = Subtract(point, spheres[flag[1]].center, normal);
                    normal = Multiply(1 / Length(normal), normal, normal);
                }
                if (flag[0] == 2)
                {
                    point = Add(origin, Multiply(closest_t, direction, ress), point);
                    float[] M = new float[3];
                    M[0] = cilinders[flag[1]].center[0];
                    M[1] = point[1];
                    M[2] = cilinders[flag[1]].center[2];
                    normal = Subtract(point, M, normal);
                    normal = Multiply(1 / Length(normal), normal, normal);
                }
                if (flag[0] == 0)
                {
                    float[] AB = new float[3];
                   
                    float[] AC = new float[3];
                    AB = Subtract(triangles[flag[1]].point2, triangles[flag[1]].point1, AB);
                    
                    AC = Subtract(triangles[flag[1]].point3, triangles[flag[1]].point1, AC);

                    normal = VecMultiply(AB, AC, normal);
                    //normal[0] = 2;
                    //normal[1] = 0;
                    //normal[2] = -2;
                    normal = Multiply(1 / Length(normal), normal, normal);
                }

                return new Tuple<Figure, float, float, float[]>(closest_figure, closest_t, sec_point, normal);
            }
            return new Tuple<Figure, float, float, float[]>(null, 0, 0, null);
        }



        public float[] TraceRay(float[] origin, float[] direction, float min_t, float max_t, float depth)
        {

            var intersection = ClosestIntersection(origin, direction, min_t, max_t);

            if (intersection.Item1 == null && intersection.Item2 == 0)
            {
                return background_color;
            }
            Figure closest_figure = intersection.Item1;
            float closest_t = intersection.Item2;

            float[] point = new float[3];
            float[] sec_point = new float[3];
            float[] res = new float[3];
            point = Add(origin, Multiply(closest_t, direction, res), point);
            sec_point = Add(origin, Multiply(intersection.Item3, direction, res), sec_point);
            float[] normal = intersection.Item4;
            

            float[] view = new float[3];
            //??
            view = Multiply(-1, direction, view);
            float lightning;
            //для блеска
            lightning = ComputeLighting(point, normal, view, closest_figure.specular);
            float[] local_color = new float[3];
            //вычисл лок цвет
            local_color = Multiply(lightning, closest_figure.color, local_color);
            if (closest_figure.reflective <= 0 || depth <= 0)
            {
                return local_color;
            }
            float[] reflected_ray = new float[3];
            reflected_ray = ReflectRay(view, normal, reflected_ray);
            float[] reflected_color = new float[3];
            //вычисл цвет отражения
            reflected_color = TraceRay(point, reflected_ray, (float)0.001, float.MaxValue, depth - 1);
            float[] r1 = new float[3];
            float[] r2 = new float[3];
            float[] r3 = new float[3];
            float[] result = new float[3];
            float[] transparent_color = new float[3];
            float[] transparent_ray = new float[3];
            float[] a = new float[3];
            float[] b = new float[3];
            float[] c = new float[3];
            float[] k = new float[3];
            k = Multiply(-1, direction, k);
            a = Multiply(DotProduct(k, normal) * (float)1.49, k, a);
            b = Multiply((float)1.49, k, b);
            c = Multiply((float)Math.Sqrt(1 + (float)1.49 * (float)1.49 * DotProduct(k, normal) * DotProduct(k, normal) - 1), normal, c);
            transparent_ray = Subtract(a, b, transparent_ray);
            transparent_ray = Subtract(transparent_ray, c, transparent_ray);




            closest_figure.flag = 1;
            transparent_color = TraceRay(sec_point, transparent_ray, 1, float.MaxValue, depth);
            closest_figure.flag = 0;
            //смеш цветов?
            r1 = Multiply(1 - closest_figure.reflective, local_color, r1);
            r2 = Multiply(closest_figure.reflective, reflected_color, r2);
            r3 = Multiply(closest_figure.transparent, transparent_color, r3);

            result = Add(r1, r2, result);
            result = Add(result, r3, result);
            return result;
        }
        public float[] CanvasToViewport(int x, int y, float[] d)
        {


            d[0] = x * viewport_size / pictureBox1.Width;
            d[1] = y * viewport_size / pictureBox1.Height;

            d[2] = projection_plane_z;
            return d;
        }

        public float[,] Multiplication_Matrix(float[,] a, float[,] b)
        {

            float[,] r = new float[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        r[i, j] += a[i, k] * b[k, j];
                    }
                }
            }
            return r;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            float[] camera_position = new float[3];
            float[,] camera_rotation = new float[,] { { 1, 0, 0},
                                                        {0, 1, 0},
                                                        {0, 0, 1}};

            camera_position[0] = (float)Convert.ToDouble(textBox1.Text);
            camera_position[1] = (float)Convert.ToDouble(textBox2.Text);
            camera_position[2] = (float)Convert.ToDouble(textBox3.Text);
            float x = (float)(((float)Convert.ToDouble(textBox6.Text) / 180) * Math.PI);
            float y = (float)(((float)Convert.ToDouble(textBox4.Text) / 180) * Math.PI);
            float z = (float)(((float)Convert.ToDouble(textBox5.Text) / 180) * Math.PI);
            float[,] MX = new float[,] { { 1, 0, 0},
                                          { 0, (float) Math.Cos(x), (float) -Math.Sin(x)},
                                          {0, (float) Math.Sin(x), (float) Math.Cos(x)}};

            float[,] MY = new float[,] { { (float) Math.Cos(y), 0, (float) Math.Sin(y)},
                                         {0, 1, 0},
                                         {(float) -Math.Sin(y), 0, (float) Math.Cos(y)}};

            float[,] MZ = new float[,] { { (float) Math.Cos(z), (float) -Math.Sin(z), 0},
                                         {(float) Math.Sin(z), (float) Math.Cos(z), 0},
                                         {0, 0, 1}};
            camera_rotation = Multiplication_Matrix(camera_rotation, MZ);
            camera_rotation = Multiplication_Matrix(camera_rotation, MY);
            camera_rotation = Multiplication_Matrix(camera_rotation, MX);

            Draw(camera_position, camera_rotation);
        }


    }
}
