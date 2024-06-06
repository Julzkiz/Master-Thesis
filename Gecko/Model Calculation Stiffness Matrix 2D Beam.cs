using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Statistics;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Gecko
{
    internal class Model_Calculation_Stiffness_Matrix_2D_Beam
    {
        public static Vector<double> FEM_CALC_Displacement_6(Class_Model model, out Matrix<double> t, out Matrix<double> K, out Matrix<double> M, out int d)
        {
            List<int> gIDs = new List<int>();
            foreach (Node node in model.nodes)
            {
                gIDs.Add(node.globalID); //all nodes and their global IDs
            }

            List<int> DOFs = gIDs.Distinct().ToList(); //no duplicates, only global ids
                                                       //Creating MathNet Matrix with d amount of DOFs per node
            d = 3;
            int Msize = DOFs.Count * d;
            M = DenseMatrix.Create(Msize, Msize, 0);
            K = DenseMatrix.Create(2 * d, 2 * d, 0);
            t = DenseMatrix.Create(2 * d, 2 * d, 0);

            List<double> L_vector = new List<double>(); //creating a vector with correct size
                                                        //Vector<double> V = Vector<double>.Build.Dense(DOFs.Count * d);
            for (int i = 0; i < Msize; i++)
            {
                L_vector.Add(0);
            }

            foreach (Load load in model.loads) //maybe just use V.Dense(DOFs.Count*d) -> vector of x length
            {
                int nodeid = model.nodes.Find((node) => node.point.DistanceTo(load.point) < 0.00003).globalID;
                L_vector[nodeid * d] = load.vector.X; //kN
                L_vector[nodeid * d + 1] = load.vector.Z;
                L_vector[nodeid * d + 2] = load.rotations[0];
            }

            int startID = 0;
            int endID = 0;

            foreach (Element beam in model.bars) //go through every bar and assign them a stiffness, start and end, and length.
            {
                double k = (beam.area * beam.E) / (beam.length * 1e3); // mm2 * kN/mm2 / mm => kN/mm 
                double k1 = (12 * beam.E * beam.Iy) / (Math.Pow(beam.length * 1e3, 3)); // kN/mm2 * mm4 / mm3 => kN/mm 
                double k2 = (6 * beam.E * beam.Iy) / (Math.Pow(beam.length * 1e3, 2)); // kN/mm2 * mm4 / mm2 => kN
                double k3 = (4 * beam.E * beam.Iy) / (beam.length * 1e3); // N/mm2 * mm4 / mm => kN*mm
                double k4 = 0.5 * k3; // kN*mm

                foreach (Node node in model.nodes) //assigning the globalID node to the start and end of bar, creating 2 indexes per node.
                {
                    if (node.point.DistanceTo(beam.axis.From) < 0.0003)
                    {
                        startID = node.globalID * d;
                    }
                    if (node.point.DistanceTo(beam.axis.To) < 0.0003)
                    {
                        endID = node.globalID * d;
                    }
                }
                //2d truss element [u1,v1,r1,u2,v2,r2]
                Matrix<double> m = DenseMatrix.OfArray(new double[,]
                //kN, kN, kNm, kN, kN, kNm
                {    { k,  0,  0,-k,  0,  0},
                 { 0, k1,-k2, 0,-k1,-k2},
                 { 0,-k2, k3, 0, k2, k4},
                 {-k,  0,  0, k,  0,  0},
                 { 0,-k1, k2, 0, k1, k2},
                 { 0,-k2, k4, 0, k2, k3}});

                //Assuming the transformation matrix is going in here
                double cos = (beam.axis.To.X - beam.axis.From.X) / (beam.length); //(m-m)/m
                double sin = (beam.axis.To.Z - beam.axis.From.Z) / (beam.length);

                t = DenseMatrix.OfArray(new double[,]
                    {{ cos, sin,  0,    0,   0,  0},
                 {-sin, cos,  0,    0,   0,  0},
                 {   0,   0,  1,    0,   0,  0},
                 {   0,   0,  0,  cos, sin,  0},
                 {   0,   0,  0, -sin, cos,  0},
                 {   0,   0,  0,    0,   0,  1}});

                Matrix<double> tT = t.Transpose();
                Matrix<double> T = tT * m * t;

                K = T;

                for (int i = 0; i < d; i++)
                    for (int j = 0; j < d; j++)
                    {
                        M[startID + i, startID + j] += K[i, j];
                        M[startID + i, endID + j] += K[i, d + j];
                        M[endID + i, startID + j] += K[d + i, j];
                        M[endID + i, endID + j] += K[d + i, d + j];
                    }
            }

            List<int> node_ints = new List<int>();
            foreach (Support support in model.supports)
            {
                int node_int = model.nodes.Find((node) => node.point.DistanceTo(support.point) < 0.00003).globalID;
                node_ints.Add(node_int);
            }

            node_ints.Sort();
            int removed = 0;
            foreach (int index in node_ints)
            {
                for (int i = 0; i < d; i++)
                {
                    if (model.supports[node_ints.IndexOf(index)].bcds[i] == 1)
                    {
                        M = M.RemoveColumn((index + index * (d - 1) + i) - removed);
                        M = M.RemoveRow((index + index * (d - 1) + i) - removed);

                        L_vector.RemoveAt((index + index * (d - 1) + i) - removed);

                        removed++;
                    }
                }
            }

            double[] array = L_vector.ToArray<double>();
            Vector<double> v = Vector<double>.Build.DenseOfEnumerable(L_vector);

            Vector<double> R = M.Solve(v); // (mm/kN) * kN => (mm)
            R = R / 1000; //m
            Matrix<double> invM = M.Inverse();

            return R;

        }

        public static List<Curve> InterpolatedGeometry_2D(Class_Model calc_model, int n, out double strain, out double stress, out double force, out double max_momenty_model, out List<int> id)
        {
            //4-dof beam element
            Vector<double> u = DenseVector.Create(6, 0);  //all displacement for all nodes
            Matrix<double> dispMatrix = DenseMatrix.OfArray(new double[n + 1, 2]);
            List<Point3d> newgeo = new List<Point3d>();
            List<Curve> displacedgeometry = new List<Curve>();
            stress = 0;
            strain = 0;
            force = 0;
            List<double> momenty_model = new List<double>();
            max_momenty_model = new double();
            double max_momenty = new double();
            id = new List<int>();

            Matrix<double> reduced_transmatrix = calc_model.transmatrix.SubMatrix(0, 3, 0, 3);

            foreach (Element beam in calc_model.bars)
            {
                double x = 0;
                double L = beam.length * 1e3; //mm

                int startid = calc_model.nodes.Find((node) => node.point.DistanceTo(beam.startnode) < 0.00003).globalID;
                int endid = calc_model.nodes.Find((node) => node.point.DistanceTo(beam.endnode) < 0.00003).globalID;

                List<Point3d> oldpoints = new List<Point3d>();
                List<Point3d> newpoints = new List<Point3d>();

                List<double> momenty = new List<double>();

                double[] t = MathNet.Numerics.Generate.LinearSpaced(n + 1, 0, 1);

                for (int i = 0; i < t.Length; i++)
                {
                    Point3d point_on_line = new Point3d();
                    point_on_line.Interpolate(beam.axis.From, beam.axis.To, t[i]);
                    oldpoints.Add(point_on_line);
                    newpoints.Add(point_on_line);
                }

                for (int i = 0; i < 3; i++)
                {
                    u[i] = calc_model.displacements[startid * 3 + i];
                    u[i + 3] = calc_model.displacements[endid * 3 + i];
                }

                u = u * calc_model.transmatrix; //to local deformations

                for (int i = 0; i < n + 1; i++)
                {
                    double N1 = 1 - x / L;
                    double N2 = x / L;
                    double N3 = 1 - 3 * Math.Pow(x, 2) / Math.Pow(L, 2) + 2 * Math.Pow(x, 3) / Math.Pow(L, 3);
                    double N4 = x - 2 * Math.Pow(x, 2) / L + Math.Pow(x, 3) / Math.Pow(L, 2);
                    double N5 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * Math.Pow(x, 3) / Math.Pow(L, 3);
                    double N6 = -Math.Pow(x, 2) / L + Math.Pow(x, 3) / Math.Pow(L, 2);

                    Matrix<double> N = DenseMatrix.OfArray(new double[,] {
                    //ux, uz, r1, ux,  uz,  r2
                    { N1,  0,  0, N2,   0,  0},    //u(x) (axial)
                    {  0, N3, -N4,  0,  N5, -N6}});  //w(x) (translation)

                    double dN3 = -6 * x / Math.Pow(L, 2) + 6 * Math.Pow(x, 2) / Math.Pow(L, 3); // 1/mm * mm 
                    double dN4 = 1 - 4 * x / L + 3 * Math.Pow(x, 2) / Math.Pow(L, 2); // unitless 
                    double dN5 = 6 * x / Math.Pow(L, 2) - 6 * Math.Pow(x, 2) / Math.Pow(L, 3); // 1/mm
                    double dN6 = -2 * x / L + 3 * Math.Pow(x, 2) / Math.Pow(L, 2); // unitless

                    Vector<double> rN = DenseVector.OfArray(new double[] {
                    //ux, uz, r1, ux,  uz,  rx 
                      0, dN3, dN4,  0,  dN5, dN6});  //rotation

                    double dN1 = -1 / L; //1/mm
                    double dN2 = 1 / L;
                    double ddN3 = -6 / Math.Pow(L, 2) + 12 * x / Math.Pow(L, 3); //1/mm2 
                    double ddN4 = -4 / L + 6 * x / Math.Pow(L, 2); // 1/mm
                    double ddN5 = 6 / Math.Pow(L, 2) - 12 * x / Math.Pow(L, 3); // 1/mm2
                    double ddN6 = -2 / L + 6 * x / Math.Pow(L, 2); //1/mm

                    Vector<double> B = DenseVector.OfArray(new double[]
                    //ux,   uz,  ry,   ux,   uz,   ry,  
                    {0, ddN3, -ddN4, 0, ddN5, -ddN6}); //(strain field)

                    dispMatrix.SetRow(i, N.Multiply(u));    //Sets row 0,1,..n+1 of the N matrix (needed for deformation in x,y and z) (i = x/L location on beam)
                    var rotation = rN * u;  //Sets strain

                    Vector<double> u_new = DenseVector.OfArray(new double[]
                    { dispMatrix[i,0], dispMatrix[i,1], rotation});

                    //Transposing from local to global
                    u_new = u_new * reduced_transmatrix.Transpose(); //

                    dispMatrix.SetRow(i, new double[] { u_new[0], u_new[1] });
                    dispMatrix = dispMatrix / 1e3; //mm -> m
                    rotation = u_new[2];

                    strain = -0.5 * beam.heigth * rN * u; //mm  
                    stress = beam.E * strain;
                    momenty.Add((-beam.E * beam.Iy * B * u) / 1e3); //kN/mm2 *mm4 * 1/m m2 * mm => kN mm => kNm
                    momenty_model.Add((-beam.E * beam.Iy * B * u) / 1e3);

                    Point3d dispp = new Point3d(dispMatrix[i, 0], 0, dispMatrix[i, 1]);
                    newgeo.Add(dispp);
                    newpoints[i] = newpoints[i] + new Point3d(dispMatrix[i, 0], 0, dispMatrix[i, 1]);

                    x += L / n;
                }

                beam.stress = stress;
                force = stress * beam.area;
                max_momenty = momenty.MaximumAbsolute();
                max_momenty_model = momenty_model.MaximumAbsolute();

                if (momenty.IndexOf(max_momenty) == -1)
                {
                    max_momenty = -max_momenty;
                }

                if (momenty_model.IndexOf(max_momenty_model) == -1)
                {
                    max_momenty_model = -max_momenty_model;
                }

                int max_momenty_index = momenty.IndexOf(max_momenty); //index of max moment (location at beam)
                Point3d max_momenty_x = oldpoints[max_momenty_index]; //point in which max moment occurs on beam.
                Curve curve = Curve.CreateControlPointCurve(newpoints, n);

                //Assigning which beam has which max moment.
                id.Add(beam.ID);
                beam.maxmoment = max_momenty;
                beam.maxmomentpoint = max_momenty_x;
                beam.deformed = curve;

                displacedgeometry.Add(curve);
            }

            return displacedgeometry;
        }

    }
}
