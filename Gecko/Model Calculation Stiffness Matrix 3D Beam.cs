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
    internal class Model_Calculation_Stiffness_Matrix_3D_Beam
    {
        public static List<double> FEM_CALC_Displacement_12(Class_Model model, out Matrix<double> t, out Matrix<double> K, out Matrix<double> M, out int d)
        {
            List<int> gIDs = new List<int>();
            foreach (Node node in model.nodes)
            {
                gIDs.Add(node.globalID); //all nodes and their global IDs
            }

            List<int> DOFs = gIDs.Distinct().ToList(); //no duplicates, only global ids
                                                       //Creating MathNet Matrix with d amount of DOFs per node
            d = 6;
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
                L_vector[nodeid * d] = load.vector.X;
                L_vector[nodeid * d + 1] = load.vector.Y;
                L_vector[nodeid * d + 2] = load.vector.Z;
                L_vector[nodeid * d + 3] = load.rotations[0];
                L_vector[nodeid * d + 4] = load.rotations[1];
                L_vector[nodeid * d + 5] = load.rotations[2];
            }

            int startID = 0;
            int endID = 0;

            foreach (Element beam in model.bars) //go through every bar and assign them a stiffness, start and end, and length.
            {
                double k = (beam.area * beam.E) / (beam.length * 1e3); // mm2 * kN/mm2 / mm => kN/mm
                double k2 = (beam.G * beam.J) / (beam.length * 1e3); // kN/mm2 * mm4 / mm => kN*mm

                double k1z = 12 * beam.E * beam.Iy / Math.Pow(beam.length * 1e3, 3); //kN/mm2 * mm4 / mm3 => kN/mm
                double k2z = 6 * beam.E * beam.Iy / Math.Pow(beam.length * 1e3, 2); //kN/mm2 * mm4 / mm2 => kN
                double k3z = 4 * beam.E * beam.Iy / (beam.length * 1e3); // kN/mm2 * mm4 / mm => kN*mm
                double k4z = 0.5 * k3z;

                double k1y = 12 * beam.E * beam.Iz / Math.Pow(beam.length * 1e3, 3);
                double k2y = 6 * beam.E * beam.Iz / Math.Pow(beam.length * 1e3, 2);
                double k3y = 4 * beam.E * beam.Iz / (beam.length * 1e3);
                double k4y = 0.5 * k3y;

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
                //3d beam element [u1,u2,u3,u4,u5,u6]
                Matrix<double> m = DenseMatrix.OfArray(new double[,]
                //ux,  uy,  uz, rx,  ry,  rz, ux,  uy,  uz, rx,  ry,  rz
                {{ k,   0,   0,  0,   0,   0, -k,   0,   0,  0,   0,   0},   //ux
                 { 0, k1y,   0,  0,   0, k2y,  0,-k1y,   0,  0,   0, k2y},   //uy
                 { 0,   0, k1z,  0,-k2z,   0,  0,   0,-k1z,  0,-k2z,   0},   //uz
                 { 0,   0,   0, k2,   0,   0,  0,   0,   0,-k2,   0,   0},   //rx 
                 { 0,   0,-k2z,  0, k3z,   0,  0,   0, k2z,  0, k4z,   0},   //ry => unit problem
                 { 0, k2y,   0,  0,   0, k3y,  0,-k2y,   0,  0,   0, k4y},   //rz
                 {-k,   0,   0,  0,   0,   0,  k,   0,  0,   0,   0,   0},   //ux
                 { 0,-k1y,   0,  0,   0,-k2y,  0, k1y,   0,  0,   0,-k2y},   //uy
                 { 0,   0,-k1z,  0, k2z,   0,  0,   0, k1z,  0, k2z,   0},   //uz 
                 { 0,   0,   0,-k2,  0,    0,  0,   0,   0, k2,   0,   0},   //rx
                 { 0,   0,-k2z,  0, k4z,   0,  0,   0, k2z,  0, k3z,   0},   //ry
                 { 0, k2y,   0,  0,   0, k4y,  0,-k2y,   0,  0,   0, k3y}}); //rz

                //Assuming the transformation matrix is going in here
                double cosx = (beam.axis.To.X - beam.axis.From.X) / beam.length;
                double cosy = (beam.axis.To.Y - beam.axis.From.Y) / beam.length;
                double cosz = (beam.axis.To.Z - beam.axis.From.Z) / beam.length;

                double cosxz = Math.Sqrt(Math.Pow(cosx, 2) + Math.Pow(cosz, 2));
                double c = Math.Cos(0);
                double s = Math.Sin(0);

                Matrix<double> RR = DenseMatrix.OfArray(new double[3, 3]);

                if ((beam.axis.To.X - beam.axis.From.X) == 0 && (beam.axis.To.Z - beam.axis.From.Z) == 0)
                {
                    RR = DenseMatrix.OfArray(new double[,] {
                    {      0, cosy,  0},
                    { cosy*c,    0,  s},
                    { cosy*s,    0,  -c}});
                }
                else
                {
                    RR = DenseMatrix.OfArray(new double[,] {
                    {                       cosx,     cosy,                        cosz},
                    {( cosx*cosy*c-cosz*s)/cosxz, -cosxz*c, ( cosy*cosz*c+cosx*s)/cosxz},
                    {( cosx*cosy*s+cosz*c)/cosxz,  cosxz*s, ( cosy*cosz*s-cosx*c)/cosxz}});
                }

                var t1 = RR.DiagonalStack(RR);
                t = t1.DiagonalStack(t1);

                beam.transmatrix = t;

                Matrix<double> tT = t.Transpose();
                K = t * m * tT;

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

            for (int i = 0; i < M.ColumnCount; i++)
                if (M[i, i] == 0)
                {
                    M[i, i] = 1;
                }

            double[] array = L_vector.ToArray<double>();
            Vector<double> v = Vector<double>.Build.DenseOfEnumerable(L_vector);

            Vector<double> R = M.Solve(v); //mm

            List<double> R12 = new List<double>();

            for (int i = 0; i < R.Count; i += 6)
            {
                // Process the first three items in the current block
                for (int j = i; j < i + 3 && j < R.Count; j++)
                {
                    R12.Add(R[j] / 1e3); //mm => m
                }

                // Skip the next three items in the current block
                if (i + 3 < R.Count)
                {
                    for (int k = i + 3; k < i + 6 && k < R.Count; k++)
                    {
                        R12.Add(R[k]); //leaving rotations out
                    }
                }
            }

            return R12; //m

        }

        public static List<Curve> InterpolatedGeometry3D(Class_Model calc_model, int n, out Vector<double> strains, out Vector<double> stresses, out List<double> forces, out double max_momenty_model, out double max_momentz_model, out List<int> id)
        {
            Vector<double> u = DenseVector.Create(12, 0);  //all displacement for all nodes
            Matrix<double> dispMatrix = DenseMatrix.OfArray(new double[n + 1, 4]);
            Matrix<double> rotMatrix = DenseMatrix.OfArray(new double[n + 1, 4]);
            List<Curve> displacedgeometry = new List<Curve>();
            strains = DenseVector.Create((n + 1) * 3, 0);
            stresses = DenseVector.Create((n + 1) * 3, 0);
            Vector<double> forces_e = DenseVector.Create((n + 1) * 3, 0);
            forces = new List<double>();

            List<double> momenty_model = new List<double>();
            List<double> momentz_model = new List<double>();

            max_momenty_model = 0;
            max_momentz_model = 0;

            double max_momenty = 0;
            double max_momentz = 0;

            id = new List<int>();


            foreach (Element beam in calc_model.bars) //need the displacement for each node of this bar. so u = [0,0,0,0,0,0,ux2,uy2,uz2,t1,t2,t3] if bar fixed - free
            {
                double x = 0;
                double L = beam.length; //m
                double z = (beam.heigth / 1e3) * 0.5; //m
                double y = (beam.width / 1e3) * 0.5; //m
                int startid = calc_model.nodes.Find((node) => node.point.DistanceTo(beam.startnode) < 0.00003).globalID;
                int endid = calc_model.nodes.Find((node) => node.point.DistanceTo(beam.endnode) < 0.00003).globalID;
                Matrix<double> reduced_transmatrix = beam.transmatrix.SubMatrix(0, 6, 0, 6);

                List<Point3d> newpoints = new List<Point3d>();
                List<Point3d> oldpoints = new List<Point3d>();
                List<double> rotations = new List<double>();

                double[] t = MathNet.Numerics.Generate.LinearSpaced(n + 1, 0, 1);

                for (int i = 0; i < t.Length; i++) //creates points on beam of even length.
                {
                    Point3d point_on_line = new Point3d();
                    point_on_line.Interpolate(beam.axis.From, beam.axis.To, t[i]);
                    newpoints.Add(point_on_line);
                    oldpoints.Add(point_on_line);
                }

                for (int i = 0; i < 6; i++) //finds start and end node, assigning its deformation values (ux1, uy1, uz1, rx1, ry1, rz1, ux2, uy2, uz2, rx2, ry2, rz2) => starting 6 and ending 6.
                {
                    u[i] = calc_model.displacements[startid * 6 + i]; //m
                    u[i + 6] = calc_model.displacements[endid * 6 + i]; //m
                }

                u = beam.transmatrix * u; //transforms from global to local deformations [12x1] x [12x12], transformation matrix of beam (can be in x, y or z-plane)

                for (int i = 0; i < n + 1; i++) //translation and rotation calculationbeam 
                {
                    double N1 = 1 - x / L;
                    double N2 = x / L;
                    double N3 = 1 - 3 * Math.Pow(x, 2) / Math.Pow(L, 2) + 2 * Math.Pow(x, 3) / Math.Pow(L, 3);
                    double N4 = x - 2 * Math.Pow(x, 2) / L + Math.Pow(x, 3) / Math.Pow(L, 2);
                    double N5 = -N3 + 1;
                    double N6 = Math.Pow(x, 3) / Math.Pow(L, 2) - Math.Pow(x, 2) / L;

                    Matrix<double> N = DenseMatrix.OfArray(new double[,] {
                //ux, uy,  uz, rx,  ry,  rz, ux,  uy, uz,  rx,  ry,  rz
                { N1, 0,    0,  0,   0,   0, N2,   0,   0,  0,   0,   0}, // == ux => m
                {  0, N3,   0,  0,   0,  N4,  0,  N5,   0,  0,   0,  N6}, // == uy => m
                {  0,  0,  N3,  0, -N4,   0,  0,   0,  N5,  0, -N6,   0}, // == uz => m
                {  0,  0,   0, N1,   0,  0,  0,   0,   0, N2,  0,   0} }); // == rot x

                    double dN1 = -1 / L;
                    double dN2 = 1 / L;
                    double dN3 = -6 * x / Math.Pow(L, 2) + 6 * Math.Pow(x, 2) / Math.Pow(L, 3);
                    double dN4 = 1 - 4 * x / L + 3 * Math.Pow(x, 2) / Math.Pow(L, 2);
                    double dN5 = -dN3;
                    double dN6 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * x / L;

                    Matrix<double> dN = DenseMatrix.OfArray(new double[,] {
                {dN1,   0,   0,   0,   0,   0, dN2,   0,   0,   0,    0,   0},
                {  0, dN3,   0,   0,   0, dN4,   0, dN5,   0,   0,    0,  dN6},
                {  0,   0, dN3,   0, -dN4,   0,   0,   0, dN5,  0, -dN6,   0},
                {  0,   0,   0, dN1,   0,   0,   0,   0,   0, dN2,    0,   0}});

                    //store the values in displacement and rotation lists
                    dispMatrix.SetRow(i, N.Multiply(u));   //Sets row 0,1,...n of the N matrix (needed for deformation in x,y and z)
                    rotMatrix.SetRow(i, dN.Multiply(u));  //Sets row 0,1,...n of the B matrix (needed to get rot_y and rot_z)

                    Vector<double> u_new = DenseVector.OfArray(new double[] //sets the new deformations in local coordinates
                    { dispMatrix[i,0], dispMatrix[i,1], dispMatrix[i,2], dispMatrix[i, 3], rotMatrix[i, 2], rotMatrix[i, 1] });

                    //Transposing from local to global
                    u_new = reduced_transmatrix.Transpose() * u_new; //transform  it back into global x, y and z coordinates

                    dispMatrix.SetRow(i, new double[] { u_new[0], u_new[1], u_new[2], u_new[3] }); //sets the displacements ux,uy,uz,theta_x globally in the correct row corresponding to location x
                    rotMatrix.SetRow(i, new double[] { rotMatrix[i, 0], u_new[5], u_new[4], rotMatrix[i, 3] });

                    List<double> rot = new List<double> { dispMatrix[i, 3], rotMatrix[i, 2], rotMatrix[i, 1] }; //rotx, roty, rotz

                    rotations.AddRange(rot);
                    newpoints[i] += new Point3d(dispMatrix[i, 0], dispMatrix[i, 1], dispMatrix[i, 2]);

                    x += L / n; //m
                }

                x = 0;
                List<double> momenty = new List<double>();
                List<double> momentz = new List<double>();

                for (int i = 0; i < n + 1; i++) //strain and stress calculation
                {
                    double dN1 = -1 / L;
                    double dN2 = 1 / L;
                    double dN3 = -6 * x / Math.Pow(L, 2) + 6 * Math.Pow(x, 2) / Math.Pow(L, 3);
                    double dN4 = 1 - 4 * x / L + 3 * Math.Pow(x, 2) / Math.Pow(L, 2);
                    double dN5 = -dN3;
                    double dN6 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * x / L;

                    Matrix<double> dN = DenseMatrix.OfArray(new double[,] {
                // ux,   uy,   uz,   rx    ry,   rz,   ux,   uy,   uz,   rx,   ry,   rz
                {dN1,   0,   0,   0,   0,   0, dN2,   0,   0,   0,   0,   0},
                {  0, dN3,   0,   0,   0, dN4,   0, dN5,   0,   0,   0,  dN6},
                {  0,   0, dN3,   0, -dN4,   0,   0,   0, dN5,   0, -dN6,   0},
                {  0,   0,   0, dN1,   0,   0,   0,   0,   0, dN2,   0,   0}});

                    double ddN1 = 0;
                    double ddN2 = 0;
                    double ddN3 = -6 / Math.Pow(L, 2) + 12 * x / Math.Pow(L, 3); // 
                    double ddN4 = -4 / L + 6 * x / Math.Pow(L, 2); // 
                    double ddN5 = 6 / Math.Pow(L, 2) - 12 * x / Math.Pow(L, 3); // 
                    double ddN6 = 6 * x / Math.Pow(L, 2) - 2 / L; // 

                    Matrix<double> B = DenseMatrix.OfArray(new double[,] {
                //ux,   uy,   uz,   rx,    ry,   rz,   ux,   uy,   uz,   rx,    ry,    rz
                { ddN1,   0,    0,    0,    0,    0,  ddN2,    0,    0,    0,     0,    0},
                {   0, ddN3,    0,    0,     0, ddN4,   0, ddN5,    0,    0,     0, ddN6},
                {   0,    0, ddN3,    0,  -ddN4,    0,    0,    0, ddN5,   0,  -ddN6,    0},
                {   0,    0,    0, ddN1,     0,    0,    0,    0,    0, ddN2,     0,    0}});


                    //store the values in displacement and rotation lists
                    var e = dN.Multiply(u);    //du_x, du_y, du_z, dtheta_x //epsilon, AXIAL STRAIN  //
                    var k = B.Multiply(u);  // 𝜺_x, 𝜺_y, 𝜺_z, 𝛾_xy, 𝛾_yz, 𝛾_zx

                    //i = 0 => sets row 0,1,2 with ex , ey, ez at x = 0
                    //i = 1 => sets row 3,4,5 with ex,ey,ex at x = L / 2
                    //each bar will have n+1 times rows with strains for each point

                    strains[i * 3] = e[0]; //ϵ_x
                    strains[i * 3 + 1] = z * k[2];
                    strains[i * 3 + 2] = -y * k[1];

                    stresses = strains * beam.E * 1e6; // kN/m2 => 𝝈_x, 𝝈_y, 𝝈_z
                    forces_e = stresses * beam.area * 1e-6; //kN

                    momenty.Add(beam.Iy * 1e-12 * beam.E * 1e6 * k[2]); //kNm
                    momentz.Add(beam.Iz * 1e-12 * beam.E * 1e6 * k[1]); //kNm

                    momenty_model.Add(beam.E * 1e6 * beam.Iy * 1e-12 * k[2]); //kNm
                    momentz_model.Add(beam.Iz * 1e-12 * beam.E * 1e6 * k[1]); //kNm

                    x += L / n;

                }

                max_momenty = momenty.MaximumAbsolute(); //kNm 
                max_momentz = momentz.MaximumAbsolute(); //kNm 

                max_momenty_model = momenty_model.MaximumAbsolute(); //kNm 
                max_momentz_model = momentz_model.MaximumAbsolute(); //kNm 

                for (int i = 0; i < forces_e.Count; i += 3) //taking only normal axial forces
                {
                    forces.Add(forces_e[i]);
                }

                if (momenty.IndexOf(max_momenty) == -1)
                {
                    max_momenty = -max_momenty;
                    beam.maxmoment = -max_momenty;
                }

                if (momentz.IndexOf(max_momentz) == -1)
                {
                    max_momentz = -max_momentz;
                }

                int max_momenty_index = momenty.IndexOf(max_momenty); //index of max moment (location at beam)
                Point3d max_momenty_x = oldpoints[max_momenty_index];

                max_momenty_x.X = Math.Round(max_momenty_x.X, 1);
                max_momenty_x.Y = Math.Round(max_momenty_x.Y, 1);
                max_momenty_x.Z = Math.Round(max_momenty_x.Z, 1);

                beam.maxmoment = max_momenty;
                beam.maxmomentpoint = max_momenty_x;
                id.Add(beam.ID);
                beam.strains = strains;
                beam.stresses = stresses;

                Curve curve = Curve.CreateControlPointCurve(newpoints, n);
                beam.deformed = curve;
                //Setting the beams displaced geometry
                displacedgeometry.Add(curve);
            }

            return displacedgeometry;

        }

    }
}
