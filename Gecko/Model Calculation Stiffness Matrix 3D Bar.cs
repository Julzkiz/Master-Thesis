using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Gecko
{
    internal class Model_Calculation_Stiffness_Matrix_3D_Bar
    {
        public static Vector<double> FEM_CALC_Displacement_4(Class_Model model, out Matrix<double> K, out Matrix<double> M, out int d)
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
            K = DenseMatrix.Create(d * 2, d * 2, 0);

            List<double> L_vector = new List<double>(); //creating a vector with correct size
            //Vector<double> V = Vector<double>.Build.Dense(DOFs.Count * d);
            for (int i = 0; i < Msize; i++)
            {
                L_vector.Add(0);
            }

            foreach (Load load in model.loads) //maybe just use V.Dense(DOFs.Count*d) -> vector of x length
            {
                int nodeid = model.nodes.Find((node) => node.point.DistanceTo(load.point) < 0.0003).globalID;
                L_vector[nodeid * d] = load.vector.X; //kN 
                L_vector[nodeid * d + 1] = load.vector.Y; //kN
                L_vector[nodeid * d + 2] = load.vector.Z; //kN
            }

            int startID = 0;
            int endID = 0;

            foreach (Element bar in model.bars) //go through every bar and assign them a stiffness, start and end, and length.
            {
                double k = (bar.area * bar.E) / (bar.length * 1e3); //mm^2 * kN/mm^2 / mm => kN/mm

                int s_nodeid = model.nodes.Find((node) => node.point.DistanceTo(bar.axis.From) < 0.0003).globalID;
                int e_nodeid = model.nodes.Find((node) => node.point.DistanceTo(bar.axis.To) < 0.0003).globalID;

                startID = s_nodeid * d;
                endID = e_nodeid * d;

                //3d truss element [u1,v1,w1,u2,v2,w2]
                Matrix<double> m = DenseMatrix.OfArray(new double[,]
                {{ k, 0, 0, -k, 0, 0},
                 { 0, 0, 0,  0, 0, 0},
                 { 0, 0, 0,  0, 0, 0},
                 {-k, 0, 0,  k, 0, 0},
                 { 0, 0, 0,  0, 0, 0},
                 { 0, 0, 0,  0, 0, 0}});

                //Assuming the transformation matrix is going in here
                double lx = (bar.axis.To.X - bar.axis.From.X) / bar.length;
                double my = (bar.axis.To.Y - bar.axis.From.Y) / bar.length;
                double nz = (bar.axis.To.Z - bar.axis.From.Z) / bar.length;

                Matrix<double> t = DenseMatrix.OfArray(new double[,] {
                {lx, nz, nz,   0,  0,  0 },
                {lx, my, nz,   0,  0,  0 },
                {lx, my, nz,   0,  0,  0 },
                { 0,  0,  0,  lx, my, nz },
                { 0,  0,  0,  lx, my, nz },
                { 0,  0,  0,  lx, my, nz }});

                Matrix<double> tT = t.Transpose();
                Matrix<double> T = tT * m * t;

                K = T;

                //now to put the global stiffness matrix K to the correct place in M.
                for (int i = 0; i < d; i++)
                {
                    for (int j = 0; j < d; j++)
                    {
                        M[startID + i, startID + j] += K[i, j];
                        M[startID + i, endID + j] += K[i, d + j];
                        M[endID + i, startID + j] += K[d + i, j];
                        M[endID + i, endID + j] += K[d + i, d + j];
                    }
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

            Vector<double> R = M.Solve(v); // mm/kN * kN => mm
            R = R / 1000; //mm -> m

            return R;
        }

    }
}
