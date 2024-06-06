using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;
using Rhino.Geometry;

namespace Gecko
{
    public class ModelAnalysis_3DBeam : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ModelAnalysis_3DBeam class.
        /// </summary>
        public ModelAnalysis_3DBeam()
          : base("3D Beam", "Nickname",
              "Description",
              "Masters", "Analysis")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "m", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Displacements", "disp", "in x and z direction", GH_ParamAccess.list);
            pManager.AddCurveParameter("Displaced geometry", "disp", "in x and z direction", GH_ParamAccess.list);
            pManager.AddNumberParameter("Max displacement", "", "in mm", GH_ParamAccess.item);
            pManager.AddGenericParameter("Calculated Model", "", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Class_Model model = new Class_Model();
            DA.GetData(0, ref model);

            List<double> R12 = Model_Calculation_Stiffness_Matrix_3D_Beam.FEM_CALC_Displacement_12(model, out Matrix<double> t, out Matrix<double> K, out Matrix<double> M, out int d); //mm

            d = d / 2;

            List<double> Rr = SkipEveryInList(R12, 3);
            List<int> node_ints = new List<int>();

            foreach (Support support in model.supports)
            {
                node_ints.Add(model.nodes.Find((node) => node.point.DistanceTo(support.point) < 0.00003).globalID);
            }

            node_ints.Sort();

            for (int i = 0; i < node_ints.Count; i++)
            {
                if (node_ints[i] * d > Rr.Count)
                {
                    Rr.Add(0);
                    Rr.Add(0);
                    Rr.Add(0);

                    R12.Add(0);
                    R12.Add(0);
                    R12.Add(0);
                    R12.Add(0);
                    R12.Add(0);
                    R12.Add(0);

                    continue;
                }

                R12.Insert(node_ints[i] * (d * 2), 0);
                R12.Insert(node_ints[i] * (d * 2) + 1, 0);
                R12.Insert(node_ints[i] * (d * 2) + 2, 0);
                R12.Insert(node_ints[i] * (d * 2) + 3, 0);
                R12.Insert(node_ints[i] * (d * 2) + 4, 0);
                R12.Insert(node_ints[i] * (d * 2) + 5, 0);

                Rr.Insert(node_ints[i] * d, 0);
                Rr.Insert(node_ints[i] * d + 1, 0);
                Rr.Insert(node_ints[i] * d + 2, 0);
            }

            List<Point3d> oldpoints = new List<Point3d>();
            List<Point3d> newpoints = new List<Point3d>();

            foreach (Node node in model.nodes)
            {
                oldpoints.Add(node.point);
            }

            oldpoints = oldpoints.Distinct().ToList();
            List<Line> curve = new List<Line>();
            List<double> distances = new List<double>();

            foreach (Element beam in model.bars)
            {
                Point3d newpoint_s = new Point3d();
                Point3d newpoint_e = new Point3d();

                int nodeid_s = model.nodes.Find((node) => node.point.DistanceTo(beam.startnode) < 0.00003).globalID;
                int nodeid_e = model.nodes.Find((node) => node.point.DistanceTo(beam.endnode) < 0.00003).globalID;

                newpoint_s.X = beam.startnode.X + Rr[nodeid_s * d];
                newpoint_s.Y = beam.startnode.Y + Rr[nodeid_s * d + 1];
                newpoint_s.Z = beam.startnode.Z + Rr[nodeid_s * d + 2];

                newpoint_e.X = beam.endnode.X + Rr[nodeid_e * d];
                newpoint_e.Y = beam.endnode.Y + Rr[nodeid_e * d + 1];
                newpoint_e.Z = beam.endnode.Z + Rr[nodeid_e * d + 2];

                Line line = new Line(newpoint_s, newpoint_e);
                newpoints.Add(newpoint_s);
                newpoints.Add(newpoint_e);
                distances.Add(beam.startnode.DistanceTo(newpoint_s));
                distances.Add(beam.endnode.DistanceTo(newpoint_e));

                curve.Add(line);
            }

            newpoints = newpoints.Distinct().ToList();
            double maxdisp = distances.MaximumAbsolute() * 100; //m => cm

            Vector<double> displacements = Vector<double>.Build.DenseOfEnumerable(R12);

            model.displacements = displacements; //m
            model.newgeometry = curve;
            model.transmatrix = t;
            model.d = d * 2;

            DA.SetDataList(0, R12);
            DA.SetDataList(1, curve);
            DA.SetData(2, maxdisp);
            DA.SetData(3, model);
        }

        static List<double> SkipEveryInList(List<double> inputList, int groupSize)
        {
            List<double> result = new List<double>();

            int currentIndex = 0;
            while (currentIndex < inputList.Count)
            {
                List<double> group = new List<double>();

                // Store the next group of elements
                for (int i = currentIndex; i < Math.Min(currentIndex + groupSize, inputList.Count); i++)
                {
                    group.Add(inputList[i]);
                }

                // Skip the next 3 indexes
                currentIndex += groupSize + 3;

                result.AddRange(group);
            }

            return result;
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        /*protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }*/

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("A4A65020-3558-4336-B4B8-264BA4955349"); }
        }
    }
}