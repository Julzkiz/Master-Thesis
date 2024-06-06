using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;
using Rhino.Geometry;

namespace Gecko
{
    public class ModelAnalysis_3DBar : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ModelAnalysis_3DBar class.
        /// </summary>
        public ModelAnalysis_3DBar()
          : base("3D Bar", "Nickname",
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
            pManager.AddNumberParameter("Max displacement", "", "in cm", GH_ParamAccess.item);
            pManager.AddGenericParameter("Calculated Model", "m", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Class_Model model = new Class_Model();
            DA.GetData(0, ref model);

            Vector<double> R4 = Model_Calculation_Stiffness_Matrix_3D_Bar.FEM_CALC_Displacement_4(model, out Matrix<double> K, out Matrix<double> M, out int d); //m
            List<double> Rr = R4.ToList();
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
                    continue;
                }

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

            foreach (Element bar in model.bars)
            {
                Point3d newpoint_s = new Point3d();
                Point3d newpoint_e = new Point3d();

                int nodeid_s = model.nodes.Find((node) => node.point.DistanceTo(bar.startnode) < 0.00003).globalID;
                int nodeid_e = model.nodes.Find((node) => node.point.DistanceTo(bar.endnode) < 0.00003).globalID;

                newpoint_s.X = bar.startnode.X + Rr[nodeid_s * d];
                newpoint_s.Y = bar.startnode.Y + Rr[nodeid_s * d + 1];
                newpoint_s.Z = bar.startnode.Z + Rr[nodeid_s * d + 2];

                newpoint_e.X = bar.endnode.X + Rr[nodeid_e * d];
                newpoint_e.Y = bar.endnode.Y + Rr[nodeid_e * d + 1];
                newpoint_e.Z = bar.endnode.Z + Rr[nodeid_e * d + 2];

                Line line = new Line(newpoint_s, newpoint_e);
                newpoints.Add(newpoint_s);
                newpoints.Add(newpoint_e);
                distances.Add(bar.startnode.DistanceTo(newpoint_s));
                distances.Add(bar.endnode.DistanceTo(newpoint_e));

                curve.Add(line);
            }

            newpoints = newpoints.Distinct().ToList();
            double maxdisp = distances.MaximumAbsolute() * 100; //m => cm

            Vector<double> displacements = Vector<double>.Build.DenseOfEnumerable(Rr);

            model.displacements = displacements;
            model.newgeometry = curve;

            DA.SetDataList(0, Rr);
            DA.SetDataList(1, curve);
            DA.SetData(2, maxdisp);
            DA.SetData(3, model);

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
            get { return new Guid("77DDD862-3FF4-414D-998B-A224421C8A4F"); }
        }
    }
}