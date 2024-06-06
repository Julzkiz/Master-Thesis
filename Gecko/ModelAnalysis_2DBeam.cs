using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;
using Rhino.Geometry;

namespace Gecko
{
    public class ModelAnalysis_2DBeam : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ModelAnalysis_2DBeam class.
        /// </summary>
        public ModelAnalysis_2DBeam()
          : base("2D Beam", "Nickname",
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

            Vector<double> R6 = Model_Calculation_Stiffness_Matrix_2D_Beam.FEM_CALC_Displacement_6(model, out Matrix<double> t, out Matrix<double> K, out Matrix<double> M, out int d); //m
            List<double> R6_int = R6.ToList();
            List<double> Rr = new List<double>();
            for (int i = 0; i < R6.ToList().Count; i++)
            {
                if ((i + 1) % 3 == 0)
                {
                    continue;
                }
                Rr.Add(R6.ToList()[i]);
            }


            List<int> node_ints = new List<int>();

            foreach (Support support in model.supports)
            {
                node_ints.Add(model.nodes.Find((node) => node.point.DistanceTo(support.point) < 0.00003).globalID);
            }

            node_ints.Sort();

            for (int i = 0; i < node_ints.Count; i++)
            {
                Rr.Insert(node_ints[i] * 2, 0);
                Rr.Insert(node_ints[i] * 2 + 1, 0);

                R6_int.Insert(node_ints[i] * 3, 0);
                R6_int.Insert(node_ints[i] * 3 + 1, 0);
                R6_int.Insert(node_ints[i] * 3 + 2, 0);
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

                newpoint_s.X = beam.startnode.X + Rr[nodeid_s * 2];
                newpoint_s.Z = beam.startnode.Z + Rr[nodeid_s * 2 + 1];

                newpoint_e.X = beam.endnode.X + Rr[nodeid_e * 2];
                newpoint_e.Z = beam.endnode.Z + Rr[nodeid_e * 2 + 1];

                Line line = new Line(newpoint_s, newpoint_e);
                newpoints.Add(newpoint_s);
                newpoints.Add(newpoint_e);
                distances.Add(beam.startnode.DistanceTo(newpoint_s));
                distances.Add(beam.endnode.DistanceTo(newpoint_e));

                curve.Add(line);
            }

            newpoints = newpoints.Distinct().ToList();
            double maxdisp = distances.MaximumAbsolute() * 100; //m => cm

            Vector<double> displacements = Vector<double>.Build.DenseOfEnumerable(R6_int);

            model.displacements = displacements * 1e3; //m => mm
            model.newgeometry = curve;
            model.transmatrix = t;
            model.d = d;

            DA.SetDataList(0, R6_int);
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
            get { return new Guid("DF42AEF3-8444-49A6-8420-C22CEAACF0A0"); }
        }
    }
}