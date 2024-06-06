using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using Rhino.Geometry;
using static System.Runtime.InteropServices.JavaScript.JSType;
using System.IO;
using System.Net;
using System.Linq;
using System.Security.Cryptography;
using Rhino.Display;
using System.Resources;
using MathNet.Numerics;
using GH_IO.Serialization;
namespace Gecko
{
    public class Force_Analysis : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Force_Analysis class.
        /// </summary>
        public Force_Analysis()
          : base("Force Analysis", "Nickname",
              "Description",
              "Masters", "Analysis")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Calculated Model", "", "Analysing the calculated model", GH_ParamAccess.item);
            pManager.AddIntegerParameter("n", "", "Discretising the elements in n-parts", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Strains", "", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("Stresses", "", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("Normal Forces", "", "kN", GH_ParamAccess.list);
            pManager.AddGenericParameter("My", "", "kN", GH_ParamAccess.list);
            pManager.AddGenericParameter("Mz", "", "kN", GH_ParamAccess.list);
            pManager.AddCurveParameter("Displaced Geometry", "", "", GH_ParamAccess.list);
            pManager.AddTextParameter("Info", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Class_Model calc_model = new Class_Model();
            DA.GetData(0, ref calc_model);

            int n = 0;
            DA.GetData(1, ref n);

            List<double> d_xdir = new List<double>();
            for (int i = 0; i < calc_model.displacements.Count; i++)
            {
                d_xdir.Add(calc_model.displacements[i]);
                i++;
            }

            List<double> strains = new List<double>();
            List<double> stresses = new List<double>();
            List<double> forces = new List<double>();
            double max_momenty = 0;
            double max_momentz = 0;

            List<string> text = new List<string>();

            //4-dof bar element
            foreach (Element bar in calc_model.bars)
            {
                double originLength = bar.axis.Length;
                double deformedLength = calc_model.newgeometry[bar.ID].Length;

                double dL = originLength - deformedLength;

                double e = dL / originLength;
                strains.Add(e);

                double s = e * bar.E;
                stresses.Add(s);

            }

            List<Curve> displacedgeometry = new List<Curve>();

            //4-dof beam element
            if (calc_model.d == 3)
            {   //update with out strain, out stress
                List<Curve> displacedgeometry_2D = Model_Calculation_Stiffness_Matrix_2D_Beam.InterpolatedGeometry_2D(calc_model, n, out double strain, out double stress, out double force, out double momenty, out List<int> ids);
                strains = new List<double> { strain };
                stresses = new List<double> { stress };
                forces = new List<double> { force };
                max_momenty = momenty;
                max_momentz = 0;

                foreach (int id in ids)
                {
                    text.Add("Max moment: " + calc_model.bars[id].maxmoment.Round(2).ToString() + " kNm at Point (" + calc_model.bars[id].maxmomentpoint.ToString() + ")");
                }

                displacedgeometry = displacedgeometry_2D;
            }

            else if (calc_model.d == 6)
            {   //update with out strain, out stress
                List<Curve> displacedgeometry_3D = Model_Calculation_Stiffness_Matrix_3D_Beam.InterpolatedGeometry3D(calc_model, n, out Vector<double> strain, out Vector<double> stress, out List<double> force, out double momenty, out double momentz, out List<int> ids);
                strains = strain.ToList();
                stresses = stress.ToList();
                forces = force;
                max_momenty = momenty;
                max_momentz = momentz;

                foreach (int id in ids)
                {
                    text.Add("Max moment: " + calc_model.bars[id].maxmoment.Round(2).ToString() + " kNm at Point (" + calc_model.bars[id].maxmomentpoint.ToString() + ")");
                }

                text = text.Distinct().ToList();
                displacedgeometry = displacedgeometry_3D;
            }

            //6-dof beam element

            DA.SetDataList(0, strains);
            DA.SetDataList(1, stresses);
            DA.SetDataList(2, forces);
            DA.SetData(3, max_momenty);
            DA.SetData(4, max_momentz);
            DA.SetDataList(5, displacedgeometry);
            DA.SetDataList(6, text);

        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("E9581D05-D593-44E6-8F7B-38B2FE9DABC3"); }
        }
    }
}