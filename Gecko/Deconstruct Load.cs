using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Gecko
{
    public class Load_Deconstruct : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Load_Deconstruct class.
        /// </summary>
        public Load_Deconstruct()
          : base("Load Deconstruct", "Nickname",
              "Description",
              "Masters", "Disassembly")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Loads", "l", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "", "", GH_ParamAccess.list);
            pManager.AddVectorParameter("Vectors", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Load> loads = new List<Load>();
            DA.GetDataList(0, loads);

            List<Point3d> points = new List<Point3d>();
            List<Vector3d> vectors = new List<Vector3d>();


            foreach (Load load in loads)
            {
                points.Add(load.point);
                vectors.Add(load.vector);
            }

            DA.SetDataList(0, points);
            DA.SetDataList(1, vectors);
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
            get { return new Guid("2FB74E7C-9286-4704-9BCD-D06E54FCB185"); }
        }
    }
}