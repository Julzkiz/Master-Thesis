using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Gecko
{
    public class Load_Component : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Load_Component class.
        /// </summary>
        public Load_Component()
          : base("Load Component", "Nickname",
              "Defines a Load",
              "Masters", "Geometry")
        {
        }
            /// <summary>
            /// Registers all the input parameters for this component.
            /// </summary>
            protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
            {
                pManager.AddTextParameter("Name", "n", "", GH_ParamAccess.item, "deadload");
                pManager.AddPointParameter("Points", "p", "", GH_ParamAccess.list);
                pManager.AddVectorParameter("Load direction", "ld", "Load in kN", GH_ParamAccess.item);
                pManager.AddNumberParameter("Rotation at point", "rt", "Rotation in kNm (x,y,z)", GH_ParamAccess.list, new List<double> { 0, 0, 0 });
            }

            /// <summary>
            /// Registers all the output parameters for this component.
            /// </summary>
            protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
            {
                pManager.AddGenericParameter("Loads", "", "", GH_ParamAccess.list);
            }

            /// <summary>
            /// This is the method that actually does the work.
            /// </summary>
            /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
            protected override void SolveInstance(IGH_DataAccess DA)
            {
                string name = "";
                DA.GetData(0, ref name);

                List<Point3d> points = new List<Point3d>();
                DA.GetDataList(1, points);

                Vector3d ldir = new Vector3d();
                DA.GetData(2, ref ldir);

                List<double> rotations = new List<double>();
                DA.GetDataList(3, rotations);

                List<Load> loads = new List<Load>();

                int i = 0;
                foreach (Point3d point in points)
                {
                    Load load = new Load();
                    load.name = name;
                    load.ID = i;
                    load.point = point;
                    load.vector = ldir; //kN
                    load.rotations = rotations;
                    i++;
                    loads.Add(load);
                }

                DA.SetDataList(0, loads);
            }

            /// <summary>
            /// Provides an Icon for the component.
            /// </summary>
            /*
            protected override System.Drawing.Bitmap Icon
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
            get { return new Guid("E6B0EAAB-D115-42CC-82DF-08046929BA36"); }
        }
    }
}