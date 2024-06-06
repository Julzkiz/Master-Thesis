using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Gecko
{
    public class Support_Component : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Support_Component class.
        /// </summary>
        public Support_Component()
          : base("Support Component", "Nickname",
              "Defining Boundary Conditions",
              "Masters", "Geometry")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Name", "n", "name of support", GH_ParamAccess.item, "Fixed support");
            pManager.AddPointParameter("Points", "p", "points", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Boundary Conditions", "bc", "list of 1s & 0s", GH_ParamAccess.list, new List<int>() { 1, 1, 1, 0, 0, 0 });
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Supports", "s", "support points", GH_ParamAccess.list);
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

            List<int> bc = new List<int>();
            DA.GetDataList(2, bc);

            List<Support> supports = new List<Support>();

            bool boundary = Convert.ToBoolean(1); //a way to convert the list of 1s and 0s to true and false statements later.

            int i = 0;
            foreach (Point3d point in points)
            {
                Support support = new Support();
                support.name = name;
                support.ID = i;
                support.point = point;
                support.bcds = bc;
                supports.Add(support);
                i++;
            }

            DA.SetDataList(0, supports);
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
    get { return new Guid("F5028BEB-5736-4C1A-827B-EDDE5588C3E6"); }
}
}
}