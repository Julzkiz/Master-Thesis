using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Gecko
{
    public class Element_Deconstruct : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Element_Deconstruct class.
        /// </summary>
        public Element_Deconstruct()
          : base("Element Deconstruct", "Nickname",
              "Description",
              "Masters", "Disassembly")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Element", "b", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Name", "n", "name of element", GH_ParamAccess.list);
            pManager.AddLineParameter("Line", "l", "line of element", GH_ParamAccess.list);
            pManager.AddNumberParameter("ID", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Element> bars = new List<Element>();
            DA.GetDataList(0, bars);

            List<String> names = new List<String>();
            List<Line> lines = new List<Line>();
            List<int> ids = new List<int>();

            foreach (Element bar in bars)
            {
                names.Add(bar.name);
                ids.Add(bar.ID);
                double l = bar.length;
                lines.Add(bar.axis);
            }

            DA.SetDataList(0, names);
            DA.SetDataList(1, lines);
            DA.SetDataList(2, ids);

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
        }
        */
        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("21026A02-EFA5-42A0-92E9-C0EEB9BB1D0B"); }
        }
    }
}