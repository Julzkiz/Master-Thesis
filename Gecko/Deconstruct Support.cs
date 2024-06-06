using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Gecko
{
    public class Support_Deconstruct : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Support_Deconstruct class.
        /// </summary>
        public Support_Deconstruct()
          : base("Support Deconstruct", "Nickname",
              "Description",
              "Masters", "Disassembly")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Supports", "s", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Name", "n", "", GH_ParamAccess.list);
            pManager.AddNumberParameter("ID", "i", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Support> supports = new List<Support>();
            DA.GetDataList(0, supports);

            List<String> names = new List<String>();
            List<int> ids = new List<int>();

            foreach (Support support in supports)
            {
                names.Add(support.name);
                ids.Add(support.ID);
            }

            DA.SetDataList(0, names);
            DA.SetDataList(1, ids);
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
            get { return new Guid("BE28C84A-295E-4956-9FCA-ECD014D61C5C"); }
        }
    }
}