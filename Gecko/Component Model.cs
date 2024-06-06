using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Gecko
{
    public class Model_Component : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Model_Component class.
        /// </summary>
        public Model_Component()
          : base("Model Component", "Nickname",
              "Assembles the Model",
              "Masters", "Geometry")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Elements", "", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("Nodes", "", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("Supports", "", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("Loads", "", "", GH_ParamAccess.list);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Element> bars = new List<Element>();
            List<Support> supports = new List<Support>();
            List<Load> loads = new List<Load>();
            List<Node> nodes = new List<Node>();

            DA.GetDataList(0, bars);
            DA.GetDataList(1, nodes);
            DA.GetDataList(2, supports);
            DA.GetDataList(3, loads);

            Class_Model model = new Class_Model();

            model.bars = bars;
            model.loads = loads;
            model.supports = supports;
            model.nodes = nodes;

            DA.SetData(0, model);
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
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        */

        public override Guid ComponentGuid
        {
            get { return new Guid("4A638F64-D571-4F86-99E8-C2403E20D4FC"); }
        }
    }
}