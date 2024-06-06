using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Gecko
{
    public class Model_Deconstruct : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Model_Deconstruct class.
        /// </summary>
        public Model_Deconstruct()
          : base("Model Deconstruct", "Nickname",
              "Description",
              "Masters", "Disassembly")
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
            pManager.AddGenericParameter("Bars", "b", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("Loads", "l", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("Supports", "s", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Class_Model model = new Class_Model();
            DA.GetData(0, ref model);

            List<Element> bars = new List<Element>();
            List<Load> loads = new List<Load>();
            List<Support> supports = new List<Support>();

            bars = model.bars;
            loads = model.loads;
            supports = model.supports;

            DA.SetDataList(0, bars);
            DA.SetDataList(1, loads);
            DA.SetDataList(2, supports);

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
            get { return new Guid("4C5E38EA-1B76-4BE1-A7F3-BD65902C6340"); }
        }
    }
}