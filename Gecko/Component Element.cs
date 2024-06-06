using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Collections;
using Rhino.Geometry;

namespace Gecko
{
    public class Element_Component : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Element_Component class.
        /// </summary>
        public Element_Component()
          : base("Element Component", "Nickname",
              "Creates Elements",
              "Masters", "Geometry")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Name", "n", "name of element", GH_ParamAccess.item, "Default Element");
            pManager.AddLineParameter("Lines", "l", "list of lines [m]", GH_ParamAccess.list);
            pManager.AddTextParameter("Section", "s", "Section of bar, wxh", GH_ParamAccess.item, "100x100");
            pManager.AddTextParameter("Material", "m", "Material of bar", GH_ParamAccess.item, "Steel");
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Elements", "b", "list of elements", GH_ParamAccess.list);
            pManager.AddGenericParameter("Nodes", "n", "nodes in model", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string name = "";
            DA.GetData(0, ref name);

            List<Element> elements = new List<Element>();

            List<Line> lines = new List<Line>();
            DA.GetDataList(1, lines);

            string section = "";
            DA.GetData(2, ref section);

            string material = "";
            DA.GetData(3, ref material);

            //Can add a different section definition eventually to include IPE's
            string[] section_split = section.Split('x');
            double sectionsize0 = Convert.ToDouble(section_split[0]); //w
            double sectionsize1 = Convert.ToDouble(section_split[1]); //h

            //Assigning ID, name and area to each bar and beam in the list
            int i = 0;
            foreach (Line line in lines)
            {
                Element bar = new Element();
                bar.name = name + section;
                bar.ID = i;
                bar.area = sectionsize0 * sectionsize1; //mm2
                bar.heigth = sectionsize1;
                bar.width = sectionsize0;
                bar.E = 210; //GPA = kN/mm2
                bar.G = bar.E / (2 * (1 + 0.3)); //kN/mm2
                bar.Iy = (sectionsize0 * Math.Pow(sectionsize1, 3)) / 12; //mm4
                bar.Iz = (sectionsize1 * Math.Pow(sectionsize0, 3)) / 12; //mm4
                bar.J = sectionsize0 * sectionsize1 * (Math.Pow(sectionsize0, 2) + Math.Pow(sectionsize1, 2)) / 12; //mm4
                bar.axis = line;
                bar.startnode = line.From;
                bar.endnode = line.To;
                bar.length = line.Length; //m
                elements.Add(bar);

                i++;
            }

            List<Node> nodes = new List<Node>();

            DataTree<Point3d> l_ID = new DataTree<Point3d>();
            List<Point3d> g_ID = new List<Point3d>();

            int j = 0;
            foreach (Element bar in elements)
            {
                Node start = new Node();
                Node end = new Node();

                int decimalplaces = 3;

                Line barline = bar.axis;
                start.point = barline.From;
                end.point = barline.To;

                //Rounding points for global indexing due to some decimal discrepancies.
                Point3d point_start = new Point3d(Math.Round(start.point.X, decimalplaces), Math.Round(start.point.Y, decimalplaces), Math.Round(start.point.Z, decimalplaces));
                Point3d point_end = new Point3d(Math.Round(end.point.X, decimalplaces), Math.Round(end.point.Y, decimalplaces), Math.Round(end.point.Z, decimalplaces));

                if (!g_ID.Contains(point_start)) //Indexing all the new points globally
                {
                    g_ID.Add(point_start);
                    start.globalID = j * 2;
                }

                if (!g_ID.Contains(point_end))
                {
                    g_ID.Add(point_end);
                    end.globalID = 1 + j * 2;
                }

                //Indexing the points locally for each bar.
                GH_Path path = new GH_Path(i);
                List<Point3d> points = new List<Point3d>{
                    start.point,
                    end.point};

                l_ID.AddRange(points, path);
                start.localID = 0;
                end.localID = 1;

                nodes.Add(start);
                nodes.Add(end);
                bar.nodes = nodes;

                j++;
            }

            //Check which node in the g_ID list this node is closest to, assign that g_Index
            foreach (Node node in nodes)
            {
                int prev_id = Point3dList.ClosestIndexInList(g_ID, node.point);
                node.globalID = prev_id;
            }

            DA.SetDataList(0, elements);
            DA.SetDataList(1, nodes);
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
        ///         */
        public override Guid ComponentGuid
{
    get { return new Guid("0579C8F6-0494-4456-89F6-1F7FC10F1D7C"); }
}
}
}