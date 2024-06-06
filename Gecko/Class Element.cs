using MathNet.Numerics.LinearAlgebra;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Gecko
{
    internal class Element
    {
        public string name;
        public int ID;
        public double area;
        public double heigth;
        public double width;
        public double E;
        public double G;
        public double J;
        public double Iy;
        public double Iz;

        public Line axis;
        public Point3d startnode;
        public Point3d endnode;
        public List<Node> nodes;
        public double length;

        public Vector<double> strains;
        public Vector<double> stresses;
        public double stress;

        public Point3d maxmomentpoint;
        public double maxmoment;
        public Matrix<double> transmatrix;
        public Curve deformed;
    }
}
