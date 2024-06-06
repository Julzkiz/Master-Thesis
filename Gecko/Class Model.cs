using MathNet.Numerics.LinearAlgebra;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Gecko
{
    internal class Class_Model
    {
        public List<Element> bars;
        public List<Load> loads;
        public List<Support> supports;
        public List<Node> nodes;
        public Vector<double> displacements;
        public List<Line> newgeometry;
        public Matrix<double> transmatrix;
        public int d;

    }
}