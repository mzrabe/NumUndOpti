package org.mzrabe.approximation;

import org.mzrabe.diffquotient.D1f;
import org.mzrabe.diffquotient.Gradient;
import org.mzrabe.lina.Function;
import org.mzrabe.lina.Gauss;
import org.mzrabe.lina.Matrix;
import org.mzrabe.lina.Vector;

public abstract class leastSquares {
	
	/** the used gradient method */
	public static Gradient grad = new D1f();

	/**
	 * @param f - the definition of the function
	 * @param x - start values of the coefficients
	 * @param y - the solution of the function f<sub>i</sub>(values[i])
	 * @param values - a list of the arguments arguments
	 * @return - the coefficients of the functions which get the best fit of the function
	 * @throws Exception
	 */
	public static double[] approx(Function f, double[] x, double[] y, double[] ... values) throws Exception {

		/* generate the matrix C */
		double[][] C = new double[values.length][x.length];
		for(int i=0;i<C.length;i++)
		{
			for(int j = 0;j<x.length;j++)
			{
				C[i][j] = grad.df(f, x, values[i], j);
			}
		}
		
//		Matrix.printMatix(C);
		
		/* calculate d */
		double[] d = new double[C.length];
		for(int i=0;i<values.length;i++)
		{
			d[i] = -f.getValue(x, values[i]) + y[i];
		}
			
		
		/* calculate matrix A for the normal from */
		double[][] A = Matrix.multi(Matrix.trans(C),C);
//		Matrix.printMatix(A);
		d = Matrix.multi(Matrix.trans(C), d);
//		Vector.print(d);
		double[] r = Gauss.getSolution(A, d,false);
		double[] rBefor;
		
		do
		{
//			Vector.print(r);
			rBefor = Vector.multiScalar(r, 1);
			x = Vector.sum(x, r);
			
			C = new double[values.length][x.length];
			for(int i=0;i<C.length;i++)
			{
				for(int j = 0;j<x.length;j++)
				{
					C[i][j] = grad.df(f, x, values[i], j);
				}
			}
			
			/* calculate d */
			d = new double[C.length];
			for(int i=0;i<values.length;i++)
			{
				d[i] = -f.getValue(x, values[i]) + y[i];
			}
			
			/* calculate matrix A for the normal from */
			A = Matrix.multi(Matrix.trans(C),C);
			d = Matrix.multi(Matrix.trans(C), d);
			r = Gauss.getSolution(A, d,false);
		}while(Vector.oneNorm(Vector.minus(r, rBefor)) >= 0.0000001);
		
		return x;
	}

}
