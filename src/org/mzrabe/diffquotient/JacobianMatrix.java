package org.mzrabe.diffquotient;

import org.mzrabe.lina.Function;

/**
 * Implementation to calculate the Jacobian matrix using the object {@link org.mzrabe.lina.Function}.
 * @author Moritz Zahn  
 */
public class JacobianMatrix
{
	/** the used gradient method */
	public static Gradient grad = new D1f();
	
	
//	/**
//	 * Get a new instance of the {@code JacobianMatrix} and specifies the used gradient 
//	 * method and the space between to point to calculate the gradient.
//	 * @param dh - the space between to points to calculate the gradient
//	 * @param grad - the gradient method
//	 */
//	public JacobianMatrix(double dh, Gradient grad)
//	{
//		this.grad = grad;
//		grad.setDh(dh);
//	}
//	
//	/**
//	 * Get a new instance of the {@code JacobianMatrix} and specifies dh.
//	 * @param dh - the space between to points to calculate the gradient
//	 */
//	public JacobianMatrix(double dh)
//	{
//		this.grad = new D1f(dh);
//	}
//	
//	/**
//	 * Get a new instance of the {@code JacobianMatrix}. Standard used gradient method is {@link D1f}.
//	 */
//	public JacobianMatrix()
//	{
//		this.grad = new D1f();
//	}
	
	/**
	 * Calculate the Jacobian matrix of the vectorial function at the given point x. 
	 * The additional parameter c are constances which will not consider by the gradient method.
	 * @param f - the vectorial function
	 * @param x - the point 
	 * @param c - constances of the function
	 * @return - the Jacobian matrix of f at x
	 */
	public static double[][] getJacobiMatrix(Function[] f, double[] x, double[] c)
	{
		double[][] jm = new double[f.length][x.length];
		
		for(int i=0;i<f.length;i++)
		{
			for(int j=0;j<x.length;j++)
			{
				try
				{
					jm[i][j] = grad.df(f[i], x, c, j);
				}
				catch (Exception e)
				{
					jm[i][j] = 0.0;
					e.printStackTrace();
				}
			}
		}
		
		return jm;
	}
	
	/**
	 * Calculate the Jacobian matrix of the vectorial function at the given point x. 
	 * @param f - the vectorial function
	 * @param x - the point 
	 * @return - the Jacobian matrix of f at x
	 */
	public static double[][] getJacobiMatrix(Function[] f, double[] x)
	{
		return getJacobiMatrix(f, x, null);
	}
}
