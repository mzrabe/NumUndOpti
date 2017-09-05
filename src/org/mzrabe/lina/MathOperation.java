package org.mzrabe.lina;

/**
 *	Simple mathematics operation. A function from Rn -> R1.
 */
public interface MathOperation
{
	
	/**
	 * Define the a function with the argument x and optional an number of coefficients s. 
	 * The parameters s will not takes into account for the calculation of gradient
	 * or the hessian matrix of the function.
	 * 
	 * @param x - the argument of the function
	 * @param c - a number of coefficients (optional)
	 * @return - the value of the function at x
	 * @throws Exception 
	 */
	public abstract double getValue(double[] x,double ... c) throws Exception;

}
