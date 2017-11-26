package org.mzrabe.lina;

import org.mzrabe.diffquotient.D1f;
import org.mzrabe.diffquotient.Gradient;
import org.mzrabe.diffquotient.HessianMatrix;

/**
 * @author Moritz Zahn
 */
public abstract class Function implements MathOperation{
	
	/**
	 * the used 
	 */
	private Gradient grad;
	private HessianMatrix Hf;
	
	/**
	 * Get a new instance of a function. As default {@link D1f} is used as gradient method and {@link HessianMatrix} for the Hessian matrix.
	 */
	public Function()
	{
		grad 	=	new D1f();
		Hf 		= 	new HessianMatrix();
	}

	/**
	 * Get the partial derivative of the function after x[n].
	 * @param x - derivation at the point x
	 * @param c - the constanc values of the function
	 * @param n - the index of the variable
	 * @return - the partial derivation of the function a the point x
	 * @throws Exception - 
	 */
	public double df( double[] x, double[] c, int... n) throws Exception {
		return grad.df(this, x, c, n);
	}
	
	/**
	 * Get the partial derivative of the function after x[0].
	 * @param x - derivation at the point x
	 * @param c - the constanc values of the function
	 * @return - the partial derivation of the function a the point x
	 * @throws Exception -
	 */
	public double df( double[] x, double ... c) throws Exception {
		return grad.df(this, x, c, 0);
	}

	/**
	 * The gradient of the function at the point x
	 * @param x - the point 
	 * @param negativ - should the gradient have the opposite direction (to the lowest descent)
	 * @return - the gradient
	 */
	public double[] grad(double[] x, boolean negativ) {
		return grad.grad(this, x, negativ);
	}
	
	/**
	 * Get the Hessian matrix of the function at the point x
	 * @param x - the point
	 * @param c - the constances of the function
	 * @return - the hessian matrix
	 */
	public double[][] getHf(double[] x, double ... c)
	{
		return Hf.getHf(this, x, c);
	}
	
	/**
	 * Set the step size to calculate the gradient.
	 * @param dh - the step size see {@link Gradient}
	 */
	public void setDh(double dh)
	{
		grad.setDh(dh);
	}
}
