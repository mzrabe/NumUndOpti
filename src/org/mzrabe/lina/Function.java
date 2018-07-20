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
	 * @throws Exception 
	 */
	public double[] grad(double[] x, boolean negativ) throws Exception {
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
	
	/**
	 * Method which calculate the solution of the equation system.
	 * @param f - equation system
	 * @param x - arguments of f
	 * @param c - constant arguments of f
	 * @return - the solution of f(x,c) as double[] array
	 * @throws Exception - if a calculation of a value f<sub>i</sub>(x,c) fails
	 *
	 */
	public static double[] getValues(Function[] f, double[] x, double ... c) throws Exception
	{
		double[] back = new double[f.length];
		for(int i=0;i<f.length;i++)
		{
			back[i] = f[i].getValue(x, c);
		}
		return back;
	}
	
	/**
	 * Get the integration of this function from the value x1 to x2. Same as {@link Function#integrate(double[], double[], double...)} with x1 as double[]{x1} and x2 as double[]{x2}.
	 * @param x1 - start parameter for the integration
	 * @param x2 - end parameter for the integration
	 * @param c - constant parameters of this function
	 * @return - solution of the integration from x1 to x2
	 * @throws Exception - if integration fails
	 */
	public double integrate(double x1, double x2, double ... c) throws Exception
	{
		return integrate(new double[] {x1}, new double[] {x2}, c);
	}
	
	/**
	 * Get the integration of this function from the value x1 to x2.
	 * @param x1 - start parameter for the integration
	 * @param x2 - end parameter for the integration
	 * @param c - constant parameters of this function
	 * @return - solution of the integration from x1 to x2
	 * @throws Exception - if integration fails
	 */
	public double integrate(double[] x1, double[] x2, double ... c) throws Exception
	{
		int n = 1;
		double[] h = new double[x1.length];
		for(int i = 0;i<h.length;i++)
		{
			h[i] = x2[i]-x1[i];
		}
		
		double d = h[0];
		
		for(int i = 1;i<h.length;i++)
		{
			d*=h[i];
		}
		
		double p = 0;
		int numP = (int) Math.pow(2,h.length);
		int[] comb = new int[h.length];
		
		// 0 0 0
		// 1 1 1
		// 000, 100, 101, 001, 010, 110, 011, 111
		// 000 001 010 011 100 101 110 111
		
		double temp[] = new double[x1.length];
		
		for(int i = 0;i<numP;i++)
		{
			for(int j = 0;j<comb.length;j++)
			{
				temp[j] = comb[j] == 0 ? x1[j] : x2[j];
			}
			
			p+=	getValue(temp, c);
			
			if(i==numP-1)
				break;
			
			if(i%2==0)
			{
				for(int j = comb.length-1;j>=0;j--)
				{
					if(comb[j] == 0)
					{
						comb[j] = 1;
						break;
					}
				}
			}
			else
			{
				for(int j = comb.length-1;j>=0;j--)
				{
					if(comb[j] == 0)
					{
						comb[j] = 1;
						break;
					}
					comb[j] = 0;
				}
			}
		}
		
		double T = d * p/(numP);
		
		for(int k = 0;k<10;k++)
		{
			double M = 0;
			for(int j = 0;j<n;j++)
			{
				M += getValue(Vector.sum(x1,Vector.multiScalar(h, j+0.5)), c);
			}
			M*=h[0];
			T=(T+M)/2.;
			for(int i=0;i<h.length;i++)
			{
				h[i]/=2.;
			}
			n*=2.;
			if(Math.abs(T-M) <= 0.01)
			{
//				System.out.println("n = "+n);
//				System.out.println("h = "+h[0]);
				break;
			}
		}
//		System.out.println(String.format("T = %f", T));
		return T;
	}
}
