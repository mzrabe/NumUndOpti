package org.mzrabe.opti;

import java.util.Arrays;

import org.mzrabe.lina.Function;
import org.mzrabe.lina.Gauss;
import org.mzrabe.lina.Vector;

/**
 * 
 * @author Moritz Zahn <mzrabe@gmail.com>
 *
 */
public class SimpleLocalNewton extends OptiAlgorithm{
	
	/**
	 * the direction vector 
	 */
	private double[] d;

	/**
	 * the local hessian matrix at the start point xs
	 */
	private double[][] localHess;
	
	/**
	 * 
	 * Initialize the optimization algorithm. This is the abstract class.
	 * 
	 * @param func
	 *            - the function where a minimum is searched
	 * @param tolerance
	 *            - the tolerance the function value
	 * @param maxNumIteration
	 *            - the maximal number of iterations
	 */
	public SimpleLocalNewton(Function func, double tolerance, int maxNumIteration) {
		super(func, tolerance, maxNumIteration);
	}
	
	/* (non-Javadoc)
	 * @see org.mzrabe.opti.OptiAlgorithm#initAlgorithms()
	 */
	@Override
	protected void initAlgorithms(double[] ... x) throws Exception {
		if(x.length > 1)
		{
			log.warn("You entered more then one start point. The " + getName() + " algorithms only needs one. "
					+ "Now the first point "+Arrays.toString(x[0])+" is used.");
		}
		
		localHess = func.getHf(x[0]);
		solution = Arrays.copyOf(x[0],x[0].length);
		d = Gauss.getSolution(localHess, func.grad(solution,true), false);
		
	}

	/* (non-Javadoc)
	 * @see org.mzrabe.opti.OptiAlgorithm#algorithms()
	 */
	@Override
	protected void algorithms() throws Exception {
		
		d = Gauss.getSolution(localHess, func.grad(solution,true), false);
		log.debug("d: " + Arrays.toString(d));
		solution = Vector.sum(solution, d);
		numberOfIterations++;
		
	}
	
	

	@Override
	protected boolean isFinish() {
		return Vector.twoNorm(d) > tolerance ? false : true;
	}
	
	@Override
	public String getName()
	{
		return "Simple Local Newton";
	}

	@Override
	public double[] getSolution()
	{
		return solution;
	}

	

}
