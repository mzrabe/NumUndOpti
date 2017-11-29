package org.mzrabe.opti;

import java.util.Arrays;

import org.mzrabe.lina.Function;
import org.mzrabe.lina.Gauss;
import org.mzrabe.lina.Vector;

public class LocalNewton extends OptiAlgorithm{

	/**
	 * the direction vector 
	 */
	private double[] d;
	
	public LocalNewton(Function func, double tolerance, int maxNumIteration) {
		super(func, tolerance, maxNumIteration);
	}
	
	/* (non-Javadoc)
	 * @see org.mzrabe.opti.OptiAlgorithm#initAlgorithms()
	 */
	@Override
	protected void initAlgorithms(double[] ... x) {
		if(x.length > 1)
		{
			log.warn("You entered more then one start point. The " + getName() + " algorithms only needs one. "
					+ "Now the first point "+Arrays.toString(x[0])+" is used.");
		}
		
		solution = Arrays.copyOf(x[0],x[0].length);
		d = Gauss.getSolution(func.getHf(solution), func.grad(solution,true), false);
		
	}

	/* (non-Javadoc)
	 * @see org.mzrabe.opti.OptiAlgorithm#algorithms()
	 */
	@Override
	protected void algorithms() {
		
		d = Gauss.getSolution(func.getHf(solution), func.grad(solution,true), false);
		log.debug("d: " + Arrays.toString(d));
		solution = Vector.sum(solution, d);
		numberOfIterations++;
		
	}
	
	

	@Override
	protected boolean isFinish() {
		return Vector.twoNorm(d) > tolerance  ? false : true;
	}
	
	@Override
	public String getName()
	{
		return "Local Newton";
	}

	@Override
	public double[] getSolution()
	{
		return solution;
	}

}
