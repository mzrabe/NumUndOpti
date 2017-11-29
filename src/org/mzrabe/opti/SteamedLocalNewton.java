package org.mzrabe.opti;

import java.util.Arrays;

import org.mzrabe.lina.Function;
import org.mzrabe.lina.Gauss;
import org.mzrabe.lina.Vector;

public class SteamedLocalNewton extends OptiAlgorithm {
	
	/**
	 * the direction vector 
	 */
	private double[] d;
	/**
	 * the step strategy for the algorithms
	 */
	private StepStrategy step;
	
	/**
	 * {@link OptiAlgorithm#OptiAlgorithm(Function, double, int)}
	 */
	public SteamedLocalNewton(Function func, double tolerance, int maxNumIteration) {
		this(func, tolerance, maxNumIteration, new ArmijoStep());
	}
	
	/**
	 * {@link OptiAlgorithm#OptiAlgorithm(Function, double, int)}
	 * @param step - the step strategy for the optimization 
	 */
	public SteamedLocalNewton(Function func, double tolerance, int maxNumIteration, StepStrategy step) {
		super(func, tolerance, maxNumIteration);
		this.step = step;
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
		log.debug("d = " +  Arrays.toString(d));
		d = Vector.multiScalar(d, step.getAlpha(func,d,solution));
		log.debug("d * alpha= " +  Arrays.toString(d));
		solution = Vector.sum(solution, d);
		log.debug("x + d * alpha= " +  Arrays.toString(solution));
		numberOfIterations++;
		
	}
	
	

	@Override
	protected boolean isFinish() {
		return Vector.twoNorm(d) > tolerance ? false : true;
	}
	
	@Override
	public String getName()
	{
		return "Steamed Local Newton";
	}

	@Override
	public double[] getSolution()
	{
		return solution;
	}

}
