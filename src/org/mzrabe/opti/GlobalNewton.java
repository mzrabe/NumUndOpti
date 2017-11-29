package org.mzrabe.opti;

import java.util.Arrays;

import org.mzrabe.lina.Function;
import org.mzrabe.lina.Gauss;
import org.mzrabe.lina.Vector;

public class GlobalNewton extends OptiAlgorithm {
	
	/**
	 * the direction vector 
	 */
	private double[] d;
	/**
	 * the step strategy for the algorithms
	 */
	private StepStrategy step;
	/**
	 * additional parameter for the algorithms
	 */
	private double delta, rho;
	
	/**
	 * {@link OptiAlgorithm#OptiAlgorithm(Function, double, int)}
	 */
	public GlobalNewton(Function func, double tolerance, int maxNumIteration) {
		this(func, tolerance, maxNumIteration, new ArmijoStep(),0.1,2.1);
	}
	
	/**
	 * {@link OptiAlgorithm#OptiAlgorithm(Function, double, int)}
	 * @param step - the step strategy for the optimization 
	 */
	public GlobalNewton(Function func, double tolerance, int maxNumIteration, StepStrategy step, double delta, double rho) {
		super(func, tolerance, maxNumIteration);
		this.step = step;
		this.delta = delta;
		this.rho = rho;
		
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
		
		solution= Arrays.copyOf(x[0],x[0].length);
		d = new double[]{1.,1.};
		
	}

	/* (non-Javadoc)
	 * @see org.mzrabe.opti.OptiAlgorithm#algorithms()
	 */
	@Override
	protected void algorithms() {
		
		d = Gauss.getSolution(func.getHf(solution), func.grad(solution,true), false);
		
		if(d == null || Vector.scalarProduct(func.grad(solution, false),d) > Math.pow(-delta*(Vector.twoNorm(d)), rho))
		{
			d = func.grad(solution, true);
			d = Vector.multiScalar(d, step.getAlpha(func, d, solution));
		}
		
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
		return "Global Newton";
	}

	@Override
	public double[] getSolution()
	{
		return solution;
	}

}
