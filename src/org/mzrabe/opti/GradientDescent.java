package org.mzrabe.opti;

import java.util.Arrays;

import org.mzrabe.lina.Function;
import org.mzrabe.lina.Vector;

public class GradientDescent extends OptiAlgorithm{

	
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
	public GradientDescent(Function func, double tolerance, int maxNumIteration) {
		this(func, tolerance, maxNumIteration, new ArmijoStep());
	}
	
	/**
	 * {@link OptiAlgorithm#OptiAlgorithm(Function, double, int)}
	 * @param step - the step strategy for the optimization 
	 */
	public GradientDescent(Function func, double tolerance, int maxNumIteration, StepStrategy step) {
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
		d = new double[]{1.,1.};
		
	}

	/* (non-Javadoc)
	 * @see org.mzrabe.opti.OptiAlgorithm#algorithms()
	 */
	@Override
	protected void algorithms() throws Exception {
		
		d = func.grad(solution, true);
		d = Vector.multiScalar(d, step.getAlpha(func, d, solution));
		solution = Vector.sum(solution, d);
		
//		numberOfIterations++;
		
	}
	
	

	@Override
	protected boolean isFinish() {
		return Vector.twoNorm(d) > tolerance  ? false : true;
	}
	
	@Override
	public String getName()
	{
		return "Gradient Descent";
	}

	@Override
	public double[] getSolution()
	{
		return solution;
	}
	
}
