package org.mzrabe.opti;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.mzrabe.lina.Function;
import org.mzrabe.lina.Vector;

public class ArmijoStep implements StepStrategy {
	
	/**
	 * Coefficient for the Armijo algorithms this value has to be
	 * between 0 and 1.
	 */
	private double sigma;
	/**
	 * Coefficient for the Armijo algorithms this value has to be
	 * between 0 and 1.
	 */
	private double beta;
	private int maxIteration;
	private static final Logger log = LogManager.getRootLogger();
	
	public ArmijoStep()
	{
		this.sigma = 0.25;
		this.beta = 0.25;
		this.maxIteration = 300;
	}
	
	@Override
	public double getAlpha(Function f, double[] d, double[] x, double ... c) {
		
		double alpha = 1;
		int numIteration = 1;
		
		double rightSide = f.getValue(Vector.sum(x, Vector.multiScalar(d, alpha)), c);
		double leftSide = (f.getValue(x, c)+sigma*alpha*Vector.scalarProduct(d, f.grad(x, false)));
		log.debug(String.format("numIteration = %d", numIteration));
		log.debug(String.format("alpha %f", alpha));
		log.debug(String.format("rigthSide %f, leftSide %f", rightSide, leftSide));
		
		while(rightSide > leftSide
				  && numIteration < maxIteration)
		{
			numIteration++;
			log.debug(String.format("numIteration = %d", numIteration));
			alpha *= beta;
			log.debug(String.format("alpha %f", alpha));
			rightSide = f.getValue(Vector.sum(x, Vector.multiScalar(d, alpha)), c);
			leftSide = (f.getValue(x, c)+sigma*alpha*Vector.scalarProduct(d, f.grad(x, false)));
			log.debug(String.format("rigthSide %f, leftSide %f", rightSide, leftSide));
			
			
			
			
		}
		log.debug("Armijo step algorithms needed " + numIteration + " iterations. Alpha = " + alpha + "\n");
		return alpha;
	}
	
	/**
	 * @return the {@link #sigma}
	 */
	public double getSigma() {
		return sigma;
	}
	
	/**
	 * @return the {@link #beta}
	 */
	public double getBeta() {
		return beta;
	}

	/**
	 * @return the {@link #maxIteration}
	 */
	public int getMaxIteration() {
		return maxIteration;
	}

	/**
	 * @param alpha the {@link #alpha} to set
	 */
	public ArmijoStep setSigma(double sigma) {
		if (sigma > 0 && sigma < 1) {
			this.sigma = sigma;
			return this;
		} else {
			throw new IllegalArgumentException("The value 'sigma' have to be greater than 0 and lower than 1. The value is " + sigma);
		}
	}

	/**
	 * @param beta the {@link #beta} to set
	 */
	public ArmijoStep setBeta(double beta) {
		if (beta > 0 && beta < 1) {
			this.beta = beta;
			return this;
		} else {
			throw new IllegalArgumentException("The value 'beta' have to be greater than 0 and lower than 1. The value is " + beta);
		}
	}

	/**
	 * @param maxIteration the {@link #maxIteration} to set
	 */
	public ArmijoStep setMaxIteration(int maxIteration) {
		if (maxIteration > 0) {
			this.maxIteration = maxIteration;
			return this;
		} else {
			throw new IllegalArgumentException("The value 'maxIteration' have to be greater than 0. The value is " + maxIteration);
		}
	}
}
