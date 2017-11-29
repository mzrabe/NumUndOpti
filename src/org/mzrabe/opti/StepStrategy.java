package org.mzrabe.opti;

import org.mzrabe.lina.Function;

public interface StepStrategy {
	
	public double getAlpha(Function f, double[] d, double[] x, double ... c);

}
