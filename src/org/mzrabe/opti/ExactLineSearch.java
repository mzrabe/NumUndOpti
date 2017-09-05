package org.mzrabe.opti;

import java.util.Arrays;

import org.mzrabe.lina.Function;
import org.mzrabe.lina.Vector;
import org.mzrabe.zeropoint.Newton;

public class ExactLineSearch implements StepStrategy {

	@Override
	public double getAlpha(Function f, double[] d, double[] x, double... c) throws Exception {
		
		double[] parameterOfF = Arrays.copyOf(x, x.length);
		
		Function line = new Function() {
			
			@Override
			public double getValue(double[] x, double... s) throws Exception {
				return f.getValue(Vector.sum(parameterOfF, Vector.multiScalar(d, x[0])), s);
			}
		};
		
		double alpha = 1;
		return Newton.getSolution(line, false, alpha);
	}

}
