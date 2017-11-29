package org.mzrabe.opti;

import java.util.Arrays;

import org.mzrabe.lina.Function;
import org.mzrabe.lina.Vector;

public class ExactLineSearch implements StepStrategy {

	@Override
	public double getAlpha(Function f, double[] d, double[] x, double... c) {
		
		double[] parameterOfF = Arrays.copyOf(x, x.length);
		
		Function line = new Function() {
			
			@Override
			public double getValue(double[] x, double... s) {
				return f.getValue(Vector.sum(parameterOfF, Vector.multiScalar(d, x[0])), s);
			}
		};
		
		double tol = 1e-10;
		int num = 0;
		double m;
		double[] alpha = {1};
		do
		{
			m = line.grad(alpha,true)[0]/line.getHf(alpha)[0][0]; 
			alpha[0] = alpha[0] + m;  
			num++;
		}
		while(Math.abs(m)>tol && num <= 100);
			
		return alpha[0];
	}

}
