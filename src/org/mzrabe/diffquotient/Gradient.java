package org.mzrabe.diffquotient;

import org.mzrabe.lina.Function;

public interface Gradient {
	
	public double df(Function func, double[] x, double[] c, int ... n) throws Exception;
	public double df(Function func, double[] x, int ... n) throws Exception;
	public double[] grad(Function func, double[] x, boolean negativ);
	public void setDh(double dh);

}
