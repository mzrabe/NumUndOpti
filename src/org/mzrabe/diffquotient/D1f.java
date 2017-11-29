package org.mzrabe.diffquotient;

import java.util.Arrays;

import org.mzrabe.lina.Function;

public class D1f implements Gradient {

	private double dh = 1e-6;
	
	public D1f(){};
	public D1f(double dh)
	{
		this.dh = dh;
	};
	
	@Override
	public void setDh(double dh)
	{
		this.dh = dh;
	}

	@Override
	public double df(Function func, double[] x, double[] c, int... n) throws Exception{
		if (n.length == 0) {
			throw new Exception("The lenght of n have to be greater or equals than 1!");
//		} else if (x.length < n.length) {
//			throw new Exception("The lenght of x have to be greater or equals than n!");
		} else if (n.length == 1) {
			double[] x1 = Arrays.copyOf(x, x.length);
			x1[n[0]] += dh;
			return (func.getValue(x1, c) - func.getValue(x, c)) / dh;
		} else {
			double[] xi = Arrays.copyOf(x, x.length);
			int[] n_new = Arrays.copyOf(n, n.length - 1);
			xi[n[n.length - 1]] += dh;
			return (df(func, xi, c, n_new) - df(func, Arrays.copyOf(x, x.length), c, n_new)) / dh;
		}
	}
	
	@Override
	public double df(Function func, double[] x, int... n) throws Exception
	{
		return df(func, x, null, n);
	}

	@Override
	public double[] grad(Function func, double[] x, boolean negativ) {
		double[] grad = new double[x.length];

		for (int i = 0; i < grad.length; i++) {
			double[] temp = Arrays.copyOf(x, x.length);
			temp[i] += dh;
			grad[i] = (func.getValue(temp) - func.getValue(x)) / dh;
		}

		if (negativ) {
			for (int i = 0; i < grad.length; i++) {
				grad[i] *= -1;
			}
		}

		return grad;
	}

}
