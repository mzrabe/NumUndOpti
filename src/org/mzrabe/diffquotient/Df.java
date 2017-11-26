package org.mzrabe.diffquotient;

import java.util.Arrays;

import org.mzrabe.lina.Function;

public abstract class Df implements Function{
	
	private double dh = 1e-3;
	
	/**
	 * @param dh the dh to set
	 */
	public void setDh(double dh)
	{
		this.dh = dh;
	}
	


	
	
	public double[] getGradient(double[] x){
		double[] grad = new double[x.length];
		
		for(int i=0;i<grad.length;i++){
			double[] temp = Arrays.copyOf(x, x.length);
			temp[i]+=dh;
			grad[i]=(function(temp)-function(x))/dh;
		}
		
		return grad;
	}
	
	public double D1f(double[] x, int... n) throws Exception{
		if(n.length == 0)
		{
			throw new Exception("The lenght of n have to be greater or equals than 1!");
		}
		else if(x.length < n.length)
		{
			throw new Exception("The lenght of x have to be greater or equals than n!");
		}
		else if(n.length == 1)
		{
			double[] x1 = Arrays.copyOf(x, x.length);
			x1[n[0]] += dh;
			return (function(x1)-function(x))/dh;
		}
		else
		{
			double[] xi = Arrays.copyOf(x, x.length);
			int[] n_new = Arrays.copyOf(n, n.length-1);
			xi[n[n.length-1]]+= dh;
			return (D1f(xi,n_new) - D1f(Arrays.copyOf(x,x.length),n_new))/dh;
		}		
	}
	

}
