package org.mzrabe.utils;

/**
 * @author ,  
 */
public class Util
{
	
	public static double[] range(double start, double end, double step)
	{
		double[] back = new double[(int) ((end-start)/step)];
		for(int i=0;i<back.length;i++)
		{
			back[i] = start + i*step;
		}
		return back;
	}

}
