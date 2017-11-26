package org.mzrabe.diffquotient;


import org.mzrabe.lina.Function;

/**
 * 
 * @author mzrabe
 * Class which calculate the Hesse matrix.
 */
public class HessianMatrix{
	
	Gradient grad;
	
	public HessianMatrix()
	{
		grad = new D1f();
	}
	
	public double[][] getHf(Function f, double[] x, double ... c)
	{
		double[][] Hf = new double[x.length][x.length];
		
		for(int i=0;i<Hf.length;i++){
			for(int j=0;j<Hf.length;j++){
				try {
//					System.out.println(String.format("Zeile %d, Spalte %d = %f", i,j,D1f(x, i,j)));
					Hf[i][j] = grad.df(f,x,c, i,j);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		
		return Hf;
	}
	
}
