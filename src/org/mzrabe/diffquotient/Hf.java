package org.mzrabe.diffquotient;

/**
 * 
 * @author mzrabe
 * Class which calculate the Hesse matrix.
 */
public abstract class Hf extends Df{
	
	public double[][] getHf(double[] x){
		double[][] Hf = new double[x.length][x.length];
		
		for(int i=0;i<Hf.length;i++){
			for(int j=0;j<Hf.length;j++){
				try {
//					System.out.println(String.format("Zeile %d, Spalte %d = %f", i,j,D1f(x, i,j)));
					Hf[i][j] = D1f(x, i,j);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		
		return Hf;
	}
	
	// Hf = {{d²f/d²x[i],d²f/dx[i]dx[j]},
	
}
