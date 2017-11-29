package approximation;

import org.mzrabe.diffquotient.Df;
import org.mzrabe.diffquotient.Hf;
import org.mzrabe.lina.Function;
import org.mzrabe.lina.Gauss;
import org.mzrabe.lina.Matrix;
import org.mzrabe.lina.Vector;

public abstract class leastSquares implements Function {

	private final double[] x;
	private final double[] y;

	public leastSquares(double[] x, double[] y) {
		this.x = x;
		this.y = y;
	}
	
	

	public double[] approx(double[] a) {
		double[] residum = new double[a.length];
		double maxResidum;

		do {

//			double[][] A = new double[a.length][a.length];
//			double[] r = new double[a.length];
//
//			for (int i = 0; i < x.length; i++) {
//
//				final double xi = x[i];
//				final double yi = y[i];
//
//				Df df = new Df() {
//
//					@Override
//					public double function(double[] x) {
//						return Math.pow(this.function(x) - yi, 2);
//					}
//				};
//				Hf hf = new Hf() {
//					@Override
//					public double function(double[] x) {
//						return Math.pow(weibull(xi, x[0], x[1]) - yi, 2);
//					}
//				};
//
//				A = Matrix.addition(A, getHf(a));
//				r = Vector.addtion(r, Vector.multiScalar(getGradient(a), -1.));
//
//			}
//
//			residum = Gauss.getSolution(A, r, false);
//			maxResidum = Vector.oneNorm(residum);
//			System.out.println(String.format("Residum %15f", Math.abs(maxResidum)));
//			Vector.addtion(a, residum);

		} while (Math.abs(1/*maxResidum*/) > 1e-8);

		return a;
	}

}
