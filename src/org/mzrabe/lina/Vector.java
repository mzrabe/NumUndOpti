package org.mzrabe.lina;

import java.util.Arrays;
import java.util.Locale;

/**
 * Mathematics operation with vectors.
 * @author Moritz Zahn, mzrabe@gmail.com  
 */
public class Vector {
	
	/** the basis vector of the x-axes */
	public static final double[] XBASIS = {1,0,0};
	/** the basis vector of the y-axes */
	public static final double[] YBASIS = {0,1,0};
	/** the basis vector of the z-axes */
	public static final double[] ZBASIS = {0,0,1};
	
	/** the negative basis vector of the x-axes */
	public static final double[] nXBASIS = {-1,0,0};
	/** the negative basis vector of the y-axes */
	public static final double[] nYBASIS = {0,-1,0};
	/** the negative basis vector of the z-axes */
	public static final double[] nZBASIS = {0,0,-1};
	
	/**
	 * Rotate the vector v counterclockwise around the x basis vector (x-axes) with the angle alpha, 
	 * around the y-axis with the angle beta and around the z-axes with the angel gamma. 
	 * @param v - the vector to rotate in Cartesian coordinates
	 * @param alpha - the angle to rotate around the x-axes in radiant
	 * @param beta - the angle to rotate around the y-axes in radiant
	 * @param gamma - the angle to rotate around the z-axes in radiant
	 * @return - the rotated vector v
	 */
	public static double[] rotate(double[] v, double alpha, double beta, double gamma)
	{
		
		if(v.length != 3)
			throw new IllegalArgumentException("The vector needs a dimesion of 3. The vector v has a dimension of " + v.length);
		
		return rotateXAxis(rotateYAxis(rotateZAxis(v, gamma), beta), alpha);
		
		
	}
	
	/**
	 * Rotate the point v counterclockwise around the given axes with the given angle.
	 * @param v - the point as vector in Cartesian coordinates 
	 * @param axes - the axes as direction vector in Cartesian coordinates which goes through the origin point
	 * @param phi - the angle in radiant
	 * @return - the rotated point
	 */
	public static double[] rotate(double[] v, double[] axes, double phi)
	{
		if(twoNorm(axes) == 0)
			throw new IllegalArgumentException("The vector [0,0,0]^T is not a direction!");
		
		double[] r = getNormVector(axes);
		
		if(r[0] == XBASIS[0] && r[1] == XBASIS[1] && r[2] == XBASIS[2])
			return rotateXAxis(v, phi);
		if(r[0] == YBASIS[0] && r[1] == YBASIS[1] && r[2] == YBASIS[2])
			return rotateYAxis(v, phi);
		if(r[0] == ZBASIS[0] && r[1] == ZBASIS[1] && r[2] == ZBASIS[2])
			return rotateZAxis(v, phi);
		if(r[0] == nXBASIS[0] && r[1] == nXBASIS[1] && r[2] == nXBASIS[2])
			return rotateXAxis(v, -phi);
		if(r[0] == nYBASIS[0] && r[1] == nYBASIS[1] && r[2] == nYBASIS[2])
			return rotateYAxis(v, -phi);
		if(r[0] == nZBASIS[0] && r[1] == nZBASIS[1] && r[2] == nZBASIS[2])
			return rotateZAxis(v, -phi);

		double[] back;
		
		double gamma = angleRightHandRule(new double[]{r[0],r[1],0},YBASIS,ZBASIS);
		double alpha = angleRightHandRule(new double[]{0,r[1],r[2]},ZBASIS,XBASIS);
		
		double[][] RxRz = Matrix.multi(Matrix.getRx(alpha), Matrix.getRz(gamma));
		
		back = Matrix.multi(Matrix.trans(RxRz), Matrix.multi(Matrix.getRz(phi), Matrix.multi(RxRz, v)));
		
		return back;
		
	}
	
	/**
	 * Rotate the point v counterclockwise around the given axes which is defined over the given point P1 and P2.
	 * Use this method if the rotation axis not goes through the origin point [0,0,0]^T.
	 * @param v - the point as vector in Cartesian coordinates 
	 * @param P1 - the start point of the rotation axis
	 * @param dir - the direction vector of the rotation axis
	 * @param phi - the angle in radiant
	 * @return - the rotated point
	 */
	public static double[] rotate(double[] v, double[] P1, double[] dir, double phi)
	{
		if(P1[0] == 0 && P1[1] == 0 && P1[2] == 0)
			/* rotation axis goes through the origin */
			return rotate(v, dir, phi);
		
		double[] back = Arrays.copyOf(v, v.length);
		
		/* rotation axis must goes to the origin */
		subtract(back, P1);
		
		double[] axes = dir;
		
		if(twoNorm(axes) == 0)
			throw new IllegalArgumentException("The vector [0,0,0]^T is not a direction!");
		
		double[] r = getNormVector(axes);
		
		if(r[0] == XBASIS[0] && r[1] == XBASIS[1] && r[2] == XBASIS[2])
			back = rotateXAxis(back, phi);
		else if(r[0] == YBASIS[0] && r[1] == YBASIS[1] && r[2] == YBASIS[2])
			back = rotateYAxis(back, phi);
		else if(r[0] == ZBASIS[0] && r[1] == ZBASIS[1] && r[2] == ZBASIS[2])
			back = rotateZAxis(back, phi);
		else if(r[0] == nXBASIS[0] && r[1] == nXBASIS[1] && r[2] == nXBASIS[2])
			back = rotateXAxis(back, -phi);
		else if(r[0] == nYBASIS[0] && r[1] == nYBASIS[1] && r[2] == nYBASIS[2])
			back = rotateYAxis(back, -phi);
		else if(r[0] == nZBASIS[0] && r[1] == nZBASIS[1] && r[2] == nZBASIS[2])
			back = rotateZAxis(back, -phi);
		else
		{
//			double alpha = Math.acos(r[2]);
//			double gamma = (r[0] == 0 && r[1] == 0) ? 0 : Math.asin(r[0]/Math.sqrt(r[0]*r[0]+r[1]*r[1]));
			
			/* rotate the rotation axes into the y-z area */
			double gamma = (r[0] == 0 && r[1] == 0) ? 0 : angleRightHandRule(new double[]{r[0],r[1],0},YBASIS,ZBASIS);
			/* rotate the rotation axes to the z axes */
			double alpha =  angleRightHandRule(new double[]{0,r[1],r[2]},ZBASIS,XBASIS);
			
			double[][] RxRz = Matrix.multi(Matrix.getRx(alpha), Matrix.getRz(gamma));
			
			back = Matrix.multi(Matrix.trans(RxRz), Matrix.multi(Matrix.getRz(phi), Matrix.multi(RxRz, back)));
		}
		
		add(back, P1);
		
		return back;
		
	}
	
	
	
	/**
	 * Get the normed vector with the length of 1
	 * @param v - the vector
	 * @return - the normed vector
	 */
	public static double[] getNormVector(double[] v)
	{
		double f = twoNorm(v);
		return f == 0 ? new double[v.length] : new double[]{v[0]/f,v[1]/f,v[2]/f};
	}
	
	/**
	 * Calculate the vector product of two vectors with die dimension 3x1.
	 * @param v - the one vector 3x1
	 * @param u - the other vector 3x1
	 * @return the vector product
	 */
	public static double[] vectorProdukt(double[] v, double[] u)
	{
		if(v == null)
			throw new NullPointerException("The vector v is null");
		
		if(u == null)
			throw new NullPointerException("The vector u is null");
		
		if(v.length != 3 || u.length != 3)
			throw new IllegalArgumentException("The vectors needs a dimesion of 3. The dim(v) = " + v.length + " and dim(u)= " + u.length);
		
		
		double[] back = new double[v.length];
		
		back[0] = v[1]*u[2] - v[2]*u[1];
		back[1] = v[2]*u[0] - v[0]*u[2];
		back[2] = v[0]*u[1] - v[1]*u[0];
		
		return back;
	}
	
	/**
	 * Rotate the vector v counterclockwise around the x basis vector (x-axes). 
	 * @param v - the vector to rotate in Cartesian coordinates
	 * @param alpha - the angle to rotate in radiant
	 * @return - the rotated vector v
	 */
	public static double[] rotateXAxis(double[] v, double alpha)
	{
		double[] back = Arrays.copyOf(v, v.length);
		
		/* nothing to rotate */
		if(alpha == 0)
		{
			return back;
		}
		
		if(v.length != 3)
			throw new IllegalArgumentException("The vector needs a dimesion of 3. The vector v has a dimension of " + v.length);
		
		double sinAlpha = Math.sin(alpha);
		double cosAlpha = Math.cos(alpha);
		
		double[][] Rx = {{1,0,0},{0,cosAlpha,-sinAlpha},{0,sinAlpha,cosAlpha}};
		
		return Matrix.multi(Rx, v);
		
		
	}
	
	/**
	 * Rotate the vector v counterclockwise around the y basis vector (y-axes). 
	 * @param v - the vector to rotate in Cartesian coordinates
	 * @param beta - the angle to rotate in radiant
	 * @return - the rotated vector v
	 */
	public static double[] rotateYAxis(double[] v, double beta)
	{
		
		double[] back = Arrays.copyOf(v, v.length);
		
		/* nothing to rotate */
		if(beta == 0)
		{
			return back;
		}
		
		
		if(v.length != 3)
			throw new IllegalArgumentException("The vector needs a dimesion of 3. The vector v has a dimension of " + v.length);
		
		double sinBeta = Math.sin(beta);
		double cosBeta = Math.cos(beta);
		
		double[][] Ry = {{cosBeta,0.,sinBeta},{0.,1.,0.},{-sinBeta,0.,cosBeta}};
		
		return Matrix.multi(Ry, v);
	}
	
	/**
	 * Rotate the vector v counterclockwise around the z basis vector (z-axes). 
	 * @param v - the vector to rotate in Cartesian coordinates
	 * @param gamma - the angle to rotate in radiant
	 * @return - the rotated vector v
	 */
	public static double[] rotateZAxis(double[] v, double gamma)
	{
		
		double[] back = Arrays.copyOf(v, v.length);
		
		/* nothing to rotate */
		if(gamma == 0)
		{
			return back;
		}
		
		if(v.length != 3)
			throw new IllegalArgumentException("The vector needs a dimesion of 3. The vector v has a dimension of " + v.length);
		
		double sinGamma = Math.sin(gamma);
		double cosGamma = Math.cos(gamma);
		
		double[][] Rz = {{cosGamma,-sinGamma,0.},{sinGamma,cosGamma,0.},{0.,0.,1.}};
		
		return Matrix.multi(Rz, v);
	}
	
	/** 
	 * Get the angle between these to vector (form 0 to 180 degree). This calculation base on the scalar product. The angle is return in radiant. 
	 * @param v - the one vector
	 * @param u - the other vector
	 * @return - the angle between the vectors in radiant
	 * 
	 */
	public static double angle(double[] v, double[] u)
	{
		
		double normV = twoNorm(v);
		double normU = twoNorm(u);
		if(normV == 0)
			throw new IllegalArgumentException("The length (Euclidean norm) of v has to be greater than zero.");
		if(normU == 0)
			throw new IllegalArgumentException("The length (Euclidean norm) of u has to be greater than zero.");
		
		double a = scalarProduct(v, u);
		
		
		return Math.acos(a/(normV * normU));
	}
	
	/** 
	 * Get the angle between these to vector (form 0 to 180 degree). 
	 * This calculation base on the vector product so you get the smallest angel to rotate the vector v in the direction of u 
	 * by using the right hand rule (the thumb is the vector product of u and v). The angle is return in radiant. If you want to rotate something with this angle you has to use 
	 * the vector product v*u as rotation axes.
	 * @param v - the one vector
	 * @param u - the other vector
	 * @return - the angle between the vectors in radiant
	 * 
	 */
	public static double angleRightHandRule(double[] v, double[] u)
	{
		
		double normV = twoNorm(v);
		double normU = twoNorm(u);
		if(normV == 0)
			throw new IllegalArgumentException("The length (Euclidean norm) of v has to be greater than zero.");
		if(normU == 0)
			throw new IllegalArgumentException("The length (Euclidean norm) of u has to be greater than zero.");
		
		double normVxU = twoNorm(vectorProdukt(v, u));
		
		
		return Math.asin(normVxU/(normV * normU));
	}
	
	/** 
	 * Get the angle between these to vector (form -180 to 180 degree). 
	 * This calculation base on the vector product so you get the smallest angel to rotate the vector v in the direction of u 
	 * by using the right hand rule around the given rotation axes. The angle is return in radiant.
	 * @param v - the one vector
	 * @param u - the other vector
	 * @param rotationAxes - the wanted rotation axes, you have to be sure that the vector product of v and u is collinear to the rotationAxes
	 * @return - the angle between the vectors in radiant, if the rotationAxes shows in the opposite direction of the vector product of u and v the angel will be negative
	 * 
	 */
	public static double angleRightHandRule(double[] v, double[] u, double[] rotationAxes)
	{
		
		double normV = twoNorm(v);
		double normU = twoNorm(u);
		if(normV == 0)
			throw new IllegalArgumentException("The length (Euclidean norm) of v has to be greater than zero.");
		if(normU == 0)
			throw new IllegalArgumentException("The length (Euclidean norm) of u has to be greater than zero.");
		
		double[] vxu = vectorProdukt(v, u);
		
		double normVxU = twoNorm(vxu);
		double negative = (vxu[0]*-1 == rotationAxes[0] && vxu[1]*-1 == rotationAxes[1] && vxu[2]*-1 == rotationAxes[2]) ? -1 : 1;
		
		return negative * Math.asin(normVxU/(normV * normU));
	}
	
	/** 
	 * Get the angle between these to vector (form -180 to 180 degree). 
	 * This calculation base on the vector product so you get the smallest angel to rotate the vector v in the direction of u 
	 * by using the right hand rule around the given rotation axes. The angle is return in radiant.
	 * @param v - the one vector
	 * @param u - the other vector
	 * @param rotationAxes - the wanted rotation axes, you have to be sure that the vector product of v and u is collinear to the rotationAxes
	 * @return - the angle between the vectors in radiant, if the rotationAxes shows in the opposite direction of the vector product of u and v the angel will be negative
	 * 
	 */
	public static double angleBetween(double[] v, double[] u)
	{
		
		//FIXME noch nich fertig
		
		double[] vxu = vectorProdukt(v, u);
		
		double normV = twoNorm(v);
		double normU = twoNorm(u);
		if(normV == 0)
			throw new IllegalArgumentException("The length (Euclidean norm) of v has to be greater than zero.");
		if(normU == 0)
			throw new IllegalArgumentException("The length (Euclidean norm) of u has to be greater than zero.");
		
		double normVxU = twoNorm(vxu);
		
		double[] rotationAxisToXYPlane = Vector.vectorProdukt(vxu, Vector.ZBASIS);
		double rotationAngleToZAxes = angleRightHandRule(vxu,Vector.ZBASIS,rotationAxisToXYPlane);
		
		double[] v_ = rotate(v, rotationAxisToXYPlane, rotationAngleToZAxes);
		double[] u_ = rotate(u, rotationAxisToXYPlane, rotationAngleToZAxes);
		
		return Math.asin(normVxU/(normV * normU));
	}
	
	

	/**
	 * Calculate the sum of the vectors
	 * @param v - the given vectors
	 * @return - the sum of all given vectors, return a new instance of double[]
	 */
	public static double[] sum(double[]... v) {
		double[] back = new double[v[0].length];
		for (double[] vector : v) {
			if (vector.length != v[0].length) {
				throw new IllegalArgumentException("The vector must have the same dimension! Current vector "
						+ vector.length + " first vector length = " + v[0].length);
			}
			for (int i = 0; i < back.length; i++) {
				back[i] += vector[i];
			}
		}
		return back;
	}
	
	/**
	 * Add the same scalar value to each component of the vector
	 * @param h - the scalar value
	 * @param v - the vector
	 * @return - the result of the addition
	 */
	public static double[] addScalar(double h,double... v){
		double[] back = new double[v.length];
		for(int i=0;i<v.length;i++){
			back[i]=v[i]+h;
		}
		return back;
	}
	
	/**
	 * Calculate the scalar multiplication of a vector with a scalar
	 * @param x - the vector
	 * @param scalar - the scalar
	 * @return - the result of the scalar multiplication - a new double[] array
	 */
	public static double[] multiScalar(double[] x, double scalar)
	{
		double[] back = Arrays.copyOf(x, x.length);
		
		if(scalar == 1.0)
			return back;
		
		for(int i=0;i<back.length;i++)
		{
			back[i]*=scalar;
		}
		return back;
	}
	
	/**
	 * The scalar product of these to vectors
	 * @param x - the one vector
	 * @param y - the other vector
	 * @return - 
	 */
	public static double scalarProduct(double[] x, double[] y)
	{
		double back = 0;
		
		for(int i=0;i<x.length;i++)
		{
			back+=x[i]*y[i];
		}
		return back;
	}
	
	/**
	 * Print the vector to System.out
	 * @param v - the vector to print
	 */
	public static void print(double... v){
		StringBuilder sb = new StringBuilder("[");
		
		for(int i=0;i<v.length-1;i++){
			sb.append(String.format(Locale.ENGLISH,"%f, ", v[i]));
		}
		sb.append(String.format(Locale.ENGLISH,"%f]^T", v[v.length-1]));
		System.out.println(sb.toString());
	}
	
	/**
	 * Print the vector to System.out
	 * @param v - the vector to print
	 */
	public static void print(int ... v){
		StringBuilder sb = new StringBuilder("[");
		
		for(int i=0;i<v.length-1;i++){
			sb.append(String.format(Locale.ENGLISH,"%d, ", v[i]));
		}
		sb.append(String.format(Locale.ENGLISH,"%d]^T", v[v.length-1]));
		System.out.println(sb.toString());
	}
	
	/**
	 * Get the vector a string (transposed)
	 * @param v - the vector to print
	 * @return - v as String
	 */
	public static String asString(double... v){
		
		if(v == null)
			return "null";
		
		StringBuilder sb = new StringBuilder("[");
		
		for(int i=0;i<v.length-1;i++){
			sb.append(String.format(Locale.ENGLISH,"%f, ", v[i]));
		}
		sb.append(String.format(Locale.ENGLISH,"%f]^T", v[v.length-1]));
		return sb.toString();
	}
	
	/**
	 * Calculate subtraction of two vectors
	 * @param a - the one vector
	 * @param b - the other vector
	 * @return - the result of the subtraction, new instance of a double[]
	 */
	public static double[] minus(double[] ... v)
	{
		double[] back = Arrays.copyOf(v[0], v[0].length);
		for (int j = 1;j<v.length;j++) {
			if (v[j].length != v[0].length) {
				throw new IllegalArgumentException("The vector must have the same dimension! Current vector "
						+ v[j].length + " first vector length = " + v[0].length);
			}
			for (int i = 0; i < back.length; i++) {
				back[i] -= v[j][i];
			}
		}
		return back;
	}
	
	/**
	 * Calculate the two norm of the vector (Euclidean norm)
	 * @param x - the vector
	 * @return - the two norm
	 */
	public static double twoNorm(double[] x){
		double twoNorm = 0;
		for(int i = 0; i<x.length;i++){
			twoNorm += Math.pow(x[i], 2);
		}
		return Math.sqrt(twoNorm);
	}
	
	/**
	 * Calculate the one norm of the vector (Absolute-value norm)
	 * @param x - the vector
	 * @return - the one norm
	 */
	public static double oneNorm(double[] x){
		double oneNorm = Math.abs(x[0]);
		
		for(int i = 1; i<x.length;i++){
			if(Math.abs(x[i])>oneNorm){
				oneNorm = Math.abs(x[i]);
			}
		}
		return Math.abs(oneNorm);
	}
	
	/**
	 * Calculate the sum of the vector components.
	 * @param x - the vector
	 * @return  - sum of the vector components
	 */
	public static double sumVector(double[] x)
	{
		double sum = 0;
		for(double d : x)
		{
			sum+=d;
		}
		return sum;
	}
	
	/**
	 * Add the vector u to the vector v. the vector v will changed
	 * @param v - the vector v which will increase, this instance will changed
	 * @param u - the vector u which will increase v
	 */
	public static void add(double[] v, double[] u)
	{
		if(v.length != u.length)
		{
			throw new IllegalArgumentException("The vectors must have the same length. The vector v = " + v.length + " and u = " + u.length +".");
		}
		
		for(int i=0;i<v.length;i++)
		{
			v[i] += u[i];
		}
	}
	
	/**
	 * Subtract the vector u from the vector v.
	 * @param v - the vector v which will decrease, this instance will changed
	 * @param u - the vector u which will decrease v
	 */
	public static void subtract(double[] v, double[] u)
	{
		if(v.length != u.length)
		{
			throw new IllegalArgumentException("The vectors must have the same length. The vector v = " + v.length + " and u = " + u.length +".");
		}
		
		for(int i=0;i<v.length;i++)
		{
			v[i] -= u[i];
		}
	}
	
	/**
	 * @param v - one vector
	 * @param u - another vector
	 * @param ablsolute - true if only the absolute values should compare
	 * @return - true if the vectors are equal, otherwise false (also if the vectors have not the same length)
	 */
	public static boolean equals(double[] v, double[] u, boolean ablsolute)
	{
		if(v.length != u.length)
			return false;
		else
			for(int i = 0;i<v.length;i++)
			{
				if(ablsolute)
				{
					if(Math.abs(v[i]) != Math.abs(u[i]))
						return false;
				}
				else
				{
					if(v[i] != u[i])
						return false;
				}
			}
		return true;
	}
	
	/**
	 * Checks if the two vector have the same direction
	 * @param v - a vector
	 * @param u - another vector 
	 * @param absolute - if flag is true, than v equals u or v equals -u
	 * @return  - true if the vectors have the same direction
	 * 
	 */
	public static boolean directonEquals(double[] v, double[] u, boolean absolute)
	{
		if(v.length != u.length)
			return false;
		else
			return equals(getNormVector(v), getNormVector(u),absolute);
	}
	
	/**
	 * Calculate the surface of the of the polygon according the Gauss's area formula.
	 * @param v - the points of the polygon ordered counterclockwise or clockwise - minimum 3 points
	 * @return - the surface of the polygon, 0 if the number of points is lower than 3 or if the points have not the dimesion 2
	 * 
	 */
	public static double getSurface(double[] ... v)
	{
		if(v.length < 3)
			return 0;
		if(v[0].length != 2)
			return 0;
		
		double s = v[v.length-1][0] * v[0][1] - v[v.length-1][1] * v[0][0];
		
		for(int i = 0; i<v.length-1;i++)
		{
			s += v[i][0] * v[i+1][1] - v[i][1] * v[i+1][0];
		}
		
		return Math.abs(s/2);
		
	}
	
	public static double[] connectVector(double[] ... v)
	{
		int lenght = 0;
		for(double[] d : v)
			lenght+=d.length;
		double[] back = new double[lenght];
		lenght = 0;
		for(double[] d : v)
		{
			for(double x : d)
			{
				back[lenght] = x;
				lenght++;
			}
		}
		return back;
	}

	/**
	 * @param c - a vector
	 * @return - the absolute vector, this mean all values of the vector convert to a absolute value
	 */
	public static double[] abs(double[] c)
	{
		double[] back = new double[c.length];
		
		for(int i = 0; i<back.length; i++)
		{
			back[i] = Math.abs(c[i]);
		}
		return back;
	}

	/**
	 * @param v - a vector
	 * @return - the minimal value of this vector. The different between {@link #oneNorm(double[])} is that negative values are take into account.
	 */
	public static double min(double[] v)
	{
		double min = v[0];
		
		for(int i = 1;i<v.length;i++)
		{
			min = Math.min(min, v[i]);
		}
		
		return min;
	}
	

}
