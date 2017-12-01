package org.mzrabe;

import static java.lang.Math.pow;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Timer;
import java.util.TimerTask;
import java.util.TreeSet;

import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.Histogram;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.style.Styler.ChartTheme;
import org.knowm.xchart.style.markers.SeriesMarkers;
import org.mzrabe.approximation.leastSquares;
import org.mzrabe.diffquotient.JacobianMatrix;
import org.mzrabe.lina.Function;
import org.mzrabe.lina.MathOperation;
import org.mzrabe.lina.Matrix;
import org.mzrabe.lina.Vector;
import org.mzrabe.methaheuristic.SimulatedAnnealingTravelSalesMan;
import org.mzrabe.methaheuristic.firefly.FireFlyAlgorithm;
import org.mzrabe.opti.DownhillSimplex;
import org.mzrabe.opti.ExactLineSearch;
import org.mzrabe.opti.GlobalNewton;
import org.mzrabe.opti.GradientDescent;
import org.mzrabe.opti.LocalNewton;
import org.mzrabe.opti.OptiAlgorithm;
import org.mzrabe.opti.PenaltyMethod;
import org.mzrabe.opti.SimpleLocalNewton;
import org.mzrabe.opti.SteamedLocalNewton;
import org.mzrabe.plot.XYContourPlot;
import org.mzrabe.plot.XYPlot;
import org.mzrabe.utils.Util;
import org.mzrabe.zeropoint.Newton;

import static java.lang.System.out;

import java.awt.Window;

public class Main {
	static ArrayList<Double> allVelocity = new ArrayList<Double>();
	static ArrayList<Double> Velocity = new ArrayList<Double>();
	static ArrayList<Double> AnzVelocity = new ArrayList<Double>();
	static double[] x = { 1, 2, 3, 4 };
	static double[] y = { 6.5, 10., 22., 59. };
	
	static Function himmelblau = new Function()
	{
		
		@Override
		public double getValue(double[] x, double... c) throws Exception
		{
			return pow(pow(x[0], 2) + x[1] - 11, 2) + pow(x[0] + pow(x[1], 2) - 7, 2);
		}
	};

	public static void main(String[] args) throws Exception {
		
		
//		newtonTest();
//		volumeCross();
//		LRTest();
//		approxTest();
		aneometer();
		contourplottest();

	}
	
	
	/**
	 * @throws Exception 
	 * 
	 */
	private static void contourplottest() throws Exception
	{
		XYContourPlot plot = new XYContourPlot("Himmelblau", "x", "y", 800, 600, new Function()
		{
			
			@Override
			public double getValue(double[] x, double... c)
			{
				// TODO Auto-generated method stub
				return pow(pow(x[0], 2) + x[1] - 11, 2) + pow(x[0] + pow(x[1], 2) - 7, 2);
			}
		}, -5, 5, 0.1, 20);
		
		plot.showChart();
		
	}


	/**
	 * @throws Exception 
	 * 
	 */
	private static void aneometer() throws Exception
	{
		double r0 = 0.022;
		double r1 = 0.052;
		double A0 = 0.024*0.002;
		double A = Math.PI/4*Math.pow(0.036, 2);
		double beta = 2./3*Math.PI;
		double rho = 1.2;
		double Mr = 0.00001;
		
		double[][] xi = {{1.42},{2.06},{2.53},{3.02},{3.55},{4.01},{4.54},{5.02},{5.55},{6.05},{6.54},{7.02},{8.05},{10},{12}};
		double[] yi = {0.133,0.203,0.233,0.255,0.272,0.279,0.292,0.300,0.301,0.304,0.301,0.305,0.308,0.316,0.320};
		
		System.out.println(xi.length + ", " + yi.length);
		
		Function f_lambda = new Function()
		{
			
			@Override
			public double getValue(double[] x, double... c)
			{
				
				return r1*(Fo(c[1], c[0], x[0], rho, A) +  Fo(c[1]+beta, c[0], x[0], rho, A) +  Fo(c[1]+2*beta, c[0],x[0], rho, A)) + r0*(F(c[1], c[0], x[0], rho, A0) +  F(c[1]+beta, c[0], x[0], rho, A0) +  F(c[1]+2*beta, c[0],x[0], rho, A0)) - c[2];
			}
		};
		
		double[] alpha = Util.range(0, Math.PI*2, 30*Math.PI/180);
//		double[] alpha = {0*Math.PI/180};

		Function approx =	new Function()
			{
				
				@Override
				public double getValue(double[] x, double... c) throws Exception
				{
					double lambda = 0;
					for(int i = 0;i<alpha.length;i++)
					{
					lambda += Newton.getSolution(f_lambda, false, 0.1, c[0], alpha[i],x[0]);
//						System.out.println(String.format("%f", lambda[i]));
					}
					return lambda/alpha.length;
				}
		};
		
		double[] sol = leastSquares.approx(approx, new double[]{Mr}, yi, xi);
		

		
		
		double[] v = Util.range(0.01, 14, 0.1);
		double[] y1 = new double[v.length];
		for(int i = 0;i<v.length;i++)
		{
			double lambda = 0;
			for(int j = 0;j<alpha.length;j++)
			{
			lambda += Newton.getSolution(f_lambda, true, 0.1, v[i], alpha[j],sol[0]);
//				System.out.println(String.format("%f", lambda[i]));
			}
			lambda = lambda/alpha.length;
			y1[i] = lambda;
		}
		
		double[] vi = new double[xi.length];
		for(int i=0;i<xi.length;i++)
		{
			vi[i] = xi[i][0];
		}
		Vector.print(sol);
		XYPlot plot = new XYPlot("blub", "x", "y", 800, 600);
		plot.chart.getStyler().setYAxisMin(0.);
		plot.addDataSet("apprx", v, y1);
		plot.addDataSet("measure", vi,yi);
		plot.showChart();
		
	}
	
	private static double Fo(double angle, double v, double lambda, double rho, double A)
	{
		double c = (v * Math.signum(Math.round(Math.cos(angle)))) - lambda * v;
		double cw = c > 0 ? 1.34 : 0.34;
		double fo = rho/2 * cw * A * Math.pow(c, 2)*Math.signum(c);
//		System.out.println(String.format("angle %f, cw %f, fo %f, v %f, c %f, lambda %f, v * Math.cos(angle) %.15f", Math.toDegrees(angle), cw, fo,v,c, lambda * v, v * Math.cos(angle)));
		return fo;
		
	}
	
	private static double F(double angle, double v, double lambda, double rho, double A)
	{
		double c = (v * Math.signum(Math.round(Math.cos(angle)))) - lambda * v;
		double cw = 1.2;
		double fo = rho/2 * cw * A * Math.pow(c, 2)*Math.signum(c);
//		System.out.println(String.format("angle %f, cw %f, fo %f, v %f, c %f, lambda %f, v * Math.cos(angle) %.15f", Math.toDegrees(angle), cw, fo,v,c, lambda * v, v * Math.cos(angle)));
		return fo;
		
	}
	


	private static void approxTest() throws Exception
	{
		Function f = new Function()
		{
			
			@Override
			public double getValue(double[] x, double... c)
			{
				return x[1]*Math.sin(x[0]*c[0]);
			}
		};
		
		double[][] x = new double[100][1];
		double[] xx = new double[100];
		for(int i=0;i<x.length;i++)
		{
			x[i][0] = (double) i/10;
			xx[i] = (double) i/10;
		}
		double[] y = new double[x.length];
		for(int i=0;i<y.length;i++)
		{
			y[i] = f.getValue(new double[]{2,1},x[i]) + (-.2 + Math.random() * .4); 
		}
	
		

		
		double[] c = leastSquares.approx(f, new double[]{2.3,1.2}, y, x);
		System.out.println("------");
		Vector.print(c);
		
		double[] y2 = new double[1000];
		double[] x2 = new double[1000];
		for(int i=0;i<1000;i++)
		{
			x2[i] = (double) i/100;
		}
		for(int i=0;i<y2.length;i++)
		{
			y2[i] = f.getValue(c, x2[i]);
		}
		
		XYPlot plot = new XYPlot("test", "x", "y", 800, 600);
		plot.addDataSet("sinRand", xx, y);
		plot.addDataSet("sinApprox", x2, y2);
		plot.showChart();
	}
	
	private static void LRTest()
	{
		double[][] mat = new double[][]
				{
						{2.1,2512,-2516},
						{-1.3,8.8,-7.6},
						{0.9,-6.2,4.6}
				};
		double[]	vec = new double[] {6.5,-5.3,2.9};
		
		Vector.print(Matrix.solveRMS(mat, vec));
		mat = new double[][]
			{
					{2.1,2512,-2516},
					{-1.3,8.8,-7.6},
					{0.9,-6.2,4.6}
			};
		vec = new double[] {6.5,-5.3,2.9};
		Vector.print(Matrix.solveRCMS(mat, vec));
		mat = new double[][]
				{
						{5, 7, 3},
						{7, 11, 2},
						{3, 2, 6}
				};
			vec = new double[] {0,0,1};
		Vector.print(Matrix.solveChorlesky(mat, vec));
	}
	
	/**
	 * Test um den Umfang einer Schnittlinie zu ermitteln. Hier wurde dies an zwei kreisen gemacht
	 */
	private static void volumeCross() throws Exception{
		
		/** radius of the cylinder 1 */
		double r1 = 4;
		/** radius of the cylinder 2 */
		double r2 = 4;
		
		
		/* parmetrisierte Schnittkurve von 2 Zylindern dessen Achse die z- und y-Achse sind */
		Function[] P = 
			{
				new Function()
				{
					
					@Override
					public double getValue(double[] x, double... c)
					{
						return r1 * Math.sin(x[0]);
					}
				},
				new Function()
				{
					
					@Override
					public double getValue(double[] x, double... c)
					{
						return r1 * Math.cos(x[0]);
					}
				},
				new Function()
				{
					
					@Override
					public double getValue(double[] x, double... c)
					{
						return r2*Math.cos(Math.asin(r1/r2 * Math.sin(x[0])));
					}
				}
			};
		
		/* the "integration range" of the curve */
		double from = 0;
		double to = 2*Math.PI;
		double step = 0.001;
		
//		System.out.println(String.format("min: %f, max: %f, step: %f", min, max, step));
		
		double dist = 0.0;
		double temp[] = {P[0].getValue(new double[]{0.0}),P[1].getValue(new double[]{0.0}),P[2].getValue(new double[]{0.0})};
		
		/* calculate the path length in the range from-to */
		for(int i = 1; i<(to-from)/step;i++)
		{
			double[] x = {from+i*step};
			double[] E = {P[0].getValue(x),P[1].getValue(x),P[2].getValue(x)};
			
			dist+= Vector.twoNorm(Vector.minus(E, temp));
			
			temp = E;
			
//			System.out.println(String.format("x=%f, (%f,%f,%f)", x[0], E[0],E[1],E[2]));
		}
		
		System.out.println(String.format("distance = %f", dist));
	}
	
	private static void newtonTest() throws Exception
	{
		Function w1 = new Function()
		{
			
			@Override
			public double getValue(double[] x, double... c)
			{
				return x[0]*x[0]*x[1]+2*x[0]-5;
			}
		};
		
		Function w2 = new Function()
		{
			
			@Override
			public double getValue(double[] x, double... c)
			{
				return x[1]*x[1]*x[0]*x[0]+x[1];
			}
		};
		
		Function[] f = {w1,w2};
		System.out.println("Jacobi Matrix");
		Matrix.printMatix(Matrix.trans(JacobianMatrix.getJacobiMatrix(f, new double[]{1,1,1})));
		System.out.println("Solution Newton");
		Vector.print(Newton.getSolution(f,true, new double[]{1,10}));
		
		
		
	}
	
	private static void optiTest() throws Exception
	{
		java.util.function.Function<List<double[]>, Double> func = new java.util.function.Function<List<double[]>, Double>()
		{
			
			@Override
			public Double apply(List<double[]> t)
			{
				double distance = 0;
				for(int i=0;i<t.size()-1;i++)
				{
					double dx = t.get(i+1)[0]-t.get(i)[0];
					double dy = t.get(i+1)[1]-t.get(i)[1];
					
					distance+= Math.sqrt(Math.pow(dx, 2)+Math.pow(dy, 2));
				}
				distance+= Math.sqrt(Math.pow(t.get(0)[0]-t.get(0)[0], 2)+Math.pow(t.get(t.size()-1)[1]-t.get(t.size()-1)[1], 2));
				return distance;
			}
		};
		
		List<double[]> points = new ArrayList<>();
//		points.add(new double[]{2.4,12.5});
//		points.add(new double[]{3.4,4.5});
//		points.add(new double[]{4.4,8.5});
//		points.add(new double[]{5.4,4.7});
//		points.add(new double[]{11.4,4.5});
//		points.add(new double[]{4,5});
//		points.add(new double[]{14,20});
//		points.add(new double[]{6,8});
//		points.add(new double[]{4,5});
//		points.add(new double[]{-1,5});
//		points.add(new double[]{10,12});
//		points.add(new double[]{8,14});
//		points.add(new double[]{3,5});
//		points.add(new double[]{-8,5});
//		points.add(new double[]{13,12});
//		points.add(new double[]{5,14});
		
		points.add(new double[]{12.306686,17.082213});
		points.add(new double[]{17.444981,13.529254});
		points.add(new double[]{5.649415,6.563093});
		points.add(new double[]{13.684960,1.156989});
		points.add(new double[]{12.822402,18.412195});
		points.add(new double[]{12.060511,13.196684});
		points.add(new double[]{16.850759,5.517385});
		points.add(new double[]{2.001419,9.249038});
		points.add(new double[]{4.346122,12.570082});
		points.add(new double[]{9.043659,15.003873});
		points.add(new double[]{13.440834,9.398337});
		points.add(new double[]{6.251404,10.627883});
		points.add(new double[]{5.916549,2.459349});
		points.add(new double[]{0.513691,0.142943});
		points.add(new double[]{19.130911,8.185860});
		points.add(new double[]{19.446356,13.819490});
		points.add(new double[]{3.131098,9.322308});
		points.add(new double[]{5.977465,7.113425});
		points.add(new double[]{6.253169,5.952993});
		points.add(new double[]{7.167164,14.303889});
		points.add(new double[]{18.541345,16.615943});
		points.add(new double[]{1.410488,12.743166});
		points.add(new double[]{16.481461,1.712561});
		points.add(new double[]{3.894367,14.508862});
		points.add(new double[]{14.828020,11.364541});
		points.add(new double[]{4.027674,6.252453});
		points.add(new double[]{1.714436,2.935946});
		points.add(new double[]{0.535872,5.167149});
		points.add(new double[]{11.771163,9.629758});
		points.add(new double[]{10.009422,9.251971});
		points.add(new double[]{18.636642,10.884600});
		points.add(new double[]{12.917801,13.427818});
		points.add(new double[]{3.031072,17.738745});
		points.add(new double[]{7.448638,19.565225});
		points.add(new double[]{15.473228,1.894936});
		points.add(new double[]{6.004097,8.112399});
		points.add(new double[]{3.755671,5.910946});
		points.add(new double[]{11.784338,1.766051});
		points.add(new double[]{14.453872,13.757201});
		points.add(new double[]{13.610517,17.506423});
		points.add(new double[]{5.029910,15.924216});
		points.add(new double[]{15.987438,12.321836});
		points.add(new double[]{18.501572,11.328981});
		points.add(new double[]{18.378642,0.076982});
		points.add(new double[]{16.182842,9.460859});
		points.add(new double[]{12.764374,12.233082});
		points.add(new double[]{15.661908,10.523940});
		points.add(new double[]{1.775505,13.506953});
		points.add(new double[]{0.967986,8.774230});
		points.add(new double[]{11.044165,5.059928});
		points.add(new double[]{19.788410,4.816433});
		
		
		
		
		/* inti 50 random points */
		
//		for(int i=0;i<=50;i++)
//		{
//			double x = Math.random() * 20;
//			double y = Math.random() * 20;
//			
//			points.add(new double[]{x,y});
//			Locale.setDefault(Locale.ENGLISH);
//			System.out.println(String.format("points.add(new double[]{%f,%f});", x,y));
//		}
		
		SimulatedAnnealingTravelSalesMan<double[],Double> travalSalesMan = new SimulatedAnnealingTravelSalesMan<double[], Double>(points,func);
		
		/**
		 * the seriesName of the yData
		 */
		String FIT_VALUE = "Fit Value";
		/**
		 * the chart to plot
		 */
//		XYChart xyChart;
//		
//		/**
//		 * the flag to get the yData an work with them
//		 */
//		SwingWrapper<XYChart> swingWrapper;
//		xyChart = new XYChartBuilder().width(1200).height(400).theme(ChartTheme.XChart).title(travalSalesMan.getClass().getSimpleName()).build();
//		
//		
//		
//		
//		/* make a number of tests */
		
//		travalSalesMan.start();
//			
//		Thread.sleep(100);

//    	xyChart.addSeries(FIT_VALUE, null, travalSalesMan.fitData);
//  		xyChart.getSeriesMap().get(FIT_VALUE).setMarker(SeriesMarkers.NONE).setLineWidth(1);
//	    /* inti SwingWrapper to show the plot */
//		swingWrapper = new SwingWrapper<XYChart>(xyChart);
//		swingWrapper.displayChart();
//		
//		 // Simulate a data feed
//	    TimerTask chartUpdaterTask = new TimerTask() {
//
//			@Override
//			public void run()
//			{
//				xyChart.updateXYSeries(FIT_VALUE, null, travalSalesMan.fitData, null);
//				swingWrapper.repaintChart();
//			}
//	    };
//
//	    Timer timer = new Timer();
//	    timer.scheduleAtFixedRate(chartUpdaterTask, 0, 10);
//
//		travalSalesMan.join();

		
		
		
//		for (int i=0;i<25000;i++)
//		{
//		List<double[]> result = travalSalesMan.findMin();
//		result.add(result.get(0));
//		List<double[]> origin = travalSalesMan.getConfiguration(travalSalesMan.startConf);
//		origin.add(origin.get(0));
		
		/* result config to string */
		
//		String configuration = Arrays.toString(travalSalesMan.findMin());
//			travalSalesMan.findMin();
//			if(travalSalesMan.fitFunc.apply(travalSalesMan.getConfiguration(best))>travalSalesMan.getOptimizedValue())
//				best = travalSalesMan.optConf;
//			
//		verteilung.add(travalSalesMan.getOptimizedValue());

			
//		if(i % 100 ==0)
//			System.out.println(i);
//		
//		
//		
////		System.out.println(String.format("ori distance : %f \t opti distance: %f", travalSalesMan.apply(origin),travalSalesMan.apply(result)));
//		}
		
//		Histogram hist = new Histogram(verteilung, 15);
//		CategoryChart chart = new CategoryChart(1200, 600);
//		chart.addSeries("100 iterations", hist.getxAxisData(), hist.getyAxisData());
//        chart.setTitle("Verteilung");
//        chart.setXAxisTitle("configuration");
//        chart.setYAxisTitle("Number");
//        SwingWrapper<CategoryChart> swingWrapper = new SwingWrapper<>(chart);
//        swingWrapper.displayChart();
//        
////		List<double[]> result = travalSalesMan.getConfiguration(travalSalesMan.optConf);
////		result.add(result.get(0));
//		List<double[]> origin = travalSalesMan.getConfiguration(travalSalesMan.optConf);
//		origin.add(origin.get(0));
////		
//		XYPlot plot = new XYPlot("Travel Sales Man" + " " + travalSalesMan.getOptimizedValue().toString(), "x", "y", 800, 600);
////		plot.addDataSet("opti", result);
//		plot.addDataSet("best",origin);
//		
//		plot.showChart();
//		
//		System.out.println(String.format("best distance : %f", travalSalesMan.fitFunc.apply(travalSalesMan.getConfiguration(travalSalesMan.optConf))));
//		
		
		
		

		Function f1 = new Function() {

			@Override
			public double getValue(double[] x, double... s) {
				return pow(x[0] - 2, 2.) + 4*  pow(x[1] - 3, 2.);
			}
		};
		
		Function f2 = new Function()
		{
			
			@Override
			public double getValue(double[] x, double... c)
			{
				return -pow(x[0],2) + pow(x[1],2);
			}
		} ;
		
		MathOperation g1 = (x, s)  ->  x[0] + 2*x[1] +1;
		MathOperation g2= (x, s)  ->  -x[0] + x[1] - 2;
		MathOperation h1 = (x, s)  ->  x[1]-1;
		MathOperation test = new MathOperation()
		{
			
			@Override
			public double getValue(double[] x, double... c)
			{
				// TODO Auto-generated method stub
				return 04+5;
			}
		};
		
		
		PenaltyMethod penaltyMethod = new PenaltyMethod(himmelblau,1e-8, 10, new DownhillSimplex(null, 1e-8, 50))
				.addRestriction(PenaltyMethod.RestrictionsType.INEQUALITY, 1e1,g1, ((r)  -> r*10. ))
				.addRestriction(PenaltyMethod.RestrictionsType.INEQUALITY, 1e1,g2, ((r)  -> r*10. ));
//				.addRestriction(PenaltyMethod.RestrictionsType.INEQUALITY, 5,(x,c) -> -x[0], ((r)  -> r+2. ))
//				.addRestriction(PenaltyMethod.RestrictionsType.INEQUALITY, 5,(x,c) -> x[0]-1., ((r)  -> r+2. ))
//				.addRestriction(PenaltyMethod.RestrictionsType.INEQUALITY, 5,(x,c) -> x[1]-1., ((r)  -> r+2. ))
//				.addRestriction(PenaltyMethod.RestrictionsType.INEQUALITY, 5,(x,c) -> -x[1]-1., ((r)  -> r+2. ));
//				.addRestriction(PenaltyMethod.RestrictionsType.INEQUALITY, 1,g2, ((r)  -> r*50));
//				.addRestriction(PenaltyMethod.RestrictionsType.EQUATION, 1,h1, ((r)  -> r+20));
		
		
//	
		penaltyMethod.findMin(new double[]{1,5});
		penaltyMethod.showChart(-5, 5, -5, 5, 100, 20);
//		
//		OptiAlgorithm.plot = new XYContourPlot("Contour", "x1", "x2", 800, 600, penaltyMethod.getFunc(), -2, 3,0,4, 0.01, 25);
//		OptiAlgorithm.plot.addDataSet("g1", Arrays.asList(new double[]{-2,0},new double[]{2,4}));
//		OptiAlgorithm.plot.addDataSet("g2", Arrays.asList(new double[]{-2,2.5},new double[]{3,0}));
//		OptiAlgorithm.plot.addDataSet("steps",penaltyMethod.getHist());
//		OptiAlgorithm.plot.showChart();
		

		
	

		// Comparator<double[]> PERSON_LASTNAME_COMPARATOR = new
		// Comparator<double[]>() {
		// @Override public int compare( double[] p1, double[] p2 ) {
		//
		//
		//
		// return p1.lastname.compareTo( p2.lastname );
		// }
		// };
		//
		// TreeSet<double[]> points = new TreeSet<double[]>();

		// for(int a=0;a<numLines;a++)
		// {
		// for(int i=0;i<X.length-1;i++)
		// {
		// for(int j=0;j<Y.length-1;j++)
		// {
		// if(ZZ[i][j] <= Z[a] && ZZ[i][j+1] > Z[a])
		// points.add(new double[]
		// {
		// Y[j]+(Y[j+1]-Y[j])/(ZZ[i][j+1]-ZZ[i][j])*(Z[a]-ZZ[i][j]),
		// X[i]
		// }
		// );
		//
		// if(ZZ[i][j] <= Z[a] && ZZ[i+i][j] > Z[a])
		// points.add(new double[]
		// {
		// X[i]+(X[i+1]-X[i])/(ZZ[i+1][j]-ZZ[i][j])*(Z[a]-ZZ[i][j]),
		// Y[j]
		// }
		// );
		// }
		// }
		// }

		// double temp;
		//
		// System.out.println("--------------------------------------------");
		// System.out.println("--------------Himmelblau--------------------");
		// System.out.println("--------------------------------------------");
		//
//		OptiAlgorithm opti = new DownhillSimplex(himmelblau, 1e-5, 1000, 1.3, 0.9);
//		 opti.findMin(new double[][]{{2,1},{3,2},{4,2}});
		//
		// opti = new SimpleLocalNewton(himmelblau, 1e-5, 1000);
		// opti.findMin(new double[]{3.,3.});
		//
		// opti = new LocalNewton(himmelblau, 1e-5, 1000);
		// opti.findMin(new double[]{3.,3.});
		//
		// opti = new SteamedLocalNewton(himmelblau, 1e-5, 1000);
		// opti.findMin(new double[]{3.,3.});
		//
		// opti = new GradientDescent(himmelblau, 1e-5, 1000,new
		// ExactLineSearch());
		// opti.findMin(new double[]{3.,3.});
		//
		// opti = new GlobalNewton(himmelblau, 1e-5, 1000);
		// opti.findMin(new double[]{3.,3.});
		//
		// System.out.println("--------------------------------------------");
		// System.out.println("--------------Rosenbrock--------------------");
		// System.out.println("--------------------------------------------");
		//
//		opti = new DownhillSimplex(rosenbrock, 1e-5, 1000, 1.3, 0.9);
//		opti.findMin(new double[][] { { 2, 1 }, { 3, 2 }, { 4, 2 } });
//		plot.addDataSet("downhill simplex", opti.getHist());
		//
//		opti = new SimpleLocalNewton(rosenbrock, 1e-5, 1000);
//		opti.findMin(new double[] { 2., 2. });
		// plot.addDataSet("simple local newton", opti.getHist());

//		opti = new LocalNewton(rosenbrock, 1e-5, 1000);
//		opti.findMin(new double[] { 2., 2. });
		// plot.addDataSet("local newton", opti.getHist());

//		opti = new SteamedLocalNewton(rosenbrock, 1e-5, 1000, new ExactLineSearch());
//		opti.findMin(new double[] { 2., 2. });
		// plot.addDataSet("steamed local newton (ELS)", opti.getHist());
//		opti = new SteamedLocalNewton(rosenbrock, 1e-5, 1000);
//		opti.findMin(new double[] { 2., 2. });
		// plot.addDataSet("steamed local newton (Armijo)", opti.getHist());

//		opti = new GradientDescent(rosenbrock, 1e-5, 1000, new ExactLineSearch());
//		opti.findMin(new double[] { 2., 2. });
		// plot.addDataSet("gradient descent (ELS)", opti.getHist());
//		opti = new GradientDescent(rosenbrock, 1e-5, 1000);
//		opti.findMin(new double[] { 2., 2. });
		// plot.addDataSet("gradient descent (Armijo)", opti.getHist());

//		opti = new GlobalNewton(rosenbrock, 1e-5, 1000);
//		opti.findMin(new double[] { 2., 2. });

		// DescentAlgorithm.localNewton(rosenbrock, x0);

		// graph.addSeries("localNewton",
		// DescentAlgorithm.solutionHistory.subList(0,
		// DescentAlgorithm.solutionHistory.size()-1));
		//
		// DescentAlgorithm.simpleLocalNewton(himmelblau,x0);
		// graph.addSeries("simpleLocalNewton",
		// DescentAlgorithm.solutionHistory);
		// DescentAlgorithm.steamedLocalNewton(himmelblau,x0);
		// DescentAlgorithm.gradientDescent(rosenbrock,x0, null);
		// graph.addSeries("gradientDescent", DescentAlgorithm.solutionHistory);
		// DescentAlgorithm.globalNewton(himmelblau, x0, null, 0.1, 2.1, false);
		// graph.addSeries("globalNewton", DescentAlgorithm.solutionHistory);

		// new XYPlot("Solution", "x1", "x2", 800, 400)
		// .addDataSet("history", DescentAlgorithm.solutionHistory)
		// .showChart();



		// File wea09 = new File("WEA09.txt");
		// try {
		// BufferedReader br = new BufferedReader(new FileReader(wea09));
		// String line = br.readLine();
		// while((line = br.readLine())!=null){
		// String[] colums = line.split("\t");
		//
		// allVelocity.add(Double.parseDouble(colums[17].replaceAll(",", ".")));
		// }
		// } catch (IOException e) {
		// e.printStackTrace();
		// }
		//
		// for(int i = 0;i<allVelocity.size();i++){
		// double d = allVelocity.get(i);
		// int idx = 0;
		// if((idx = Velocity.indexOf(d)) != -1)
		// {
		// double temp = AnzVelocity.get(idx) +1.;
		// AnzVelocity.set(idx, temp);
		// //System.out.println(String.format("Index %d , Anzahl %f",
		// idx,temp));
		// }else{
		// Velocity.add(d);
		// AnzVelocity.add(1.0);
		// }
		// }
		// double sum = 0;
		// double sum2 = 0;
		// x = new double[AnzVelocity.size()];
		// y = new double[AnzVelocity.size()];
		// for(int i=0;i<AnzVelocity.size();i++){
		// x[i] = Velocity.get(i);
		// // es muss hier durch die Breite der "bins" geteilt werden
		// // diese beträgt hier 0.1 da die messgenauigkeit 0.1 beträgt
		// y[i] = AnzVelocity.get(i)/allVelocity.size()/0.1;
		// sum+=AnzVelocity.get(i);
		// sum2+=y[i];
		// }
		// System.out.println(String.format("%f = %f", (double)
		// allVelocity.size(),sum));
		// System.out.println(String.format("Summe 2 = %f", sum2));
		//
		//
		//
		// // Testwerte erzeugen ...
		//
		// double[] a = {2,2};
		//
		//
		//// y = new double[x.length];
		//// for(int i=0;i<y.length;i++)
		//// {
		//// y[i] = a[0] * Math.sin(a[1]*x[i])+Math.pow(-1, i)*0.1;
		////
		//// }
		//// x = new double[] {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.,
		//// 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.};
		//// y = new double[]{1.36011718e-02, 2.74876829e-02, 6.69405924e-02,
		//// 1.00135061e-01, 1.33120280e-01, 1.40234739e-01,
		//// 1.41566322e-01, 1.34527954e-01, 1.09779528e-01,
		//// 6.62177329e-02, 3.23004052e-02, 1.80905096e-02,
		//// 8.75040423e-03, 3.67136525e-03, 2.11151059e-03,
		//// 9.32108277e-04, 3.42407122e-04, 1.14135707e-04,
		//// 3.80452358e-05, 1.90226179e-05, 0.00000000e+00,
		//// 1.90226179e-05};
		////
		//// sum2 = 0;
		////
		//// for(double y : y){
		//// sum2+=y;
		//// }
		//// System.out.println(String.format("Summe 2 = %f", sum2));
		//
		// System.out.println(Arrays.toString(x));
		// System.out.println(Arrays.toString(y));
		// System.out.println(x.length);
		// System.out.println(y.length);
		//
		// //Start den Algoritmus
		//
		// a = new double[] {4.2,1.8};
		// double[] residum = new double[a.length];
		// double maxResidum;
		// do{
		//
		// double[][] A = new double[a.length][a.length];
		// double[] r = new double[a.length];
		//
		// for(int i=0;i<x.length;i++)
		// {
		//
		// final double xi = x[i];
		// final double yi = y[i];
		// Df df = new Df() {
		//
		// @Override
		// public double function(double[] x)
		// {
		// return Math.pow( weibull(xi, x[0], x[1])-yi,2);
		// }
		// };
		// Hf hf = new Hf()
		// {
		// @Override
		// public double function(double[] x)
		// {
		// return Math.pow( weibull(xi, x[0], x[1])-yi,2);
		// }
		// };
		//
		// A = Matrix.addition(A, hf.getHf(a));
		// r = Vector.addtion(r, Vector.multiScalar(df.getGradient(a), -1.));
		//
		// }
		//
		// residum = Gauss.getSolution(A, r, false);
		// maxResidum = Vector.oneNorm(residum);
		// System.out.println(String.format("Residum %15f",
		// Math.abs(maxResidum)));
		// Vector.addtion(a, residum);
		//
		//
		// }while(Math.abs(maxResidum)>1e-8);
		//
		//
		// System.out.println(String.format("A = %f k = %f", a[0],a[1]));
		//
		//
		//
		//
		//
		//// double[] x = {4,2,3};
		//// double[] gradf;
		//// gradf = Df1.getGradient(x);
		//// try {
		//// System.out.println(Df1.D1f(x,0,0,0));
		//// } catch (Exception e) {
		//// e.printStackTrace();
		//// }
		////
		////
		//// for (int i=0;i<temp.length;i++){
		//// System.out.println(Arrays.toString(temp[i])+"\n");
		//// }
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "Main []";
	}

	public static double weibull(double v, double A, double k) {
		return k / A * Math.pow(v / A, k - 1) * Math.exp(-Math.pow(v / A, k));
	}

}
