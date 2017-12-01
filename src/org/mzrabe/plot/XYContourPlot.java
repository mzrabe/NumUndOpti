package org.mzrabe.plot;

import org.mzrabe.lina.Function;
import org.mzrabe.lina.MathOperation;

public class XYContourPlot extends XYPlot
{
	
	MathOperation func;
	private int numberContourLines;
	
	private double maxContourValue = Double.MAX_VALUE;
	private double minContourValue = -Double.MAX_VALUE;
	
	public XYContourPlot(String title, String xLable, String yLable, int width, int height, MathOperation func, double XYmin, double XYmax, double step, int LNum) throws Exception
	{
		super(title, xLable, yLable, width, height);
		this.func = func;
		this.Xmin = XYmin;
		this.Xmax = XYmax;
		this.Ymin = XYmin;
		this.Ymax = XYmax;
		this.Xstep = step;
		this.Ystep = step;
		this.numberContourLines = LNum;

		initContour();

	}

	public XYContourPlot(String title, String xLable, String yLable, int width, int height, MathOperation func, double Xmin, double Xmax, double Ymin, double Ymax, double step, int LNum) throws Exception
	{
		super(title, xLable, yLable, width, height);
		this.func = func;
		this.Xmin = Xmin;
		this.Xmax = Xmax;
		this.Ymin = Ymin;
		this.Ymax = Ymax;
		this.Xstep = step;
		this.Ystep = step;
		this.numberContourLines = LNum;

		initContour();

	}
	
	public XYContourPlot(String title, String xLable, String yLable, int width, int height, MathOperation func, double Xmin, double Xmax, double Ymin, double Ymax, int stepPoints, int LNum) throws Exception
	{
		super(title, xLable, yLable, width, height);
		this.func = func;
		this.Xmin = Xmin;
		this.Xmax = Xmax;
		this.Ymin = Ymin;
		this.Ymax = Ymax;
		this.Xstep = Math.abs(Xmin-Xmax) / stepPoints;
		this.Ystep = Math.abs(Ymin-Ymax) / stepPoints;
		this.numberContourLines = LNum;

		initContour();

	}
	
	public XYContourPlot(String title, String xLable, String yLable, int width, int height, MathOperation func, double Xmin, double Xmax, double Ymin, double Ymax, int stepPoints, int LNum, double LMax) throws Exception
	{
		super(title, xLable, yLable, width, height);
		this.func = func;
		this.Xmin = Xmin;
		this.Xmax = Xmax;
		this.Ymin = Ymin;
		this.Ymax = Ymax;
		this.Xstep = Math.abs(Xmin-Xmax) / stepPoints;
		this.Ystep = Math.abs(Ymin-Ymax) / stepPoints;
		this.numberContourLines = LNum;
		this.setMaxContourValue(LMax);

		initContour();

	}
	
	/**
	 * Get a new instance of a xy-plot with contour lines. 
	 * 
	 * @param title - the title of the plot
	 * @param xLable - description of the x label
	 * @param yLable - desctription of the y label
	 * @param width - the width of the chart in pixel
	 * @param height - the height of the chart in pixel
	 * @param LNum - number of contour lines
	 * @param X - the x values
	 * @param Y - the y values
	 * @param Z - the z values as 2D array <br> [<br>[f(x[0],y[0]), ... ,f(x[0],y[-1])]<br> , ... , <br>[f(x[-1],y[0]), ... , f(x[-1],y[-1])]<br>]  
	 * @throws IllegalAccessException
	 */
	public XYContourPlot(String title, String xLable, String yLable, int width, int height, int LNum, double[] X, double[] Y, double[][] Z) throws IllegalAccessException
	{
		super(title, xLable, yLable, width, height);
		this.numberContourLines = LNum;
		
		double[] xMinMax = findMinMax(X);
		double[] yMinMax = findMinMax(Y);
		
		this.Xmin = xMinMax[0];
		this.Xmax = xMinMax[1];
		this.Ymin = yMinMax[0];
		this.Ymax = yMinMax[1];
		
		addContour(X, Y, Z,LNum, getMinContourValue(), getMaxContourValue());
	}
	
	private double[] findMinMax(double[] values)
	{
		double min = Double.MAX_VALUE;
		double max = Double.MIN_VALUE;
		
		for(double value : values)
		{
			if(min > value)
				min = value;
			if(max < value)
				max = value;
		}
		
		return  new double[]{min,max};
	}
	
	public void setContourLineValues(double[] contourLines)
	{
		chart.getContourData().setContourLineValues(contourLines);
	}
	
	public void updateContour() throws Exception
	{
		initContour();
	}

	private void initContour() throws Exception
	{
		double Xrange = Xmax - Xmin;
		double Yrange = Ymax - Ymin;
		double residue = 0;

		if ((residue = (Xrange / Xstep) % 1) != 0)
		{
			Xstep = Xrange / (Math.round(Xrange / Xstep + 0.5));
			residue = Xrange % Xstep;
			log.info("corrected the Xstep to the value: " + Xstep);
			// Xmin += residue/2;
			// Xmax -= residue/2;
			// Xrange -= residue;

		}

		log.debug("Xresidue: " + residue);

		if ((residue = (Yrange / Ystep) % 1) != 0)
		{
			Ystep = Yrange / (Math.round(Yrange / Ystep + 0.5));
			residue = Yrange % Ystep;
			log.info("corrected the Ystep to the value: " + Ystep);
			// Ymin += residue/2;
			// Ymax -= residue/2;
			// Yrange -= residue;
		}

		log.debug("Yresidue: " + residue);

		log.debug("Xnum: " + (Math.abs(Xmax - Xmin) / Xstep + 1));
		log.debug("Ynum: " + (Math.abs(Ymax - Ymin) / Ystep + 1));

		int XNum = (int) (Math.abs(Xmax - Xmin) / Xstep + 1);
		int YNum = (int) (Math.abs(Ymax - Ymin) / Ystep + 1);

		double[][] Z = new double[XNum][YNum];
		double[] X = new double[XNum];
		double[] Y = new double[YNum];

		for (int i = 0; i < XNum; i++)
		{
			X[i] = Xmin + i * Xstep;
		}

		for (int i = 0; i < YNum; i++)
		{
			Y[i] = Ymin + i * Ystep;
		}
		
		double validMaxValue = Double.MIN_VALUE;

		for (int i = 0; i < Z.length; i++)
		{
			for (int j = 0; j < Z[i].length; j++)
			{
//				System.out.println(String.format("------------[%f,%f]-------------", X[i],Y[j]));
				double z = func.getValue(new double[] { X[i], Y[j] });
				
				if(Double.isNaN(z))
					z = 0;
				if(Double.isInfinite(z))
					z = Double.MAX_VALUE;
				
				Z[i][j] = z;
			}

		}
		
		addContour(X, Y, Z,numberContourLines, getMinContourValue(), getMaxContourValue());
		
//		plotContour(X, Y, Z);
	}
	
	/**
	 * @param X - the x values
	 * @param Y - the y values
	 * @param ZZ - the z values
	 * @param numberOfLines - the number of contour lines
	 * @param minContourValue - minimal contour value
	 * @param maxContourValue - maximum contour value
	 * @throws IllegalAccessException - if the minimal contour line value is greater than the maximum contour line value
	 */
    public void addContour(double[] X, double[] Y, double[][] ZZ, Integer numberOfLines, double minContourValue, double maxContourValue) throws IllegalAccessException
    {
    	chart.setContour(X, Y, ZZ,numberOfLines, minContourValue, maxContourValue);
    	
    }

//	private void plotContour(double[] X, double[] Y, double[][] Z)
//	{
//		double Zmin = Double.MAX_VALUE;
//		double Zmax = Double.MIN_VALUE;
//
//		for (int i = 0; i < Z.length; i++)
//		{
//			for (int j = 0; j < Z[i].length; j++)
//			{
//				if (Zmax < Z[i][j])
//					Zmax = Z[i][j];
//				if (Zmin > Z[i][j])
//					Zmin = Z[i][j];
//			}
//
//		}
//
//		 out.printf("Minimum: %f \t Maximum: %f \n", Zmin,Zmax);
//
//		double[] ZLine = new double[LNum];
//
//		for (int i = 0; i < LNum; i++)
//		{
//			ZLine[i] = Zmin + (i + 1) * (double) (Zmax - Zmin) / (LNum);
//			// out.printf("line %d has value %f \n",i,ZLine[i]);
//		}
//
//		int num = 0;
//		for (int a = 0; a < ZLine.length; a++)
//		{
//			Color color = new Color(0, (int) (255 * (1 - (ZLine[a] - Zmin) / (Zmax - Zmin))), (int) (255 * ((ZLine[a] - Zmin) / (Zmax - Zmin))));
//			for (int i = 0; i < X.length - 1; i++)
//			{
//				for (int j = 0; j < Y.length - 1; j++)
//				{
//					List<double[]> points = new ArrayList<>();
//					double lp, rp, tp, bp;
//					lp = (ZLine[a] - Z[i][j]) / (Z[i][j + 1] - Z[i][j]);
//					rp = (ZLine[a] - Z[i + 1][j]) / (Z[i + 1][j + 1] - Z[i + 1][j]);
//					tp = (ZLine[a] - Z[i][j + 1]) / (Z[i + 1][j + 1] - Z[i][j + 1]);
//					bp = (ZLine[a] - Z[i][j]) / (Z[i + 1][j] - Z[i][j]);
//
//					// left
//					if (0 <= lp && lp <= 1)
//					{
//						points.add(new double[] { X[i], Y[j] + (Y[j + 1] - Y[j]) / (Z[i][j + 1] - Z[i][j]) * (ZLine[a] - Z[i][j]) });
//					}
//					// bottom
//					if (0 < bp && bp < 1)
//					{
//						points.add(new double[] { X[i] + (X[i + 1] - X[i]) / (Z[i + 1][j] - Z[i][j]) * (ZLine[a] - Z[i][j]), Y[j] });
//					}
//					// top
//					if (0 < tp && tp < 1)
//					{
//						points.add(new double[] { X[i] + (X[i + 1] - X[i]) / (Z[i + 1][j + 1] - Z[i][j + 1]) * (ZLine[a] - Z[i][j + 1]), Y[j + 1] });
//					}
//					// right
//					if (0 <= rp && rp <= 1)
//					{
//						points.add(new double[] { X[i + 1], Y[j] + (Y[j + 1] - Y[j]) / (Z[i + 1][j + 1] - Z[i + 1][j]) * (ZLine[a] - Z[i + 1][j]) });
//					}
//
//					// out.printf("Fond points in square %d,%d :
//					// %d\n",i,j,points.size());
//					if (points.size() == 2 /* && i == 0 && j == 1 */)
//					{
//
//						// out.printf("point p1 %f,%f \t p2 %f,%f\n",
//						// points.get(0)[0],points.get(0)[1],points.get(1)[0],points.get(1)[1]);
//						// plot.addDataSet("line " + num, points);
//						chart.addSeries("line " + num, new double[] { points.get(0)[0], points.get(1)[0] }, new double[] { points.get(0)[1], points.get(1)[1] });
//						chart.getSeriesMap().get("line " + num).setShowInLegend(false);
//						chart.getSeriesMap().get("line " + num).setMarker(new None());
//						chart.getSeriesMap().get("line " + num).setLineColor(color);
//						chart.getSeriesMap().get("line " + num).setLineWidth(1);
//						chart.getSeriesMap().get("line " + num).setLineStyle(SeriesLines.SOLID);
//					}
//					num++;
//
//				}
//			}
//		}
//	}
	
	/**
	 * Set the number of the contour Lines and update the ContourData of the chart.
	 * @throws Exception 
	 */
	public void setNumberOfContourLines( int numberContourLines ) throws Exception
	{
		this.numberContourLines = numberContourLines;
		initContour();
	}
	
	/**
	 * Get the maximum contour value.
	 * @return - maximum contour value
	 */
	public double getMaxContourValue()
	{
		return maxContourValue;
	}
	/**
	 * Set a maximum contour value to visualize a interested region better
	 * @param maxContourValue - the maximal value
	 */
	public void setMaxContourValue(double maxContourValue)
	{
		this.maxContourValue = maxContourValue;
	}
	
	/**
	 * Get the minimum contour value.
	 * @return - minimum contour value
	 */
	public double getMinContourValue()
	{
		return minContourValue;
	}
	/**
	 * Set a minimum contour value to visualize a interested region better
	 * @param minContourValue - the maximal value
	 */
	public void setMinContourValue(double minContourValue)
	{
		this.minContourValue = minContourValue;
	}

}
