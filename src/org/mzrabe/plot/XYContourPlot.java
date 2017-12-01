package org.mzrabe.plot;

import org.mzrabe.lina.Function;

public class XYContourPlot extends XYPlot
{

	Function func;
	private int numberContourLines;

	public XYContourPlot(String title, String xLable, String yLable, int width, int height, Function func, double XYmin, double XYmax, double step, int LNum)
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

	public XYContourPlot(String title, String xLable, String yLable, int width, int height, Function func, double Xmin, double Xmax, double Ymin, double Ymax, double step, int LNum)
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
	
	public XYContourPlot(String title, String xLable, String yLable, int width, int height, int LNum, double[] X, double[] Y, double[][] Z)
	{
		super(title, xLable, yLable, width, height);
		this.numberContourLines = LNum;
		
		double[] xMinMax = findMinMax(X);
		double[] yMinMax = findMinMax(Y);
		
		this.Xmin = xMinMax[0];
		this.Xmax = xMinMax[1];
		this.Ymin = yMinMax[0];
		this.Ymax = yMinMax[1];
		
		addContour(X, Y, Z,LNum);
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
	
	public void updateContour()
	{
		initContour();
	}

	private void initContour()
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

		log.debug("Xnum: " + ((Xmax - Xmin) / Xstep + 1));
		log.debug("Ynum: " + ((Ymax - Ymin) / Xstep + 1));

		int XNum = (int) ((Xmax - Xmin) / Xstep + 1);
		int YNum = (int) ((Ymax - Ymin) / Xstep + 1);

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

		for (int i = 0; i < Z.length; i++)
		{
			for (int j = 0; j < Z[i].length; j++)
			{
				Z[i][j] = func.getValue(new double[] { X[i], Y[j] });
			}

		}
		
		addContour(X, Y, Z,numberContourLines);
		
//		plotContour(X, Y, Z);
	}
	
    public void addContour(double[] X, double[] Y, double[][] ZZ, Integer numberOfLines)
    {
    	chart.setContour(X, Y, ZZ,numberOfLines);
    	
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
	 */
	public void setNumberOfContourLines( int numberContourLines )
	{
		this.numberContourLines = numberContourLines;
		initContour();
	}

}
