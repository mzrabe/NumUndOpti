package org.mzrabe.plot;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.knowm.xchart.BitmapEncoder.BitmapFormat;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.style.Styler.ChartTheme;
import org.mzrabe.lina.Function;

public class XYPlot {
	
	public XYChart chart;
	public SwingWrapper<XYChart> swingWrapper;
	double Xmin, Xmax, Ymin, Ymax, Xstep, Ystep;
	protected static final Logger log = LogManager.getRootLogger();
    
    public XYPlot(String title, String xLable, String yLable, int width, int height)
    {
        chart = new XYChartBuilder().width(width).height(height).theme(ChartTheme.XChart).build();
        chart.setTitle(title);
        chart.setXAxisTitle(xLable);
        chart.setYAxisTitle(yLable);
        swingWrapper = new SwingWrapper<XYChart>(chart);
    }
    
    public XYPlot addDataSet(String title, double[] X, double[] Y)
    {
    	chart.addSeries(title, X, Y);
    	return this;
    }
    public XYPlot addDataSet(String title,List<double[]> points)
    {
    	int length = points.size();
    	double[] X = new double[length];
    	double[] Y = new double[length];
    	
    	for(int i= 0;i<length;i++)
    	{
    		//TODO checken ob das double array.size() == 2 ist
    		X[i] = points.get(i)[0];
    		Y[i] = points.get(i)[1];
    	}
    	
    	chart.addSeries(title, X, Y);
		chart.getStyler().setXAxisMin(Xmin);
	    chart.getStyler().setXAxisMax(Xmax);
	    chart.getStyler().setYAxisMin(Ymin);
	    chart.getStyler().setYAxisMax(Ymax);
    	return this;
    }
    
    public XYPlot updateDataSet(String title,List<double[]> points)
    {
    	int length = points.size();
    	double[] X = new double[length];
    	double[] Y = new double[length];
    	
    	for(int i= 0;i<length;i++)
    	{
    		//TODO checken ob das double array.size() == 2 ist
    		X[i] = points.get(i)[0];
    		Y[i] = points.get(i)[1];
    	}
    	
    	chart.updateXYSeries(title, X, Y,null);
		chart.getStyler().setXAxisMin(Xmin);
	    chart.getStyler().setXAxisMax(Xmax);
	    chart.getStyler().setYAxisMin(Ymin);
	    chart.getStyler().setYAxisMax(Ymax);
    	return this;
    }
    
    public XYPlot addFunction(String title, Function f, double from_x, double to_x)
    {
    	double num = 10;
    	double step = (to_x-from_x)/ (num-1);
    	
    	double[] X = new double[(int) num];
    	double[] Y = new double[(int) num];
    	
    	for(int i=0;i<num;i++)
    	{
    		X[i] = from_x + i * step;
    		Y[i] = f.getValue(new double[]{X[i]});
    	}
    	
    	addDataSet(title, X, Y);
    	return this;
    	
    }
    
    
    
//    public void addContour(double[] X, double[] Y, double[][] ZZ, Integer numberOfLines)
//    {
//    	chart.setContour(X, Y, ZZ,numberOfLines);
//    	
//    }
    
    public void showChart()
    {
        // Show it
        
        swingWrapper.displayChart();
    }
    
    /**
     * Saves the chart as the given image type.
     * @param imageType - image type like png, jpg ect. from the enumeration {@code BitmapFormat}
     * @param path - the directory to 
     * @return
     */
    public boolean saveAs(BitmapFormat imageType, String path)
    {
    	//BitmapEncoder.saveBitmap(chart, "./Sample_Chart", BitmapFormat.PNG);
    	return false;
    }
 
    

}
