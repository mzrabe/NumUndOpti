package org.mzrabe.opti;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Locale;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.mzrabe.lina.Function;
import org.mzrabe.lina.Vector;
import org.mzrabe.plot.XYContourPlot;
import org.mzrabe.plot.XYPlot;

/**
 * 
 * @author Moritz Zahn <mzrabe@gmail.com>
 *
 */
public abstract class OptiAlgorithm {

	
	protected static final Logger log = LogManager.getLogger("OptiAlgo");
	/**
	 * the tolerance for algorithm
	 */
	protected double tolerance;
	/**
	 * the allowed maximal iterations
	 */
	protected int maxNumIterations;
	/**
	 * the number of iterations to find the minimum
	 */
	protected int numberOfIterations = 0;
	/**
	 * the function where the minimum is searched
	 */
	protected Function func;
	/**
	 * boolean if all steps of the iteration should save
	 */
	protected boolean savesSteps = true;
	/**
	 * the list which contains all iterations points, if the value savesSteps is true
	 */
	protected ArrayList<double[]> hist = new ArrayList<>();
	/**
	 * the solution of the algorithms
	 */
	protected double[] solution;
	
	public XYContourPlot plot; 
	
	/**
	 * empty Constructor
	 */
	public OptiAlgorithm(){}
	
	/**
	 * Initialize the optimization algorithm. This is the abstract class.
	 * 
	 * @param tolerance
	 *            - the tolerance the function value
	 * @param maxNumIteration
	 *            - the maximal number of iterations
	 * @param points
	 *            - the start points for the downhill simplex algorithm, the
	 *            number of points must be greater by 1 then the dimension of
	 *            the points.
	 */
	public OptiAlgorithm(double tolerance, int maxNumIteration)
	{
		this(null, tolerance, maxNumIteration);
	}
	
	/**
	 * Initialize the optimization algorithm. This is the abstract class.
	 * 
	 * @param func
	 *            - the function where a minimum is searched
	 * @param tolerance
	 *            - the tolerance the function value
	 * @param maxNumIteration
	 *            - the maximal number of iterations
	 * @param points
	 *            - the start points for the downhill simplex algorithm, the
	 *            number of points must be greater by 1 then the dimension of
	 *            the points.
	 */
	public OptiAlgorithm(Function func,double tolerance, int maxNumIteration)
	{
		this.func = func;
		setTolerance(tolerance);
		setMaxNumIterations(maxNumIteration);
	}
	
	/**
	 * Init the all needed parameter for the optimization algorithms
	 * @param x - the start point(s)
	 * @throws Exception 
	 */
	protected abstract void initAlgorithms(double[] ... x) throws Exception;
	/**
	 * Define the optimization algorithms which will run until the method {@link #isFinish()} returns true.
	 * @throws Exception 
	 */
	protected abstract void algorithms() throws Exception;
	
	/**
	 * Find the minimum of the given Function
	 * @param x 
	 * @return the solution of the optimization algorithms as point of type double[]
	 * @throws Exception 
	 */
	public double[] findMin(double[] ... x) throws Exception {
		hist.clear();
		log.info("Start "+getName()+" optimization ...");
		Date startTime = new Date();
		numberOfIterations = 0;
		initAlgorithms(x);
		addPointsToHistory(getSolution());
//		if(OptiAlgorithm.plot.chart.getSeriesMap().containsKey(getName()) == false)
//			OptiAlgorithm.plot.addDataSet(getName(), getHist());
		
		while(!isFinish() && numberOfIterations < maxNumIterations)
		{
			numberOfIterations++;
			log.debug("------------- numberOfIterations " + numberOfIterations +" -------------");
			algorithms();
			addPointsToHistory(getSolution());
//			OptiAlgorithm.plot.updateDataSet(getName(), getHist());
//			OptiAlgorithm.plot.updateContour();
//			Thread.sleep(100);
//			OptiAlgorithm.plot.swingWrapper.repaintChart();
		}
		
		
		Date endTime = new Date();
		if(numberOfIterations >= maxNumIterations)
			log.warn("The maximal number of iterations was reached.");
		log.info(getName() + " algorithms found minimum after "+numberOfIterations+" iterations at " + Arrays.toString(getSolution()));
		long time = (endTime.getTime()-startTime.getTime());
		
		
		long ms = time % 1000;
		time = (time - ms)/1000;
		long s = time % 60; 
		time = (time - s)/60;
		long min = time % 60;
		time = (time - min)/60;
		long h = time;
		
		log.info(String.format(Locale.ENGLISH, "Elapsed time %d h, %d min, %d s and %d ms (%d ms)", h,min,s,ms,(endTime.getTime()-startTime.getTime())));
		
		
		return getSolution();
	}
	
	/**
	 * Checks if the tolerance of the maximal number of iterations is reached.
	 * @return true if the tolerance of the maximal number of iterations is reached otherwise false
	 * @throws Exception 
	 */
	protected abstract boolean isFinish() throws Exception;
	/**
	 * Get the Solution
	 * @throws Exception 
	 */
	public abstract double[] getSolution() throws Exception;
	
	/**
	 * Get the name of the algorithms
	 */
	public abstract String getName();
	
	
	/**
	 * @return the {@link #savesSteps}
	 */
	public boolean isSavesSteps() {
		return savesSteps;
	}

	/**
	 * @return the {@link #tolerance}
	 */
	public double getTolerance() {
		return tolerance;
	}

	/**
	 * @return the {@link #maxNumIterations}
	 */
	public int getMaxNumIterations() {
		return maxNumIterations;
	}

	/**
	 * @return the {@link #numberOfIterations}
	 */
	public int getNumberOfIterations() {
		return numberOfIterations;
	}

	/**
	 * @return the {@link #func}
	 */
	public Function getFunc() {
		return func;
	}

	/**
	 * @return the {@link #hist}
	 */
	public ArrayList<double[]> getHist() {
		return hist;
	}
	
	/**
	 * @param savesSteps the {@link #savesSteps} to set
	 */
	public void setSavesSteps(boolean savesSteps) {
		this.savesSteps = savesSteps;
	}

	/**
	 * @param tolerance the {@link #tolerance} to set
	 */
	public void setTolerance(double tolerance) {
		if (tolerance > 0) {
			this.tolerance = tolerance;
		} else {
			throw new IllegalArgumentException("The value 'tolerance' have to be greater than 0. The value is " + tolerance);
		}
	}

	/**
	 * @param maxNumIterations the {@link #maxNumIterations} to set
	 */
	public void setMaxNumIterations(int maxNumIterations) {
		if (maxNumIterations > 0) {
			this.maxNumIterations = maxNumIterations;
		} else {
			throw new IllegalArgumentException("The value 'maxNumIterations' have to be greater than 0. The value is " + maxNumIterations);
		}
	}
	
	/**
	 * makes a hard copy of the array of points
	 */
	protected void addPointsToHistory(double[] point)
	{
		if(isSavesSteps())
		{
			hist.add(Arrays.copyOf(point, point.length));
		}
	}

	public void setFunc(Function func)
	{
		this.func = func;
		
	}
	
	public void showChart(double minX, double maxX, double minY, double maxY, int step, int lines)
	{
		
		try
		{
			if(getSolution().length == 2)
			{
				plot = new XYContourPlot("Solution Plot - " +  getName() + ", x= " + Vector.asString(getSolution()), "x1", "x2", 800, 600, getFunc(), minX, maxX, minY, maxY, step, lines);
				plot.addDataSet("solution steps", getHist());
				plot.addDataSet("solution", Arrays.asList(getSolution()));
				plot.showChart();
			}
			else
			{
				
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
}
