package org.mzrabe.methaheuristic.firefly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.mzrabe.lina.Vector;

public class FlySwarm
{
	/**
	 * the static squad counter to get the id of a new squad
	 */
	public static int nextSwarmID = 0;
	/**
	 * the id of this swarm
	 */
	public int swarmID;
	/**
	 * the index of the best fly (leader)
	 */
	int leaderID;
	/**
	 * the members of flies in the swarm
	 */
	HashSet<Integer> members = new HashSet<>();
	/**
	 * if the flies swarm is already blocked 
	 */
	boolean isLocked = false;
	/**
	 * 
	 */
	public List<double[]> path = new ArrayList<>();
	
	public FlySwarm()
	{
		/* set the id and increment the global id counter */
		swarmID = FlySwarm.nextSwarmID ++;
		isLocked = true;
	}
	
	
	
	public double[] getAverage(List<Fly> flies)
	{
		double[] s = new double[flies.get(0).position.length];
		/* init with zero */
		Arrays.fill(s, 0.0);
		
		/* calculate the average */
		for (Integer i : members)
		{
			s = Vector.sum(s, flies.get(i).position);
		}
		s = Vector.multiScalar(s, (double) 1. / members.size());
		
		return s;
	}
	
	public void updatePath(List<Fly> flies)
	{
		path.add(getAverage(flies));
	}
	
	@Override
	public String toString()
	{
		return "swarmID= " + swarmID + ", leader= " + leaderID + ", number of members= " + members.size() + ", members= " + members.toString()
		+ ", lastSwarmStep= " + lastSwarmStep() + ", position= " + (path.isEmpty() ? "null" : Arrays.toString(path.get(path.size()-1)));
	}


	/**
	 * Set the leader id and add the id of the fly to the member set
	 * @param flyID - the id of the fly which should the leader
	 */
	public void setLeader(int flyID)
	{
		this.leaderID = flyID;
		if(members.add(flyID) == false)
		{
			System.out.println("The leader of the group " + swarmID + " changed. New leader is " + flyID);
		}
	}



	
	public boolean isLocked()
	{
		return isLocked;
	}
	public void lockSwarm()
	{
		this.isLocked = true;
	}
	public void unlockSwarm()
	{
		this.isLocked = false;
	}
	public Double lastSwarmStep()
	{
		if(path.size() < 2)
			return null;
		return Vector.twoNorm(Vector.minus(path.get(path.size()-2),path.get(path.size()-1)));
	}
	
}
