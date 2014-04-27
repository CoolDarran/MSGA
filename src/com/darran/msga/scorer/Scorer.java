package com.darran.msga.scorer;

import com.darran.msga.population.Sequence;

/**
 * 分数计算接口
 * @author danran
 * Created on 2014年4月22日
 */
public interface Scorer
{
	/**
	   Returns a score for the 2 sequences, as they are aligned.
	 */
	public float compareSequences(Sequence s1, Sequence s2);
}
