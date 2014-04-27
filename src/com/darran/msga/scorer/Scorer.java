package com.darran.msga.scorer;

import com.darran.msga.population.Sequence;

/**
 * ��������ӿ�
 * @author danran
 * Created on 2014��4��22��
 */
public interface Scorer
{
	/**
	   Returns a score for the 2 sequences, as they are aligned.
	 */
	public float compareSequences(Sequence s1, Sequence s2);
}
