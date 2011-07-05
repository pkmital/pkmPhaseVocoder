/*
 *  pkmPhaseVocoder.h
 *  memoryMosaic
 *
 *  Created by Mr. Magoo on 7/5/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#pragma once

#include <Accelerate/Accelerate.h>

class pkmPhaseVocoder
{
public:
	pkmPhaseVocoder(int fS = 512)
	{
		frameSize = fS;
		
		prev_phase_frame = (float *)malloc(sizeof(float) * frameSize);
		phase_diff = (float *)malloc(sizeof(float) * frameSize);
		
		// expected phase difference
		expected_phase_diff = (float *)malloc(sizeof(float) * frameSize);
		float a = 0;
		float b = 1;
		vDSP_vramp(&a, &b, expected_phase_diff, 1, frameSize);
		float expct = -2.0f*M_PI*1.0f/(float)frameSize;
		vDSP_vsmul(expected_phase_diff, 1, &expct, expected_phase_diff, 1, frameSize);
		
		memset(prev_phase_frame, 0, sizeof(float)*frameSize);
		memset(phase_diff, 0, sizeof(float)*frameSize);
		
	}
	
	void correctPhase(float *current_phase_frame)
	{
		// phase vocoder:
		
		// get phase diff
		vDSP_vsub(current_phase_frame, 1, prev_phase_frame, 1, phase_diff, 1, frameSize);
		// store prev phase diff
		cblas_scopy(frameSize, current_phase_frame, 1, prev_phase_frame, 1);
		
		// subtract expected phase difference
		vDSP_vadd(phase_diff, 1, expected_phase_diff, 1, phase_diff, 1, frameSize);
		
		// wrap to +/- PI and accumulate phase diff for reconstruction
		int i = 1;
		current_phase_frame[0] = phase_diff[0];
		while (i < frameSize) {
			float *tmp = phase_diff + i;
			int qpd = *tmp/M_PI;
			if (qpd >= 0) 
				qpd += qpd&1;
			else 
				qpd -= qpd&1;
			*tmp -= M_PI*(float)qpd;
			
			current_phase_frame[i] = phase_diff[i] + current_phase_frame[i-1];
			
			i++;
		}
		
	}
	
	int frameSize;
	float *expected_phase_diff, *prev_phase_frame, *phase_diff;

};