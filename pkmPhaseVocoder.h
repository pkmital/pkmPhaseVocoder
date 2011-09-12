/*
 *  pkmPhaseVocoder.h
 *  Phase Vocoder using the Accelerate.framework
 
 *  Created by Parag K. Mital - http://pkmital.com 
 *  Contact: parag@pkmital.com
 *
 *  Copyright 2011 Parag K. Mital. All rights reserved.
 * 
 *	Permission is hereby granted, free of charge, to any person
 *	obtaining a copy of this software and associated documentation
 *	files (the "Software"), to deal in the Software without
 *	restriction, including without limitation the rights to use,
 *	copy, modify, merge, publish, distribute, sublicense, and/or sell
 *	copies of the Software, and to permit persons to whom the
 *	Software is furnished to do so, subject to the following
 *	conditions:
 *	
 *	The above copyright notice and this permission notice shall be
 *	included in all copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,	
 *	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *	OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *	WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *	OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#pragma once

#include <Accelerate/Accelerate.h>
#include "pkmFFT.h"

#ifndef M_2_PI
#define M_2_PI 6.283185307179586
#endif

class pkmPhaseVocoder
{
public:
	pkmPhaseVocoder(int audioFrameSize = 512)
	{
		binSize = audioFrameSize / 2;
		
		magnitude_frame = (float *)malloc(sizeof(float) * binSize);
		phase_frame = (float *)malloc(sizeof(float) * binSize);
		prev_phase_frame = (float *)malloc(sizeof(float) * binSize);
		phase_diff = (float *)malloc(sizeof(float) * binSize);
		
		// expected phase difference
		expected_phase_diff = (float *)malloc(sizeof(float) * binSize);
		float a = 0;
		float b = 1;
		vDSP_vramp(&a, &b, expected_phase_diff, 1, binSize);
		float expct = M_2_PI/(float)binSize;
		vDSP_vsmul(expected_phase_diff, 1, &expct, expected_phase_diff, 1, binSize);
		
		memset(prev_phase_frame, 0, sizeof(float)*binSize);
		memset(phase_diff, 0, sizeof(float)*binSize);
		
		fft = new pkmFFT(audioFrameSize);
		
	}
	
	~pkmPhaseVocoder()
	{
		free(magnitude_frame);
		free(phase_frame);
		free(prev_phase_frame);
		free(phase_diff);
		free(expected_phase_diff);
		delete fft;
	}
	
	void correctPhaseInPlace(float *sample_data)
	{
		fft->forward(0, sample_data, magnitude_frame, phase_frame, false);
		correctPhase(phase_frame);
		fft->inverse(0, sample_data, magnitude_frame, phase_frame, false);
	}
	
	void correctPhaseOutOfPlace(float *sample_data_in, float *sample_data_out)
	{
		fft->forward(0, sample_data_in, magnitude_frame, phase_frame, false);
		correctPhase(phase_frame);
		fft->inverse(0, sample_data_out, magnitude_frame, phase_frame, false);
	}
	
	void correctPhase(float *current_phase_frame)
	{	
		// get phase diff
		vDSP_vsub(current_phase_frame, 1, prev_phase_frame, 1, phase_diff, 1, binSize);
		// store prev phase diff
		cblas_scopy(binSize, current_phase_frame, 1, prev_phase_frame, 1);
		
		// subtract expected phase difference
		vDSP_vadd(phase_diff, 1, expected_phase_diff, 1, phase_diff, 1, binSize);
		
		// wrap to +/- PI and accumulate phase diff for reconstruction
		int i = 1;
		current_phase_frame[0] = phase_diff[0];
		while (i < binSize) {
			
			float *tmp = phase_diff + i;
			*tmp = *tmp - M_2_PI * roundf(*tmp/(M_2_PI));
			/*
			int qpd = *tmp/M_PI;
			if (qpd >= 0) 
				qpd += qpd&1;
			else 
				qpd -= qpd&1;
			*tmp -= M_PI*(float)qpd;
			*/
			
			current_phase_frame[i] = *tmp + current_phase_frame[i];
			i++;
		}
		
	}
	
	pkmFFT *fft;
	int binSize;
	float *expected_phase_diff, *prev_phase_frame, *phase_diff, *phase_frame, *magnitude_frame;

};