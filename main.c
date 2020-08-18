#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CHANNEL_NUM 32
#define SAMPLE_LENGTH 400000
#define SAMPLE_RATE 1000000
#define SIGNAL_FREQ 25000
#define ECHO_POINT 1000
#define SND_SPD 340
#define PI 3.1415926

#define mic_d 0.03
#define r_min 35
#define r_max 40
#define r_step 0.1
#define t_width 20*PI/180
#define t_step 0.5*PI/180
#define grid_h -0.1
#define trans_x 0
#define trans_y 0
#define trans_z 0
#define r_size 50
#define t_size 40
#define g_size 2000

void CoordinateInit(); 
void DataProcess();

float AD_Sample[CHANNEL_NUM][SAMPLE_LENGTH];
float grid[g_size][3];
float mic[CHANNEL_NUM][3];
int CutIndex[g_size][CHANNEL_NUM];
float FRE_Amp[g_size];
float sum[ECHO_POINT];

int main(void) {
	int i, j;
	//static float data[CHANNEL_NUM][SAMPLE_LENGTH];

	FILE * fp;
	
	fp = fopen ("E:\\chengxiang.txt", "r");
	for(i=0; i<SAMPLE_LENGTH; i++)
	{
		for(j=0; j<CHANNEL_NUM-1; j++)
		{
			fscanf(fp, "%E\t", &AD_Sample[j][i]);
		}
		fscanf(fp, "%E\r\n", &AD_Sample[CHANNEL_NUM-1][i]);
	}
	fclose(fp);
	
	CoordinateInit();
	DataProcess();
	
	fp = fopen ("Amp.txt", "w+");
	for(j = 0; j < 50; j++)
	{
		for(i = 0; i < 39; i++)
		{
			fprintf(fp, "%f\t", FRE_Amp[j*40 + i]);
		}
		fprintf(fp, "%f\r\n", FRE_Amp[j*40 + 39]);
	}
	
//	for(i = 0; i < g_size; i++) 
//	{
//		for(j = 0; j < CHANNEL_NUM; j++)
//		fprintf(fp, "%d\t", CutIndex[i][j]);
//		fprintf(fp, "\n");
//	}

	fclose(fp);
	
	return 0;
}

void CoordinateInit()
{
	int i, j, k;

	/* Mic Init */
	for(i = 0; i < CHANNEL_NUM; i++)
	{
		mic[i][0]=(i - CHANNEL_NUM/2)*mic_d;
		mic[i][1]=0.0556;
		mic[i][2]=0.0035;
	}
	
	/* Grid Init */
	for(i = 0; i < r_size; i++)
	{
		for(j = 0; j < t_size; j++)
		{
			grid[t_size*i + j][0]=(r_min + r_step*i) * sin((0.5*j - 10)*PI/180);
			grid[t_size*i + j][1]=(r_min + r_step*i) * cos((0.5*j - 10)*PI/180);
			grid[t_size*i + j][2]=-0.1;
		}
	}
	
	/* CutIndexCal */
	for(i = 0; i < CHANNEL_NUM; i++)
	{
		for(j = 0; j < g_size; j++)
		{
			CutIndex[j][i] = (sqrt( pow((grid[j][0] - mic[i][0]), 2)  +\
								    pow((grid[j][1] - mic[i][1]), 2)  +\
								    pow((grid[j][2] - mic[i][2]), 2)) +\
							  sqrt( pow((grid[j][0] - trans_x  ), 2)  +\
								    pow((grid[j][1] - trans_y  ), 2)  +\
								    pow((grid[j][2] - trans_z  ), 2)))/\
							 SND_SPD*SAMPLE_RATE;
		}
	}
}

void DataProcess()
{
	float rel_L, img_L, rel_M, img_M, rel_H, img_H;
	int i, k, j;

	for(i = 0; i < g_size; i++)
	{
		for(j = 0; j < ECHO_POINT; j++) sum[j] = 0;
		rel_L = img_L = rel_M = img_M = rel_H = img_H = 0.0;

		for(j = 0; j < ECHO_POINT; j++)
			for(k = 0; k < CHANNEL_NUM; k++)
				sum[j] += AD_Sample[k][CutIndex[i][k] + j];
		
		for(j = 0; j < ECHO_POINT; j++)
		{
			//rel_L += sum[j] * cos(2*PI*(SIGNAL_FREQ - 1)/SAMPLE_RATE * j);
			rel_M += sum[j] * cos(2*PI*(SIGNAL_FREQ    )/SAMPLE_RATE * j);
			//rel_H += sum[j] * cos(2*PI*(SIGNAL_FREQ + 1)/SAMPLE_RATE * j);
			//img_L -= sum[j] * sin(2*PI*(SIGNAL_FREQ - 1)/SAMPLE_RATE * j);
			img_M -= sum[j] * sin(2*PI*(SIGNAL_FREQ    )/SAMPLE_RATE * j);
			//img_H -= sum[j] * sin(2*PI*(SIGNAL_FREQ + 1)/SAMPLE_RATE * j);
		}
		//FRE_Amp[i] = sqrt(pow(rel_L, 2) + pow(img_L, 2)) +\
		//	sqrt(pow(rel_M, 2) + pow(img_M, 2)) +\
		//	sqrt(pow(rel_H, 2) + pow(img_H, 2));
		FRE_Amp[i] = sqrt(pow(rel_M, 2) + pow(img_M, 2));
	}
}
