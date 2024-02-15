/*
Copyright (c)
Slowik, A. (Ed.). (2020). Swarm Intelligence Algorithms: A Tutorial (1st ed.). CRC Press. https://doi.org/10.1201/9780429422614
*/


#include <iostream>
#include <math.h>  
#include <ctime>  
#include <limits>  
#include <algorithm>
#include <vector>

using namespace std;


double r() { return (double)(rand() % RAND_MAX) / RAND_MAX; }

double FUN_TSP(vector<int>& ant_path, int ant0, int nE, vector<vector<double>>& Distance)
{
	vector<int> ant(nE + 1);
	for (int i = 0; i < nE; i++)
	{
		ant[i] = ant_path[i];
	}
	ant[nE] = ant0;

	double Dist = 0;
	for (int t = 0; t < nE; t++)
	{
		Dist = Dist + Distance[ant[t + 1]][ant[t]];
	}
	return Dist;
}



int roulette(vector<double>& fitness, int nE)
{
	double partsum = 0, sumfit = 0; int j = 0;

	for(int i = 0; i < nE; i++)
	{
		sumfit = sumfit + fitness[i]; 
	}

	double rn = r() * sumfit;
	
	while (j < nE)
	{
		partsum = partsum + fitness[j];

		if(partsum >= rn)  
		{
			return j; 
		}

		j++; 
	}

	if(j > nE - 1) 
	{ 
		j = nE - 1; 
	}

	return j;
}


int mainAntColony()
{
	cout << "Ant Colony Optimization algorithm demo" << endl;

	srand((int) time(NULL));
	int Max_Iteration = 500;
	int ant_pop = 50;
	double Q = 1;

	//Number of points:
	int nE = 10;
	
	//Input points:
	vector<int> X{1, 4, 10, 30, 11, 23, 28, 30, 15, 11};
	vector<int> Y{42, 36, 10, 46, 30, 44, 8, 38, 19, 4};

	vector<vector<double>> Distance(nE, vector<double>(nE));
	vector<vector<double>> ETA(nE, vector<double>(nE));
	vector<vector<double>> TAU(nE, vector<double>(nE));
	vector<vector<int>> ant(ant_pop, vector<int>(nE));
	vector<double> fitness(ant_pop);
	double best_fit;
	double alpha = 1, beta = 0.01f, rho = 0.02f;
	vector<double> NUM(nE);
	vector<double> Prob(nE);

	vector<int> ant0(nE);
	vector<int> best_ant(nE);
	vector<int> ANT0(nE + 1);
	int END, ix;

	int Iteration = 0;
	for(int i = 0; i < nE - 1; i++)
	{
		for(int j = 0; j < nE; j++)
		{
			Distance[i][j] = sqrt(pow(X[i] - X[j], 2) + pow(Y[i] - Y[j], 2));
			ETA[i][j] = 1 / Distance[i][j];
			if(i == j) 
			{ 
				ETA[i][j] = 0; 
			}
			TAU[i][j] = 0.01 * r();
		}
	}
	
	for(int i = 0; i < ant_pop; i++)
	{
		fitness[i] = std::numeric_limits<double>::max(); 
	}

	best_fit = *std::min_element(fitness.begin(), fitness.end());

	//Initialization
	while (Iteration < Max_Iteration)
	{
	Iteration = Iteration + 1;
	
	for(int k = 0; k < ant_pop; k++)
	{
		for(int s = 0; s < nE; s++)
		{ 
			ant0[s] = -1;
		}

		END=0;
		ant0[END] = rand() % nE;

		for(int l = 1; l < nE; l++)
		{
			for (int s = 0; s < nE; s++)
			{
				if(ant0 [s] > -1)
				{
					END = s;
				}
			}

			ix = ant0[END];
			double DEN = 0;

			for(int s = 0; s < nE; s++)
			{
				NUM[s] = TAU[s][ix] * alpha * ETA[s][ix] * beta;
				DEN = DEN + NUM[s];	
			}

			for(int s = 0; s < nE; s++)
			{
				Prob[s] = NUM[s] / DEN;
			}

			for(int s = 0; s <= END; s++)
			{
				if(ant0[s] > -1)
				{
					Prob[ant0[s]] = 0;
				}
			}
		
			ant0[l] = roulette(Prob, nE);
		}

		fitness[k] = FUN_TSP(ant0, ant0[0], nE, Distance);

		if(fitness[k] < best_fit)
		{
			best_fit= fitness[k];
			for(int s = 0; s < nE; s++)
			{
				best_ant [s] = ant0[s];
			}
		
		}

		for(int s = 0; s < nE; s++)
		{
			ant[k][s] = ant0[s];
			ant0[s] = -1;
		}
	}

	//Update Pheromone
	for(int k = 0; k < ant_pop; k++)
	{
		for(int s = 0; s < nE; s++)
		{
			ANT0[s] = ant[k][s];
		}
		ANT0[nE] = ant[k][0];
		for(int i = 0; i < nE; i++)
		{
			TAU[ANT0[i]][ANT0[i + 1]] = TAU[ANT0[i]][ANT0[i + 1]] + (Q / fitness[k]);
		}	
	}

	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < nE; j++)
		{
			TAU[i][j] = (1 - rho) * TAU[i][j];
		}
	}
	
	cout << best_fit << endl;
	
	}

	//Printing best solution
	cout<<"Best ant:["<<best_ant[0];
	for (int i = 1; i < nE; i++)
	{
		cout << ", " << best_ant[i];
	}
	
	cout<<"] "<<endl;
	cout<<"Best fitness: "<<best_fit <<endl;

	getchar();
	
	
	return 0;
}







