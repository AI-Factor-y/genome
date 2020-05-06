/*
genome.h

a genetic algorithm library made for c++ by abhinav-p

the function to visualise neural network requires latest version of sfml graphics library for c++
uncomment the code section at the end for visualisation functions if sfml is available on you system

Copyright(C) Abhinav-p
*/

#pragma once

#ifndef _GENOME_
#define _GENOME_
/* Uncomment only if sfml library is available and working */
//#include <SFML/Graphics.hpp>   
#include<iostream>
#include<vector>
#include<ctime>
#include<cstdlib>

using namespace std;

const int populationSize = 100;
class dna {
private:
	int random_acc = 1000;
	int input;
	int output;
public:
	int gene_size;
	float mutation_rate;


	vector<float> gene; // gene carrying the weights


	int fitness;
	dna(int inp, int out, float mutation = 0.01) {
		input = inp;
		output = out;
		gene_size = inp * out + out;
		mutation_rate = mutation;
		float random_entry;
		gene = {};
		for (int i = 0; i < gene_size; i++) {
			random_entry = (float)(rand() % (random_acc + 1)) / (float)random_acc;
			gene.push_back(random_entry);
		}


	}

	void init(int input_size, int output_size, float mutate) {
		gene_size = input_size * output_size + output_size;
		mutation_rate = mutate;
	}

	dna crossover(dna pair) {
		//cout << gene_size << endl;
		int midpoint = rand() % (gene_size - 0 + 1);
		dna child(input, output, mutation_rate);
		//child.init(input,output, mutation_rate);  // can be used to reinitialise, not requires in current scenario
		child.gene = {};
		for (int i = 0; i < gene_size; i++) {
			if (i < midpoint) {
				child.gene.push_back(gene[i]);
			}
			else {
				child.gene.push_back(pair.gene[i]);
			}
		}
		return child;

	}

	void mutation() {
		double prob;
		float random_entry;
		for (int i = 0; i < gene_size; i++) {
			prob = (double)(rand() % (random_acc + 1)) / random_acc;

			if (prob < mutation_rate) {
				random_entry = (float)(rand() % (random_acc + 1)) / (float)random_acc;
				gene[i] = random_entry;

			}
		}
	}

	void create_random_species() {

		float random_entry;
		for (int i = 0; i < gene_size; i++) {



			random_entry = (float)(rand() % (random_acc + 1)) / (float)random_acc;
			gene[i] = random_entry;


		}
	}



};
//
class genetic {

	vector<double> inputs;
	int input_size;
	int output_size;
	float mutation_rate;
	int population_size = 100;
	int generation = 1;
	vector<double> curr_input;



public:

	vector<dna> population;

	genetic(int no_inputs, int no_outputs, float mutation = 0.01, long int population_s = populationSize) {

		//seeding the time for random operations
		srand(time(0));
		//create the dna for all the species of the population
		for (int elem = 0; elem < population_s; elem++) {
			population.push_back(dna(no_inputs, no_outputs, mutation));
		}
		input_size = no_inputs;
		output_size = no_outputs;
		mutation_rate = mutation;
		population_size = population_s;

		for (int n = 0; n < input_size; n++) {
			curr_input.push_back(0.0);
		}
		/*for (int i = 0; i < population_s; i++) {
			population[i].init(input_size,output_size, mutation_rate);
		}*/

	}

	vector<float> get_genes(int index) {

		return population[index].gene;
	}

	int start_mating() {
		generation += 1;
		// creating a mating pool for species to reproduce
		vector<dna> matingpool = {};

		for (int i = 0; i < population_size; i++) {
			population[i].fitness = (population[i].fitness % 101);  // normalising the fitness values for the mating pool
			int fitness_per = (int)(population[i].fitness);  // *100 is removed
			//cout << endl << fitness_per << endl;
			for (int j = 0; j < fitness_per; j++) {
				matingpool.push_back(population[i]);
			}

		}

		// reproducing species
		if (matingpool.size() != 0) {            //check whether the mating pool is empty or not
			for (int i = 0; i < population_size; i++) {
				int choice_1 = rand() % (matingpool.size());
				int choice_2 = rand() % (matingpool.size());

				dna parent_1 = matingpool[choice_1];
				dna parent_2 = matingpool[choice_2];

				dna child(input_size, output_size, mutation_rate);
				child = parent_1.crossover(parent_2);
				population[i] = child;

			}
		}
		else {  // if its empty then create a new random population from start
			cout << endl << "warning : no improvement in population , try changing the fitness function";
			for (int i = 0; i < population_size; i++) {
				population[i].create_random_species();
			}
		}

		for (int i = 0; i < population_size; i++) {
			population[i].mutation();
		}

		return generation;
	}

	double  output(vector<double> input, int agent) {

		vector<double> result_output;
		curr_input = input;
		vector<float> gene_data = population[agent].gene;
		double sum = 0;
		/*for (int i = 0; i < input_size-1; i++) {
			sum += input[i] * gene_data[i];
		}
		sum += gene_data[gene_data.size() - 1];  */      // adding the bias--
		int data_pos;
		int k;
		for (int i = 0; i < output_size; i++) {
			sum = 0;
			data_pos = i * (gene_data.size() / output_size);
			sum += gene_data[data_pos];
			k = 0;
			for (int j = 1; j <= input_size; j++) {
				sum += input[k] * gene_data[j + data_pos];
				k += 1;
			}
			result_output.push_back(sum);
		}
		int max_result_index = max_element(result_output.begin(), result_output.end()) - result_output.begin();

		//sum = (double)((int)sum % 101) / 100;  // making the result value between 0 and 1

		return max_result_index;
	}


	void set_fitness(vector<int> fit_arr) {

		for (int i = 0; i < fit_arr.size(); i++) {
			population[i].fitness = fit_arr[i];
		}

	}

	void set_popuation_size(int size) {
		population_size = size;
	}

	void set_mutation_rate(float mutation) {
		mutation_rate = mutation;
	}
	void decay_mutation(float decay_rate) {
		mutation_rate *= decay_rate;
	}
	vector<float> get_weights(int agent) {

		return population[agent].gene;
	}

	int best_fit_agent() {
		int maxfit = population[0].fitness;
		int maxfit_index = 0;
		for (int i = 0; i < population_size; i++) {

			if (population[i].fitness > maxfit) {
				maxfit = population[i].fitness;
				maxfit_index = i;
			}

		}

		return maxfit_index;
	}

	void show_network() {
		int elem = best_fit_agent();
		cout << "\npopulation size : " << population_size << "\t\t mutation rate : " << mutation_rate * 100 << "%";;
		cout << "\nshowing details of the best species of the generation";
		cout << "\n-------------------------------------------------------------------------------------";
		cout << endl << "| output |" << "     | inputs |" << " \t\t\t" << "| bias |" << "\t\t\t" << "| weights |";
		cout << "\n_____________________________________________________________________________________";
		int data_pos;

		for (int i = 0; i < output_size; i++) {
			cout << "\n   " << i << "\t\t";
			for (int k = 0; k < curr_input.size(); k++) {
				cout << curr_input[k] << " ";
			}
			cout << "\t\t  ";
			data_pos = i * (get_weights(elem).size() / output_size);
			cout << get_weights(elem)[data_pos] << "\t\t  ";

			for (int j = 1; j <= input_size; j++) {
				cout << get_weights(elem)[j + data_pos] << " ";

			}
		}
		cout << "\n" << "=====================================================================================";

	}

	int get_generation() {
		return generation;
	}

	vector<vector<float>> all_weights() {
		vector<vector<float>> all;  // 2d array to store weights in  a ordered fashion

		for (int i = 0; i < output_size; i++) {
			all.push_back(get_weights(i));
		}
		return all;
	}

	vector<double> get_inputs() {
		return curr_input;
	}



};
double range_map(double var, double xi, double xf, double Xi, double Xf) {

	double percent = (var - xi) / (xf - xi);

	double new_var = Xi + (Xf - Xi) * percent;

	return new_var;

}
//NOTE
//==================================
// Functions for neural network visualisation 
// works on sfml graphics library of c++
// uncomment if you have sfml graphics library installed in your system
/*

vector<sf::RectangleShape> draw_network_lines(int output_no, int input_no, vector<vector<float>> out_w, vector<float> positions = { 400,100,500,130,50,50 }) {
	vector<sf::CircleShape> outline_circles;

	vector<float> inp_nodes;
	const float PI = 3.14159265;
	vector<float> out_nodes;
	vector<sf::CircleShape> inp_img(input_no);
	vector<sf::CircleShape> out_img(output_no);
	for (int i = 0; i < input_no; i++) {
		inp_nodes.push_back(positions[1] + positions[4] * i);
	}
	for (int i = 0; i < output_no; i++) {
		out_nodes.push_back(positions[3] + positions[5] * i);
	}
	for (int i = 0; i < input_no; i++) {
		inp_img[i].setRadius(10);
		inp_img[i].setFillColor(sf::Color::Blue);
		inp_img[i].setPosition(positions[0], inp_nodes[i]);
		outline_circles.push_back(inp_img[i]);

	}
	for (int i = 0; i < output_no; i++) {
		out_img[i].setRadius(10);
		out_img[i].setFillColor(sf::Color::Blue);
		out_img[i].setPosition(positions[2], out_nodes[i]);
		outline_circles.push_back(out_img[i]);
	}

	float dist;
	float angle, slope, rad, weight;
	int color;
	vector<sf::RectangleShape> net_lines;

	for (int i = 0; i < output_no; i++) {
		for (int j = 0; j < input_no; j++) {

			dist = pow((positions[0] - positions[2]) * (positions[0] - positions[2]) + pow((inp_nodes[j] - out_nodes[i]), 2), 0.5);
			weight = out_w[i][j + 1];
			weight = range_map(weight, 0, 1, 1, 5);
			color = range_map(weight, 0, 1, 60, 255);
			sf::RectangleShape line({ dist,weight });
			line.setFillColor(sf::Color(color, 100, 100));
			line.setPosition({ positions[2],out_nodes[i] });
			slope = ((inp_nodes[j] - out_nodes[i]) / (positions[0] - positions[2]));
			rad = atan(slope);
			angle = rad * 180 / PI;

			line.rotate(angle - 180);
			net_lines.push_back(line);
		}
	}


	return net_lines;

}
vector<sf::CircleShape> draw_network(int output_no, int input_no, vector<vector<float>> out_w, vector<float> positions = { 400,100,500,130,50,50 }) {
	vector<sf::CircleShape> outline_circles;

	vector<float> inp_nodes;
	const float PI = 3.14159265;
	vector<float> out_nodes;
	vector<sf::CircleShape> inp_img(input_no);
	vector<sf::CircleShape> out_img(output_no);
	for (int i = 0; i < input_no; i++) {
		inp_nodes.push_back(positions[1] + positions[4] * i);
	}
	for (int i = 0; i < output_no; i++) {
		out_nodes.push_back(positions[3] + positions[5] * i);
	}
	for (int i = 0; i < input_no; i++) {
		inp_img[i].setRadius(10);
		inp_img[i].setFillColor(sf::Color::Blue);
		inp_img[i].setPosition(positions[0] - 10, inp_nodes[i] - 10);
		outline_circles.push_back(inp_img[i]);

	}
	for (int i = 0; i < output_no; i++) {
		out_img[i].setRadius(10);
		out_img[i].setFillColor(sf::Color::Blue);
		out_img[i].setPosition(positions[2] - 10, out_nodes[i] - 10);
		outline_circles.push_back(out_img[i]);
	}

	float dist;
	float angle, slope, rad, weight;
	int color;
	vector<sf::RectangleShape> net_lines;




	return outline_circles;

}


*/

#endif

